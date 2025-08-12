# run_all.R — Minimal reproducible pipeline (exports: model_results + uncertainty + lab importance)
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(lubridate)
  library(survival);  library(ranger)
  library(DBI);       library(duckdb)
})

set.seed(42)

# ---------------------- Config ----------------------
config <- list(
  sample_size = 1000,   # NA = full run; set to an integer (e.g., 1000) for a quick test
  min_age = 18,
  lab_window_h = 24,
  rf_min_coverage_pct = 5,         # labs seen in ≥ X% for RF
  rf_max_n = 200000,               # cap for RF training
  n_features_list = c(10,15,20,25,30,35,40,50,60,70),
  feature_floors   = c(0,5,10,15,20),
  n_random_sets    = 20,
  include_observations = TRUE,     # vital signs as additional "tests" (sensitivity)
  
  # ---- file paths (edit these) ----
  ed_paths = list(
    edstays    = "edstays.csv.gz",
    vitalsigns = "vitalsign.csv.gz",
    triage     = "triage.csv.gz",
    medrecon   = "medrecon.csv.gz"
  ),
  mimic_paths = list(
    admissions = "../admissions.csv.gz",
    patients = "../patients.csv.gz", 
    labevents = "../labevents.csv.gz",
    d_labitems = "../d_labitems.csv.gz"
  ),
  out_dir = "outputs"
)

dir.create(file.path(config$out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", sprintf(...), "\n")
parse_datetime <- function(x){ if (inherits(x, c("POSIXct","Date"))) return(x); ymd_hms(x, quiet=TRUE) %||% ymd(x, quiet=TRUE) }
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------- Load cohort ----------------------
load_cohort <- function() {
  log("Loading ED stays + patients/admissions")
  ed <- read_csv(config$ed_paths$edstays, show_col_types = FALSE) %>%
    mutate(across(any_of(c("intime","outtime")), parse_datetime))
  pats <- read_csv(config$mimic_paths$patients, show_col_types = FALSE) %>%
    mutate(dod = parse_datetime(dod))
  adm  <- read_csv(config$mimic_paths$admissions, show_col_types = FALSE) %>%
    mutate(across(any_of(c("admittime","dischtime","deathtime")), parse_datetime))
  
  cohort <- ed %>%
    inner_join(pats %>% select(subject_id, anchor_age, anchor_year), by="subject_id") %>%
    mutate(age_at_ed = anchor_age + (year(intime) - anchor_year)) %>%
    filter(age_at_ed >= config$min_age, age_at_ed <= 120)
  
  list(edstays = cohort, patients = pats, admissions = adm)
}

# ---------------------- Supplementary ED data ----------------------
load_supplementary <- function(stay_ids) {
  safe <- function(path){
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || !"stay_id" %in% names(df)) return(NULL)
    # parse likely time cols
    time_cols <- names(df)[grepl("time|date", names(df), ignore.case = TRUE)]
    for (tc in time_cols) if (is.character(df[[tc]])) df[[tc]] <- parse_datetime(df[[tc]])
    df %>% filter(stay_id %in% stay_ids)
  }
  list(vitals = safe(config$ed_paths$vitalsigns),
       triage = safe(config$ed_paths$triage),
       medrecon = safe(config$ed_paths$medrecon))
}

# ---------------------- Labs (0–24h) via DuckDB ----------------------
extract_labs_window <- function(cohort){
  log("Extracting labs within 0–%dh via DuckDB", config$lab_window_h)
  dbdir <- tempfile("duckdb_filab")
  con <- dbConnect(duckdb::duckdb(), dbdir = dbdir)
  on.exit({ try(dbDisconnect(con, shutdown=TRUE), silent=TRUE); unlink(dbdir, recursive=TRUE) }, add=TRUE)
  
  dbWriteTable(con, "cohort_stays",
               cohort$edstays %>% select(subject_id, stay_id, intime), overwrite=TRUE)
  
  # load labevents using DuckDB's CSV reader
  lp <- normalizePath(config$mimic_paths$labevents, winslash = "/", mustWork = TRUE)
  
  dbExecute(
    con,
    "CREATE TABLE labevents AS
     SELECT * FROM read_csv_auto(?, sample_size=10000, ignore_errors=true)",
    params = list(lp)
  )
  
  qry <- sprintf("
    SELECT l.subject_id, l.hadm_id, l.itemid, l.charttime, l.valuenum, l.flag,
           c.stay_id, c.intime,
           (EXTRACT(EPOCH FROM (l.charttime - c.intime))/3600.0) AS hours_from_ed
    FROM labevents l
    JOIN cohort_stays c USING(subject_id)
    WHERE (EXTRACT(EPOCH FROM (l.charttime - c.intime))/3600.0) BETWEEN 0 AND %d", config$lab_window_h)
  
  labs <- dbGetQuery(con, qry) %>% as_tibble() %>% mutate(charttime = parse_datetime(charttime))
  
  labs_first <- labs %>%
    group_by(stay_id, itemid) %>%
    slice_min(hours_from_ed, n = 1, with_ties = FALSE) %>% ungroup()
  
  lab_map <- if (file.exists(config$mimic_paths$d_labitems)) {
    read_csv(config$mimic_paths$d_labitems, show_col_types = FALSE) %>%
      select(itemid, label, fluid, category)
  } else tibble(itemid = integer(), label=character(), fluid=character(), category=character())
  
  freq <- labs_first %>%
    group_by(itemid) %>%
    summarise(
      n_patients   = n_distinct(stay_id),
      pct_patients = 100 * n_patients / n_distinct(labs_first$stay_id),
      n_abnormal   = sum(flag == "abnormal", na.rm = TRUE),
      pct_abnormal = 100 * n_abnormal / n(),
      mean_value   = mean(valuenum, na.rm = TRUE),
      .groups = "drop"
    ) %>% left_join(lab_map, by="itemid") %>% arrange(desc(n_patients))
  
  list(processed_labs = labs_first, lab_frequency = freq)
}

# ---------------------- Add observations as "tests" (optional) ----------------------
obs_as_tests <- function(supp, cohort) {
  v <- if (!is.null(supp$vitals) && nrow(supp$vitals)>0) supp$vitals else supp$triage
  if (is.null(v) || nrow(v)==0) return(NULL)
  f2c <- function(t){ ifelse(is.na(t), NA_real_, ifelse(t > 50, (t-32)*5/9, t)) }
  first_non_na <- function(x){ i <- which(!is.na(x))[1]; if (is.na(i)) NA else x[i] }
  
  time_col <- intersect(names(v), c("charttime","time","recordedtime","measured_time"))[1] %||% names(v)[1]
  first <- v %>% arrange(stay_id, .data[[time_col]]) %>%
    group_by(stay_id) %>% summarise(across(everything(), first_non_na), .groups="drop") %>%
    mutate(temp_c = f2c(temperature))
  
  flags <- first %>%
    mutate(
      resprate_abn  = if_else(!is.na(resprate)  & (resprate  <12 | resprate  >20), 1L, 0L, missing = NA_integer_),
      o2sat_abn     = if_else(!is.na(o2sat)     & (o2sat     <95),                1L, 0L, missing = NA_integer_),
      temp_abn      = if_else(!is.na(temp_c)    & (temp_c <36.1 | temp_c >38.0),  1L, 0L, missing = NA_integer_),
      sbp_abn       = if_else(!is.na(sbp)       & (sbp <90 | sbp >140),           1L, 0L, missing = NA_integer_),
      hr_abn        = if_else(!is.na(heartrate) & (heartrate <60 | heartrate>100),1L, 0L, missing = NA_integer_),
      dbp_abn       = if_else(!is.na(dbp)       & (dbp <60 | dbp >90),            1L, 0L, missing = NA_integer_)
    )
  
  idmap <- cohort$edstays %>% select(stay_id, subject_id, intime)
  make_rows <- function(df, itemid, value_col, flag_col){
    df %>% filter(!is.na(.data[[value_col]])) %>%
      select(stay_id, valuenum = all_of(value_col), flag = all_of(flag_col)) %>%
      left_join(idmap, by="stay_id") %>%
      mutate(itemid = itemid,
             flag   = if_else(flag==1L, "abnormal", "normal"),
             charttime = intime, hours_from_ed = 0)
  }
  
  out <- bind_rows(
    make_rows(flags, 900001L, "resprate",  "resprate_abn"),
    make_rows(flags, 900002L, "o2sat",     "o2sat_abn"),
    make_rows(flags, 900003L, "temp_c",    "temp_abn"),
    make_rows(flags, 900004L, "sbp",       "sbp_abn"),
    make_rows(flags, 900005L, "heartrate", "hr_abn"),
    make_rows(flags, 900006L, "dbp",       "dbp_abn")
  )
  
  freq <- out %>% group_by(itemid) %>%
    summarise(
      n_patients   = n_distinct(stay_id),
      pct_patients = 100 * n_patients / n_distinct(out$stay_id),
      n_abnormal   = sum(flag == "abnormal", na.rm = TRUE),
      pct_abnormal = 100 * n_abnormal / n(),
      mean_value   = mean(valuenum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      itemid==900001 ~ "Respiratory rate", itemid==900002 ~ "SpO2",
      itemid==900003 ~ "Temperature",      itemid==900004 ~ "Systolic BP",
      itemid==900005 ~ "Heart rate",       itemid==900006 ~ "Diastolic BP",
      TRUE ~ "Observation"
    ), fluid = "vital sign", category = "Observations")
  
  list(processed_labs = out, lab_frequency = freq)
}

merge_labs_obs <- function(labs, obs){
  if (is.null(obs)) return(labs)
  list(
    processed_labs = bind_rows(labs$processed_labs, obs$processed_labs),
    lab_frequency  = bind_rows(labs$lab_frequency,  obs$lab_frequency) %>% arrange(desc(n_patients))
  )
}

# ---------------------- Outcomes ----------------------
define_outcomes <- function(cohort){
  base <- cohort$edstays %>%
    select(subject_id, stay_id, hadm_id, intime, outtime, disposition, age_at_ed, gender)
  
  hosp <- base %>%
    filter(!is.na(hadm_id)) %>%
    left_join(cohort$admissions %>% select(subject_id, hadm_id, admittime, dischtime, hospital_expire_flag),
              by = c("subject_id","hadm_id")) %>%
    mutate(dischtime = if_else(is.na(dischtime), outtime, dischtime),
           los_days = as.numeric(dischtime - admittime, units="days"),
           los_event = if_else(hospital_expire_flag==0, 1L, 0L),
           in_hosp_death = as.integer(hospital_expire_flag))
  
  mort <- base %>%
    left_join(cohort$patients %>% select(subject_id, dod), by="subject_id") %>%
    mutate(ttd = as.numeric(dod - intime, units="days"),
           death_1yr_event = if_else(!is.na(ttd) & ttd>=0 & ttd <= 365.25, 1L, 0L),
           time_to_event_1yr = if_else(death_1yr_event==1L, ttd, 365.25))
  
  readm <- base %>%
    filter(!is.na(hadm_id)) %>%
    left_join(cohort$patients %>% select(subject_id, dod), by="subject_id") %>%
    left_join(
      cohort$admissions %>% select(subject_id, hadm_id, admittime, dischtime) %>%
        arrange(subject_id, admittime) %>% group_by(subject_id) %>% mutate(next_admittime = lead(admittime)) %>% ungroup(),
      by = c("subject_id","hadm_id")
    ) %>%
    mutate(days_to_readm = as.numeric(next_admittime - dischtime, units="days"),
           ttd_from_dc   = as.numeric(dod - dischtime, units="days"),
           readm_30d_event = as.integer(!is.na(days_to_readm) & days_to_readm>=0 & days_to_readm<=30),
           t_readm_30d     = if_else(readm_30d_event==1L, days_to_readm,
                                     if_else(!is.na(ttd_from_dc) & ttd_from_dc>=0 & ttd_from_dc<=30, ttd_from_dc, 30)),
           readm_1yr_event = as.integer(!is.na(days_to_readm) & days_to_readm>=0 & days_to_readm<=365.25),
           t_readm_1yr     = if_else(readm_1yr_event==1L, days_to_readm,
                                     if_else(!is.na(ttd_from_dc) & ttd_from_dc>=0 & ttd_from_dc<=365.25, ttd_from_dc, 365.25)))
  
  dispo <- base %>% mutate(
    disposition_category = case_when(
      str_detect(disposition, regex("HOME|DISCHARGE", TRUE)) ~ "HOME",
      str_detect(disposition, regex("ADMIT", TRUE))          ~ "ADMITTED",
      str_detect(disposition, regex("TRANSFER", TRUE))       ~ "TRANSFER",
      str_detect(disposition, regex("EXPIRED|DIED|DEATH", TRUE)) ~ "EXPIRED",
      str_detect(disposition, regex("AMA", TRUE))            ~ "LEFT_AMA",
      str_detect(disposition, regex("LWBS", TRUE))           ~ "LEFT_WITHOUT_BEING_SEEN",
      TRUE ~ "OTHER"
    ),
    disposition_admitted = disposition_category=="ADMITTED")
  
  base %>%
    left_join(hosp %>% select(stay_id, los_days, los_event, in_hosp_death), by="stay_id") %>%
    left_join(mort %>% select(stay_id, death_1yr_event, time_to_event_1yr), by="stay_id") %>%
    left_join(readm %>% select(stay_id, readm_30d_event, t_readm_30d, readm_1yr_event, t_readm_1yr), by="stay_id") %>%
    left_join(dispo %>% select(stay_id, disposition_category, disposition_admitted), by="stay_id")
}

# ---------------------- NEWS2 ----------------------
calc_news2 <- function(supp){
  v <- if (!is.null(supp$vitals) && nrow(supp$vitals)>0) supp$vitals else supp$triage
  if (is.null(v) || nrow(v)==0) return(tibble(stay_id=integer(), news2_score=integer()))
  f2c <- function(t) ifelse(is.na(t), NA_real_, ifelse(t>50, (t-32)*5/9, t))
  time_col <- intersect(names(v), c("charttime","time","recordedtime","measured_time"))[1] %||% names(v)[1]
  v1 <- v %>% group_by(stay_id) %>% arrange(.data[[time_col]]) %>% slice_head(n=1) %>%
    ungroup() %>% mutate(temp_c = f2c(temperature),
                         conc_alert = TRUE,  # default Alert if not documented
                         supp_o2 = if ("o2flow" %in% names(.)) !is.na(o2flow) & o2flow>0 else FALSE)
  score_fun <- function(rr, spo2, supp_o2, temp, sbp, hr, alert){
    miss <- 0; s <- 0
    if (is.na(rr)) miss <- miss+1 else s <- s + ifelse(rr<=8,3, ifelse(rr<=11,1, ifelse(rr<=20,0, ifelse(rr<=24,2,3))))
    if (is.na(spo2)) miss <- miss+1 else s <- s + ifelse(spo2<=91,3, ifelse(spo2<=93,2, ifelse(spo2<=95,1,0)))
    s <- s + ifelse(isTRUE(supp_o2), 2, 0)
    if (is.na(temp)) miss <- miss+1 else s <- s + ifelse(temp<=35,3, ifelse(temp<=36,1, ifelse(temp<=38,0, ifelse(temp<=39,1,2))))
    if (is.na(sbp)) miss <- miss+1 else s <- s + ifelse(sbp<=90,3, ifelse(sbp<=100,2, ifelse(sbp<=110,1, ifelse(sbp<=219,0,3))))
    if (is.na(hr))  miss <- miss+1 else s <- s + ifelse(hr<=40,3, ifelse(hr<=50,1, ifelse(hr<=90,0, ifelse(hr<=110,1, ifelse(hr<=130,2,3)))))
    s <- s + ifelse(isTRUE(alert), 0, 3)
    if (miss > 2) NA_integer_ else as.integer(s)
  }
  v1 %>% mutate(news2_score = pmap_int(list(resprate,o2sat,supp_o2,temp_c,sbp,heartrate,conc_alert), score_fun)) %>%
    select(stay_id, news2_score)
}

# ---------------------- FI-Lab: importance, sets, scores ----------------------
calc_lab_importance <- function(labs, outcomes){
  eligible <- labs$lab_frequency %>% filter(pct_patients >= config$rf_min_coverage_pct) %>% pull(itemid)
  dt <- as.data.table(labs$processed_labs)[itemid %in% eligible]
  dt <- unique(dt, by=c("stay_id","itemid"))
  wide <- dcast(dt, stay_id ~ paste0("lab_", itemid), value.var="valuenum", fill=NA_real_)
  base <- outcomes %>% group_by(subject_id) %>% arrange(intime) %>% slice_head(n=1) %>% ungroup() %>%
    select(stay_id, death_1yr_event) %>% filter(!is.na(death_1yr_event)) %>%
    mutate(death_1yr_event = factor(death_1yr_event, levels=c(0,1)))
  rf_df <- base %>% inner_join(wide, by="stay_id")
  lab_cols <- grep("^lab_", names(rf_df), value=TRUE)
  if (nrow(rf_df) == 0 || length(lab_cols) == 0) {
    return(labs$lab_frequency %>%
             mutate(importance_score = pct_abnormal * log1p(n_patients),
                    included_in_rf = FALSE) %>% arrange(desc(importance_score)))
  }
  if (nrow(rf_df) > config$rf_max_n) rf_df <- rf_df %>% slice_sample(n = config$rf_max_n)
  w <- prop.table(table(rf_df$death_1yr_event)); wt <- ifelse(rf_df$death_1yr_event=="1", 1/w["1"], 1/w["0"])
  fit <- ranger(death_1yr_event ~ ., data = rf_df %>% select(-stay_id),
                importance = "permutation", case.weights = wt, num.trees = 500,
                num.threads = max(1, parallel::detectCores()-1), verbose = FALSE)
  imp <- enframe(fit$variable.importance, name="feature", value="importance_score") %>%
    filter(str_starts(feature,"lab_")) %>% mutate(itemid = as.integer(sub("lab_","",feature))) %>%
    select(itemid, importance_score)
  labs$lab_frequency %>%
    left_join(imp, by="itemid") %>%
    mutate(included_in_rf = itemid %in% eligible,
           importance_score = case_when(
             !included_in_rf ~ 0,
             is.na(importance_score) ~ 0,
             TRUE ~ importance_score * sqrt(pct_patients/100)
           )) %>% arrange(desc(importance_score))
}

make_feature_sets <- function(lab_imp){
  common <- lab_imp %>% arrange(desc(n_patients)) %>% pull(itemid)
  important <- lab_imp %>% arrange(desc(importance_score)) %>% pull(itemid)
  sets <- list()
  for (k in config$n_features_list) if (length(common) >= k) {
    sets[[paste0("Common_", k)]]     <- list(itemids = common[1:k],    strategy="Common",    n_features=k)
    sets[[paste0("Importance_", k)]] <- list(itemids = important[1:k], strategy="Importance", n_features=k)
    for (i in seq_len(config$n_random_sets)) {
      sets[[paste0("Random_", k, "_", i)]] <- list(itemids = sample(common, k), strategy="Random", n_features=k, random_id=i)
    }
  }
  sets
}

calc_filab_scores <- function(labs, feature_sets){
  base <- labs$processed_labs
  variants <- imap_dfr(feature_sets, function(info, name){
    base %>% filter(itemid %in% info$itemids) %>%
      group_by(stay_id) %>%
      summarise(set_name = name,
                strategy = info$strategy,
                n_features_requested = info$n_features,
                random_id = info$random_id %||% NA_integer_,
                n_features_used = n_distinct(itemid),
                n_abnormal = sum(flag=="abnormal", na.rm=TRUE),
                fi_lab_score = n_abnormal / n_features_used,
                .groups="drop")
  })
  bind_rows(lapply(config$feature_floors, function(floor){
    variants %>% filter(n_features_used >= floor, n_features_requested >= floor) %>%
      mutate(feature_floor = floor)
  }))
}

# ---------------------- Survival models ----------------------
fit_models <- function(outcomes, filab_scores, news2_scores){
  outs <- list(
    list(name="1-Year Mortality", time="time_to_event_1yr", event="death_1yr_event", first_stay=TRUE),
    list(name="In-Hospital Death", time="los_days", event="in_hosp_death", first_stay=TRUE),
    list(name="Length of Stay", time="los_days", event="los_event"),
    list(name="30-Day Readmission", time="t_readm_30d", event="readm_30d_event"),
    list(name="1-Year Readmission", time="t_readm_1yr", event="readm_1yr_event"),
    list(name="ED Admission", time="const1", event="disposition_admitted", binary=TRUE)
  )
  res <- list()
  for (o in outs){
    dat <- outcomes
    if (!is.null(o$first_stay) && o$first_stay) dat <- dat %>% group_by(subject_id) %>% arrange(intime) %>% slice_head(n=1) %>% ungroup()
    dat <- dat %>%
      left_join(news2_scores, by="stay_id") %>%
      mutate(time_var  = if (!is.null(o$binary) && o$binary) 1 else as.numeric(.data[[o$time]]),
             event_var = as.integer(.data[[o$event]]),
             age_s     = as.numeric(scale(age_at_ed)[,1]),
             news2_s   = if (sum(!is.na(news2_score))>1) as.numeric(scale(news2_score)[,1]) else NA_real_,
             gender_f  = factor(gender)) %>%
      filter(!is.na(time_var), !is.na(event_var), time_var > 0)
    
    cfgs <- filab_scores %>% filter(strategy %in% c("Common","Importance")) %>%
      distinct(strategy, n_features_requested, feature_floor)
    
    for (i in seq_len(nrow(cfgs))){
      cfg <- cfgs[i,]
      sc  <- filab_scores %>%
        filter(strategy==cfg$strategy,
               n_features_requested==cfg$n_features_requested,
               feature_floor==cfg$feature_floor) %>%
        select(stay_id, fi_lab_score)
      
      md <- dat %>% inner_join(sc, by="stay_id") %>%
        mutate(filab_s = as.numeric(scale(fi_lab_score)[,1])) %>%
        filter(complete.cases(time_var, event_var, filab_s, age_s, gender_f))
      
      if (nrow(md) < 30 || sum(md$event_var, na.rm=TRUE) < 10) next
      
      m1 <- coxph(Surv(time_var, event_var) ~ age_s + gender_f, data=md)
      m3 <- coxph(Surv(time_var, event_var) ~ filab_s + age_s + gender_f, data=md)
      
      m2 <- m4 <- NULL
      if (sum(!is.na(md$news2_s)) > 20){
        md2 <- md %>% filter(!is.na(news2_s))
        if (nrow(md2) > 20){
          m2 <- coxph(Surv(time_var, event_var) ~ news2_s + age_s + gender_f, data=md2)
          m4 <- coxph(Surv(time_var, event_var) ~ filab_s + news2_s + age_s + gender_f, data=md2)
        }
      }
      
      get_hr <- function(model, term){
        if (is.null(model) || !(term %in% rownames(coef(summary(model))))) return(c(NA, NA, NA, NA))
        sm <- summary(model)
        hr  <- exp(coef(model)[term])
        ci  <- exp(confint(model)[term, ])
        p   <- sm$coefficients[term, "Pr(>|z|)"]
        c(hr, ci[1], ci[2], p)
      }
      
      key <- tibble(
        outcome = o$name,
        strategy = cfg$strategy,
        n_features = cfg$n_features_requested,
        feature_floor = cfg$feature_floor,
        n_obs = nrow(md),
        n_events = sum(md$event_var, na.rm=TRUE),
        c_index_demographics = unname(summary(m1)$concordance[1]),
        c_index_news2 = if (!is.null(m2)) unname(summary(m2)$concordance[1]) else NA_real_,
        c_index_filab = unname(summary(m3)$concordance[1]),
        c_index_combined = if (!is.null(m4)) unname(summary(m4)$concordance[1]) else NA_real_,
        aic_demographics = AIC(m1),
        aic_news2 = if (!is.null(m2)) AIC(m2) else NA_real_,
        aic_filab = AIC(m3),
        aic_combined = if (!is.null(m4)) AIC(m4) else NA_real_
      )
      
      hr_filab    <- get_hr(m3, "filab_s")
      hr_news2    <- get_hr(m2, "news2_s")
      hr_filab_c  <- get_hr(m4, "filab_s")
      hr_news2_c  <- get_hr(m4, "news2_s")
      hr_age_d    <- get_hr(m1, "age_s")
      hr_age_n    <- get_hr(m2, "age_s"); hr_age_f <- get_hr(m3, "age_s"); hr_age_c <- get_hr(m4, "age_s")
      
      res[[length(res)+1]] <- bind_cols(
        key,
        tibble(
          hr_filab = hr_filab[1], hr_filab_lower = hr_filab[2], hr_filab_upper = hr_filab[3], p_filab = hr_filab[4],
          hr_news2 = hr_news2[1], hr_news2_lower = hr_news2[2], hr_news2_upper = hr_news2[3], p_news2 = hr_news2[4],
          hr_filab_in_combined = hr_filab_c[1], hr_filab_in_combined_lower = hr_filab_c[2], hr_filab_in_combined_upper = hr_filab_c[3], p_filab_in_combined = hr_filab_c[4],
          hr_news2_in_combined = hr_news2_c[1], hr_news2_in_combined_lower = hr_news2_c[2], hr_news2_in_combined_upper = hr_news2_c[3], p_news2_in_combined = hr_news2_c[4],
          hr_age_demographics = hr_age_d[1], hr_age_demographics_lower = hr_age_d[2], hr_age_demographics_upper = hr_age_d[3], p_age_demographics = hr_age_d[4],
          hr_age_news2 = hr_age_n[1], hr_age_news2_lower = hr_age_n[2], hr_age_news2_upper = hr_age_n[3], p_age_news2 = hr_age_n[4],
          hr_age_filab = hr_age_f[1], hr_age_filab_lower = hr_age_f[2], hr_age_filab_upper = hr_age_f[3], p_age_filab = hr_age_f[4],
          hr_age_combined = hr_age_c[1], hr_age_combined_lower = hr_age_c[2], hr_age_combined_upper = hr_age_c[3], p_age_combined = hr_age_c[4]
        )
      )
    }
  }
  if (length(res)==0) return(tibble())
  bind_rows(res)
}

# ---------------------- Uncertainty summary (random sets) ----------------------
uncertainty_summary <- function(filab_scores){
  filab_scores %>%
    filter(strategy=="Random") %>%
    group_by(stay_id, n_features_requested, feature_floor) %>%
    summarise(n_random_sets = n(),
              mean_score = mean(fi_lab_score, na.rm=TRUE),
              sd_score   = sd(fi_lab_score, na.rm=TRUE),
              cv_score   = sd_score/mean_score,
              min_score  = min(fi_lab_score, na.rm=TRUE),
              max_score  = max(fi_lab_score, na.rm=TRUE),
              .groups="drop") %>%
    group_by(n_features_requested, feature_floor) %>%
    summarise(n_patients = n(),
              median_cv = median(cv_score, na.rm=TRUE),
              q75_cv    = quantile(cv_score, 0.75, na.rm=TRUE),
              median_range = median(max_score - min_score, na.rm=TRUE),
              .groups="drop")
}

# ---------------------- Main ----------------------
main <- function(){
  log("Start")
  cohort <- load_cohort()
  if (!is.na(config$sample_size) && config$sample_size > 0) {
    set.seed(42)  # reproducible subset
    n_avail <- nrow(cohort$edstays)
    n_take  <- min(as.integer(config$sample_size), n_avail)
    cohort$edstays <- cohort$edstays %>% slice_sample(n = n_take)
    log("Sample run: kept %d of %d ED stays", n_take, n_avail)
  }
  supp   <- load_supplementary(cohort$edstays$stay_id)
  labs   <- extract_labs_window(cohort)
  
  if (config$include_observations){
    obs <- obs_as_tests(supp, cohort)
    if (!is.null(obs)) labs <- merge_labs_obs(labs, obs)
  }
  
  # Filter to stays with ≥1 lab
  keep_stays <- labs$processed_labs %>% distinct(stay_id) %>% pull()
  cohort$edstays <- cohort$edstays %>% filter(stay_id %in% keep_stays)
  
  outcomes <- define_outcomes(cohort)
  news2    <- calc_news2(supp)
  lab_imp  <- calc_lab_importance(labs, outcomes)
  fsets    <- make_feature_sets(lab_imp)
  filab    <- calc_filab_scores(labs, fsets)
  models   <- fit_models(outcomes, filab, news2)
  uncert   <- uncertainty_summary(filab)
  
  # ------------ Exports (ONLY the essentials, edit as needed) ------------
  write_csv(models, file.path(config$out_dir, "tables", "model_results.csv"))
  
  log("Done")
}

main()
