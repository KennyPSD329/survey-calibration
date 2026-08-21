# ==============================================================================
# 01_prepare_data_application_data.R
#
# Reads all raw NHANES (2001-2006) and NHIS (2003-2004) public-use files,
# processes covariates and mortality linkage, and writes .rds cache files to
# data/processed/.  Subsequent scripts (02, 03) load the cache and never
# re-read raw data.
#
# Working directory must be set to  application/  before sourcing.
# Expected raw data layout:
#   data/raw/nhanes/99/   data/raw/nhanes/99x/  data/raw/nhanes/99p/
#   data/raw/nhanes/01/   data/raw/nhanes/01x/  data/raw/nhanes/01p/
#   data/raw/nhanes/03/   data/raw/nhanes/03x/  data/raw/nhanes/03p/
#   data/raw/nhanes/05/   data/raw/nhanes/05x/  data/raw/nhanes/05p/
#   data/raw/mortality/   (NHANES and NHIS .dat files)
#   data/raw/nhis/        (nhis_2003_2004.csv)
# ==============================================================================

rm(list = ls())
library(haven)
library(stringr)
library(readr)
library(dplyr)

if (!dir.exists("data/raw/nhanes")) {
  stop("Set working directory to data_application/ before running this script.")
}

# ------------------------------------------------------------------------------
# Helper: read all non-.dat XPT files in one cycle folder, merge on SEQN
# ------------------------------------------------------------------------------
read_nhanes_cycle <- function(cycle_path) {
  fnames <- list.files(path = cycle_path)
  fnames <- fnames[!str_detect(fnames, "\\.dat$")]
  dat <- NULL
  for (f in fnames) {
    tmp <- read_xpt(file.path(cycle_path, f))
    dat <- if (is.null(dat)) tmp else merge(tmp, dat, by = "SEQN", all = TRUE)
  }
  dat
}

# ==============================================================================
# 1. Load NHANES XPT files (main, DXA, physical activity)
# ==============================================================================
cat("Loading NHANES XPT files...\n")
nhanes00  <- read_nhanes_cycle("data/raw/nhanes/99")
nhanes00x <- read_nhanes_cycle("data/raw/nhanes/99x")
nhanes00p <- read_nhanes_cycle("data/raw/nhanes/99p")
nhanes01  <- read_nhanes_cycle("data/raw/nhanes/01")
nhanes01x <- read_nhanes_cycle("data/raw/nhanes/01x")
nhanes01p <- read_nhanes_cycle("data/raw/nhanes/01p")
nhanes03  <- read_nhanes_cycle("data/raw/nhanes/03")
nhanes03x <- read_nhanes_cycle("data/raw/nhanes/03x")
nhanes03p <- read_nhanes_cycle("data/raw/nhanes/03p")
nhanes05  <- read_nhanes_cycle("data/raw/nhanes/05")
nhanes05x <- read_nhanes_cycle("data/raw/nhanes/05x")
nhanes05p <- read_nhanes_cycle("data/raw/nhanes/05p")

# ==============================================================================
# 2. DXA processing: keep valid whole-body scans, average per SEQN
# ==============================================================================
cat("Processing DXA data...\n")
vars_dxa <- c("DXDTOFAT", "DXXTRFAT", "DXDTRPF", "DXDTOPF")

proc_dxa <- function(df, cycle05 = FALSE) {
  if (cycle05) {
    df <- df %>% filter(DXAEXSTS == 1 & DXITOTBN == 0 & DXITOTST == 0)
  } else {
    df <- df %>% filter(DXAEXSTS == 1 & DXITOT == 0)
  }
  df %>%
    group_by(SEQN) %>%
    summarise(across(all_of(vars_dxa), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

nhanes00x <- proc_dxa(nhanes00x)
nhanes01x <- proc_dxa(nhanes01x)
nhanes03x <- proc_dxa(nhanes03x)
nhanes05x <- proc_dxa(nhanes05x, cycle05 = TRUE)

nhanes00 <- merge(nhanes00, nhanes00x, by = "SEQN", all.x = TRUE)
nhanes01 <- merge(nhanes01, nhanes01x, by = "SEQN", all.x = TRUE)
nhanes03 <- merge(nhanes03, nhanes03x, by = "SEQN", all.x = TRUE)
nhanes05 <- merge(nhanes05, nhanes05x, by = "SEQN", all.x = TRUE)
rm(nhanes00x, nhanes01x, nhanes03x, nhanes05x)

# ==============================================================================
# 3. Physical activity processing: compute weekly MET score
# ==============================================================================
cat("Processing physical activity data...\n")
proc_pa <- function(df) {
  df$PADTIMES[is.na(df$PADTIMES)] <- 0
  df$PADDURAT[is.na(df$PADDURAT)] <- 0
  df$PADMETS [is.na(df$PADMETS)]  <- 0
  df %>%
    mutate(MET_value = PADTIMES * PADDURAT * PADMETS) %>%
    group_by(SEQN) %>%
    summarise(total_MET_wk = sum(MET_value, na.rm = TRUE) / 4.3, .groups = "drop")
}

nhanes00p <- proc_pa(nhanes00p)
nhanes01p <- proc_pa(nhanes01p)
nhanes03p <- proc_pa(nhanes03p)
nhanes05p <- proc_pa(nhanes05p)

nhanes00 <- merge(nhanes00, nhanes00p, by = "SEQN", all.x = TRUE)
nhanes01 <- merge(nhanes01, nhanes01p, by = "SEQN", all.x = TRUE)
nhanes03 <- merge(nhanes03, nhanes03p, by = "SEQN", all.x = TRUE)
nhanes05 <- merge(nhanes05, nhanes05p, by = "SEQN", all.x = TRUE)
rm(nhanes00p, nhanes01p, nhanes03p, nhanes05p)

# ==============================================================================
# 4. NHANES mortality linkage
# ==============================================================================
cat("Reading NHANES mortality files...\n")
read_nhanes_mort_dat <- function(path) {
  read_fwf(
    file      = path,
    col_types = "iiiiiiii",
    fwf_cols(
      seqn         = c(1, 6),
      eligstat     = c(15, 15),
      mortstat     = c(16, 16),
      ucod_leading = c(17, 19),
      diabetes     = c(20, 20),
      hyperten     = c(21, 21),
      permth_int   = c(43, 45),
      permth_exm   = c(46, 48)
    ),
    na = c("", ".")
  )
}

mort00 <- read_nhanes_mort_dat("data/raw/mortality/NHANES_1999_2000_MORT_2019_PUBLIC.dat")
mort01 <- read_nhanes_mort_dat("data/raw/mortality/NHANES_2001_2002_MORT_2019_PUBLIC.dat")
mort03 <- read_nhanes_mort_dat("data/raw/mortality/NHANES_2003_2004_MORT_2019_PUBLIC.dat")
mort05 <- read_nhanes_mort_dat("data/raw/mortality/NHANES_2005_2006_MORT_2019_PUBLIC.dat")

# Rename seqn -> SEQN for merging
for (nm in c("mort00", "mort01", "mort03", "mort05")) {
  df <- get(nm)
  names(df)[names(df) == "seqn"] <- "SEQN"
  assign(nm, df)
}

nhanes_mort00 <- merge(nhanes00, mort00, by = "SEQN")
nhanes_mort01 <- merge(nhanes01, mort01, by = "SEQN")
nhanes_mort03 <- merge(nhanes03, mort03, by = "SEQN")
nhanes_mort05 <- merge(nhanes05, mort05, by = "SEQN")
rm(nhanes00, nhanes01, nhanes03, nhanes05, mort00, mort01, mort03, mort05)

# ==============================================================================
# 5. Align variable names across cycles, assign cycle and source labels
# ==============================================================================
nhanes_mort00$ALQ101 <- nhanes_mort00$ALQ100
nhanes_mort01$ALQ101 <- nhanes_mort01$ALD100
nhanes_mort05$SMQ640 <- nhanes_mort05$SMD641
nhanes_mort05$SMQ650 <- nhanes_mort05$SMD650

nhanes_mort00$cycle <- 1L; nhanes_mort00$source <- "external"
nhanes_mort01$cycle <- 2L; nhanes_mort01$source <- "external"
nhanes_mort03$cycle <- 3L; nhanes_mort03$source <- "internal"
nhanes_mort05$cycle <- 4L; nhanes_mort05$source <- "external"

dsgn.vars  <- c("SEQN", "SDMVPSU", "SDMVSTRA", "SDDSRVYR", "WTINT2YR", "WTMEC2YR")
Demgph.vars <- c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC3", "DMDEDUC2", "DMDMARTL")
bdy.vars    <- c("BMXBMI", "BMXWAIST", "BMXHT", "BMXWT", "WHD010", "WHD020")
fat.vars    <- c("DXDTOFAT", "DXXTRFAT", "DXDTRPF", "DXDTOPF")
smk.vars    <- c("SMQ020", "SMD030", "SMQ040", "SMQ050Q", "SMQ640", "SMQ620", "SMQ650")
Alch.vars   <- c("ALQ101", "ALQ110", "ALQ120Q", "ALQ120U", "ALQ130", "ALQ140Q", "ALQ140U", "ALQ150")
pa.vars     <- c("total_MET_wk")
mo.vars     <- c("mortstat", "permth_int", "permth_exm", "eligstat")
cmn.vars    <- c(dsgn.vars, Demgph.vars, bdy.vars, fat.vars, smk.vars, Alch.vars, pa.vars, mo.vars)

select_cmn <- function(df) df[, c(cmn.vars[cmn.vars %in% names(df)], "cycle", "source")]

nhanes_mort00 <- select_cmn(nhanes_mort00)
nhanes_mort01 <- select_cmn(nhanes_mort01)
nhanes_mort03 <- select_cmn(nhanes_mort03)
nhanes_mort05 <- select_cmn(nhanes_mort05)

# Main analysis uses cycles 01+03+05 (excludes 1999-2000)
nhanes_raw <- rbind(nhanes_mort01, nhanes_mort05, nhanes_mort03)

# Save 1999-2000 separately for EDA summary table
nhanes_mort00_eda <- nhanes_mort00
rm(nhanes_mort00, nhanes_mort01, nhanes_mort03, nhanes_mort05)

# ==============================================================================
# 6. Recode NHANES covariates
# ==============================================================================
cat("Recoding NHANES covariates...\n")

# Alcohol ---
nhanes_raw$alchl <- 2 - nhanes_raw$ALQ101
nhanes_raw$alchl[nhanes_raw$ALQ110 == 1] <- 1
nhanes_raw$alchl[nhanes_raw$ALQ110 == 9] <- 0

dys <- nhanes_raw$ALQ120U
dys[(nhanes_raw$ALQ101 == 2) & (nhanes_raw$ALQ110 == 2)] <- 0
dys[nhanes_raw$ALQ120Q == 0] <- 1
dys[nhanes_raw$ALQ120U == 1] <- 7
dys[nhanes_raw$ALQ120U == 2] <- 365.25 / 12
dys[nhanes_raw$ALQ120U == 3] <- 365.25
nhanes_raw$ach_freq <- nhanes_raw$ALQ120Q / dys
nhanes_raw$ach_freq[nhanes_raw$alchl == 0] <- 0

nhanes_raw$ALQ130[nhanes_raw$ALQ130 > 82] <- NA
nhanes_raw$ALQ130[nhanes_raw$ALQ130 > 15] <- 15
nhanes_raw$ALQ130[(nhanes_raw$ALQ101 == 2) & (nhanes_raw$ALQ110 == 2)] <- 0
nhanes_raw$drk_pdy <- nhanes_raw$ach_freq * nhanes_raw$ALQ130

nhanes_raw$alc3 <- 1
nhanes_raw$alc3[((nhanes_raw$ALQ101 == 1) | (nhanes_raw$ALQ110 == 1)) &
                (nhanes_raw$ALQ120Q <= 365) & (nhanes_raw$ALQ120Q > 0)] <- 3
nhanes_raw$alc3[((nhanes_raw$ALQ101 == 1) | (nhanes_raw$ALQ110 == 1)) &
                (nhanes_raw$ALQ120Q == 0)] <- 2
nhanes_raw$drk_pdy[nhanes_raw$alc3 == 2] <- 0
nhanes_raw$drk_pdy[(nhanes_raw$alc3 == 3) & is.na(nhanes_raw$drk_pdy)] <-
  mean(nhanes_raw$drk_pdy[(nhanes_raw$alc3 == 3) & !is.na(nhanes_raw$drk_pdy)])
nhanes_raw$alc2 <- nhanes_raw$alc3
nhanes_raw$alc2[((nhanes_raw$drk_pdy * 365.25) < 12) | (nhanes_raw$alc3 == 2)] <- 1
nhanes_raw$alc2[nhanes_raw$alc2 == 3] <- 2
nhanes_raw$alc2 <- nhanes_raw$alc2 - 1

nhanes_raw$hvy.drk <- 1
nhanes_raw$hvy.drk[(nhanes_raw$drk_pdy > 1) & (nhanes_raw$RIAGENDR == 2)] <- 2
nhanes_raw$hvy.drk[(nhanes_raw$drk_pdy > 2) & (nhanes_raw$RIAGENDR == 1)] <- 2
nhanes_raw$hvy.drk[nhanes_raw$alc2 == 0] <- 0

# Smoking ---
nhanes_raw$smk3 <- 0
nhanes_raw$smk3[(nhanes_raw$SMQ020 == 1) & (nhanes_raw$SMQ040 %in% c(3, 7))] <- 1
nhanes_raw$smk3[(nhanes_raw$SMQ020 == 1) & (nhanes_raw$SMQ040 %in% c(1, 2))] <- 2
nhanes_raw$smk3[is.na(nhanes_raw$SMQ020) | nhanes_raw$SMQ020 %in% c(7, 9)] <- 0
nhanes_raw$smk3[!is.na(nhanes_raw$SMQ050Q)] <- 1
nhanes_raw$smk3[!is.na(nhanes_raw$SMQ640) | !is.na(nhanes_raw$SMQ650)] <- 2

# Mortality ---
nhanes_elig <- subset(nhanes_raw, eligstat == 1)
nhanes_elig$mortality <- ifelse(
  nhanes_elig$mortstat == 1 &
  (nhanes_elig$permth_exm <= 120 | is.na(nhanes_elig$permth_exm)),
  1, 0
)

# Age ---
nhanes_proc <- nhanes_elig %>% filter(RIDAGEYR >= 18 & RIDAGEYR <= 69)
nhanes_proc$age_c6 <- cut(nhanes_proc$RIDAGEYR, breaks = c(18, seq(30, 80, 10)),
                           include.lowest = TRUE)
nhanes_proc$age_match_mort <- cut(nhanes_proc$RIDAGEYR,
                                   breaks = c(18, 25, 35, 45, 55, 70),
                                   right = TRUE, include.lowest = TRUE)
nhanes_proc <- nhanes_proc %>% filter(RIDAGEYR >= 40)

# DXA constraint (requires DXA to be available) ---
nhanes_proc <- nhanes_proc %>% filter(!is.na(DXDTOFAT))

# Multi-cycle weight adjustments ---
# WTINT2YR4: 3-cycle combined design (01+03+05)
# WTINT2YR3: 2-cycle external only (01+05)
nhanes_proc$WTINT2YR4 <- nhanes_proc$WTINT2YR / 3
nhanes_proc$WTINT2YR3 <- nhanes_proc$WTINT2YR / 2

# Gender ---
nhanes_proc$RIAGENDR <- factor(nhanes_proc$RIAGENDR,
                                levels = c(1, 2), labels = c("Male", "Female"))

# Education ---
nhanes_proc$DMDEDUC2[nhanes_proc$DMDEDUC2 %in% c(7, 9)] <- 4
nhanes_proc$DMDEDUC3[nhanes_proc$DMDEDUC3 %in% c(7, 9)] <- 13
nhanes_proc$educ5 <- nhanes_proc$DMDEDUC2
nhanes_proc$educ5[is.na(nhanes_proc$DMDEDUC2) &
                  nhanes_proc$DMDEDUC3 %in% c(0:8, 55, 66)] <- 1
nhanes_proc$educ5[is.na(nhanes_proc$DMDEDUC2) &
                  nhanes_proc$DMDEDUC3 %in% 9:12] <- 2
nhanes_proc$educ5[is.na(nhanes_proc$DMDEDUC2) &
                  nhanes_proc$DMDEDUC3 %in% c(13, 14)] <- 3
nhanes_proc$educ5[is.na(nhanes_proc$DMDEDUC2) &
                  nhanes_proc$DMDEDUC3 == 15] <- 4

# Race ---
nhanes_proc$race4 <- nhanes_proc$RIDRETH1
nhanes_proc$race4[nhanes_proc$race4 == 1] <- 2   # merge Mexican-American into Hispanic
nhanes_proc$race4 <- as.factor(nhanes_proc$race4)
levels(nhanes_proc$race4) <- c("Hispanic", "NH-White", "NH-Black", "NH-Other")
nhanes_proc$race4 <- factor(nhanes_proc$race4,
                             levels = c("NH-White", "NH-Black", "Hispanic", "NH-Other"))

# BMI (self-reported; fall back to measured if self-report missing) ---
nhanes_proc$WHD020[nhanes_proc$WHD020 > 1000] <- NA
nhanes_proc$WHD010[nhanes_proc$WHD010 > 1000] <- NA
nhanes_proc$slfBMI <- nhanes_proc$WHD020 * 0.45359237 / (nhanes_proc$WHD010 * 0.0254)^2
nhanes_proc$slfBMI[is.na(nhanes_proc$slfBMI) & !is.na(nhanes_proc$BMXBMI)] <-
  nhanes_proc$BMXBMI[is.na(nhanes_proc$slfBMI) & !is.na(nhanes_proc$BMXBMI)]
nhanes_proc <- nhanes_proc %>% filter(!is.na(slfBMI) & slfBMI >= 18.5)

# Physical activity ---
nhanes_proc$total_MET_wk[is.na(nhanes_proc$total_MET_wk)] <- 0
nhanes_proc$phys.all <- nhanes_proc$total_MET_wk /
  sum(nhanes_proc$total_MET_wk) * nrow(nhanes_proc)

# Height ---
nhanes_proc$height <- nhanes_proc$WHD010 * 2.54
nhanes_proc$height[is.na(nhanes_proc$height) & !is.na(nhanes_proc$BMXHT)] <-
  nhanes_proc$BMXHT[is.na(nhanes_proc$height) & !is.na(nhanes_proc$BMXHT)]
nhanes_proc <- nhanes_proc %>% filter(height >= 120 & height <= 220)

# Select final analysis columns ---
keep_cols <- c("SDMVPSU", "SDMVSTRA", "WTINT2YR", "WTINT2YR4", "WTINT2YR3",
               "RIDAGEYR", "RIAGENDR", "race4", "hvy.drk", "smk3", "phys.all",
               "mortality", "educ5", "DXDTOFAT", "DXXTRFAT", "DXDTRPF", "DXDTOPF",
               "age_match_mort", "slfBMI", "height", "source", "cycle")
nhanes_analysis <- nhanes_proc[, keep_cols]

rm(nhanes_raw, nhanes_elig, nhanes_proc)

cat(sprintf("NHANES analysis: n = %d  (internal = %d, external = %d)\n",
    nrow(nhanes_analysis),
    sum(nhanes_analysis$source == "internal"),
    sum(nhanes_analysis$source == "external")))

# ==============================================================================
# 7. Read and process NHIS 2003-2004
# ==============================================================================
cat("Loading NHIS data...\n")
nhis2003 <- read.csv("data/raw/nhis/nhis_2003_2004.csv", stringsAsFactors = FALSE)

read_nhis_mort_dat <- function(path) {
  read_fwf(
    file      = path,
    col_types = "ciiiiiiidd",
    fwf_cols(
      publicid     = c(1, 14),
      eligstat     = c(15, 15),
      mortstat     = c(16, 16),
      ucod_leading = c(17, 19),
      diabetes     = c(20, 20),
      hyperten     = c(21, 21),
      dodqtr       = c(22, 22),
      dodyear      = c(23, 26),
      wgt_new      = c(27, 34),
      sa_wgt_new   = c(35, 42)
    ),
    na = c("", ".")
  )
}

nhis_mort03 <- read_nhis_mort_dat("data/raw/mortality/NHIS_2003_MORT_2019_PUBLIC.dat")
nhis_mort04 <- read_nhis_mort_dat("data/raw/mortality/NHIS_2004_MORT_2019_PUBLIC.dat")

nhis2003$PUBLICID <- sprintf("%04d%06d%02d%02d",
                              nhis2003$YEAR, nhis2003$HHX, nhis2003$FMX, nhis2003$PX)
nhis03 <- nhis2003[nhis2003$YEAR == 2003, ]
nhis04 <- nhis2003[nhis2003$YEAR == 2004, ]

nhis_03_mort <- merge(nhis03, nhis_mort03, by.x = "PUBLICID", by.y = "publicid")
nhis_04_mort <- merge(nhis04, nhis_mort04, by.x = "PUBLICID", by.y = "publicid")
nhis_raw     <- rbind(nhis_03_mort, nhis_04_mort)
rm(nhis03, nhis04, nhis_03_mort, nhis_04_mort, nhis_mort03, nhis_mort04, nhis2003)

# Mortality ---
nhis_elig <- subset(nhis_raw, eligstat == 1)
nhis_elig$mortality <- ifelse(
  (nhis_elig$mortstat == 1 & nhis_elig$dodyear <= 2013 & nhis_elig$YEAR == 2003) |
  (nhis_elig$mortstat == 1 & nhis_elig$dodyear <= 2014 & nhis_elig$YEAR == 2004),
  1, 0
)

nhis_elig$WTINT2YR <- nhis_elig$SAMPWEIGHT / 2   # 2-year pooled weight

nhis_elig <- nhis_elig %>%
  filter(AGE >= 18 & AGE <= 69,
         WEIGHT < 300, HEIGHT < 77,
         !(SEX == 2 & PREGNANTNOW %in% c(2, 7, 8, 9)))

# Age ---
nhis_elig <- nhis_elig %>% filter(!is.na(AGE))
nhis_elig$age_c6 <- cut(nhis_elig$AGE,
                         breaks = c(18, seq(30, 80, 10)), include.lowest = TRUE)
nhis_elig$age_match_mort <- cut(nhis_elig$AGE,
                                 breaks = c(18, 25, 35, 45, 55, 70),
                                 right = TRUE, include.lowest = TRUE)
nhis_elig$RIDAGEYR <- nhis_elig$AGE
nhis_elig <- nhis_elig %>% filter(RIDAGEYR >= 40)

# Education ---
nhis_elig$educ5 <- NA
nhis_elig$educ5[nhis_elig$EDUC %in% c(101:111, 996)] <- 1
nhis_elig$educ5[nhis_elig$EDUC %in% c(112:116, 100)] <- 2
nhis_elig$educ5[nhis_elig$EDUC %in% c(200:202)]       <- 3
nhis_elig$educ5[nhis_elig$EDUC %in% c(300:303)]       <- 4
nhis_elig$educ5[nhis_elig$EDUC %in% c(500:530, 400)]  <- 5
nhis_elig <- nhis_elig %>% filter(!is.na(educ5))

# Race ---
nhis_elig$race4 <- NA
nhis_elig$race4[nhis_elig$HISPYN == 2] <- "Hispanic"
nhis_elig$race4[nhis_elig$HISPYN == 1 & nhis_elig$RACENEW == 100] <- "NH-White"
nhis_elig$race4[nhis_elig$HISPYN == 1 & nhis_elig$RACENEW == 200] <- "NH-Black"
nhis_elig$race4[nhis_elig$HISPYN == 1 & nhis_elig$RACENEW >= 300 &
                nhis_elig$RACENEW < 997] <- "NH-Other"
nhis_elig$race4 <- factor(nhis_elig$race4,
                           levels = c("NH-White", "NH-Black", "Hispanic", "NH-Other"))
nhis_elig <- nhis_elig %>% filter(!is.na(race4))

# Gender ---
nhis_elig$RIAGENDR <- factor(nhis_elig$SEX,
                              levels = c(1, 2), labels = c("Male", "Female"))

# BMI ---
nhis_elig$WEIGHT[nhis_elig$WEIGHT %in% c(996, 997, 998, 999)] <- NA
nhis_elig$WEIGHT[nhis_elig$WEIGHT < 50 | nhis_elig$WEIGHT > 299] <- NA
nhis_elig$HEIGHT[nhis_elig$HEIGHT %in% c(95, 96, 97, 98, 99)] <- NA
nhis_elig$HEIGHT[nhis_elig$HEIGHT < 36 | nhis_elig$HEIGHT > 90] <- NA
nhis_elig$height   <- nhis_elig$HEIGHT * 2.54
nhis_elig$slfBMI   <- (nhis_elig$WEIGHT * 0.45359237) / ((nhis_elig$HEIGHT * 0.0254)^2)
nhis_elig <- nhis_elig %>% filter(!is.na(slfBMI) & slfBMI >= 18.5)

# Physical activity ---
na_special <- c(995, 996, 997, 998, 999, 9997, 9998, 9999)
nhis_elig$MOD10FNO[nhis_elig$MOD10FNO %in% na_special] <- NA
nhis_elig$MOD10DNO[nhis_elig$MOD10DNO %in% na_special] <- NA
nhis_elig$MOD10FTP[nhis_elig$MOD10FTP %in% c(97, 98, 99)] <- 0
nhis_elig$MOD10DTP[nhis_elig$MOD10DTP %in% c(7, 8, 9)]   <- 0
nhis_elig$VIG10FNO[nhis_elig$VIG10FNO %in% na_special] <- NA
nhis_elig$VIG10DNO[nhis_elig$VIG10DNO %in% na_special] <- NA
nhis_elig$VIG10FTP[nhis_elig$VIG10FTP %in% c(97, 98, 99)] <- 0
nhis_elig$VIG10DTP[nhis_elig$VIG10DTP %in% c(7, 8, 9)]   <- 0

freq_per_week <- function(no, tp) {
  freq <- rep(NA_real_, length(no))
  freq[tp == 2] <- no[tp == 2]             # per day (treated as per week here)
  freq[tp == 3] <- no[tp == 3] / 1         # per week
  freq[tp == 4] <- no[tp == 4] / 4.3       # per month
  freq[tp == 5] <- no[tp == 5] / 52        # per year
  freq
}
dur_in_min <- function(no, tp) {
  d <- rep(NA_real_, length(no))
  d[tp == 1] <- no[tp == 1]
  d[tp == 2] <- no[tp == 2] * 60
  d
}

nhis_elig$mod_MET <- freq_per_week(nhis_elig$MOD10FNO, nhis_elig$MOD10FTP) *
                     dur_in_min(nhis_elig$MOD10DNO, nhis_elig$MOD10DTP) * 4.0
nhis_elig$vig_MET <- freq_per_week(nhis_elig$VIG10FNO, nhis_elig$VIG10FTP) *
                     dur_in_min(nhis_elig$VIG10DNO, nhis_elig$VIG10DTP) * 8.0
nhis_elig$mod_MET[is.na(nhis_elig$mod_MET)] <- 0
nhis_elig$vig_MET[is.na(nhis_elig$vig_MET)] <- 0
nhis_elig$LTPA_MET <- nhis_elig$mod_MET + nhis_elig$vig_MET
nhis_elig$phys.all <- nhis_elig$LTPA_MET /
  sum(nhis_elig$LTPA_MET, na.rm = TRUE) * nrow(nhis_elig)

# Smoking ---
nhis_elig$smk3 <- 0
nhis_elig$smk3[(nhis_elig$SMOKEV == 2) & (nhis_elig$SMOKFREQNOW %in% c(1, 7))] <- 1
nhis_elig$smk3[(nhis_elig$SMOKEV == 2) & (nhis_elig$SMOKFREQNOW %in% c(2, 3))] <- 2
nhis_elig$smk3[nhis_elig$SMOKEV %in% c(7, 9)] <- 0
nhis_elig$smk3[!nhis_elig$QUITNO %in% c(0, 997, 999)] <- 1
nhis_elig$smk3[!nhis_elig$CIGDAYMO %in% c(0, 96, 97, 98)] <- 2

# Alcohol ---
nhis_elig$alchl <- nhis_elig$ALC1YR - 1
nhis_elig$alchl[nhis_elig$alchl %in% c(6, 8)] <- 0
nhis_elig$alchl[nhis_elig$ALCLIFE == 2] <- 1
nhis_elig$alchl[nhis_elig$ALCLIFE %in% c(7, 9)] <- 0

nhis_elig$ALCDAYSYR[nhis_elig$ALCDAYSYR %in% c(996, 997, 999)] <- 0
nhis_elig$ach_freq <- nhis_elig$ALCDAYSYR / 365.25
nhis_elig$ALCAMT[nhis_elig$ALCAMT %in% c(96, 97, 98, 99)] <- 0
nhis_elig$ALCAMT[(nhis_elig$ALC1YR == 1) & (nhis_elig$ALCLIFE == 1)] <- 0
nhis_elig$drk_pdy <- nhis_elig$ach_freq * nhis_elig$ALCAMT

nhis_elig$alc3 <- 1
nhis_elig$alc3[((nhis_elig$ALC1YR == 2) | (nhis_elig$ALCLIFE == 2)) &
               (nhis_elig$ALCDAYSYR <= 365) & (nhis_elig$ALCDAYSYR > 0)] <- 3
nhis_elig$alc3[((nhis_elig$ALC1YR == 2) | (nhis_elig$ALCLIFE == 2)) &
               (nhis_elig$ALCDAYSYR == 0)] <- 2
nhis_elig$drk_pdy[nhis_elig$alc3 == 2] <- 0
nhis_elig$drk_pdy[(nhis_elig$alc3 == 3) & is.na(nhis_elig$drk_pdy)] <-
  mean(nhis_elig$drk_pdy[(nhis_elig$alc3 == 3) & !is.na(nhis_elig$drk_pdy)])
nhis_elig$alc2 <- nhis_elig$alc3
nhis_elig$alc2[((nhis_elig$drk_pdy * 365.25) < 12) | (nhis_elig$alc3 == 2)] <- 1
nhis_elig$alc2[nhis_elig$alc2 == 3] <- 2
nhis_elig$alc2 <- nhis_elig$alc2 - 1

nhis_elig$hvy.drk <- 1
nhis_elig$hvy.drk[(nhis_elig$drk_pdy > 1) & (nhis_elig$RIAGENDR == "Female")] <- 2
nhis_elig$hvy.drk[(nhis_elig$drk_pdy > 2) & (nhis_elig$RIAGENDR == "Male")]   <- 2
nhis_elig$hvy.drk[nhis_elig$alc2 == 0] <- 0

# DXA: not collected in NHIS ---
nhis_elig$DXDTOFAT <- NA_real_
nhis_elig$DXXTRFAT <- NA_real_
nhis_elig$DXDTRPF  <- NA_real_
nhis_elig$DXDTOPF  <- NA_real_

nhis_elig$SDMVPSU  <- nhis_elig$PSU
nhis_elig$SDMVSTRA <- nhis_elig$STRATA
nhis_elig$source   <- "nhis"

nhis_keep <- c("SDMVPSU", "SDMVSTRA", "WTINT2YR", "RIDAGEYR", "RIAGENDR", "race4",
               "hvy.drk", "smk3", "phys.all", "mortality", "educ5",
               "DXDTOFAT", "DXXTRFAT", "DXDTRPF", "DXDTOPF",
               "age_match_mort", "slfBMI", "height", "source")
nhis_analysis <- nhis_elig[, nhis_keep]

rm(nhis_raw, nhis_elig)
cat(sprintf("NHIS analysis: n = %d\n", nrow(nhis_analysis)))

# ==============================================================================
# 8. Save processed data (main analysis: cycles 01+03+05)
# ==============================================================================
saveRDS(nhanes_analysis, "data/processed/nhanes_analysis.rds")
saveRDS(nhis_analysis,   "data/processed/nhis_analysis.rds")

# ==============================================================================
# 9. Build 4-cycle EDA dataset (includes 1999-2000; for diagnostic figures only)
#    The calibration analysis uses nhanes_analysis.rds (cycles 01+03+05 only).
# ==============================================================================
cat("Building 4-cycle EDA dataset for diagnostic figures...\n")

eda00 <- nhanes_mort00_eda   # select_cmn() already applied; has ALQ101, total_MET_wk

# Alcohol (same recoding as applied to nhanes_raw above)
eda00$alchl <- 2 - eda00$ALQ101
eda00$alchl[eda00$ALQ110 == 1] <- 1
eda00$alchl[eda00$ALQ110 == 9] <- 0
dys00 <- eda00$ALQ120U
dys00[(eda00$ALQ101 == 2) & (eda00$ALQ110 == 2)] <- 0
dys00[eda00$ALQ120Q == 0]    <- 1
dys00[eda00$ALQ120U == 1]    <- 7
dys00[eda00$ALQ120U == 2]    <- 365.25 / 12
dys00[eda00$ALQ120U == 3]    <- 365.25
eda00$ach_freq <- eda00$ALQ120Q / dys00
eda00$ach_freq[eda00$alchl == 0] <- 0
eda00$ALQ130[eda00$ALQ130 > 82] <- NA
eda00$ALQ130[eda00$ALQ130 > 15] <- 15
eda00$ALQ130[(eda00$ALQ101 == 2) & (eda00$ALQ110 == 2)] <- 0
eda00$drk_pdy <- eda00$ach_freq * eda00$ALQ130
eda00$alc3 <- 1
eda00$alc3[((eda00$ALQ101 == 1) | (eda00$ALQ110 == 1)) &
            (eda00$ALQ120Q <= 365) & (eda00$ALQ120Q > 0)] <- 3
eda00$alc3[((eda00$ALQ101 == 1) | (eda00$ALQ110 == 1)) &
            (eda00$ALQ120Q == 0)] <- 2
eda00$drk_pdy[eda00$alc3 == 2] <- 0
eda00$drk_pdy[(eda00$alc3 == 3) & is.na(eda00$drk_pdy)] <-
  mean(eda00$drk_pdy[(eda00$alc3 == 3) & !is.na(eda00$drk_pdy)])
eda00$alc2 <- eda00$alc3
eda00$alc2[((eda00$drk_pdy * 365.25) < 12) | (eda00$alc3 == 2)] <- 1
eda00$alc2[eda00$alc2 == 3] <- 2
eda00$alc2 <- eda00$alc2 - 1
eda00$hvy.drk <- 1
eda00$hvy.drk[(eda00$drk_pdy > 1) & (eda00$RIAGENDR == 2)] <- 2
eda00$hvy.drk[(eda00$drk_pdy > 2) & (eda00$RIAGENDR == 1)] <- 2
eda00$hvy.drk[eda00$alc2 == 0] <- 0

# Smoking
eda00$smk3 <- 0
eda00$smk3[(eda00$SMQ020 == 1) & (eda00$SMQ040 %in% c(3, 7))] <- 1
eda00$smk3[(eda00$SMQ020 == 1) & (eda00$SMQ040 %in% c(1, 2))] <- 2
eda00$smk3[is.na(eda00$SMQ020) | eda00$SMQ020 %in% c(7, 9)] <- 0
eda00$smk3[!is.na(eda00$SMQ050Q)] <- 1
eda00$smk3[!is.na(eda00$SMQ640) | !is.na(eda00$SMQ650)] <- 2

# Mortality, age, DXA, BMI, height filters
eda00 <- subset(eda00, eligstat == 1)
eda00$mortality <- ifelse(
  eda00$mortstat == 1 & (eda00$permth_exm <= 120 | is.na(eda00$permth_exm)), 1, 0
)
eda00 <- eda00 %>% filter(RIDAGEYR >= 40 & RIDAGEYR <= 69, !is.na(DXDTOFAT))
eda00$WHD020[eda00$WHD020 > 1000] <- NA
eda00$WHD010[eda00$WHD010 > 1000] <- NA
eda00$slfBMI <- eda00$WHD020 * 0.45359237 / (eda00$WHD010 * 0.0254)^2
eda00$slfBMI[is.na(eda00$slfBMI) & !is.na(eda00$BMXBMI)] <-
  eda00$BMXBMI[is.na(eda00$slfBMI) & !is.na(eda00$BMXBMI)]
eda00 <- eda00 %>% filter(!is.na(slfBMI) & slfBMI >= 18.5)
eda00$height <- eda00$WHD010 * 2.54
eda00$height[is.na(eda00$height) & !is.na(eda00$BMXHT)] <-
  eda00$BMXHT[is.na(eda00$height) & !is.na(eda00$BMXHT)]
eda00 <- eda00 %>% filter(height >= 120 & height <= 220)

# Gender and race (as factors matching nhanes_analysis encoding)
eda00$RIAGENDR <- factor(eda00$RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
eda00$race4 <- eda00$RIDRETH1
eda00$race4[eda00$race4 == 1] <- 2
eda00$race4 <- as.factor(eda00$race4)
levels(eda00$race4) <- c("Hispanic", "NH-White", "NH-Black", "NH-Other")
eda00$race4 <- factor(eda00$race4, levels = c("NH-White", "NH-Black", "Hispanic", "NH-Other"))

# Education
eda00$DMDEDUC2[eda00$DMDEDUC2 %in% c(7, 9)] <- 4
eda00$DMDEDUC3[eda00$DMDEDUC3 %in% c(7, 9)] <- 13
eda00$educ5 <- eda00$DMDEDUC2
eda00$educ5[is.na(eda00$DMDEDUC2) & eda00$DMDEDUC3 %in% c(0:8, 55, 66)] <- 1
eda00$educ5[is.na(eda00$DMDEDUC2) & eda00$DMDEDUC3 %in% 9:12] <- 2
eda00$educ5[is.na(eda00$DMDEDUC2) & eda00$DMDEDUC3 %in% c(13, 14)] <- 3
eda00$educ5[is.na(eda00$DMDEDUC2) & eda00$DMDEDUC3 == 15] <- 4

# phys.all normalized within cycle 00
eda00$total_MET_wk[is.na(eda00$total_MET_wk)] <- 0
eda00$phys.all <- eda00$total_MET_wk / sum(eda00$total_MET_wk) * nrow(eda00)

# WTINT2YR4 / WTINT2YR3 are NA for cycle 00 (not used in 3- or 2-cycle analysis)
eda00$WTINT2YR4 <- NA_real_
eda00$WTINT2YR3 <- NA_real_
eda00$age_match_mort <- cut(eda00$RIDAGEYR, breaks = c(18, 25, 35, 45, 55, 70),
                            right = TRUE, include.lowest = TRUE)

# Align column sets and combine (convert factors to character for safe rbind)
to_chr <- function(df) {
  df$RIAGENDR <- as.character(df$RIAGENDR)
  df$race4    <- as.character(df$race4)
  df$age_match_mort <- as.character(df$age_match_mort)
  df
}
eda00_sel        <- to_chr(eda00[, keep_cols[keep_cols %in% names(eda00)]])
nhanes_analysis2 <- to_chr(nhanes_analysis)
eda00_sel[, setdiff(keep_cols, names(eda00_sel))] <- NA_real_
eda00_sel <- eda00_sel[, keep_cols]

nhanes_4cycles <- rbind(eda00_sel, nhanes_analysis2)
nhanes_4cycles$RIAGENDR <- factor(nhanes_4cycles$RIAGENDR,
                                  levels = c("Male", "Female"))
nhanes_4cycles$race4 <- factor(nhanes_4cycles$race4,
                               levels = c("NH-White", "NH-Black", "Hispanic", "NH-Other"))

cat(sprintf("4-cycle EDA dataset: n = %d  (cycle 1: %d)\n",
    nrow(nhanes_4cycles), sum(nhanes_4cycles$cycle == 1, na.rm = TRUE)))

saveRDS(nhanes_4cycles, "data/processed/nhanes_4cycles.rds")

rm(nhanes_mort00_eda, eda00, eda00_sel, nhanes_analysis2, dys00, nhanes_4cycles)

cat("\nSaved:\n")
cat("  data/processed/nhanes_analysis.rds\n")
cat("  data/processed/nhis_analysis.rds\n")
cat("  data/processed/nhanes_4cycles.rds\n")
cat("Script 01 complete.\n")
