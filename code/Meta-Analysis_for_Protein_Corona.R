############################################################
# Random Effects Meta-Analysis of Protein Corona Detection
# Author: Alexa Canchola
# Advisor: Wei-Chun Chou
# Date: Jul 22, 2025
############################################################
#### 1.  Libraries  --------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

pacman::p_load(
  readxl, dplyr, tidyr, stringr, metafor,
  ggplot2, scales, knitr, kableExtra
)
#### 2.  Import & minimal cleaning  ----------------------------------------
data_path <- "PC-DB_for_Meta-Analysis - v2.xlsx"
raw <- read_excel(data_path, sheet = 1)

df <- raw %>%
  transmute(
    study_id = `Study ID`,
    np_id    = `NP Entry ID`,
    material = `NP Sub-Type`,
    size_nm  = as.numeric(`Size (nm)`),
    zeta_mV  = as.numeric(`Zeta Potential (mV)`),
    
    # protein intensity columns (keep raw names from sheet)
    APOE,
    APOB100 = APOB,
    C3      = CO3,
    CLU     = CLUS
  )
# convert intensities → 0 / 1 detection flags
prot_cols <- c("APOE", "APOB100", "C3", "CLU")
detectify  <- function(v) as.integer(!is.na(v) & v != 0 & v != "ND")
df <- df %>% mutate(across(all_of(prot_cols), detectify))

# derive material / size / charge strata
df <- df %>%
  mutate(
    mat_grp = case_when(
      str_detect(material, regex("silica", TRUE))        ~ "Silica",
      str_detect(material, regex("metal|oxide", TRUE))   ~ "Metal/Metal-oxide",
      TRUE                                               ~ "Lipid / Polymer"
    ),
    size_bin   = if_else(size_nm < 100, "<100 nm", "≥100 nm"),
    charge_bin = cut(
      zeta_mV,
      breaks  = c(-Inf, -20, 0, Inf),
      labels  = c("≤−20 mV", "−20–0 mV", ">0 mV")
    )
  )
#### 4.  Helper: random-effects pooled proportion  -------------------------
meta_prop_pool <- function(data) {
  agg <- data %>%
    group_by(study_id) %>%
    summarise(x = sum(detect, na.rm = TRUE),
              n = n(),
              .groups = "drop") %>%
    filter(n > 0)
  
  if (nrow(agg) < 2) return(NULL)
  
  esc <- escalc(measure = "PLO", xi = x, ni = n, data = agg)
  fit <- rma(yi, vi, method = "REML", data = esc)
  
  tibble(
    k       = fit$k,
    N       = sum(agg$n),
    prop    = transf.ilogit(fit$b),
    ci_low  = transf.ilogit(fit$ci.lb),
    ci_high = transf.ilogit(fit$ci.ub),
    I2      = fit$I2
  )
}
#### 5.  Define analysis strata  -------------------------------------------
strata_def <- tribble(
  ~mat_grp,              ~size_bin,  ~label,                       ~is_ref,
  "Silica",              "<100 nm",  "Silica <100 nm",             TRUE,
  "Silica",              "≥100 nm",  "Silica ≥100 nm",             FALSE,
  "Metal/Metal-oxide",   "<100 nm",  "Metal/Metal-oxide <100 nm",  FALSE,
  "Metal/Metal-oxide",   "≥100 nm",  "Metal/Metal-oxide ≥100 nm",  FALSE,
  "Lipid / Polymer",     "<100 nm",  "Lipid / Polymer <100 nm",    FALSE,
  "Lipid / Polymer",     "≥100 nm",  "Lipid / Polymer ≥100 nm",    FALSE
)
#### 6.  Loop over proteins to build main table  ---------------------------
results <- purrr::map(prot_cols, function(prot) {
  part <- df_long %>% filter(protein == prot)
  
  strata_def %>%
    rowwise() %>%
    mutate(
      pool = list(
        meta_prop_pool(
          part %>%
            filter(mat_grp  == .env$mat_grp,
                   size_bin == .env$size_bin)
        )
      )
    ) %>%
    unnest(pool) %>%
    mutate(protein = prot)
}) %>% bind_rows() %>% relocate(protein)
#### 7.  Fixed meta-regression contrasts  ----------------------------------
get_contrast_p <- function(df_long, prot, ref_label, strata_df) {
  
  # keep only rows for this protein
  d <- df_long %>% filter(protein == prot)
  
  # tag each NP row with its stratum label
  d_tagged <- purrr::pmap_dfr(
    strata_df,
    function(mat_grp, size_bin, label, ...) {
      d %>%
        filter(mat_grp == !!mat_grp, size_bin == !!size_bin) %>%
        mutate(stratum = label)
    }
  )
  
  if (nrow(d_tagged) == 0)
    return(tibble(label = strata_df$label, p_val = NA_real_))
  
  # aggregate study × stratum
  agg <- d_tagged %>%
    group_by(study_id, stratum) %>%
    summarise(x = sum(detect, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n > 0)
  
  if (length(unique(agg$stratum)) < 2)
    return(tibble(label = strata_df$label, p_val = NA_real_))
  
  # meta-regression
  agg$stratum <- factor(agg$stratum, levels = strata_df$label)
  agg$stratum <- relevel(agg$stratum, ref = ref_label)
  esc <- escalc(measure = "PLO", xi = x, ni = n, data = agg)
  fit <- rma(yi, vi, mods = ~ stratum, data = esc, method = "REML")
  
  tibble(
    label = strata_df$label[-1],          # exclude intercept
    p_val = summary(fit)$pval[-1]
  )
}

# compute p-values for each protein
pvals <- purrr::imap_dfr(
  split(results, results$protein),
  function(tbl, prot) {
    ref_lbl <- tbl$label[tbl$is_ref][1]
    get_contrast_p(df_long, prot, ref_lbl, strata_def) %>%
      mutate(protein = prot)
  }
)

# merge & format p column
tableX <- results %>%
  left_join(pvals, by = c("protein", "label")) %>%
  mutate(`p vs ref*` = case_when(
    is_ref        ~ "—ref—",
    is.na(p_val)  ~ "n.e.",
    p_val < 0.001 ~ "<0.001",
    TRUE          ~ formatC(p_val, digits = 3, format = "f")
  )) %>%
  select(protein, label, k, N,
         prop, ci_low, ci_high, I2, `p vs ref*`)
#### 8.  Display table  ----------------------------------------------------
kable(tableX, digits = 2)
write.csv(tableX, "meta_analysis_results.csv", row.names = FALSE)
# 1) how many entries survive each step
n_raw <- nrow(df)                         # 598

n_complete <- df %>%
  filter(!is.na(size_nm), !is.na(material)) %>% nrow()

n_in_strata <- df %>%
  filter(!is.na(size_nm), !is.na(material)) %>%
  mutate(mat_grp = case_when(
    str_detect(material, regex("silica", TRUE))            ~ "Silica",
    str_detect(material, regex("metal|oxide", TRUE))       ~ "Metal/Metal-oxide",
    TRUE                                                   ~ "Lipid / Polymer"
  )) %>%
  filter(mat_grp != "Other") %>%  nrow()

cat("Raw file:", n_raw,
    "\nAfter complete vars:", n_complete,
    "\nAfter strata filters:", n_in_strata, "\n")
## 8 rows missing size or material
df %>% filter(is.na(size_nm) | is.na(material))