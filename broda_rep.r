library(vars)
library(ggplot2)
library(zoo)

df <- read.csv("broda_panel_1973_2023_UPDATED_v4.csv")
df <- df[order(df$countrycode, df$year), ]



df$dlng <- ave(df$dlng, df$countrycode, FUN = function(x) {
  if (sum(!is.na(x)) >= 2) na.approx(x, na.rm = FALSE) else x
})
df$dlng <- ave(df$dlng, df$countrycode, FUN = function(x) {
  if (any(is.na(x)) && any(!is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
df$dlng[is.na(df$dlng)] <- mean(df$dlng, na.rm = TRUE)

df$b_regime <- NA_character_
df$b_regime[df$broda_regime == "peg"]   <- "peg"
df$b_regime[df$broda_regime == "float"] <- "float"

#Lags

df$b_lag1  <- ave(df$b_regime, df$countrycode, FUN = function(x) c(NA, head(x, -1)))
df$b_lag2  <- ave(df$b_regime, df$countrycode, FUN = function(x) c(NA, NA, head(x, -2)))
df$b_lead1 <- ave(df$b_regime, df$countrycode, FUN = function(x) c(tail(x, -1), NA))
df$b_lead2 <- ave(df$b_regime, df$countrycode, FUN = function(x) c(tail(x, -2), NA, NA))

#Stability restriction 
                  
df$b_stable <- !is.na(df$b_regime) &
  !is.na(df$b_lag1) & !is.na(df$b_lag2) &
  !is.na(df$b_lead1) & !is.na(df$b_lead2) &
  df$b_regime == df$b_lag1  & df$b_regime == df$b_lag2 &
  df$b_regime == df$b_lead1 & df$b_regime == df$b_lead2

#CFA dummies

cfa_countries <- c("BFA","COG","CAF","CMR","CIV","TCD","GAB","SEN","NER","MLI")
df$cfa1994 <- as.integer(df$countrycode %in% cfa_countries & df$year == 1994)

df <- df[!is.na(df$b_stable) & df$b_stable, ]
df <- df[df$b_regime %in% c("peg","float"), ]
df <- df[df$year <= 1996, ]
df <- df[abs(df$dlny) < 0.20, ]

needed <- c("dln_tot_best","dlny","dln_rer_merged","d_ln_cpi",
            "trade_openness","dlng","cfa1994","countrycode","b_regime")
df <- df[complete.cases(df[, needed]), ]

cat("Estimation sample:", range(df$year), "\n")

peg_data <- df[df$b_regime == "peg",   ]
flt_data <- df[df$b_regime == "float", ]

cat("Peg obs:", nrow(peg_data), " Float obs:", nrow(flt_data), "\n")

# core function: demean, fit VAR, return cumulative IRF paths (all responses)
fit_irf_all <- function(d, has_cfa) {
  vars_dm <- c("dln_tot_best","dlny","dln_rer_merged","d_ln_cpi","trade_openness","dlng")
  for (v in vars_dm)
    d[[v]] <- d[[v]] - ave(d[[v]], d$countrycode, FUN = function(x) mean(x, na.rm = TRUE))
  d <- d[complete.cases(d[, vars_dm]), ]
  if (nrow(d) < 30) return(NULL)

  endog <- as.matrix(d[, c("dln_tot_best","dlny","dln_rer_merged","d_ln_cpi")])
  exog  <- if (has_cfa && sum(d$cfa1994, na.rm = TRUE) > 0)
    as.matrix(data.frame(trade_openness = d$trade_openness,
                         dlng = d$dlng, cfa1994 = d$cfa1994))
  else
    as.matrix(data.frame(trade_openness = d$trade_openness, dlng = d$dlng))

  v <- tryCatch(VAR(endog, p = 4, type = "const", exogen = exog), error = function(e) NULL)
  if (is.null(v)) return(NULL)

  responses <- c("dlny","dln_rer_merged","d_ln_cpi")
  out <- list()
  for (r in responses) {
    ir <- tryCatch(
      irf(v, impulse = "dln_tot_best", response = r,
          n.ahead = 10, boot = FALSE, cumulative = TRUE),
      error = function(e) NULL
    )
    out[[r]] <- if (!is.null(ir)) as.numeric(-ir$irf$dln_tot_best[, r] * 100) else rep(NA, 11)
  }
  out
}

resample_countries <- function(d) {
  cc <- unique(d$countrycode)
  drawn <- sample(cc, length(cc), replace = TRUE)
  do.call(rbind, lapply(seq_along(drawn), function(i) {
    sub <- d[d$countrycode == drawn[i], ]
    sub$countrycode <- paste0(drawn[i], "_", i)
    sub
  }))
}

# point estimates
pt_peg <- fit_irf_all(peg_data, has_cfa = TRUE)
pt_flt <- fit_irf_all(flt_data, has_cfa = FALSE)

# cluster bootstrap (200 small)
B <- 200
responses <- c("dlny","dln_rer_merged","d_ln_cpi")

boot_peg  <- lapply(responses, function(r) matrix(NA, B, 11))
boot_flt  <- lapply(responses, function(r) matrix(NA, B, 11))
names(boot_peg) <- responses
names(boot_flt) <- responses

set.seed(42)
for (b in 1:B) {
  if (b %% 50 == 0) cat("bootstrap", b, "/", B, "\n")
  bp <- fit_irf_all(resample_countries(peg_data), has_cfa = TRUE)
  bf <- fit_irf_all(resample_countries(flt_data), has_cfa = FALSE)
  for (r in responses) {
    if (!is.null(bp)) boot_peg[[r]][b, ] <- bp[[r]]
    if (!is.null(bf)) boot_flt[[r]][b, ] <- bf[[r]]
  }
}

# build plot 
irf_df <- data.frame()
labels <- c(dlny = "Real GDP", dln_rer_merged = "Real Exchange Rate", d_ln_cpi = "CPI")

for (r in responses) {
  lo_peg <- apply(boot_peg[[r]], 2, quantile, 0.05, na.rm = TRUE)
  hi_peg <- apply(boot_peg[[r]], 2, quantile, 0.95, na.rm = TRUE)
  lo_flt <- apply(boot_flt[[r]], 2, quantile, 0.05, na.rm = TRUE)
  hi_flt <- apply(boot_flt[[r]], 2, quantile, 0.95, na.rm = TRUE)

  irf_df <- rbind(irf_df,
    data.frame(horizon = 0:10, point = pt_peg[[r]],
               lo = lo_peg, hi = hi_peg,
               variable = labels[r], regime = "Fixed"),
    data.frame(horizon = 0:10, point = pt_flt[[r]],
               lo = lo_flt, hi = hi_flt,
               variable = labels[r], regime = "Flexible")
  )
}

irf_df$variable <- factor(irf_df$variable, levels = c("Real GDP","Real Exchange Rate","CPI"))
irf_df$regime   <- factor(irf_df$regime,   levels = c("Fixed","Flexible"))

p <- ggplot(irf_df, aes(x = horizon)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(aes(y = point), linewidth = 0.8) +
  geom_point(aes(y = point), size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(variable ~ regime, scales = "free_y") +
  labs(x = "Years after shock", y = "% change",
       title = "Broda Year Range: Responses to a 1-SD Permanent Fall in Terms of Trade",
       subtitle = "Panel VAR, p=4, within-regime country-demeaned, CFA1994 control, 90% cluster-bootstrap CI") +
  theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"))


ggsave("broda_new.png", p, width = 7, height = 8, dpi = 150)
cat("Done.\n")

#inference tables
RESPONSE <- "d_ln_cpi"

pt_diff  <- pt_peg[[RESPONSE]] - pt_flt[[RESPONSE]]
boot_diff <- boot_peg[[RESPONSE]] - boot_flt[[RESPONSE]]

ci_peg  <- list(lo = apply(boot_peg[[RESPONSE]],  2, quantile, 0.05, na.rm = TRUE),
                hi = apply(boot_peg[[RESPONSE]],  2, quantile, 0.95, na.rm = TRUE))
ci_flt  <- list(lo = apply(boot_flt[[RESPONSE]],  2, quantile, 0.05, na.rm = TRUE),
                hi = apply(boot_flt[[RESPONSE]],  2, quantile, 0.95, na.rm = TRUE))
ci_diff <- list(lo = apply(boot_diff, 2, quantile, 0.05, na.rm = TRUE),
                hi = apply(boot_diff, 2, quantile, 0.95, na.rm = TRUE))

fmt <- function(pt, ci) data.frame(
  horizon = 0:10,
  irf     = round(pt, 2),
  lo_90   = round(ci$lo, 2),
  hi_90   = round(ci$hi, 2),
  sig     = ifelse(ci$lo > 0 | ci$hi < 0, "*", "")
)

cat("\nResponse:", RESPONSE, "| B =", B, "| * = 90% CI excludes zero\n")
cat("\nPEG\n");               print(fmt(pt_peg[[RESPONSE]],  ci_peg),  row.names = FALSE)
cat("\nFLOAT\n");             print(fmt(pt_flt[[RESPONSE]],  ci_flt),  row.names = FALSE)
cat("\nDIFF (peg - float)\n"); print(fmt(pt_diff, ci_diff), row.names = FALSE)


library(knitr)

HORIZONS <- c(impact = 1, short_run = 3, long_run = 11)

labels <- c(dlny = "Real GDP",
            dln_rer_merged = "Real Exchange Rate",
            d_ln_cpi = "CPI")

fmt_cell <- function(pt_val, boot_mat, h_idx) {
  pt  <- round(pt_val[h_idx], 2)
  lo  <- round(quantile(boot_mat[, h_idx], 0.05, na.rm = TRUE), 2)
  hi  <- round(quantile(boot_mat[, h_idx], 0.95, na.rm = TRUE), 2)
  sig <- if (!is.na(lo) && !is.na(hi) && (lo > 0 | hi < 0)) "*" else ""
  sprintf("%.2f%s", pt, sig)
}

rows <- list()
for (r in responses) {
  lab <- labels[r]
  rows[[length(rows)+1]] <- c(
    Variable   = lab,
    Regime     = "Fixed",
    Impact     = fmt_cell(pt_peg[[r]], boot_peg[[r]], HORIZONS["impact"]),
    Short_Run  = fmt_cell(pt_peg[[r]], boot_peg[[r]], HORIZONS["short_run"]),
    Long_Run   = fmt_cell(pt_peg[[r]], boot_peg[[r]], HORIZONS["long_run"])
  )
  rows[[length(rows)+1]] <- c(
    Variable   = "",
    Regime     = "Flexible",
    Impact     = fmt_cell(pt_flt[[r]], boot_flt[[r]], HORIZONS["impact"]),
    Short_Run  = fmt_cell(pt_flt[[r]], boot_flt[[r]], HORIZONS["short_run"]),
    Long_Run   = fmt_cell(pt_flt[[r]], boot_flt[[r]], HORIZONS["long_run"])
  )
}

tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(tbl) <- c("Variable","Regime","Impact (h=0)","Short Run (h=2)","Long Run (h=10)")

cat("\n")
print(kable(tbl, align = "llrrr",
            caption = sprintf(
              "Responses to a 1-SD Permanent Fall in ToT | peg: %d obs, float: %d obs | * = 90%% CI excludes zero",
              nrow(peg_data), nrow(flt_data))))

write.csv(tbl, "broda_table_rep.csv", row.names = FALSE)
cat("Table saved to broda_table_rep.csv\n")
