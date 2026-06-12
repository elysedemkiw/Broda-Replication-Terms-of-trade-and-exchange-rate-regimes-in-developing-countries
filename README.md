# Terms of Trade and Exchange Rate Regimes, Broda (2004) Replication & Extension

A panel-VAR replication of Christian Broda, *"Terms of trade and exchange rate regimes in developing countries"* (Journal of International Economics 63, 2004, 31–58), plus three extensions that push the sample to 2023 and test what financialization, capital-account openness, and inflation targeting do to Broda's "floats insulate" result.

Author: Elyse Demkiw · Columbia University

---

## What this repo does

Broda's claim, in one line: after a negative terms-of-trade (ToT) shock, **pegs take the hit through output and floats take it through the exchange rate**. I reproduce that on a modernized dataset for his 1973–1996 window, then extend the panel to 2023 to ask whether the mechanism still holds once you account for things Broda's sample mostly predates — high capital mobility, foreign-currency balance-sheet exposure, and inflation-targeting central banks.

Everything is estimated as a **regime-stratified panel VAR(4)** on log-differenced series, with the ToT shock ordered first (treated as exogenous to the small open economy), country fixed effects removed by within-regime demeaning, and 90% confidence intervals from a country-block bootstrap.

---

## Repo layout

```
.
├── broda_panel_1973_2023_UPDATED_v4.csv   # master panel: 73 countries, 1973–2023, 61 cols
│
├── broda_rep.r                  # BASELINE: Broda window (1973–1996), peg vs float IRFs
├── broda_filter_it.r            # EXTENSION: drop every inflation-targeting country, 1975–2015
├── broda_post_no_manufacturing.r# EXTENSION: 1997–2020, drop manufacturing-transition countries
├── ka_fx_explore.r              # descriptive trends in kaopen & FX exposure (no VAR)
│
├── broda_rep.png                # baseline figure
├── broda_fin_split_kaopen.png   # within-regime median split on capital-account openness
└── broda_fin_split_fx_exposure.png # within-regime median split on FX exposure (IT excluded)
```

Run any script from the folder containing the CSV. Each `broda_*.r` is self-contained: it filters the panel, fits the VARs, bootstraps, writes a PNG and a `broda_table_*.csv`

```r
# dependencies
install.packages(c("vars", "ggplot2", "zoo", "dplyr", "knitr", "ggbeeswarm"))
source("broda_rep.r")
```

---

## The VAR specification (all scripts share this core)

| Choice | What I do | Why |
|---|---|---|
| Endogenous block | `dln_tot_best`, `dlny`, `dln_rer_merged`, `d_ln_cpi` | Broda's four variables: ToT, real GDP, real exchange rate, CPI |
| Ordering / ID | ToT ordered **first** | small-open-economy assumption — the country can't move world prices, so ToT is exogenous to its own GDP/RER/CPI |
| Lags | `p = 4` | Broda's choice |
| Transformation | first log-differences, then **within-country demeaned inside each regime** | series are I(1); differencing kills level fixed effects, demeaning removes country-specific drift in the differenced series |
| Exogenous controls | `trade_openness`, `dlng` (gov't spending growth), `cfa1994` | Broda's significant controls + the CFA devaluation dummy (see below) |
| Regimes | estimated **separately** for peg and float | this is Broda's "Method A" (his footnote 19 robustness), not the pooled-with-interactions baseline — cleaner to read, same story |
| IRFs | `cumulative = TRUE`, then `× −1 × 100` | cumulate because the VAR is in differences; flip sign so the plot shows the response to a **fall** in ToT; scale to % |
| Inference | country-block bootstrap, `B = 200`, 90% percentile CI (5th/95th), `set.seed(42)` | resamples whole countries with replacement (renamed to avoid collisions) so CIs respect within-country correlation; this is the cluster-robust analogue of Broda's parametric Monte Carlo |

A cell is only estimated if it has ≥ 30 observations after filtering; otherwise the function returns `NULL` and the bootstrap skips it.

---

## Sample filters

Applied in every estimation script, in this order:

1. **Regime = peg or float only.** `intermediate` (managed/crawling) and `exclude` (freely-falling, dual-market) are dropped.
2. **Regime stability.** A country-year survives only if its regime is identical and non-missing at **t−2, t−1, t, t+1, t+2** (`b_stable`). This is the ±2-year window — it drops transitional years around a regime switch so the IRFs aren't contaminated by the switch itself.
3. **Outlier trim.** `abs(dlny) < 0.20` drops the worst hyperinflation/collapse country-years (≥ 20% real GDP swings), which otherwise blow up the RER block.
4. **Complete cases** on the four endogenous variables plus the controls.

The estimation window is the only thing that changes across scripts: `≤ 1996` for the baseline, `1975–2015` for the IT filter, `1997–2020` for the manufacturing-exclusion extension.

---

## Variable construction & sources

The panel was assembled from modern vintages of every series — Broda built his in 2004, and 20 years of revisions and better methods mean nothing here is his exact data. Where my series and his diverge, that's expected and noted.

### Terms of trade — `tot_best` → `dln_tot_best` (the shock)
- **Sources:** World Bank WDI **net barter terms of trade** index, spliced/backchained with the IMF **commodity terms-of-trade** (CTOT, Gruss–Kebhaj).
- **Backchaining:** WDI is the preferred series but doesn't cover every country back to 1973, and CTOT *ends in 1995*. So `tot_best` is a chained composite — I extend the preferred series backward using the growth rates of the alternative where the preferred one is missing, rather than level-splicing (which would import the alternative's base-year level). The result is one continuous `ln_tot_best` per country; `dln_tot_best` is its first difference.
- **Caveat:** CTOT fits commodity exporters well and commodity-light countries poorly. Robustness was checked on the 20 lowest-R² countries, and the Granger test (below) is the main defense of treating ToT as exogenous.

### Real GDP — `rgdpna` → `dlny`
- **Source:** Penn World Table 10.01, `rgdpna` (real GDP at constant 2017 national prices).
- `dlny = Δ ln(rgdpna)`. Annual; a few minor pre-1980 gaps.

### CPI — `cpi_index` → `d_ln_cpi`
- **Sources:** World Bank WDI headline CPI, spliced with IMF WEO (Oct 2025 vintage) to fill 1970s gaps and reconcile base years.
- `d_ln_cpi = Δ ln(CPI)`. This is a log-difference, **not** the WDI inflation-rate percentage — important, they are not interchangeable.

### Real exchange rate — `rer_merged` → `dln_rer_merged`
- **Headline source used in the VAR:** Darvas / Bruegel **narrow (65-country) CPI-deflated REER** (`rer_merged_src = darvas65` for ~97% of rows).
- **Why not Broda's series:** Broda's original EREER IX (IMF) is unavailable. The panel also carries IMF REER, CEPII, the broad Bruegel REER, and a bilateral-vs-USD fallback (`bilateral_rer`) for comparison, but the merged VAR series is Darvas-narrow because it has the best coverage.
- **Big caveat:** IMF REER and bilateral RER are poorly correlated in Δln for a lot of the sample (worse in the short run than the long run). RER is the variable I trust least and the one where my numbers depart most from Broda's. Individual RER responses are often insignificant at 10%; the **peg–float gap** is the robust object, not the level.

### Exchange-rate regime — `broda_regime`
- **Source:** Ilzetzki–Reinhart–Rogoff (IRR) **coarse** de facto classification (`irr_coarse`), mapped:
  - `1 → peg`, `2 → intermediate`, `3 → float`
  - `4, 5, 6 → exclude` (freely-falling, dual/parallel-market, managed — these are noise, not regime information)
- **Robustness alternative:** Levy-Yeyati–Sturzenegger (`lys_class`) is carried in the panel for a de-facto cross-check. Broda's own primary classification is Ghosh et al. (de jure + frequent-adjuster); IRR is purely de facto and will disagree with it on 30–40% of country-years, which is exactly the LYS robustness check Broda himself ran.
- **Note:** because IRR is de facto and crisis-sensitive, it tosses a lot of 1980s Latin America into "freely falling." That shrinks the float sample relative to Broda but is arguably the right call.

### Controls
- **Trade openness** — `trade_openness`: WB WDI, (exports + imports)/GDP. ~76% coverage.
- **Government spending** — `dlng`: growth of real gov't consumption. **Interpolated** within country (`na.approx`), then country-mean-filled, then global-mean-filled for any residual gaps, so the control never drops an otherwise-valid country-year. One of only two controls Broda found systematically significant (openness is the other).
- **CFA 1994 dummy** — `cfa1994`: 1 for the 10 CFA-franc-zone countries (BFA, COG, CAF, CMR, CIV, TCD, GAB, SEN, NER, MLI) in **1994 only**, the year of the 50% CFA devaluation. Included as an exogenous regressor in the **peg** VAR so that one massive common devaluation doesn't masquerade as a generic peg response. (The float VAR has no CFA-1994 observations, so the dummy is only added when `sum(cfa1994) > 0`.)

### Extension variables
- **Capital-account openness** — `kaopen` / `ka_open`: Chinn–Ito KAOPEN, de jure, the workhorse measure. `ka_open` is the 0–1 normalized version.
- **FX exposure** — `fx_exposure = ewn_debt_liab / ewn_gdp_usd`: balance-sheet vulnerability from the Lane–Milesi-Ferretti **External Wealth of Nations** dataset. `net_fx_exposure` subtracts FX reserves; `net_ext_liab` is `−net_IIP/GDP`. Conceptually this is the BLS/BGJS "net FX position as % of GDP" idea (debt liabilities − reserves), built from EWN components.
- **Inflation targeting** — `it_ever`, `it_active`, `it_adoption_year`, plus central-bank mandate flags (`mandate_price_stability`, `mandate_dual_growth`): hand-coded from adoption dates. 26 countries in the panel ever adopt IT; far fewer survive into the stable peg/float estimation sample (most IT adopters classify as managed floats and get dropped as `intermediate`).

---

## Extensions

### 1. Extend to 2023 + drop the manufacturing transitions — `broda_post_no_manufacturing.r`
The full 1973–2023 panel produces an **anomaly**: GDP still behaves like Broda, but floats show a *massive real appreciation* (wrong sign) and CPI hovers at zero. Two suspects, both tested:

- **Manufacturing transition.** A handful of countries stopped being commodity price-takers and became manufacturing exporters with real market power, which breaks the small-open-economy / exogenous-ToT assumption. This script **drops CHN, KOR, TUR, POL, MEX** and re-estimates on **1997–2020**.
- This is the structural-break check behind the deck's "rerun the market-share test excluding CHN, POL, MEX, KOR, TUR" line.

### 2. Strip out inflation targeters — `broda_filter_it.r`
Filters `it_ever != 1` (every country that *ever* adopted IT is removed, not just its IT years), window **1975–2015**. The motivation: an IT central bank mechanically anchors CPI near its target, so a negative ToT shock can't show up in inflation the way Broda's model predicts. Excluding IT countries **reverses the sign of the RER response** in the extended sample — evidence that the modern anomaly is partly an artifact of monetary regime, not exchange-rate regime.

Open hypothesis left in the deck: if the tradables share α is high, depreciation raises CPI rather than stabilizing non-tradable prices, so high-tradables IT central banks should undershoot their target *less* after a ToT shock.

### 3. Capital-account openness & FX exposure splits — `broda_fin_split_kaopen.png`, `broda_fin_split_fx_exposure.png`
The financialization extension: Broda assumes free capital mobility, but in his window many of these countries had closed capital accounts. Greater borrowing access should let a country smooth a *temporary* ToT shock, muting the output response — but high FX-denominated debt should *amplify* the cost of the depreciation that's supposed to do the buffering (balance-sheet / "fear of floating" channel).

- **Design:** within each regime, **median-split** countries into high vs low on the financialization variable, then fit four separate VARs (peg-low, peg-high, float-low, float-high) and overlay the IRFs. If financialization neutralizes Friedman's mechanism, the high-exposure float IRF should look more peg-like.
- **kaopen split** (`broda_fin_split_kaopen.png`): within-regime median split on Chinn–Ito, 1973–1996 window, 90% bootstrap CI.
- **fx_exposure split** (`broda_fin_split_fx_exposure.png`): median split on `fx_exposure`, **IT countries excluded**, ±1-year stability, 90% cluster-bootstrap CI.
- **Honesty note (the thing to not fool yourself about):** the kaopen and fx_exposure heterogeneity is driven largely by the *same* ~12-country overlap cluster (many CFA + a few large dollarized EMs). So this is **not** a clean "modern financialization broke Broda" finding — it's "a specific group of closed-and-dollarized countries drives both splits." To actually isolate a *modern* threshold effect you'd need a pre/post-1996 temporal split of the FX-exposure result, which the current spec does not do. `ka_open` and `fx_exposure` are only weakly correlated, so disagreement between the two splits is informative, not a bug.

---

## Notes on the code (gotchas)

- **Sign convention.** IRFs are multiplied by `−1`. Every plotted line is the response to a **negative** ToT shock. A line going down in the GDP panel = output falling. Don't double-flip it when reading.
- **Cumulative IRFs are non-negotiable.** The VAR is in differences; the level response is the cumulative sum. Reporting the non-cumulative IRF would be wrong.
- **Bootstrap is by country, not by observation.** `resample_countries()` draws whole countries with replacement and renames them (`CODE_i`) so a country drawn twice doesn't collapse into one demeaning group. This is what makes the CIs cluster-robust.
- **`dlng` interpolation runs before filtering** so gov't-spending gaps don't silently shrink the sample. If you move it after the filters you'll lose country-years.
- **`B = 200`** is fine for the shape of the CIs but is light for tail precision; bump it to 1000+ for anything going into a final table.
- **The CFA dummy only attaches to the peg VAR** (guarded by `sum(cfa1994) > 0`). If you re-window to a period without 1994, it drops automatically.
- **Significance stars** (`*`) mean the 90% bootstrap CI excludes zero at that horizon. The headline inference object is the **peg − float difference** (`DIFF`), not either regime alone.

---

## Replication scorecard (baseline, 1973–1996)

| Variable | Regime | Mine | Broda |
|---|---|---|---|
| Real GDP | Fixed | −1.54% | −1.9% |
| Real GDP | Flexible | −0.3% | −0.2% |
| RER | Fixed | 0.3% | 1.3% |
| RER | Flexible | 1.59% | 5.1% |
| CPI | Fixed | −1.05% (stronger than Broda) | −0.8% |
| CPI | Flexible | +0.95% | ~2% |

GDP replicates well and is significant short-run; CPI is significant short and long run and actually beats Broda's fit; RER is individually insignificant at 10% but the peg–float gap holds (≈0.33× the size of Broda's RER gap). The qualitative Friedman result — **pegs pay in output, floats pay in the exchange rate** — reproduces.

---

## Known limitations

- RER is the weak link: source disagreement in Δln means the level responses are noisy and shouldn't be over-read.
- The extension splits are concentrated in a small overlap cluster, so causal claims about "financialization" are not yet identified — they're suggestive heterogeneity.
- IRR-coarse vs Broda's Ghosh classification reshuffles ~a third of country-years; this is a deliberate de-facto robustness choice, but it does make the samples non-identical to his.
- `B = 200` bootstrap draws — adequate for figures, light for final tables.
