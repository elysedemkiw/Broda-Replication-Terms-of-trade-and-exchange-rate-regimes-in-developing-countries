# Terms of Trade and Exchange Rate Regimes, Broda (2004) Replication and Extension

Panel VAR replication of Broda (2004) on whether floating exchange rates insulate output from terms of trade shocks, extended from 1973 to 2023, plus three extensions testing financialization, capital account openness, and inflation targeting

Author: Elyse Demkiw, Columbia University

Slides recommended for context: [Elyse_Demkiw__ToT_and_Exchange_Rate_Regimes_Presentation.pdf](https://github.com/user-attachments/files/28878193/Elyse_Demkiw__ToT_and_Exchange_Rate_Regimes_Presentation.2.pdf)

## What it does

Reproduces Broda's result that after a negative terms of trade shock pegs take the hit through output and floats take it through the exchange rate, then asks whether that still holds once you add capital mobility, FX balance sheet exposure, and inflation targeting

Method is a regime stratified panel VAR(4) on log differenced series, terms of trade ordered first as exogenous, country fixed effects demeaned within regime, 90% CIs from a country block bootstrap

## Repo layout

```
broda_panel_1973_2023_UPDATED_v4.csv   master panel, 73 countries, 1973 to 2023, 61 cols
broda_rep.r                            baseline, 1973 to 1996, peg vs float IRFs
broda_filter_it.r                      drop inflation targeters, 1975 to 2015
broda_post_no_manufacturing.r          1997 to 2020, drop manufacturing transitions
ka_fx_explore.r                        descriptive kaopen and FX trends, no VAR
broda_rep.png                          baseline figure
broda_fin_split_kaopen.png             median split on capital openness
broda_fin_split_fx_exposure.png        median split on FX exposure, IT excluded
```

Run any script from the folder holding the CSV, each one filters the panel, fits the VARs, bootstraps, and writes a PNG plus a broda_table csv

```r
install.packages(c("vars", "ggplot2", "zoo", "dplyr", "knitr", "ggbeeswarm"))
source("broda_rep.r")
```

## VAR spec

Variables dln_tot_best, dlny, dln_rer_merged, d_ln_cpi with terms of trade ordered first, p=4 lags, peg and float estimated separately, IRFs cumulated then flipped to show a fall in terms of trade scaled to percent, controls trade_openness, dlng, cfa1994, inference by country block bootstrap B=200 seed 42, a cell needs 30 or more obs or it is skipped

## Sample filters

Peg or float only, regime stable across a plus or minus 2 year window, abs(dlny) below 0.20 to trim collapses, complete cases on all four variables plus controls

## Variables and sources

* `dln_tot_best` the terms of trade shock and the exogenous driver everything responds to, source WDI net barter terms of trade backchained with IMF CTOT
* `dlny` real GDP growth, source PWT 10.01 rgdpna
* `d_ln_cpi` inflation as a log difference of CPI, source WDI spliced with IMF WEO
* `dln_rer_merged` change in the real exchange rate and the noisiest variable so read the peg float gap not the level, source Darvas Bruegel narrow REER
* `broda_regime` peg or float classification, source Ilzetzki Reinhart Rogoff coarse de facto, 1 is peg, 3 is float, the rest excluded
* `trade_openness` exports plus imports over GDP, source WDI
* `dlng` real gov spending growth, source WDI gov consumption
* `cfa1994` dummy for the 1994 CFA devaluation, attached to the peg VAR only, hand coded
* `kaopen` and `ka_open` capital account openness, source Chinn Ito KAOPEN
* `fx_exposure` FX balance sheet vulnerability, source Lane Milesi-Ferretti External Wealth of Nations
* `it_ever` and `it_active` inflation targeting flags, hand coded from adoption dates

Note we use log differences so the model tracks short run changes rather than levels, and cumulating the IRFs recovers the level response

## Extensions

1 drop manufacturing transitions CHN KOR TUR POL MEX on 1997 to 2020, fixes the wrong sign float appreciation in the full panel
2 strip every inflation targeter on 1975 to 2015, which flips the sign of the RER response
3 within regime median splits on capital openness and FX exposure, suggestive heterogeneity driven by a small closed and dollarized country cluster, not yet a clean financialization finding

## Scorecard, 1973 to 1996, mine vs Broda

real GDP fixed −1.54 vs −1.9, real GDP flexible −0.3 vs −0.2, RER fixed 0.3 vs 1.3, RER flexible 1.59 vs 5.1, CPI fixed −1.05 vs −0.8, CPI flexible +0.95 vs about 2, all percent

GDP and CPI replicate well, RER is individually insignificant but the peg float gap holds, so the Friedman result reproduces

## Limitations

RER is noisy from source disagreement, the extension splits ride a small overlap cluster so causal claims are not identified, IRR vs Broda's Ghosh classification reshuffles about a third of country years, and B=200 is fine for figures but light for final tables


## Known limitations

- RER is the weak link: source disagreement in Δln means the level responses are noisy and shouldn't be over-read.
- The extension splits are concentrated in a small overlap cluster, so causal claims about "financialization" are not yet identified — they're suggestive heterogeneity.
- IRR-coarse vs Broda's Ghosh classification reshuffles ~a third of country-years; this is a deliberate de-facto robustness choice, but it does make the samples non-identical to his.
- `B = 200` bootstrap draws — adequate for figures, light for final tables.
