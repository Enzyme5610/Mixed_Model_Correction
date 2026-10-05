# Mixed Model Correction

Browser-based tool for testing treatment effects on physiology or qPCR data
with a linear mixed model, accounting for Line and Batch variability.
Includes pairwise comparisons and downloadable plots. No installation needed.
Data is local to your platform (PC/Mac/Linux).

**App:** https://enzyme5610.github.io/Mixed_Model_Correction/

## How to use

1. Upload a CSV (format below).
2. Check the parameter list. Numeric columns are preselected; unselect any
   that aren't parameters (e.g. age, culture days, well number).
3. Choose the pairwise adjustment (Tukey or Bonferroni).
4. Click **Run models**, then download results, pairwise tables and plots.

The first load takes about 30 seconds while R loads in the browser.

## CSV format

One row per measurement, with columns named exactly `Tx`, `Line` and
`Batch` (any position), plus one column per parameter:

```
Tx,Line,Batch,GRIA1,GAD1
Control,L1,B1,4.21,6.10
Control,L1,B2,4.35,6.02
AD,L5,B1,5.02,7.44
```

- **Line:** one label per actual cell line.
- **Batch:** numbered within each line (L1's first batch is B1).
- **qPCR:** enter ΔCt values. Average technical replicates or enter each as
  its own row.

## Model

`parameter ~ Tx + (1 | Line/Batch)`, fit with `lmerTest::lmer`, with Tx
tested by a Type II F test using Kenward-Roger degrees of freedom. With one
Line the model uses `(1 | Batch)`; with one Batch, or one row per line and
batch, `(1 | Line)`. Pairwise comparisons use `emmeans`.

If lines share batches (all lines run together in each batch, e.g. on the
same qPCR plate), answer **"Yes"** to **"Do lines share batches?"** to use
the crossed model
`parameter ~ Tx + (1 | Line) + (1 | Batch)`. Batch labels must then mean the
same run for every line.

## Plots

Dots, bars or violins showing each sample and the model mean ± 95% CI. The
Y axis can show values as entered, relative to a reference group (linear
data, e.g. physiology), or as fold change 2^-ΔΔCt (ΔCt data, qPCR).
Statistics always use the values as entered.

## Credits

- Original R script: **Dr. Luis Gustavo Hernandez Carballo**
- Shiny app and visualizations: **Prachetas Jai Patel**

If you use this tool in a publication, poster or presentation, please
acknowledge both authors (see **Cite this repository** on GitHub).

Released under the [MIT License](LICENSE).

## References

- Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects
  models using lme4. *Journal of Statistical Software* 67(1):1–48.
- Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests
  in linear mixed effects models. *Journal of Statistical Software*
  82(13):1–26.
- Halekoh U, Højsgaard S (2014). A Kenward-Roger approximation and parametric
  bootstrap methods for tests in linear mixed models – the R package
  pbkrtest. *Journal of Statistical Software* 59(9):1–32.
- Kenward MG, Roger JH (1997). Small sample inference for fixed effects from
  restricted maximum likelihood. *Biometrics* 53(3):983–997.
- Lenth RV. emmeans: Estimated Marginal Means, aka Least-Squares Means.
  R package.
- Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data
  using real-time quantitative PCR and the 2^-ΔΔCT method. *Methods*
  25(4):402–408.
