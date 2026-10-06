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

If lines share batches (multiple lines run per batch, e.g. on the
same qPCR plate), answer **"Yes"** to **"Do lines share batches?"** to use
the crossed model
`parameter ~ Tx + (1 | Line) + (1 | Batch)`. Batch labels must then mean the
same run for every line.

## Plots

Dots, bars or violins showing each sample, with error bars as the model's
95% CI or SE, or the raw SEM or SD (caps, direction and mean marker are
adjustable). The Y axis can show values as entered, relative to a reference group (linear
data, e.g. physiology), or as fold change 2^-ΔΔCt (ΔCt data, qPCR).
Statistics always use the values as entered.

**Batch-adjusted values** (optional) subtract each batch's estimated shift
from the plotted values, using the model's random-effect estimates (BLUPs).
With the crossed model this removes run-wide (e.g. qPCR plate) shifts; with
the nested model, each line's batch-to-batch deviation. This is for display
only: statistics are unchanged, and adjusted values shouldn't be re-tested.

A reference (control) group sets the comparison direction in tables and plots.
Significance can be shown as p-values or stars (ns, *, **, ***, ****), and
dots can be colored by Line or Batch.

**All parameters in one figure** shows the selected parameters side by side on
one shared Y axis, with the groups next to each other for each parameter.

## Credits

- Original R script: **Dr. Luis Gustavo Hernandez Carballo**
- Shiny app and visualizations: **Prachetas Jai Patel**

If you use this tool in a publication, poster or presentation, please
acknowledge both authors (see **Cite this repository** on GitHub).

Released under the [MIT License](LICENSE).

## Acknowledgments

Loading screen animation: [loading-bar](https://github.com/loadingio/loading-bar)
by loading.io (MIT License, © 2017 loading.io).

## References

### Methods

- Aarts E, Verhage M, Veenvliet JV, Dolan CV, van der Sluis S (2014). A
  solution to dependency: using multilevel analysis to accommodate nested
  data. *Nature Neuroscience* 17(4):491–496.
- Kenward MG, Roger JH (1997). Small sample inference for fixed effects from
  restricted maximum likelihood. *Biometrics* 53(3):983–997.
- Schielzeth H, Nakagawa S (2013). Nested by design: model fitting and
  interpretation in a mixed model era. *Methods in Ecology and Evolution*
  4(1):14–24.
- Bolker B, et al. GLMM FAQ: Nested or crossed?
  https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#nested-or-crossed
- Tukey JW (1953). The problem of multiple comparisons. Unpublished
  manuscript, reprinted in *The Collected Works of John W. Tukey*, Vol. VIII
  (1994). Chapman & Hall.
- Kramer CY (1956). Extension of multiple range tests to group means with
  unequal numbers of replications. *Biometrics* 12(3):307–310.
- Dunn OJ (1961). Multiple comparisons among means. *Journal of the American
  Statistical Association* 56(293):52–64.
- Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data
  using real-time quantitative PCR and the 2^-ΔΔCT method. *Methods*
  25(4):402–408.

### Batch adjustment (plots)

- Robinson GK (1991). That BLUP is a good thing: the estimation of random
  effects. *Statistical Science* 6(1):15–32.
- Ruijter JM, Thygesen HH, Schoneveld OJLM, Das AT, Berkhout B, Lamers WH
  (2006). Factor correction as a tool to eliminate between-session variation
  in replicate experiments: application to molecular biology and
  retrovirology. *Retrovirology* 3:2.
- Hellemans J, Mortier G, De Paepe A, Speleman F, Vandesompele J (2007).
  qBase relative quantification framework and software for management and
  automated analysis of real-time quantitative PCR data. *Genome Biology*
  8(2):R19.
- Ruijter JM, Ruiz Villalba A, Hellemans J, Untergasser A, van den Hoff MJB
  (2015). Removal of between-run variation in a multi-plate qPCR experiment.
  *Biomolecular Detection and Quantification* 5:10–14.
- Johnson WE, Li C, Rabinovic A (2007). Adjusting batch effects in
  microarray expression data using empirical Bayes methods. *Biostatistics*
  8(1):118–127.
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma
  powers differential expression analyses for RNA-sequencing and microarray
  studies. *Nucleic Acids Research* 43(7):e47.

### Software

- R Core Team (2025). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  https://www.R-project.org/
- Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects
  models using lme4. *Journal of Statistical Software* 67(1):1–48.
- Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests
  in linear mixed effects models. *Journal of Statistical Software*
  82(13):1–26.
- Halekoh U, Højsgaard S (2014). A Kenward-Roger approximation and parametric
  bootstrap methods for tests in linear mixed models – the R package
  pbkrtest. *Journal of Statistical Software* 59(9):1–32.
- Lenth R, Piaskowski J (2026). emmeans: Estimated marginal means, aka
  least-squares means. R package version 2.0.4.
  doi:10.32614/CRAN.package.emmeans
- Chang W, Cheng J, Allaire JJ, Sievert C, Schloerke B, Aden-Buie G, Xie Y,
  Allen J, McPherson J, Dipert A, Borges B (2026). shiny: Web application
  framework for R. R package version 1.14.0. doi:10.32614/CRAN.package.shiny
- Sievert C, Cheng J, Aden-Buie G (2026). bslib: Custom 'Bootstrap' 'Sass'
  themes for 'shiny' and 'rmarkdown'. R package version 0.12.0.
  doi:10.32614/CRAN.package.bslib
- Schloerke B, Chang W, Stagg G, Aden-Buie G (2026). shinylive: Run 'shiny'
  applications in the browser. R package version 0.5.0.
  doi:10.32614/CRAN.package.shinylive
- Stagg GW, Lionel H, et al. (2023). webR: The statistical language R
  compiled to WebAssembly via Emscripten. https://github.com/r-wasm/webr
