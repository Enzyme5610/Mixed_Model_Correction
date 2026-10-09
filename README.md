# Mixed Model Correction

Nested model by default; the **Advanced designs** switch (sidebar) adds crossed
and other designs ([Advanced designs](#advanced-designs)).

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

- **One file = one experiment** (same comparison, same measurement type).
- **Line:** one label per actual cell line, e.g. its ID (C20300M). Same label
  across groups for the same line.
- **Batch:** your main batch source, a unit that can shift all its cells
  together: the culture batch if you track it, otherwise the coverslip or
  plate. Number each line's batches B1, B2; different lines can all be B1.
  A second source (coverslips within culture batches, or a shared run) needs
  Advanced designs.
- **qPCR:** enter ΔCt values. Average technical replicates or enter each as
  its own row.

## Model

`parameter ~ Tx + (1 | Line/Batch)` (nested; most ephys, some qPCR), fit with
`lmerTest::lmer`, with Tx tested by a Type II F test using Kenward-Roger
degrees of freedom. With one Line the model uses `(1 | Batch)`; with one
Batch, or no line with repeat batches, `(1 | Line)`. Pairwise comparisons use
`emmeans`.

**FDR q-values** (Benjamini-Hochberg) for panels of many similar parameters
(e.g. gene panels). Adjusts Tx p-values across all parameters in a run.

**Possible outliers:** cells more than 3 SD from the model's prediction for
their group, line and batch (scaled residual). Circled in red and listed
under ANOVA results. Nothing is removed. Optional residual checks (normal
Q-Q, residuals vs fitted) test the model's assumptions on the residuals.

**Skewed data** (a few values much higher than the rest, e.g. event
frequency): log10-transform values before upload. No 0 or negative values.

**Batches shared by lines** (one qPCR plate with several lines) are corrected
per line; the shared plate shift isn't separated. The crossed model is under
[Advanced designs](#advanced-designs).

The results table lists each parameter's random effects (fallbacks can differ
between parameters).

## Plots

Dots, bars, box plots or violins showing each sample, with error bars (95% CI or SE, or 
the raw SEM or SD; none by default for box plots). The Y axis can show values as entered, relative to a 
reference group (linear data), or as fold change 2^-ΔΔCt (ΔCt data, qPCR). 
Statistics always use the values as entered.

**Batch-adjusted values** (optional) subtract each batch's estimated shift
from the plotted values, using the model's random-effect estimates (BLUPs).
With the nested model: each line's batch-to-batch deviation. With the crossed
model (Advanced): each run's shared shift and each line's shift within it.
For display only. Statistics are unchanged, adjusted values shouldn't be re-tested.

A reference (control) group sets the comparison direction in tables and plots.
Significance can be shown as p-values or stars (ns, *, **, ***, ****).

Groups appear in CSV order (reference first) and can be reordered by
dragging; each group's color can be picked. Dots can be colored by Line,
Batch or Group and shaped by Line or Batch, with the color and shape of
each level selectable.

Plots download as PNG, TIFF, JPEG, PDF, SVG or EMF. EMF is editable in PowerPoint.

**All parameters in one figure** shows the selected parameters side by side on
one shared Y axis, with the groups next to each other for each parameter.

**Split and combined** (Plots and Pairwise tabs): groups named like `Ctrl_M,
Ctrl_F, KO_M, KO_F` are split at the last `_`. Ctrl vs KO is fit within each
subset (M, F) and on all cells (Combined), each its own model. Every group
needs every subset. Does not test whether subsets differ (Advanced designs
adds that test).

## Advanced designs

The **Advanced designs** switch (sidebar) adds the following.

**Crossed batch design** (some ephys, most qPCR): runs that held several lines
and were repeated on other days (qPCR plate, or several lines recorded the
same day). Number the shared runs B1, B2 (or use the plate ID or date); the
same label must mean the same run for every line. Lines don't need to be in
every run; also when each line was on one plate only. Model:
`parameter ~ Tx + (1 | Line) + (1 | Batch) + (1 | Line:Batch)`, the nested
model plus a batch shift shared by lines. `(1 | Line:Batch)` is left out with
one row per line per batch, or one batch per line. Without a shared batch it
falls back to nested.

**Optional columns**, used by exact name:

- **Pair:** matched lines (e.g. parental and corrected clone). Adds `(1 | Pair)`.
- **Coverslip:** coverslips within each culture batch (Batch = culture
  batch). Adds `(1 | Line:Batch:Coverslip)`.
- **Run:** a second batch source shared by lines (plate, recording day) when
  Batch is the culture batch. Adds `(1 | Run) + (1 | Line:Run)`.

Each term is added only when the data can estimate it.

**Group × Subset interaction** (Pairwise tab, Split and combined view): does
the group difference change between subsets? With names like `Ctrl_veh,
Ctrl_drug, KO_veh, KO_drug` this is genotype × treatment. Needs 2+ lines per
subset.

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

- Aarts E, Verhage M, Veenvliet JV, Dolan CV, van der Sluis S (2014). A
  solution to dependency: using multilevel analysis to accommodate nested
  data. *Nature Neuroscience* 17(4):491–496.
- Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects
  models using lme4. *Journal of Statistical Software* 67(1):1–48.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. *Journal of the Royal
  Statistical Society: Series B* 57(1):289–300.
- Bland JM, Altman DG (1996). Statistics notes: Transforming data. *BMJ*
  312(7033):770.
- Bolker B, et al. GLMM FAQ: Nested or crossed?
  https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#nested-or-crossed
- Chang W, Cheng J, Allaire JJ, Sievert C, Schloerke B, Aden-Buie G, Xie Y,
  Allen J, McPherson J, Dipert A, Borges B (2026). shiny: Web application
  framework for R. R package version 1.14.0. doi:10.32614/CRAN.package.shiny
- Dunn OJ (1961). Multiple comparisons among means. *Journal of the American
  Statistical Association* 56(293):52–64.
- Halekoh U, Højsgaard S (2014). A Kenward-Roger approximation and parametric
  bootstrap methods for tests in linear mixed models – the R package
  pbkrtest. *Journal of Statistical Software* 59(9):1–32.
- Johnson P (2026). devEMF: EMF graphics output device. R package version 4.6.
  doi:10.32614/CRAN.package.devEMF
- Kenward MG, Roger JH (1997). Small sample inference for fixed effects from
  restricted maximum likelihood. *Biometrics* 53(3):983–997.
- Kramer CY (1956). Extension of multiple range tests to group means with
  unequal numbers of replications. *Biometrics* 12(3):307–310.
- Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests
  in linear mixed effects models. *Journal of Statistical Software*
  82(13):1–26.
- Lenth R, Piaskowski J (2026). emmeans: Estimated marginal means, aka
  least-squares means. R package version 2.0.4.
  doi:10.32614/CRAN.package.emmeans
- Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data
  using real-time quantitative PCR and the 2^-ΔΔCT method. *Methods*
  25(4):402–408.
- R Core Team (2025). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  https://www.R-project.org/
- Robinson GK (1991). That BLUP is a good thing: the estimation of random
  effects. *Statistical Science* 6(1):15–32.
- Schielzeth H, Nakagawa S (2013). Nested by design: model fitting and
  interpretation in a mixed model era. *Methods in Ecology and Evolution*
  4(1):14–24.
- Schloerke B, Chang W, Stagg G, Aden-Buie G (2026). shinylive: Run 'shiny'
  applications in the browser. R package version 0.5.0.
  doi:10.32614/CRAN.package.shinylive
- Sievert C, Cheng J, Aden-Buie G (2026). bslib: Custom 'Bootstrap' 'Sass'
  themes for 'shiny' and 'rmarkdown'. R package version 0.12.0.
  doi:10.32614/CRAN.package.bslib
- Stagg GW, Lionel H, et al. (2023). webR: The statistical language R
  compiled to WebAssembly via Emscripten. https://github.com/r-wasm/webr
- Tukey JW (1953). The problem of multiple comparisons. Unpublished
  manuscript, reprinted in *The Collected Works of John W. Tukey*, Vol. VIII
  (1994). Chapman & Hall.
