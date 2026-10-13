# Mixed Model Correction

Browser-based tool for testing treatment effects on physiology or qPCR data
with a linear mixed model, accounting for Line and Batch variability.
Includes pairwise comparisons and downloadable plots. No installation needed.
Data is local to your platform (PC/Mac/Linux).

The regular version uses the nested model. The **Advanced designs** switch
(sidebar) adds crossed batches, matched lines, line-specific treatment effects
and each-vs-reference comparisons
([Advanced designs](#advanced-designs)).

**App:** https://enzyme5610.github.io/Mixed_Model_Correction/

## How to use

1. Upload a CSV (format below).
2. Check the parameter list. Numeric columns are preselected; unselect any
   that aren't parameters (e.g. age, culture days, well number).
3. Choose the pairwise adjustment (Bonferroni or Tukey).
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
- **qPCR:** enter ΔCt values. Average technical replicates or enter each as
  its own row. The ANOVA results and Pairwise tabs have a one-click switch to
  fold change (2^-ΔΔCt).

**Batch column:** groups cells that share more than their line, such as cells
grown or treated together. What to enter:

| Your experiment | Batch column |
|---|---|
| Each line in one group (e.g. Control lines vs AD lines) | Culture batch (differentiation round) or coverslip |
| Same line in several groups (e.g. vehicle and drug) | The unit that got one treatment, usually the coverslip or well |
| Several lines shared plates or recording days | The shared plate or day ([crossed design](#advanced-designs)) |
| No batch recorded | `B1` on every row (model: `parameter ~ Tx + (1 \| Line)`) |

Number each line's batches B1, B2. Another line's B1 is a different batch,
except in the crossed design, where B1 is one plate or day shared by lines.

## Model

**Nested** (the regular version): `parameter ~ Tx + (1 | Line/Batch)`, fit with
`lmerTest::lmer` (REML), with Tx tested by a Type II F test using the
Kenward-Roger method (small-sample corrected standard errors and degrees of
freedom). With one Line the model uses `(1 | Batch)`; with one Batch, or no
line with repeat batches, `(1 | Line)`; with one of each, the parameter is
skipped. Lines and batches are counted per parameter, where it was measured. Pairwise comparisons use `emmeans`;
their confidence intervals use the same adjustment. Exact calls:
[Methods](#methods).

**Model warnings** (e.g. no convergence, or a singular fit: a variance
estimated as zero, so that term has no effect) are listed under ANOVA results.

**FDR q-values** (Benjamini-Hochberg) for panels of many similar parameters
(e.g. gene panels). Adjusts the Tx p-values (not pairwise) across all
parameters in a run.

**Possible outliers:** cells more than 3 SD from the model's prediction for
their group, line and batch (scaled residual). Circled in red and listed
under ANOVA results. Nothing is removed. Optional residual checks (normal
Q-Q, residuals vs fitted) test the model's assumptions on the residuals.

**Skewed data** (a few values much higher than the rest, e.g. event
frequency): log10-transform values before upload. No 0 or negative values.

**Several lines on one plate or day:** the nested model corrects each line's
part of the plate as its own batch (crossed model:
[Advanced designs](#advanced-designs)).

The results table lists each parameter's random effects (fallbacks can differ
between parameters).

## Plots

Dots, bars, box plots or violins showing each sample, with error bars (95% CI or SE, or 
the raw SEM or SD; none by default for box plots). The Y axis can show values as entered, relative to a 
reference group (linear data), or as fold change 2^-ΔΔCt (ΔCt data, qPCR). 
Statistics always use the values as entered.

**Batch-adjusted values** (optional) subtract each batch's estimated shift
from the plotted values, using the model's random-effect estimates (BLUPs):
each line's batch-to-batch deviation. For display only. Statistics are unchanged, adjusted values shouldn't be re-tested.
Error bars then use the model (95% CI or SE); SEM and SD would shrink with the
adjustment, so they return when it is turned off.

Comparisons read as the group plotted later minus the one plotted earlier. The
reference (control) group plots first, so its comparisons read other −
reference; reordering the groups flips a comparison. P-values don't change.
Significance can be shown as p-values or stars (ns, *, **, ***, ****).

Groups appear in CSV order (reference first) and can be reordered by
dragging; each group's color can be picked and its name changed (plots only;
tables keep the original names). Dots can be colored by Line,
Batch or Group and shaped by Line or Batch, with the color and shape of
each level selectable.

Plots download as PNG, TIFF, JPEG, PDF, SVG or EMF. EMF is editable in PowerPoint.

**All parameters in one figure** shows the selected parameters side by side on
one shared Y axis, with the groups next to each other for each parameter.

**Split and combined** (Plots and Pairwise tabs): groups named like `Ctrl_M,
Ctrl_F, KO_M, KO_F` are split at the last `_`. Ctrl vs KO is fit within each
subset (M, F) and on all cells (Combined), each its own model. Every group
needs every subset. The pooled groups (Ctrl, KO) get their own colors and
names.

Significant in M but not in F doesn't mean M and F differ. The **Group ×
Subset interaction** (Pairwise tab, Split and combined view) tests whether the
group difference changes between subsets. With names like `Ctrl_veh,
Ctrl_drug, KO_veh, KO_drug` this is genotype × treatment. Needs 2+ lines per
subset. Model: `parameter ~ Group * Subset + <random effects>`, plus
`(1 | Line:Group)` or `(1 | Line:Subset)` when that factor varies within
lines (each line's own response), with a Kenward-Roger F test of Group × Subset.
With 2 groups and 2 subsets the table also gives the change: the group
difference in the later subset minus the earlier one (e.g. `(KO - Ctrl in drug)
- (KO - Ctrl in veh)`), with a 95% CI. Groups read later minus earlier in plot
order, as in the other tables.

In both figures the parameters or subsets on the x axis can be reordered and
renamed the same way (e.g. Combined to All).

## Advanced designs

Turn on **Advanced designs** in the sidebar. Everything above still applies;
Advanced adds the following. The sidebar's example file then switches to one
with every Advanced feature: healthy lines and mutant clones made from them
(matched lines), each with vehicle and drug, on plates shared by all lines.

**Crossed batch design** (some ephys, most qPCR): runs that held several lines
and were repeated on other days (qPCR plate, or several lines recorded the
same day). Number the shared runs B1, B2 (or use the plate ID or date); the
same label must mean the same run for every line. A line can be in one run or
several. Model:
`parameter ~ Tx + (1 | Line) + (1 | Batch) + (1 | Line:Batch)`, the nested
model plus a batch shift shared by lines. `(1 | Line:Batch)` is left out with
one row per line per batch, or one batch per line. Without a shared batch it
falls back to nested. Batch-adjusted values remove each run's shared shift and
each line's shift within it.

**Matched lines** (sidebar checkbox, shown when some lines aren't in the
reference group): lines come in matched pairs (e.g. a healthy line and a
mutant clone made from it). Each line is matched to its partner in the reference group. Comparing
within each pair removes donor differences, so real effects are easier to
detect. Adds
`(1 | Pair)`.

**Allow line-specific treatment effects** (sidebar checkbox, shown when the
same line is in several groups, e.g. vehicle and drug): adds `(1 | Line:Tx)`.
Each line can respond differently; the test then asks whether the effect holds
across lines. Needs several lines.

**Both together:** matched lines, each in several groups (healthy line and
mutant clone, each with vehicle and drug): `(1 | Pair) + (1 | Line:Tx)`.
Background is removed and each line can respond differently. Matched lines and
line-specific effects work with nested and crossed.

**Comparisons:** all pairs (Bonferroni or Tukey) or each group vs the reference
(Dunnett, emmeans' `dunnettx` approximation). Same model, no refit; tables and
plot brackets show only the tested pairs. With each vs reference, the
reference group decides which comparisons are tested.

## Methods

What runs for each parameter, with links to each package's documentation:

| Step | R call | Docs |
|---|---|---|
| Fit | `lmerTest::lmer(parameter ~ Tx + <random effects>, REML = TRUE)` | [lme4](https://CRAN.R-project.org/package=lme4), [lmerTest](https://CRAN.R-project.org/package=lmerTest) |
| Tx test | `anova(fit, type = "II", ddf = "Kenward-Roger")` | [lmerTest](https://CRAN.R-project.org/package=lmerTest), [pbkrtest](https://CRAN.R-project.org/package=pbkrtest) |
| Model means, SE, 95% CI | `emmeans(fit, "Tx", lmer.df = "kenward-roger")` | [emmeans: models](https://rvlenth.github.io/emmeans/articles/sophisticated.html) |
| All pairs | `pairs(emm, adjust = "bonferroni")` or `"tukey"` | [emmeans: comparisons](https://rvlenth.github.io/emmeans/articles/comparisons.html) |
| Each vs reference (Advanced) | `contrast(emm, "trt.vs.ctrl", adjust = "dunnettx")` | [emmeans: comparisons](https://rvlenth.github.io/emmeans/articles/comparisons.html) |
| Group × Subset | `anova(lmer(parameter ~ Group * Subset + <random effects> + (1 \| Line:Group) or (1 \| Line:Subset)), type = "II", ddf = "Kenward-Roger")` | [lmerTest](https://CRAN.R-project.org/package=lmerTest) |
| Group × Subset change | `contrast(emmeans(fit, ~ Group * Subset), interaction = "revpairwise")` | [emmeans: interactions](https://rvlenth.github.io/emmeans/articles/interactions.html) |
| FDR q-values | `p.adjust(p, "BH")` across parameters | [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html) |
| Possible outliers | `residuals(fit, type = "pearson", scaled = TRUE)` beyond ±3 | [lme4](https://CRAN.R-project.org/package=lme4) |
| Batch-adjusted values | `lme4::ranef(fit)` (BLUPs) | [lme4](https://CRAN.R-project.org/package=lme4) |

Kenward-Roger corrects both the standard errors and the degrees of freedom for
small samples. Pairwise confidence intervals use the same adjustment as their
p-values. `dunnettx` is emmeans' approximation to Dunnett's test.

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
  doi:10.18637/jss.v067.i01
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
  doi:10.18637/jss.v059.i09
- Johnson P (2026). devEMF: EMF graphics output device. R package version 4.6.
  doi:10.32614/CRAN.package.devEMF
- Kenward MG, Roger JH (1997). Small sample inference for fixed effects from
  restricted maximum likelihood. *Biometrics* 53(3):983–997.
- Kramer CY (1956). Extension of multiple range tests to group means with
  unequal numbers of replications. *Biometrics* 12(3):307–310.
- Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests
  in linear mixed effects models. *Journal of Statistical Software*
  82(13):1–26. doi:10.18637/jss.v082.i13
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
