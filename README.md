# Methodology for assessing damage to calcareous rocky habitats caused by military activities

R scripts accompanying the manuscript of the same title by Iuliia Vasheniak and Yakiv Didukh. The repository provides a visual comparison of seven habitat types and an analysis of the sensitivity of four habitat assessment indices to changes in expert scores.

## Files

| File | Purpose |
| --- | --- |
| [Four-dimensional bubble plot of habitat assessment parameters.R](Four-dimensional%20bubble%20plot%20of%20habitat%20assessment%20parameters.R) | Produces the four-dimensional bubble plot corresponding to Figure 3. |
| [Reliability assessment of the proposed methodology.R](Reliability%20assessment%20of%20the%20proposed%20methodology.R) | Calculates baseline indices, leave-one-criterion-out (LOCO) sensitivity, rank stability, and Monte Carlo sensitivity; produces the panels corresponding to Figure 4. |

Both scripts contain their input data and run independently. No external spreadsheet or data file is required. The bubble plot uses rounded index values; the sensitivity script contains all 161 primary scores (seven habitats × 23 criteria) and calculates indices without intermediate rounding. If primary scores are revised, update the sensitivity script's `expected` baseline checks and the rounded values in the bubble-plot script accordingly.

## Requirements and running the scripts

Install R; RStudio is optional. Required packages are `ggplot2`, `ggrepel`, `viridisLite`, and `patchwork`.

1. Download the repository using **Code → Download ZIP**, then extract the archive.
2. Set the R working directory to the extracted folder containing the scripts. Use `getwd()` to check the current directory.
3. Install the packages once in the R console:

   ```r
   install.packages(
     c("ggplot2", "ggrepel", "viridisLite", "patchwork"),
     repos = "https://cloud.r-project.org"
   )
   ```

4. In RStudio, open either script and click **Source**, or run these commands in the R console:

   ```r
   source("Four-dimensional bubble plot of habitat assessment parameters.R", encoding = "UTF-8")
   source("Reliability assessment of the proposed methodology.R", encoding = "UTF-8")
   ```

Alternatively, after installing the packages, run these commands in a terminal from the repository folder:

```sh
Rscript "Four-dimensional bubble plot of habitat assessment parameters.R"
Rscript "Reliability assessment of the proposed methodology.R"
```

Internet access and a writable R package library are required for package installation. The scripts also attempt to install missing packages automatically; the bubble-plot script may prompt for a CRAN mirror if one has not been configured. Using the installation command above avoids that prompt. R and package versions are not pinned; the sensitivity script records the versions used in `R_sessionInfo.txt`.

Paths for generated files are relative to the current R working directory, returned by `getwd()`. Re-running a script overwrites its generated files with the same names.

## Indices and primary data

Each index is the equally weighted arithmetic mean of its criteria on the scale 0.1, 0.3, 0.6, and 1.0.

| Index | Criteria | Interpretation of higher values |
| --- | --- | --- |
| Z | z1–z10 | Better habitat conservation status. |
| T | t1–t6 | Greater conservation value. |
| V | v1–v3 | Greater expected recovery difficulty and longer recovery time. |
| ST | st1–st4 | More favourable combined spatiotemporal characteristics. |

The habitat types are R15.1, R15.2, R15.3, R15.4, R15.5, R16b.2, and R13.2. The primary scores correspond to the source workbook `Для валідації методики.xlsx`, sheet `Аркуш1`, used for the analysis dated 30 October 2026. The workbook is identified for provenance and is not required to run the scripts. The manuscript provides the ecological definitions of the individual criteria. The scripts analyse the four indices separately and do not calculate a combined priority index.

## Figure 3

The bubble plot uses the two-decimal index values reported in Table 5. Z is plotted on a reversed horizontal axis and T on the vertical axis. Colour represents V; bubble size represents ST. The dashed lines at 0.3 and 0.6 are visual guides and do not define independently validated priority classes. Bubble sizes are rescaled to the observed ST range; apparent area differences should be read with the legend and do not directly express proportional differences in ecological impact.

The bubble-plot script does not set a label-placement seed, so label positions may vary between runs. Font appearance can also vary across graphics environments.

Output in the current working directory: `Priority_matrix.tiff` (**8 × 8 inches, 600 dpi, LZW compression**).

## LOCO and rank stability

For each index, one criterion is removed at a time and the mean of the remaining criteria is recalculated. Relative change is `100 * abs(modified - baseline) / baseline`. Mean and maximum changes are summarized across the seven habitat types for each removed criterion.

Spearman's rho and Kendall's tau-b compare baseline and perturbed habitat orderings within each index. Tied values receive average ranks. The exported LOCO ranks assign lower rank numbers to lower Z and ST values and to higher T and V values. Primary scores are converted to integers before summation to avoid artificial separation of ties caused by floating-point arithmetic. Calculations use unrounded index values.

## Monte Carlo scenarios

Each primary score is perturbed independently. The original score is retained with probability `p_base`; otherwise it moves by one adjacent level of the four-level scale.

| Original score | Allowed outcomes | Probabilities |
| --- | --- | --- |
| 0.1 | 0.1, 0.3 | `p_base`, `1 - p_base` |
| 0.3 | 0.1, 0.3, 0.6 | `(1 - p_base)/2`, `p_base`, `(1 - p_base)/2` |
| 0.6 | 0.3, 0.6, 1.0 | `(1 - p_base)/2`, `p_base`, `(1 - p_base)/2` |
| 1.0 | 0.6, 1.0 | `1 - p_base`, `p_base` |

There are 10,000 simulations per scenario. The primary scenario uses `p_base = 0.60`; additional scenarios use 0.80 and 0.90. The random seed is `20261030`. Common uniform random numbers are used across scenarios. These probabilities specify sensitivity scenarios and are not empirically calibrated expert error rates.

Unequal distances between score levels and one-sided boundary transitions mean that the perturbations do not preserve the baseline expectation. For example, when `p_base = 0.60`, a score of 1.0 has an expected simulated value of 0.84.

The script reports means, medians, standard deviations, the 2.5th and 97.5th percentiles of simulated index values, theoretical means under the transition rules, and shifts from baseline. Quantiles use R's `type = 7` convention. The percentile limits define conditional simulation intervals, not empirical confidence intervals. Baseline inclusion in an interval is not evidence of accuracy or external validity. The script also summarizes Spearman and Kendall correlations across simulation runs. Correlations are undefined when all seven simulated values are identical; such cases are counted and omitted only when summarizing the rank correlations.

## Sensitivity outputs

The sensitivity script creates `Habitat_sensitivity_20261030/` under the current working directory:

- `00_primary_criteria.csv`: primary scores.
- `01_baseline_indices_exact.csv`: unrounded baseline indices.
- `02_LOCO_full.csv`: each habitat and criterion removal.
- `03_LOCO_summary.csv`: changes and rank correlations by removed criterion.
- `04_MC_intervals_all_scenarios.csv`: index distributions under all scenarios.
- `05_MC_intervals_primary.csv`: index distributions for `p_base = 0.60`.
- `06_MC_rank_stability.csv`: rank-correlation summaries and undefined-case counts.
- `MC_primary_samples.rds`: simulated index values for the primary scenario.
- `analysis_settings.txt` and `R_sessionInfo.txt`: analysis settings and software versions.
- `Sensitivity_combined.tiff`: combined figure, 10 × 11 inches, 600 dpi, LZW compression.
- `Sensitivity_combined.png`: combined figure, 10 × 11 inches, 200 dpi.
- `A_LOCO.tiff`: bars show mean relative changes; circles show maximum changes across habitats.
- `B_Rank_stability.tiff`: Spearman and Kendall correlations following criterion removal.
- `C_Monte_Carlo.tiff`: simulated means, baseline values, and percentile intervals for the primary scenario.

Individual TIFF panels are **10 × 3.8 inches, 600 dpi, with LZW compression**. In panel C, filled points show simulation **means**, crosses show baseline values, and lines show the 2.5th–97.5th percentile intervals. Medians are included in the CSV tables. The RDS file can be opened with `readRDS()`.

## Interpretation

The analysis assesses sensitivity conditional on the supplied expert scores, equal criterion weights, and specified adjacent-level perturbations. It does not establish empirical accuracy of the expert scores or independently validate the chosen numerical scale and weights. ST contains five identical primary profiles, which limits detailed ordering of these habitats.

## Software documentation

- [R rank correlations](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)
- [ggrepel label placement](https://ggrepel.slowkow.com/reference/geom_text_repel.html)
- [ggplot2 size scales](https://ggplot2.tidyverse.org/reference/scale_size.html)
