# hip-bearing-meta-regression

Data and R code for a registry-based random-effects meta-regression comparing revision risk across total hip replacement bearing surfaces, with metal-on-highly-cross-linked polyethylene as the reference.

## Files

| File | Contents |
|---|---|
| hip_meta_data.csv | 29 hazard ratios extracted from four registry sources (NJR, NARA, Dutch LROI, German EPRD) |
| hip_meta_analysis.R | Main model, sensitivity analyses, and figures |

## Data columns

| Column | Meaning |
|---|---|
| registry, source | Registry dataset and publication |
| label, contrast | Implant group and the bearing comparison (compared vs. reference) |
| hr, lo, hi | Hazard ratio and 95% confidence interval as reported |
| t | Follow-up time in years at which the estimate applies |
| hr2, lo2, hi2 | NJR 2-year estimates used in a sensitivity analysis |
| MoCPE, CoCPE, CoXLPE, CoC, MoM | Contrast coding: +1 for the compared bearing, -1 for the reference bearing |
| adjusted | Whether the source reported an adjusted hazard ratio |

## How to run

Install R, then run install.packages(c("metafor", "clubSandwich")) in R. Put both files in the same folder, set that folder as the working directory, and run source("hip_meta_analysis.R"). The script prints the main model and sensitivity analyses and saves forest.png and bubble.png.

## Sources

Whitehouse et al., PLOS Medicine, 2024 (NJR); Mikkelsen et al., Acta Orthopaedica, 2023 (NARA); Peters et al., Acta Orthopaedica, 2018 (Dutch Arthroplasty Register); Elliott et al., Acta Orthopaedica, 2026 (German Arthroplasty Registry).

## License

Code: MIT (see LICENSE). Data (hip_meta_data.csv): CC BY 4.0. The hazard ratios were extracted from the published sources above, which should also be cited.
