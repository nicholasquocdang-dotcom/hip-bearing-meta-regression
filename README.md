# hip-bearing-meta-regression

Data and R code for a registry-based random-effects meta-regression comparing revision risk across total hip replacement bearing surfaces, with metal-on-highly-cross-linked polyethylene as the reference.

## Files

| File | Contents |
|---|---|
| hip_meta_data.csv | 38 hazard ratios extracted from ten registry sources (NJR, Nordic registries, Dutch LROI, German EPRD, US Kaiser Permanente, US AJRR, US HealthEast) |
| hip_meta_analysis.R | Main model, sensitivity analyses, and forest plot |
| sensitivity_results.csv | Full results of all sensitivity analyses |
| PubMed_screening_log.xlsx | Search string, screening decision for every PubMed record, and full-text decisions |

## Data columns

| Column | Meaning |
|---|---|
| registry, source | Registry dataset and publication |
| label, contrast | Implant group and the bearing comparison (compared vs. reference) |
| hr, lo, hi | Hazard ratio and 95% confidence interval as reported |
| t | Follow-up time in years at which the estimate applies |
| t_source | How the follow-up time was obtained (reported, derived, or approximated) |
| hr2, lo2, hi2 | NJR 2-year estimates used in a sensitivity analysis |
| MoCPE, CoCPE, CoXLPE, CoC, MoM | Contrast coding: +1 for the compared bearing, -1 for the reference bearing |
| adjusted | Whether the source reported an adjusted hazard ratio |

## How to run

Install R, then run install.packages(c("metafor", "clubSandwich")) in R. Put hip_meta_data.csv and hip_meta_analysis.R in the same folder, set that folder as the working directory, and run source("hip_meta_analysis.R"). The script prints the main model and sensitivity analyses and saves forest.png.

## Sources

Whitehouse et al., PLOS Medicine, 2024 (NJR); Mikkelsen et al., Acta Orthopaedica, 2023 (NARA); Pakarinen et al., Journal of Arthroplasty, 2025 (NARA); Kadar et al., Clinical Orthopaedics and Related Research, 2012 (Norwegian Arthroplasty Register); Peters et al., Acta Orthopaedica, 2018 (Dutch Arthroplasty Register); Elliott et al., Acta Orthopaedica, 2026 (German Arthroplasty Registry); Paxton et al., Clinical Orthopaedics and Related Research, 2015 (Kaiser Permanente); Cafri et al., Clinical Orthopaedics and Related Research, 2017 (Kaiser Permanente); Reddy et al., Journal of Arthroplasty, 2026 (AJRR); Huang et al., Clinical Orthopaedics and Related Research, 2013 (HealthEast).

## License

Code: MIT (see LICENSE). Data (hip_meta_data.csv): CC BY 4.0. The hazard ratios were extracted from the published sources above, which should also be cited.
