# Multicenter study analysis

Analysis code used for the manuscript, "Performance of Gram Stain Characteristics for Predicting Culture Growth in Lower Respiratory Specimens Across Multiple Institutions". The primary workflow is an R [`targets`](https://books.ropensci.org/targets/) pipeline defined in `_targets.R`; supporting plotting and table functions are in `R/`.

## Data availability

Participant-level source data and derived analysis data are not included in this repository, but can be shared upon reasonable request to the corresponding author, Carey-Ann Burnham <carey-ann@pattern.bio>. 

## Reproducing the analysis

1. Open the project in R (version recorded in `renv.lock`) and restore packages with `renv::restore()`.
2. Set `MULTICENTER_STUDY_DATA_DIR` to the directory containing the derived data files.
3. Optionally set `MULTICENTER_STUDY_FIGURE_DIR` to the desired output directory. If it is not set, output is written to `figures/multicenter_study` relative to the project root.
4. Run `targets::tar_make()`.

For example, from an R session:

```r
Sys.setenv(
  MULTICENTER_STUDY_DATA_DIR = "/path/to/approved/data",
  MULTICENTER_STUDY_FIGURE_DIR = "figures/multicenter_study"
)
targets::tar_make()
```

The exploratory notebooks and abstract scripts are retained for provenance. They may require additional source files and local path configuration.
