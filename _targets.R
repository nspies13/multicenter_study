library(targets)
library(tarchetypes) 
library(tidyverse)
library(gt)
library(gtsummary)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("spec", "yardstick")

tar_option_set(packages = c("tidyverse", "gt", "gtsummary", "ggpubr", "glue"), error = "continue")
tar_source()

load_data <- list(
  
  ### Load the Data
  tar_target(gram_stain_raw, read_csv(file.path(data_path, "gram_stain_general_input.csv")) |> filter(LRT_REF_SQ_EPI_AMT != "Abundant") |> mutate(LRT_REF_PLS_AMT = factor(LRT_REF_PLS_AMT, levels = c("None", "Rare", "Few", "Moderate", "Abundant"), ordered = T), across(matches("RBC|SQ"), ~factor(., levels = c("None", "Few", "Abundant"), ordered = T)))),
  tar_target(gram_stain_results_input, read_rds(file.path(data_path, "gram_stain_results_input.rds")) |> filter(LRT_REF_SQ_EPI_AMT != "Abundant")),
  tar_target(gram_culture_joined, read_rds(file.path(data_path, "culture_gram_results_joined.rds")) |> filter(`Epithelial Cells` != "Abundant")),
  tar_target(gram_culture_joined_by_subject, read_rds(file.path(data_path, "culture_gram_join_by_subject.rds")) |> filter(`Epithelial Cells` != "Abundant")),
  tar_target(culture_results_joined_by_organism, makeCultureInput(gram_culture_joined)), 
  tar_target(culture_results_joined_by_subject, makeCultureInputBySubject(culture_results_joined_by_organism)),
  
  ### Nest Data for Mapping
  tar_target(gram_culture_joined_nested, gram_culture_joined |> group_by(SiteID, LRT_REF_GS_SPECTYPE) |> nest()),
  tar_target(gram_culture_joined_combined_list, bind_rows(gram_culture_joined_nested, tibble(SiteID = c("All", "All", "All"), LRT_REF_GS_SPECTYPE = c("All", "BAL", "ETA"), data = list(gram_culture_joined, gram_culture_joined |> filter(LRT_REF_GS_SPECTYPE == "BAL"), gram_culture_joined |> filter(LRT_REF_GS_SPECTYPE == "ETA"))))),
  tar_target(gram_culture_joined_by_subject_nested, gram_culture_joined_by_subject |> group_by(Site, Procedure) |> nest()),
  tar_target(gram_culture_joined_by_subject_combined_list, bind_rows(gram_culture_joined_nested, tibble(Site = c("All", "All", "All"), Procedure = c("All", "BAL", "ETA"), data = list(gram_culture_joined_by_subject, gram_culture_joined_by_subject |> filter(Procedure == "BAL"), gram_culture_joined_by_subject |> filter(Procedure == "ETA"))))),
  
  ### Make Modeling Inputs
  tar_target(significant_growth_model_input, gram_culture_joined_by_subject |> mutate(target = `Any Significant Growth`) |> select(-`Any Significant Growth`, -`Any Concordant Growth`)),
  tar_target(concordant_growth_model_input, gram_culture_joined_by_subject |> mutate(target = `Any Concordant Growth`) |> select(-`Any Significant Growth`, -`Any Concordant Growth`))
  
)

paper_figures <- list(
  
  # Make Demographics Table
  tar_target(table_1_demographic_overview, makeDemographicsTable(gram_culture_joined_by_subject)),
  
  # Make Descriptive Statistic Figures
  tar_target(table_2a_gram_stain_metadata, makeGramStainMetadataTables(gram_culture_joined_by_subject |> ungroup() |> select(-c(SubjectID, Age, Sex, Location, `Unique Organisms Grown`, `Any Significant Growth`, `Predominant Morphology`, `Predominant Morphology (w/ Other)`, `Any Concordant Growth`)))),
  tar_target(table_2a_gram_stain_results, makeGramStainMorphologyTables(gram_culture_joined |> ungroup() |> select(SubjectID, Site = SiteID, Morphology = LRT_REF_GS_MORPH, Procedure = LRT_REF_GS_SPECTYPE, `Organism Abdunance` = LRT_REF_GS_ORG_AMT) |> select(-SubjectID))),
  tar_target(supp_table_culture_results_table, makeCultureResultsTables(input_by_organism = culture_results_joined_by_organism, input_by_subject = culture_results_joined_by_subject)),
  tar_target(fig2b_culture_significant_growth_counts, makeSignificantGrowthCountFigures(gram_culture_joined)),
  tar_target(fig2c_culture_concordant_growth_counts, makeConcordantGrowthCountFigures(gram_culture_joined)),
  
  # Make Predictors of Significant Growth Figures
  tar_target(fig3_predicting_significant_growth_figures, 
             makeModelFigures(input = gram_culture_joined_by_subject, 
                              target = "Any Significant Growth", procedure = "All", 
                              out_path = paste0(figure_path, "final_figures/any_significant_growth_all"))),
  tar_target(supp_fig_predicting_significant_growth_figures_eta, 
             makeModelFigures(input = gram_culture_joined_by_subject, 
                              target = "Any Significant Growth", procedure = "ETA", 
                              out_path = paste0(figure_path, "final_figures/any_significant_growth_ETA"))),
  tar_target(supp_fig_predicting_significant_growth_figures_bal, 
             makeModelFigures(input = gram_culture_joined_by_subject, 
                              target = "Any Significant Growth", procedure = "BAL", 
                              out_path = paste0(figure_path, "final_figures/any_significant_growth_BAL"))),
  
  # Make Predictors of Concordant Growth Figures
  tar_target(fig4_predicting_concordant_growth_figures, 
             makeModelFigures(input = gram_culture_joined_by_subject, 
                              target = "Any Concordant Growth", procedure = "All", 
                              out_path = paste0(figure_path, "final_figures/any_concordant_growth_all"))),
  tar_target(supp_fig_predicting_concordant_growth_figures_eta, 
             makeModelFigures(input = gram_culture_joined_by_subject |> mutate(`Organism Abundance` = ifelse(`Organism Abundance` == "None", NA, `Organism Abundance`)), 
                              target = "Any Concordant Growth", procedure = "ETA", 
                              out_path = paste0(figure_path, "final_figures/any_concordant_growth_ETA"))),
  tar_target(supp_fig_predicting_concordant_growth_figures_bal, 
             makeModelFigures(input = gram_culture_joined_by_subject |> mutate(`Organism Abundance` = ifelse(`Organism Abundance` == "None", NA, `Organism Abundance`)), 
                              target = "Any Concordant Growth", procedure = "BAL", 
                              out_path = paste0(figure_path, "final_figures/any_concordant_growth_BAL")))
  
)
  
other_figures <- list(
  
  ### Make Organism Abundance Figures
  tar_target(gram_morphology_abundance_by_predominant_organism_figure, makeGramStainAbundanceByMorphologyFigure(gram_culture_joined_by_subject)),
  tar_target(morphology_by_significant_growth_figure, makeSignificantGrowthByMorphologyFigures(gram_culture_joined))

)

## Run Pipeline
list(
  load_data,
  paper_figures
  ) 
