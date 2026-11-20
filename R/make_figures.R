makeCultureInput <- function(input = tar_read(gram_culture_joined)){
  
  output <- 
    input |>
      select(
        SubjectID, Procedure = LRT_REF_GS_SPECTYPE, Site = SiteID, Organism = MALDI_ID_COLLAPSED, Organism_full_name = MALDI_ID_GROUPED, `Colony Count` = LRT_REF_RPT_BIN) |> 
      mutate(
        Organism_full_name = str_squish(str_to_upper(Organism_full_name)),
        `Organism Gram Stain` = 
          case_when(
            grepl("STAPH|STREP|ENTEROCOC|CORYN|BACILLUS", Organism_full_name) ~ "Positive",
            grepl("KLEBS|ESCHERI|ENTEROBAC|SERRAT|PROTEU|CITROBA|HAEMOP|MORAX|ACHROM|BORDE|PSEUDOM|STENO|ACINETO|CHRYSEO", Organism_full_name) ~ "Negative",
            grepl("YEAST|CANDIDA|MOLD", Organism_full_name) ~ "Yeast/Mold",
          ),
        `Growth Significance` =
          case_when(
            Procedure == "BAL" &
              `Colony Count` %in% c("10,000-100,000", ">100,000") ~ "Significant",
            Procedure == "BAL" &
              `Colony Count` %in% c("<1,000", "1,000-10,000") ~ "Insignificant",
            Procedure == "BAL" &
              `Colony Count` %in% c(NA) ~ "No Growth",
            Procedure == "ETA" &
              `Colony Count` %in% c(">100,000") ~ "Significant",
            Procedure == "ETA" &
              `Colony Count` %in% c("<1,000", "1,000-10,000", "10,000-100,000") ~ "Insignificant",
            Procedure == "ETA" &
              `Colony Count` %in% c(NA) ~ "No Growth"),
        `Growth Significance` = factor(`Growth Significance`, levels = c("No Growth", "Insignificant", "Significant"))) |> 
    distinct()
  
  output
  
}

makeCultureInputBySubject <- function(input = tar_read(culture_results_joined_by_organism)){
  
  output_long <- 
    input |> 
    group_by(SubjectID) |>
    mutate(
      `No Growth` = all(is.na(Organism_full_name), na.rm = T),
      `Any Gram-Negative` = any(`Organism Gram Stain` == "Negative", na.rm = T),
      `Any Gram-Positive` = any(`Organism Gram Stain` == "Positive", na.rm = T),
      `Any Yeast/Mold` = any(`Organism Gram Stain` == "Yeast/Mold", na.rm = T),
      `Any S. aureus` = any(Organism_full_name == "STAPHYLOCOCCUS AUREUS", na.rm = T),
      `Any P. aeruginosa` = any(Organism_full_name == "PSEUDOMONAS AERUGINOSA", na.rm = T),
      `Any K. pneumoniae` = any(Organism_full_name == "KLEBSIELLA PNEUMONIAE", na.rm = T),
      `Significant Gram-Negative` = any(`Organism Gram Stain` == "Negative" & `Growth Significance` == "Significant", na.rm = T),
      `Significant Gram-Positive` = any(`Organism Gram Stain` == "Positive" & `Growth Significance` == "Significant", na.rm = T),
      `Significant Yeast` = any(`Organism Gram Stain` == "Yeast/Mold" & `Growth Significance` == "Significant", na.rm = T),
      `Significant S. aureus` = any(Organism_full_name == "STAPHYLOCOCCUS AUREUS" & `Growth Significance` == "Significant", na.rm = T),
      `Significant P. aeruginosa` = any(Organism_full_name == "PSEUDOMONAS AERUGINOSA" & `Growth Significance` == "Significant", na.rm = T),
      `Significant K. pneumoniae` = any(Organism_full_name == "KLEBSIELLA PNEUMONIAE" & `Growth Significance` == "Significant", na.rm = T))
  
  output_counts <-  
    output_long |> select(Organism_full_name, `Growth Significance`) |> distinct() |> count(`Growth Significance`, Organism_full_name) |> 
      group_by(SubjectID, `Growth Significance`) |> 
      summarise(`Unique Organisms` = sum(n)) |> 
      group_by(SubjectID) |>
      mutate(`Unique Organisms` = ifelse(`Growth Significance` == "No Growth", 0, `Unique Organisms`)) |>
      mutate(`Unique Organisms: ALL` = sum(`Unique Organisms`)) |> 
      pivot_wider(names_from = `Growth Significance`, values_from = `Unique Organisms`, names_prefix = "Unique Organisms: ", values_fill = 0) |> 
      select(-`Unique Organisms: No Growth`)
  
  output <- 
    left_join(output_long, output_counts) |> 
      select(-c(`Organism`:`Growth Significance`)) |> 
      distinct() |> 
      ungroup()
  
  output <- output |> left_join(read_csv("../../../data/culture_contam.csv")) |> mutate(Contaminated = ifelse(is.na(Contaminated), F, Contaminated))
  
  output

}

makeContaminatedCultureTables <- function(input = tar_read(culture_results_joined_by_subject)){
  


  
}

makeDemographicsTable <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  gt_demographics <-
    input |> 
      ungroup() |> 
      select(Site, Age, Sex, Location, Procedure) |> 
      tbl_summary(by = Site, missing_text = "None") |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p() |>  
      add_q() |>
      bold_p(q = T) |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
      as_gt() |> 
      gt::tab_header(title = "Cohort Demographics") |>
      theme_pattern()
  gtsave(gt_demographics, filename = paste0(figure_path, "final_figures/demographics_table.html"))
  gtsave(gt_demographics, filename = paste0(figure_path, "final_figures/demographics_table.png"), zoom = 10)
  gtsave(gt_demographics, filename = paste0(figure_path, "final_figures/demographics_table.docx"))
  
  gt_demographics
  
}

makeGramStainMetadataTables <- function(input = tar_read(gram_culture_joined_by_subject) |> ungroup() |> select(-SubjectID, -`Unique Organisms Grown`, -`Predominant Morphology`, -`No Organism Seen`)){
  
  gt_gram_stain_general_by_site <- 
    input |> 
      select(-Procedure, "Predominant Morphology" = `Predominant Morphology (w/ Other)`, -`Epithelial Cells`, -Age, -Sex, -c(`Gram-Positive Cocci`:`Any Concordant Growth`, `Any Significant Growth`)) |> 
      tbl_summary(by = Site, sort = list(`Predominant Morphology` = "frequency"), missing = "no") |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      modify_footnote(c(all_stat_cols(), p.value) ~ NA) |> 
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_site_per_specimen_table.html"))
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_site_per_specimen_table.pdf"))
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_site_per_specimen_table.docx"))
  
  gt_gram_stain_general_by_procedure <- 
    input |>
      select(-Site, "Predominant Morphology" = `Predominant Morphology (w/ Other)`, -`Epithelial Cells`, -Age, -Sex, -c(`Gram-Positive Cocci`:`Any Concordant Growth`, `Any Significant Growth`)) |> 
      tbl_summary(by = Procedure, sort = list(`Predominant Morphology (w/ Other)` = "frequency"), missing = "no") |>
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |> 
      add_overall() |>  
      modify_footnote(c(all_stat_cols(), p.value) ~ NA) |> 
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") 
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_procedure_per_specimen_table.html"))
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_procedure_per_specimen_table.pdf"))
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "final_figures/general_gram_stain_overview_by_procedure_per_specimen_table.docx"))
  
  gt_gram_stain_combo <- 
    tbl_merge(list(gt_gram_stain_general_by_procedure, gt_gram_stain_general_by_site), tab_spanner = c("By Procedure", "By Site")) |> as_gt() |> theme_pattern() |>
    tab_header(title = "Gram Stain Results", subtitle = "By Site and Collection Procedure")
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "final_figures/general_gram_stain_overview_combo_table.html"))
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "final_figures/general_gram_stain_overview_combo_table.png"), zoom = 6)
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "final_figures/general_gram_stain_overview_combo_table.pdf"))
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "final_figures/general_gram_stain_overview_combo_table.docx"))
  
  list("gram_stain_combined" = gt_gram_stain_combo, "gram_stain_by_site" = gt_gram_stain_general_by_site, "gram_stain_by_procedure" = gt_gram_stain_general_by_procedure)
  
}

makeGramStainMorphologyTables <- function(gram_stain_results_input = tar_read(gram_culture_joined) |> ungroup() |> select(SubjectID, Site = SiteID, LRT_REF_GS_MORPH, Procedure = LRT_REF_GS_SPECTYPE, LRT_REF_GS_ORG_AMT) |> select(-SubjectID)){

  gt_organism_count_by_site <- 
    gram_stain_results_input |> 
    tbl_summary(by = Site, missing = "no") |> 
    italicize_levels() |>
    bold_labels() |> 
    add_overall() |> 
    add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |> 
    bold_p()
  gt::gtsave(gt_organism_count_by_site |> as_gt() |> tab_header(title = "Gram Stain Morphologies", subtitle = "By Site, All Morphologies"), 
             filename = paste0(figure_path, "final_figures/gram_stain_organism_count_by_site_table.html"))
  gt::gtsave(gt_organism_count_by_site |> as_gt() |> tab_header(title = "Gram Stain Morphologies", subtitle = "By Site, All Morphologies"), 
             filename = paste0(figure_path, "final_figures/gram_stain_organism_count_by_site_table.docx"))
  
  gt_organism_count_by_procedure <- 
    gram_stain_results_input |> 
    tbl_summary(by = Procedure, missing = "no") |> 
    italicize_levels() |>
    bold_labels() |> 
    add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |> 
    bold_p() |> 
    add_overall() |> 
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_organism_count_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Morphologies", subtitle = "By Collection Procedure"), 
             filename = paste0(figure_path, "final_figures/gram_stain_organism_count_by_procedure_table.html"))
  gt::gtsave(gt_organism_count_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Morphologies", subtitle = "By Collection Procedure"), 
             filename = paste0(figure_path, "final_figures/gram_stain_organism_count_by_procedure_table.docx"))
  
  gt_gram_stain_organism_combo <- 
    tbl_merge(list(gt_organism_count_by_site, gt_organism_count_by_procedure)) |> as_gt() |> 
    tab_header(title = "Gram Stain Morphologies", subtitle = "By Site and Procedure")
  gt::gtsave(gt_gram_stain_organism_combo, filename = paste0(figure_path, "final_figures/gram_stain_organism_combo_table.html"))
  gt::gtsave(gt_gram_stain_organism_combo, filename = paste0(figure_path, "final_figures/gram_stain_organism_combo_table.pdf"))
  gt::gtsave(gt_gram_stain_organism_combo, filename = paste0(figure_path, "final_figures/gram_stain_organism_combo_table.docx"))
  
  list("gram_stain_organisms_combined" = gt_gram_stain_organism_combo, "gram_stain_organisms_by_site" = gt_organism_count_by_site, "gram_stain_organisms_by_procedure" = gt_organism_count_by_procedure)
  
}

makeGramStainCountFigures <- function(input = tar_read(gram_culture_joined_by_subject) |> ungroup() |> select(-`Predominant Morphology`, "Predominant Morphology" = `Predominant Morphology (w/ Other)`)){
 
  features_to_plot <- c("Site", "Procedure", "Sex", "Location", 
                        "Epithelial Cells", "Polymorphonucleocytes", 
                        "Organism Abundance", "Predominant Morphology")
  
  plots <- map(features_to_plot, ~input[.x] |> count(!!sym(.x)) |> drop_na() |>  
                 ggplot(aes(y = fct_rev(factor(!!sym(.x))), x = n)) + 
                  geom_col(fill = "#63666aff", color = "black") +
                  geom_text(aes(x = n + 2, label = n), hjust = 0, color = "black", fontface = "bold") + 
                  geom_text(aes(x = 1, label = fct_rev(factor(!!sym(.x)))), hjust = 0, color = "white", fontface = "bold") + 
                  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
                  labs(title = .x) + 
                  coord_cartesian(clip = "off") + 
                  theme(axis.text.y.left = element_blank(), axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank(), axis.line = element_blank(), axis.text.x.bottom = element_blank()))
  
  gg_combo <- ggpubr::ggarrange(plotlist = plots, ncol = 2, nrow = 4)
  ggsave(paste0(figure_path, "final_figures/gram_stain_results_by_feature.png"), gg_combo, width = 10.5, height = 9.6, units = "in", dpi = 600)
  ggsave(paste0(figure_path, "final_figures/gram_stain_results_by_feature.svg"), gg_combo, width = 10.5, height = 9.6, units = "in", dpi = 600)
  
  gg_input <- 
    input |>
      select(-Age, -`Unique Organisms Grown`, -`Predominant Morphology`, -where(is.logical), "Predominant Morphology" = `Predominant Morphology (w/ Other)`) |> 
      mutate(
        across(where(is.character) | where(is.logical) | where(is.numeric), ~factor(.x)),
        across(where(is.ordered), ~factor(.x, ordered = F))) |> 
      pivot_longer(cols = -SubjectID, names_to = "feature", values_to = "value") |>
      mutate(feature = factor(feature, levels = c("Site", "Procedure", "Sex", "Location", 
                                                  "Epithelial Cells", "Red Blood Cells", "Polymorphonucleocytes", 
                                                  "Organism Abundance", "Predominant Morphology", "Unique Morphologies Seen"))) |> 
      ungroup() |> 
      count(feature, value)

  
    
  gg_input |> 
    ggplot(aes(y = value, x = n)) +
      geom_col(position = "dodge") + 
      facet_wrap(~feature, scales = "free_y", ncol = 2)
  

  }

makeGramStainAbundanceByMorphologyFigure <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  top_5_morphs <- 
    input |> ungroup() |> count(`Predominant Morphology (w/ Other)`, sort = T) |> drop_na() |> head(5) |> pull(`Predominant Morphology (w/ Other)`)
  
  gg_abundance_by_organism_input <- 
    input |> 
    filter(`Predominant Morphology (w/ Other)` %in% top_5_morphs) |> 
    mutate(Abundance_bins = factor(`Organism Abundance`, levels = c("Rare", "Few", "Moderate", "Abundant"), labels = c("Rare/Few", "Rare/Few", "Moderate/Abundant", "Moderate/Abundant")),
           `Predominant Morphology (w/ Other)` = factor(`Predominant Morphology (w/ Other)`, levels = top_5_morphs))
  
  gg_abundance_by_organism <- 
    ggplot(data = gg_abundance_by_organism_input, aes(y = fct_rev(`Predominant Morphology (w/ Other)`), fill = fct_rev(Abundance_bins))) + 
    geom_bar(position = "fill", color = "black") +
    stat_count(data = gg_abundance_by_organism_input |> drop_na(`Predominant Morphology (w/ Other)`), aes(label = paste0("n = ", after_stat(count)), y = fct_rev(`Predominant Morphology (w/ Other)`)), x = 0.01, hjust = 0, vjust = 0.5, geom = "text", inherit.aes = F, color = "black", fontface = "bold", size = 4) + 
    scico::scale_fill_scico_d(palette = "bam", begin = 1, end = 0.5) +
    guides(fill = guide_legend(title = "Abundance", label.position = "top", label.vjust = -1, keyheight = unit(0.1, "in"), keywidth = unit(0.4, "in"), label.theme = element_text(size = 8, face = "bold"), direction = "horizontal", reverse = T)) +
    scale_x_continuous(name = "Proportion", limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1), labels = c("0.25", "0.50", "0.75", "1"), expand = c(0.01,0)) + 
    ggtitle("Abundance Proportions By Predominant Organism", subtitle = "Reference Lab, All Sites, All Procedures") +
    theme(legend.position = c(1, 1.05), legend.justification = "right", legend.title = element_blank(), axis.title.y.left = element_blank(), axis.text.y.left = element_text(face = "bold", color = "black"), plot.title = element_text(size = 16, face = "bold.italic", margin = margin(0, 0, 4, 0), hjust = 0), plot.subtitle = element_text(size = 14, face = "plain", hjust = 0, margin = margin(0, 0, 2, 0)))
  ggsave(paste0(figure_path, "final_figures/abundance_proportions_by_organism.svg"), gg_abundance_by_organism, width = 8.5, height = 4, units = "in")
  ggsave(paste0(figure_path, "final_figures/abundance_proportions_by_organism.png"), gg_abundance_by_organism, width = 8.5, height = 4, units = "in", dpi = 600)
  
  gg_abundance_by_organism
  
}

makeCultureResultTables <- function(input_by_organism = tar_read(culture_results_joined_by_organism), input_by_subject = tar_read(culture_results_joined_by_subject)){
  
  gt_culture_organisms_by_site <- 
    input_by_organism |> 
      tbl_summary(by = Site, missing = "no", include = c("Procedure", "Site", "Organism", "Growth Significance")) |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_culture_organisms_by_site |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Site, All Sources, One Row Per Organism"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_organism_per_site_table.html"))
  gt::gtsave(gt_culture_organisms_by_site |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Site, All Sources, One Row Per Organism"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_organism_per_site_table.docx"))

  gt_culture_organisms_by_procedure <- 
    input_by_organism |> 
      tbl_summary(by = Procedure, missing = "no", include = c("Procedure", "Site", "Organism", "Growth Significance")) |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_culture_organisms_by_procedure |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Procedure, All Sites, One Row Per Organism"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_organism_per_procedure_table.html"))
  gt::gtsave(gt_culture_organisms_by_procedure |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Procedure, All Sites, One Row Per Organism"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_organism_per_procedure_table.docx"))
  
  gt_culture_organisms_by_significance <- 
    input_by_organism |> 
    filter(`Growth Significance` != "No Growth") |> 
    mutate(`Growth Significance` = factor(`Growth Significance`, levels = c("Insignificant", "Significant"))) |> 
    tbl_summary(by = `Growth Significance`, missing = "no", include = c("Organism"), percent = "row") |> 
    italicize_levels() |>
    bold_labels() |> 
    add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
    bold_p() |>
    add_overall(statistic = ~"{n}") |>  
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |> 
    as_gt() |>
    tab_header(title = "Culture Results (MALDI-ToF) by Significance", subtitle = '"Significant" = >10k CFU for BAL, >100k CFU for ETA') |> 
    theme_pattern() 
  gt::gtsave(gt_culture_organisms_by_significance, 
             filename = paste0(figure_path, "final_figures/culture_results_by_significance_table.html"))
  gt::gtsave(gt_culture_organisms_by_significance |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Procedure, All Sites, One Row Per Organism"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_organism_per_procedure_table.docx"))
  
  gt_culture_subjects_by_site <-
    input_by_subject |> 
      select(-SubjectID, -Procedure) |> 
      tbl_summary(by = Site, missing = "no") |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_culture_subjects_by_site |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Site, All Procedures, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_subject_per_site_table.html"))
  gt::gtsave(gt_culture_subjects_by_site |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Site, All Procedures, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_subject_per_site_table.docx"))
  
  gt_culture_subjects_by_procedure <-
    input_by_subject |> 
    select(-SubjectID, -Site) |> 
    tbl_summary(by = Procedure, missing = "no") |> 
    italicize_levels() |>
    bold_labels() |> 
    add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
    bold_p() |>
    add_overall() |>  
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_culture_subjects_by_procedure |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_subject_per_procedure_table.html"))
  gt::gtsave(gt_culture_subjects_by_procedure |> as_gt() |> tab_header(title = "Culture Results", subtitle = "By Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_culture_results_by_subject_per_procedure_table.docx"))
  
  

}

makeCultureResultTreeMap <- function(input = tar_read(culture_results_joined_by_organism)){
  
  library(treemapify)
  
  input |> 
    count(Organism, `Organism Gram Stain`) |>
    drop_na(`Organism Gram Stain`) |>
    ggplot(aes(area = n, fill = `Organism Gram Stain`, label = Organism, subgroup = `Organism Gram Stain`)) +
      geom_treemap(layout = "squarified") + 
      geom_treemap_text(colour = "black", place = "centre", reflow = TRUE)
    
  
}

makeMixedMicroFigures <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  gg_mixed_microorganisms <- 
    input |> 
      filter(`Predominant Morphology (w/ Other)` == "Mixed microorganism") |> 
      ggplot(aes(y = `Unique Organisms Grown`, fill = `Organism Abundance`)) +
      geom_bar() + 
      theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 10), plot.title = element_text(size = 12, face = "bold.italic"), plot.subtitle = element_text(size = 10, face = "plain")) + 
      labs(title = "Predominant Organism Counts From Mixed Gram Stains", subtitle = "Reference Lab, By Site, All Procedures", x = "Proportion of Mixed Gram Stains", y = "Site")
  ggsave(paste0(figure_path, "final_figures/mixed_microorganism_predominant_counts_by_site.svg"), gg_mixed_microorganisms, width = 8.5, height = 4, units = "in")
  ggsave(paste0(figure_path, "final_figures/mixed_microorganism_predominant_counts_by_site.png"), gg_mixed_microorganisms, width = 8.5, height = 4, units = "in", dpi = 600)
  
  gg_mixed_microorganisms
  
  input |> 
    filter(`Predominant Morphology (w/ Other)` == "Mixed microorganism" & `Organism Abundance` %in% c("Rare", "Few", "None")) |> 
    count(MALDI_ID_COLLAPSED, sort = T)
  
}

makeSignificantGrowthCountBarchart <- function(input = tar_read(culture_results_joined_by_organism)){
  
  gg_bar_input_sig <- 
    input |> 
      filter(!is.na(Organism) & Organism != "Other" & Organism != "Yeast" & `Growth Significance` == "Significant") |>
      count(Organism, `Organism Gram Stain`, sort = T)
  
  gg_bar_input_all <- 
    input |> 
    filter(!is.na(Organism) & Organism != "Other" & Organism != "Yeast") |>
    count(Organism, `Organism Gram Stain`, sort = T)
  
  gg_bar <- 
    ggplot() +
      geom_col(data = gg_bar_input_all, aes(x = n, y = fct_rev(Organism), fill = `Organism Gram Stain`), alpha = 0.25) +
      geom_col(data = gg_bar_input_sig, aes(x = n, y = fct_rev(Organism), fill = `Organism Gram Stain`)) +
      geom_text(data = gg_bar_input_all[1:2, ], aes(x = n * 0.75, y = fct_rev(Organism)), label = "All", fontface = "bold", alpha = 0.75) + 
      geom_text(data = gg_bar_input_sig[1:2, ], aes(x = n * 0.5, y = fct_rev(Organism)), label = "Significant", color = "white", fontface = "bold") + 
      geom_text(data = gg_bar_input_all, aes(x = n + 1, y = fct_rev(Organism), label = n), hjust = 0) + 
      geom_text(data = gg_bar_input_sig, aes(x = n * 0.95, y = fct_rev(Organism), label = n), hjust = 1, color = "white") + 
      scico::scale_fill_scico_d(palette = "acton", begin = 0.6, end = 0) +
      theme(legend.position = "none", axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank(), axis.text.y.left = element_text(face = "italic"))
  
  ggsave(paste0(figure_path, "final_figures/significant_growth_counts_by_organism.svg"), gg_bar, width = 10, height = 8, bg = "white")
  ggsave(paste0(figure_path, "final_figures/significant_growth_counts_by_organism.png"), gg_bar, width = 10, height = 8, dpi = 600, bg = "white")
  
}

makeSignificantGrowthCountTreemap <- function(input = tar_read(culture_results_joined_by_organism)){

  library(treemapify)
  
  order <- 
    input |> 
      mutate(Organism = ifelse(Organism %in% c("Other", "Yeast"), paste0(" ", Organism), as.character(Organism))) |>
      count(`Growth Significance`, Organism, sort = T) |>  
      filter(`Growth Significance` == "Significant") |> 
      pluck("Organism")
  
  split <- ifelse(!grepl("Yeast|Other", as.character(input$Organism)), str_split(as.character(input$Organism), " "), as.character(input$Organism))
  genus <- ifelse(!grepl("Yeast|Other", as.character(input$Organism)) & !is.na(input$Organism), paste0(split |> map_chr(1) |> str_sub(1, 1), ". "), "")
  species <- map(split, ~ ifelse(length(.x) > 1, pluck(.x,2), ""))
  labels = str_c(genus, species)
  input$label = labels |> as_tibble() |> mutate(label = ifelse(grepl("Other|Yeast", input$Organism), as.character(input$Organism), value)) |> pluck("label")
  
  gg_significant_growth_counts_input <- 
    input |> 
      filter(!(is.na(`Organism Gram Stain`) & `Growth Significance` != "No Growth")) |> 
      select(label, `Growth Significance`, `Organism Gram Stain`, Procedure) |> 
      group_by(`Organism Gram Stain`) |>
      mutate(label = fct_lump_n(label, 10)) |>
      ungroup() |>
      mutate(label = 
               case_when(
                 label == "Other" & `Organism Gram Stain` == "Positive" ~ "Other Gram-Positive",
                 label == "Other" & `Organism Gram Stain` == "Negative" ~ "Other Gram-Negative",
                 `Organism Gram Stain` == "Yeast/Mold" ~ "Yeast/Mold",
                 is.na(`Organism Gram Stain`) ~ "No Growth",
                 T ~ label
               ), 
             label = str_remove_all(label, "-calcoaceticus"),
             `Organism Gram Stain` = ifelse(is.na(`Organism Gram Stain`), "No Growth", `Organism Gram Stain`)) |> 
      count(label, `Organism Gram Stain`, sort = T) |> 
      distinct() |>
      mutate(prop_of_all = n/sum(n)) |>
      group_by(label) |>
      mutate(total_organism = sum(n),
             prop_significant = n/total_organism) |>
      ungroup()
  
  gg_culture_results_treemap <-  
    gg_significant_growth_counts_input |> 
      mutate(`Organism Gram Stain` = factor(`Organism Gram Stain`, levels = c("Positive", "Yeast/Mold", "Negative", "No Growth"))) |>
      ggplot(aes(area = total_organism, fill = `Organism Gram Stain`, label = label, subgroup = `Organism Gram Stain`)) +
        geom_treemap(layout = "srow", start = "topleft") + 
        geom_treemap_text(layout = "srow", colour = "white", place = "topleft", fontface = "bold.italic", size = 24, start = "topleft") +
        geom_treemap_text(aes(label = total_organism), layout = "srow", colour = "white", place = "bottomright", size = 36, start = "topleft") +
        geom_treemap_subgroup_border(layout = "srow", colour = "white", size = 5, start = "topleft") +
        scico::scale_fill_scico_d(palette = "acton", begin = 0, end = 0.8) +
        theme(legend.position = "none")
  ggsave(paste0(figure_path, "final_figures/culture_growth_tile_map_all.svg"), gg_culture_results_treemap, width = 15, height = 9, dpi = 600)
  ggsave(paste0(figure_path, "final_figures/culture_growth_tile_map_all.png"), gg_culture_results_treemap, width = 15, height = 9, dpi = 600)
  
  
  gg_significant_growth_counts <-   
    ggplot(gg_significant_growth_counts_input, aes(y = Organism, x = n, fill = `Organism Gram Stain`, alpha = `Growth Significance`)) + 
      geom_col(linewidth = 0) + 
      geom_text(data = gg_significant_growth_counts_input |> filter(`Growth Significance` == "Significant"), aes(label = n), nudge_x = -1, hjust = 1, fontface = "bold.italic", color = "white") + 
      geom_text(data = gg_significant_growth_counts_input |> filter(`Growth Significance` == "Significant"), aes(x = total, label = total), alpha = 0.25, nudge_x = 1, hjust = 0, fontface = "italic") + 
      annotate("text", label = "Significant", fontface = "italic", x = c(26, 24), y = c(16, 15), hjust = 0.5, color = "white") + 
      annotate("text", label = "All", fontface = "italic", x = c(78, 94), y = c(16, 15), hjust = 0.5, color = "gray65") + 
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) + 
      scale_y_discrete() +
      scico::scale_fill_scico_d(palette = "acton", begin = 0.6, end = 0, na.value = "gray75") +
      scale_alpha_manual(values = c(0.2, 1)) +
      coord_cartesian(clip = "off") +
      ggtitle("Organism Counts in Culture", subtitle = "Significant: >10,000 CFU for BAL, >100,000 CFU for ETA") + 
      theme(legend.position = "none", 
            axis.title.y.left = element_blank(), axis.text.y.left = element_text(size = 10, face = "italic", color = "black"), axis.line.y.left = element_blank(),
            axis.title.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank())
  
  gg_significant_growth_counts_input_by_procedure <- 
    input |> 
      filter(!is.na(Organism)) |> 
      select(Organism, `Growth Significance`, Procedure) |> 
      count(Organism, `Growth Significance`, Procedure, sort = T) |>
      distinct() |>
      pivot_wider(names_from = `Growth Significance`, values_from = n, values_fill = 0) |>
      mutate(Total = `Significant` + `Insignificant`, 
             `Significant` = glue::glue("{`Significant`} ({round(`Significant`/Total*100, 0)}%)"))
  
  gt_culture_counts_by_procedure <- 
    left_join(
    gg_significant_growth_counts_input_by_procedure |> filter(Procedure == "ETA") |> select(-Procedure, -Insignificant),
    gg_significant_growth_counts_input_by_procedure |> filter(Procedure == "BAL") |> select(-Procedure, -Insignificant),
    by = "Organism", suffix = c("_ETA", "_BAL")
  ) |> 
    filter(Organism != "Yeast" & Organism != "Other") |>
    # replace NA with 0 is all colums
    replace_na(list(Significant_ETA = "0", Significant_BAL = "0", Total_ETA = 0, Total_BAL = 0)) |>
    arrange(desc(Total_ETA + Total_BAL)) |>
    gt() |> 
    tab_spanner(label = md("***ETA***"), columns = 2:3) |>
    tab_spanner(label = md("***BAL***"), columns = 4:5) |>
    cols_label(Organism = "Organism", Total_ETA = "Total", Significant_ETA = "Significant", Total_BAL = "Total", Significant_BAL = "Significant") |>
    #right justify Organism column text
    tab_style(style = cell_text(align = "right", style = "italic"), locations = cells_body(columns = "Organism")) |>
    #bold and italicize column labels
    tab_style(style = cell_text(weight = "bold", align = "right"), locations = cells_column_labels())
  gtsave(gt_culture_counts_by_procedure, paste0(figure_path, "final_figures/culture_counts_by_procedure_table.png"))
  gtsave(gt_culture_counts_by_procedure, paste0(figure_path, "final_figures/culture_counts_by_procedure_table.html"))
  gtsave(gt_culture_counts_by_procedure, paste0(figure_path, "final_figures/culture_counts_by_procedure_table.docx"))
  
  
  gg_significant_growth_counts_by_procedure <-   
    ggplot(gg_significant_growth_counts_input_by_procedure, aes(y = Organism, x = n, fill = `Organism Gram Stain`, alpha = `Growth Significance`)) + 
    geom_col(linewidth = 0) + 
    geom_text(data = gg_significant_growth_counts_input_by_procedure |> filter(`Growth Significance` == "Significant"), aes(label = n), nudge_x = -1, hjust = 1, fontface = "bold.italic", color = "white") + 
    geom_text(data = gg_significant_growth_counts_input_by_procedure |> filter(`Growth Significance` == "Significant"), aes(x = total, label = total), alpha = 0.25, nudge_x = 1, hjust = 0, fontface = "italic") + 
    facet_wrap(~fct_rev(Procedure), ncol = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_y_discrete() +
    scico::scale_fill_scico_d(palette = "acton", begin = 0.6, end = 0, na.value = "gray75") +
    scale_alpha_manual(values = c(0.25, 1)) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none", 
          axis.title.y.left = element_blank(), axis.text.y.left = element_text(size = 10, face = "italic", color = "black"), axis.line.y.left = element_blank(),
          axis.title.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank())

  gg_significant_growth_combo <- gg_significant_growth_counts + patchwork::inset_element(gg_significant_growth_counts_by_procedure, left = 0.3, bottom = 0, right = 1, top = 0.7)
  ggsave(paste0(figure_path, "final_figures/significant_growth_counts_by_organism.svg"), gg_significant_growth_combo, width = 16, height = 9.5)
  ggsave(paste0(figure_path, "final_figures/significant_growth_counts_by_organism.png"), gg_significant_growth_combo, width = 16, height = 9.5, dpi = 600)
  
  gg_significant_growth_combo
  
}

makeSignificantGrowthByMorphologyFigures <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  gg_significant_growth_by_predominant_organism_abundance <- 
    input |> 
    filter(!is.na(`Predominant Morphology (w/ Other)`)) |>
    mutate(Abundance_bins = factor(`Organism Abundance`, levels = c("Rare", "Few", "Moderate", "Abundant"), labels = c("Rare/Few", "Rare/Few", "Moderate/Abundant", "Moderate/Abundant")),
           `Predominant Morphology (w/ Other)` = factor(fct_other(`Predominant Morphology (w/ Other)`, 
                             keep = c("Gram-positive cocci", "Gram-negative bacilli", "Mixed microorganism", "Other")),
                          levels = c("Other", "Gram-negative bacilli", "Mixed microorganism", "Gram-positive cocci"))) |>
    ggplot(aes(y = `Organism Abundance`, fill = fct_rev(`Predominant Morphology (w/ Other)`), alpha = `Any Significant Growth`)) + 
    geom_bar() +
    facet_wrap(~fct_rev(`Predominant Morphology (w/ Other)`)) +
    guides(alpha = guide_legend(title = "", reverse = T, direction = "horizontal"), fill = guide_none()) +
    scico::scale_fill_scico_d(palette = "acton", begin = 0, end = 0.5) +
    scale_alpha_manual(values = c(0.5, 1)) + 
   # scale_y_discrete(labels = c("None/Insignificant", "Significant Growth")) +
    scale_x_continuous(name = "Count") +
    coord_cartesian(clip = "off") +
    ggtitle("Significant Growth by Gram Stain Morphology and Abundance") +
    theme(legend.position = "none", legend.text.position = "top", axis.title.y.left = element_blank(), legend.background = element_blank(), legend.box = element_blank())
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_predominant_organism_abundance.svg"), gg_significant_growth_by_predominant_organism_abundance, width = 8, height = 6)
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_predominant_organism_abundance.png"), gg_significant_growth_by_predominant_organism_abundance, width = 8, height = 6, dpi = 600)
  
  gg_significant_growth_by_predominant_organism_abundance
  
  gg_concordant_growth_by_predominant_organism_abundance <- 
    input |> 
    filter(!is.na(`Predominant Morphology (w/ Other)`)) |>
    mutate(Abundance_bins = factor(`Organism Abundance`, levels = c("Rare", "Few", "Moderate", "Abundant"), labels = c("Rare/Few", "Rare/Few", "Moderate/Abundant", "Moderate/Abundant")),
           `Predominant Morphology (w/ Other)` = factor(fct_other(`Predominant Morphology (w/ Other)`, 
                                                                  keep = c("Gram-positive cocci", "Gram-negative bacilli", "Mixed microorganism", "Other")),
                                                        levels = c("Other", "Gram-negative bacilli", "Mixed microorganism", "Gram-positive cocci"))) |>
    ggplot(aes(y = `Organism Abundance`, fill = fct_rev(`Predominant Morphology (w/ Other)`), alpha = `Any Concordant Growth`)) + 
    geom_bar() +
    facet_wrap(~fct_rev(`Predominant Morphology (w/ Other)`)) +
    guides(alpha = guide_legend(title = "", reverse = T, direction = "horizontal"), fill = guide_none()) +
    scico::scale_fill_scico_d(palette = "acton", begin = 0, end = 0.5) +
    scale_alpha_manual(values = c(0.5, 1)) + 
    # scale_y_discrete(labels = c("None/Insignificant", "Significant Growth")) +
    scale_x_continuous(name = "Count") +
    coord_cartesian(clip = "off") +
    ggtitle("Concordant Growth by Gram Stain Morphology and Abundance") +
    theme(legend.position = "none", legend.text.position = "top", axis.title.y.left = element_blank(), legend.background = element_blank(), legend.box = element_blank())
  ggsave(paste0(figure_path, "rachel_asm_concordant_growth_by_predominant_organism_abundance.svg"), gg_concordant_growth_by_predominant_organism_abundance, width = 8, height = 6)
  ggsave(paste0(figure_path, "rachel_asm_concordant_growth_by_predominant_organism_abundance.png"), gg_concordant_growth_by_predominant_organism_abundance, width = 8, height = 6, dpi = 600)
  
  list(gg_significant_growth_by_predominant_organism_abundance, gg_concordant_growth_by_predominant_organism_abundance)
  
}

makeConcordantGrowthCountFigures <- function(gram_culture_joined = tar_read(gram_culture_joined)){
  
  gg_input <- 
    gram_culture_joined |> 
    group_by(SubjectID, MALDI_ID_GROUPED) |>
    mutate(concordant_growth_organism = any(gram_equals_culture, na.rm = T)) |> 
    mutate(MALDI_ID_COLLAPSED = ifelse(MALDI_ID_COLLAPSED %in% c("Other", "Yeast"), paste0(" ", MALDI_ID_COLLAPSED), as.character(MALDI_ID_COLLAPSED))) |> 
    ungroup()
  
  order <- 
    gg_input |> 
    count(concordant_growth_organism, MALDI_ID_COLLAPSED, sort = T) |>  
    filter(concordant_growth_organism) |> 
    pluck("MALDI_ID_COLLAPSED")
  
  split <- str_split(order, " ")
  genus <- paste0(split |> map_chr(1) |> str_sub(1, 1), ". ")
  species <- split |> map_chr(2)
  labels = str_c(genus, species)
  
  gg_concordant_growth_counts_input <- 
    gg_input |> 
    filter(!is.na(MALDI_ID_COLLAPSED)) |> 
    select(MALDI_ID_COLLAPSED, concordant_growth_organism, LRT_REF_GS_SPECTYPE) |> 
    count(MALDI_ID_COLLAPSED, concordant_growth_organism, sort = T) |> 
    left_join(gram_culture_joined |> select(MALDI_ID_COLLAPSED, culture_gram) |> distinct()) |> 
    mutate(culture_gram = ifelse(MALDI_ID_COLLAPSED %in% c("Other", "Yeast"), NA, culture_gram),
           culture_gram = ifelse(MALDI_ID_COLLAPSED == "Streptococcus pneumoniae", "Positive", culture_gram)) |> 
    distinct() |>
    group_by(MALDI_ID_COLLAPSED) |>
    mutate(MALDI_ID_COLLAPSED = fct_rev(factor(MALDI_ID_COLLAPSED, levels = order, labels = labels)),
           total = sum(n),
           prop = n/total) |> 
    ungroup() |> 
    drop_na()
  
  
  gg_concordant_growth_counts <-   
    ggplot(gg_concordant_growth_counts_input, aes(y = MALDI_ID_COLLAPSED, x = n, fill = culture_gram, alpha = concordant_growth_organism)) + 
    geom_col(linewidth = 0) + 
    geom_text(data = gg_concordant_growth_counts_input |> filter(concordant_growth_organism), aes(label = n), nudge_x = -1, hjust = 1, fontface = "bold.italic", color = "white") + 
    geom_text(data = gg_concordant_growth_counts_input |> filter(!concordant_growth_organism), aes(x = total, label = total), alpha = 0.25, nudge_x = 1, hjust = 0, fontface = "italic") + 
    annotate("text", label = "Concordant", fontface = "italic", x = c(55.5, 33), y = c(16, 15), hjust = 0.5, color = "white") + 
    annotate("text", label = "All", fontface = "italic", x = c(122.5, 83.5), y = c(16, 15), hjust = 0.5, color = "gray65") + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) + 
    scale_y_discrete() +
    scico::scale_fill_scico_d(palette = "acton", begin = 0.6, end = 0, na.value = "gray75") +
    scale_alpha_manual(values = c(0.2, 1)) +
    coord_cartesian(clip = "off") +
    ggtitle("Organism Counts in Culture by Concordance with Gram Stain") + 
    theme(legend.position = "none", 
          axis.title.y.left = element_blank(), axis.text.y.left = element_text(size = 10, face = "italic", color = "black"), axis.line.y.left = element_blank(),
          axis.title.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank())
  
  
  gg_concordant_growth_counts_input_by_procedure <- 
    gg_input |> 
    filter(!is.na(MALDI_ID_COLLAPSED)) |> 
    select(MALDI_ID_COLLAPSED, concordant_growth_organism, LRT_REF_GS_SPECTYPE) |> 
    count(MALDI_ID_COLLAPSED, concordant_growth_organism, LRT_REF_GS_SPECTYPE, sort = T) |> 
    left_join(gram_culture_joined |> select(MALDI_ID_COLLAPSED, culture_gram) |> distinct()) |> 
    mutate(culture_gram = ifelse(MALDI_ID_COLLAPSED %in% c("Other", "Yeast"), NA, culture_gram),
           culture_gram = ifelse(MALDI_ID_COLLAPSED == "Streptococcus pneumoniae", "Positive", culture_gram)) |> 
    distinct() |>
    group_by(MALDI_ID_COLLAPSED, LRT_REF_GS_SPECTYPE) |>
    mutate(MALDI_ID_COLLAPSED = fct_rev(factor(MALDI_ID_COLLAPSED, levels = order, labels = labels)),
           total = sum(n),
           prop = n/total) |> 
    ungroup() |> 
    drop_na()
  
  gg_concordant_growth_counts_by_procedure <-   
    ggplot(gg_concordant_growth_counts_input_by_procedure, aes(y = MALDI_ID_COLLAPSED, x = n, fill = culture_gram, alpha = concordant_growth_organism)) + 
    geom_col(linewidth = 0) + 
    geom_text(data = gg_concordant_growth_counts_input_by_procedure |> filter(concordant_growth_organism), aes(label = n), nudge_x = -1, hjust = 1, fontface = "bold.italic", color = "white") + 
    geom_text(data = gg_concordant_growth_counts_input_by_procedure |> filter(!concordant_growth_organism), aes(x = total, label = total), alpha = 0.25, nudge_x = 1, hjust = 0, fontface = "italic") + 
    facet_wrap(~fct_rev(LRT_REF_GS_SPECTYPE), ncol = 2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_y_discrete() +
    scico::scale_fill_scico_d(palette = "acton", begin = 0.6, end = 0, na.value = "gray75") +
    scale_alpha_manual(values = c(0.25, 1)) +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none", 
          axis.title.y.left = element_blank(), axis.text.y.left = element_text(size = 10, face = "italic", color = "black"), axis.line.y.left = element_blank(),
          axis.title.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank())

  gg_concordant_combo <- gg_concordant_growth_counts + patchwork::inset_element(gg_concordant_growth_counts_by_procedure, left = 0.3, bottom = 0, right = 1, top = 0.7, align_to = "plot")
  ggsave(paste0(figure_path, "final_figures/concordant_growth_counts_by_organism_combined.svg"), gg_concordant_combo, width = 16, height = 9.5)
  ggsave(paste0(figure_path, "final_figures/concordant_growth_counts_by_organism_combined.png"), gg_concordant_combo, width = 16, height = 9.5, dpi = 600)
  
  gg_concordant_combo
  
}


makeConcordantGrowthByMorphologyFigures <- function(gram_culture_joined = tar_read(gram_culture_joined)){
  
  gram_culture_joined |> 
    group_by(SubjectID) |>
    mutate(any_concordant = any(gram_equals_culture, na.rm = T)) |>
    select(SubjectID, LRT_REF_GS_MORPH, LRT_REF_GS_ORG_AMT, any_concordant) |> 
    filter(LRT_REF_GS_MORPH %in% c("Gram-positive cocci", "Gram-negative bacilli", "Mixed microorganism")) |> 
    mutate(Abundance_bins = factor(LRT_REF_GS_ORG_AMT, levels = c("Rare", "Few", "Moderate", "Abundant"), labels = c("Rare/Few", "Rare/Few", "Moderate/Abundant", "Moderate/Abundant")),
           Morph = factor(fct_other(LRT_REF_GS_MORPH, 
                                    keep = c("Gram-positive cocci", "Gram-negative bacilli", "Mixed microorganism", "Other")),
                          levels = c("Other", "Gram-negative bacilli", "Mixed microorganism", "Gram-positive cocci"))) |>
    ggplot(aes(y = any_concordant, fill = fct_rev(Morph), alpha = Abundance_bins)) + 
    geom_bar() +
    facet_wrap(~fct_rev(Morph)) +
    guides(alpha = guide_legend(title = "", reverse = T, direction = "horizontal"), fill = guide_none()) +
    scico::scale_fill_scico_d(palette = "acton", begin = 0, end = 0.5) +
    scale_alpha_manual(values = c(0.5, 1)) + 
    scale_y_discrete(labels = c("No/Inconcordant Growth", "Concordant Growth")) +
    scale_x_continuous(name = "Count") +
    coord_cartesian(clip = "off") +
    ggtitle("Concordant Growth by Gram Stain Morphology and Abundance") +
    theme(legend.position = "none", legend.text.position = "top", axis.title.y.left = element_blank(), legend.background = element_blank(), legend.box = element_blank())
  ggsave(paste0(figure_path, "final_figures/concordant_growth_by_morphology_and_abundance.svg"), width = 8, height = 3)
  ggsave(paste0(figure_path, "final_figures/concordant_growth_by_morphology_and_abundance.png"), width = 8, height = 3, dpi = 600)
  
}

makeModelFigures <- function(input = tar_read(gram_culture_joined_by_subject), target = "Any Significant Growth", procedure = "All", out_path = paste0(figure_path, "final_figures/any_significant_growth_all")){
  
  library(tidymodels)
  library(geomtextpath)
  library(gridExtra)
  library(patchwork)
  library(grid)
  library(gtsummary)
  
  input$target <- factor(input[[target]], levels = c(FALSE, TRUE))
  model_input <- 
    input |> ungroup() |> select(-matches("Any|Predom|Grow|Subje|No Org")) |> 
      mutate(
        across(c(Location, Sex), ~fct_infreq(.x)),
        across(is.ordered, ~factor(as.character(.x), levels = c("None", "Rare", "Few", "Moderate", "Abundant"), ordered = F))) 
  
  if(procedure != "All"){ model_input <- model_input |> filter(Procedure == procedure) |> select(-Procedure)}
  
  ### Make UV Tables
  # gt_uv <- 
  #   model_input |> 
  #     tbl_uvregression(y = target, method = glm, method.args = list(family = binomial), exponentiate = T) |>
  #       bold_labels() |> 
  #       italicize_levels() |> 
  #       add_n(location = "level") |> 
  #       add_q() |> 
  #       bold_p(q = T)
  # gt::gtsave(as_gt(gt_uv), paste0(out_path, "_uv_table.html"))
  # gt::gtsave(as_gt(gt_uv), paste0(out_path, "_uv_table.docx"))
  
  ### Make MV Tables
  mv_model <- glm(target ~ ., data = model_input, family = binomial)
  gt_mv <- 
    tbl_regression(mv_model, exponentiate = T) |>
      bold_labels() |> 
      italicize_levels() |> 
      add_n(location = "level") |> 
      add_q() |> 
      bold_p(q= T)
  gt::gtsave(as_gt(gt_mv), paste0(out_path, "_mv_table.html"))
  gt::gtsave(as_gt(gt_mv), paste0(out_path, "_mv_table.docx"))
  
  ### Convert binned quantifications to numerics
  roc_inputs <- 
    model_input |> 
      select(`Epithelial Cells`:`Unique Morphologies Seen`, target) |> 
      mutate(across(-target, ~as.numeric(.x)))
  
  ### Calculate auROCs
  auc_rocs <- 
    map(names(roc_inputs |> select(-target)), ~roc_auc(roc_inputs, .x, truth = target, event_level = "first") |> mutate(Feature = .x)) |> 
    bind_rows()
  
  ### Flip negatively correlated features
  features_to_flip <- auc_rocs |> filter(.estimate < 0.5) |> pluck("Feature")
  features_to_keep <- auc_rocs |> filter(.estimate >= 0.5) |> pluck("Feature")
  auc_rocs <- auc_rocs |> mutate(.estimate = round(if_else(.estimate < 0.5, 1 - .estimate, .estimate), digits = 2))
  
  ### Make ROC curves
  roc_curves_pos <- map(features_to_flip, ~roc_curve(roc_inputs, .x, truth = target, event_level = "second") |> mutate(Feature = .x)) |> bind_rows()
  roc_curves_neg <- map(features_to_keep, ~roc_curve(roc_inputs, .x, truth = target, event_level = "first") |> mutate(Feature = .x)) |> bind_rows()
  roc_curves <- bind_rows(roc_curves_pos, roc_curves_neg)
  
  ### Add metadata for plot
  feature_levels <- c("Gram-Positive Cocci", "Yeast/Fungus", "Mixed Microorganisms", "Gram-Negative Bacilli", "Other Morphology", 
                      "Organism Abundance", "Unique Morphologies Seen", "Polymorphonucleocytes", "Epithelial Cells", "Red Blood Cells")
  colors <- c(scico::scico(5, palette = "acton", begin = 0, end = 0.8), "black", "grey33", scico::scico(3, palette = "bilbao", begin = 0.9, end = 0.2))
  roc_curves <- roc_curves |> left_join(auc_rocs) |> mutate(Feature = factor(Feature, levels = feature_levels), 
                                              Type = factor(ifelse(Feature %in% c("Unique Morphologies Seen", "Epithelial Cells", "Polymorphonucleocytes", "Organism Abundance", "Red Blood Cells"), 
                                                                   "Metadata", "Morphology")),
                                              Label = paste0(Feature, "\n(AUC = ", .estimate, ")"))
  legend_table <- tibble(color = colors, Feature = feature_levels)
  roc_curves <- roc_curves |> left_join(legend_table)
  
  ### Make AUC tables
  auc_table_metadata <- 
    tableGrob(
      auc_rocs |> 
        arrange(desc(.estimate)) |> 
        filter(Feature %in% c("Unique Morphologies Seen", "Epithelial Cells", "Polymorphonucleocytes", "Organism Abundance", "Red Blood Cells")) |> 
        select(AUC = .estimate, Feature), 
      rows = NULL,
      theme = ttheme_minimal(base_size = 8,
                             core=list(fg_params=list(hjust=0, x=0.01), bg_params=list(fill=NA), padding = unit(c(1.5, 1.5), "mm")), 
                             colhead=list(fg_params=list(hjust=0, x=0.01), bg_params=list(fill=NA), padding = unit(c(1.5, 1.5), "mm"))))
  
  auc_table_morphology <- 
    tableGrob(
      auc_rocs |> 
        arrange(desc(.estimate)) |> 
        filter(!Feature %in% c("Unique Morphologies Seen", "Epithelial Cells", "Polymorphonucleocytes", "Organism Abundance", "Red Blood Cells")) |> 
        select(AUC = .estimate, Feature), 
      rows = NULL,
      theme = ttheme_minimal(base_size = 8, 
                             core=list(fg_params=list(hjust=0, x=0.01), bg_params=list(fill=NA), padding = unit(c(1.5, 1.5), "mm")), 
                             colhead=list(fg_params=list(hjust=0, x=0.01), bg_params=list(fill=NA), padding = unit(c(1.5, 1.5), "mm"))))
    
  
  ### Plot ROC curves
  gg_roc_metadata <- 
    ggplot(roc_curves |> filter(Type == "Metadata"), aes(1-specificity, sensitivity, color = color)) + 
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray75", alpha = 0.5) + 
      geom_textline(aes(label = Feature), linewidth = 1.5, fontface = 'bold', show.legend = T, hjust = 0.4, size = 3) +
      scale_color_identity() + 
      xlab("1 - Specificity") + ylab("Sensitivity") +
      theme()

  gg_roc_morphology <- 
    ggplot(roc_curves |> filter(Type == "Morphology"), aes(1-specificity, sensitivity, color = color)) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray75", alpha = 0.5) + 
    geom_textline(aes(label = Feature, vjust = Feature, hjust = Feature), linewidth = 1.5, fontface = 'bold', show.legend = T, size = 3) +
    scale_color_identity() +
    scale_vjust_manual(values = c(0.5, 0.5, 0.5, 2, 2)) + 
    scale_hjust_manual(values = c(0.5, 0.3, 0.65, 0.7, 0.3)) + 
    xlab("1 - Specificity") + ylab("Sensitivity") +
    theme()

  gg_roc_all <- 
    ggpubr::ggarrange(
      gg_roc_metadata + inset_element(auc_table_metadata, clip = F, 0.6, 0.3, 0.9, 0.05), 
      gg_roc_morphology + inset_element(auc_table_morphology, clip = F, 0.6, 0.3, 0.9, 0.05), 
    ncol = 2)
  ggsave(paste0(out_path, "_roc_curves.svg"), gg_roc_all, width = 8.5, height = 4, units = "in")
  ggsave(paste0(out_path, "_roc_curves.png"), gg_roc_all, width = 8.5, height = 4, units = "in", dpi = 600)
  
  forest_plot <- 
    forestmodel::forest_model(mv_model, 
                              exponentiate = T, 
                              factor_separate_line = T,
                              format_options = forestmodel::forest_model_format_options(color = "black", shape = 15, text_size = 4, point_size = 2, banded = TRUE))
  ggsave(paste0(out_path, "_forest_plot.svg"), forest_plot, width = 8.5, height = 8, units = "in")
  ggsave(paste0(out_path, "_forest_plot.png"), forest_plot, width = 8.5, height = 8, units = "in", dpi = 600)
  
  thresholds <-  roc_inputs |> summarise(across(everything(), ~list(0:length(unique(.)))))

  library(probably)
  
  predictive_values <- 
    roc_inputs |>
      select(`Epithelial Cells`:`Organism Abundance`, target) |>
      pivot_longer(cols = -c(target)) |> 
      group_by(name) %>% 
      reframe(threshold_perf(., 
        estimate = value, 
        truth = target,
        metrics = metric_set(ppv, npv),
        thresholds = thresholds[[cur_group_id()]][[1]],
        event_level = "second")) |> 
      select(-.estimator) |>
      filter(.threshold > 0) |> 
      distinct() |> 
      mutate(.threshold = factor(.threshold, labels = c("None", "Rare", "Few", "Moderate", "Abundant")))

  predictive_values_plot <- 
    predictive_values |> 
      ggplot(aes(x = .estimate, y = name, fill = fct_rev(.threshold)), group = fct_rev(.threshold)) + 
      geom_col(position = position_dodge(), alpha = 0.9) + 
      geom_text(aes(label = round(.estimate, 2), x = .estimate + 0.01, color = fct_rev(.threshold)), position = position_dodge(0.9), hjust = 0, fontface = "bold") + 
      geom_text(data = predictive_values |> filter(.metric == "npv"), aes(label = .threshold, x = 0.05, color = fct_rev(.threshold),), position = position_dodge(0.9), hjust = 0, fontface = "bold.italic") + 
      facet_wrap(~.metric, labeller = as_labeller(c(npv = "NPV", ppv = "PPV"))) +
      scico::scale_fill_scico_d(palette = "lipari", begin = 0.1, end = 0.9) +
      scico::scale_color_scico_d(palette = "lipari", begin = 0, end = 0.8) +
      ylab("") + 
      scale_x_continuous(name = "Value", expand = expansion(mult = c(0, 0.05))) + 
      coord_cartesian(clip = "off") +
      theme(legend.position = "none", axis.text.y.left = element_text(face = "bold", color = "black"))
  ggsave(paste0(out_path, "_predictive_values.svg"), predictive_values_plot, width = 8.5, height = 6)
  ggsave(paste0(out_path, "_predictive_values.png"), predictive_values_plot, width = 8.5, height = 6, dpi = 600)
  
  ##### Make Machine Learning Models #####
  train <- model_input |> mutate(across(!where(is.numeric), ~factor(.x)))
  
  ### Make cross-validation folds
  cv <- vfold_cv(train, v = 5, strata = target)
  
  ### Choose performance metrics
  binary_metrics <- metric_set(accuracy, sensitivity, specificity, ppv, npv, mcc)
  continuous_metrics <- metric_set(roc_auc, pr_auc)
  all_metrics <- metric_set(roc_auc, pr_auc, sens, spec, ppv, npv, mcc)
  
  ### Build rule-fit model
  rf_model <- rand_forest() |> set_engine("ranger") |> set_mode("classification")
  rf_cv_set <- rf_model |> fit_resamples(target ~ ., resamples = cv, metrics = all_metrics, 
                                                control = control_resamples(save_pred = T, save_workflow = T, event_level = "second"))
  rf_fit <- rf_model |> fit(target ~ ., data = train)
  
  ### Plot CV performance
  gg_rf_cv_metrics <- 
    rf_cv_set |> 
    collect_metrics() |> 
    transmute(Metric = .metric, mean, upper = mean + std_err, lower = mean - std_err) |> 
    ggplot(aes(mean, Metric)) +
    geom_col(alpha = 0.5) +
    geom_segment(aes(x = lower, xend = upper, yend = Metric), linewidth = 1) + 
    geom_point(aes(x = upper), shape = "|", size = 4) +
    geom_point(aes(x = lower), shape = "|", size = 4) + 
    geom_text(aes(label = paste0("[", round(lower, 2), " - ", round(upper, 2), "]"), x = upper + 0.01), hjust = 0) +
    scale_x_continuous(name = "Value", limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) + 
    scale_y_discrete(limits = c("mcc", "npv", "ppv", "spec", "sens", "pr_auc", "roc_auc"), 
                     labels = c("MCC", "NPV", "PPV", "Spec", "Sens", "auPR", "auROC")) +
    coord_cartesian(clip = "off") + 
    ggtitle("Cross-Validation Performance: Decision-Tree Model", subtitle = "Predicting Significant Culture Growth from Gram Stain Results") +
    theme(axis.line.y.left = element_blank())
  ggsave(paste0(out_path, "_rf_cv_metrics.svg"), gg_rf_cv_metrics, width = 13, height = 5)
  ggsave(paste0(out_path, "_rf_cv_metrics.png"), gg_rf_cv_metrics, width = 13, height = 5, dpi = 600)
  
  ### Make Pred Prob Distribtuions
  if(procedure == "All"){
    gg_input <- rf_cv_set |> collect_predictions() |> left_join(train |> transmute(.row = row_number(), Site, Procedure, `Organism Abundance`, target))
  }
  
  if(procedure != "All"){
    gg_input <- rf_cv_set |> collect_predictions() |> left_join(train |> transmute(.row = row_number(), Site, `Organism Abundance`, target))
  }
  
  ### Plot ROC curves
  gg_rf_roc <- 
    gg_input |>
      group_by(id) |>
      roc_curve(truth = target, .pred_TRUE, event_level = "second") |> 
      ggplot(aes(x = 1 - specificity, y = sensitivity, group = id, color = id)) +
      geom_path() + 
      scale_color_grey() +
      ggtitle("ROC Curve: Machine Learning Model", subtitle = "Predicting Significant Culture Growth from Gram Stain Results and Clinical Data") + 
      theme(legend.position = c(0.8, 0.4))
  ggsave(paste0(out_path, "_rf_roc_curves.svg"), gg_rf_roc, width = 5, height = 5)
  ggsave(paste0(out_path, "_rf_roc_curves.png"), gg_rf_roc, width = 5, height = 5, dpi = 600)
  
  ### Plot Distributions
  gg_rf_histograms <- 
    gg_input |>
      ggplot(aes(.pred_TRUE, fill = target)) +
      geom_histogram(breaks = seq(0, 1, by = 0.05), linewidth = 1, color = "black", closed = "left") + 
      stat_bin(aes(y = after_stat(count), label = ifelse(after_stat(count) < 4, "", after_stat(count)), color = target), position = position_stack(vjust = 0.5), breaks = seq(0, 1, by = 0.05), fontface = "bold", geom="text", show.legend = F) + 
      scale_fill_manual(values = c("gray80", "#534C7A")) +
      scale_color_manual(values = c("black", "white")) +
      scale_x_continuous(name = "Predicted Probability", breaks = seq(0, 1, by = 0.25)) +
      ggtitle("Predicted Probability Distributions: Machine Learning Model", subtitle = "Predicting Significant Culture Growth from Gram Stain Results and Clinical Data") + 
      theme(axis.line.y.left = element_blank(), axis.title.y.left = element_blank(), legend.position = c("inside"), 
            legend.position.inside = c(0.75, 0.75), legend.title = element_blank())
  ggsave(paste0(out_path, "_rf_cv_histogram.svg"), gg_rf_histograms, width = 8, height = 5)
  ggsave(paste0(out_path, "_rf_cv_histogram.png"), gg_rf_histograms, width = 8, height = 5, dpi = 600)

  ### Site-Wise Performance
  gg_rf_roc_by_site <- 
    gg_input |>
      group_by(Site) |>
      roc_curve(truth = target, .pred_TRUE, event_level = "second") |> 
      ggplot(aes(x = 1 - specificity, y = sensitivity, group = Site, color = Site, linetype = fct_rev(Site))) +
      geom_path() + 
      scico::scale_color_scico_d(palette = "lipari", begin = 0.1, end = 0.9) + 
      ggtitle("ROC Curve: Machine Learning Model", subtitle = "Site-Wise Performance") + 
      theme(legend.position = c(0.8, 0.4))
  ggsave(paste0(out_path, "_rf_roc_curves_by_site.svg"), gg_rf_roc_by_site, width = 5, height = 5)
  ggsave(paste0(out_path, "_rf_roc_curves_by_site.png"), gg_rf_roc_by_site, width = 5, height = 5, dpi = 600)
  
  gg_rf_histograms_by_site <- 
    gg_input |>
      ggplot(aes(.pred_TRUE, fill = target)) +
      geom_histogram(breaks = seq(0, 1, by = 0.05), linewidth = 1, color = "black", closed = "left") + 
      stat_bin(aes(y = after_stat(count), label = ifelse(after_stat(count) < 4, "", after_stat(count)), color = target), position = position_stack(vjust = 0.5), breaks = seq(0, 1, by = 0.05), fontface = "bold", geom="text", show.legend = F) + 
      facet_wrap(~Site) + 
      scale_fill_manual(values = c("gray80", "#534C7A")) +
      scale_color_manual(values = c("black", "white")) +
      scale_x_continuous(name = "Predicted Probability", breaks = seq(0, 1, by = 0.25)) +
      ggtitle("Predicted Probability Distributions: Machine Learning Model", subtitle = "Predicting Significant Culture Growth from Gram Stain Results and Clinical Data") + 
      theme(axis.line.y.left = element_blank(), axis.title.y.left = element_blank(), legend.position = c("inside"), 
            legend.position.inside = c(0.75, 0.75), legend.title = element_blank())
  ggsave(paste0(out_path, "_rf_cv_histogram_by_site.svg"), gg_rf_histograms_by_site, width = 8, height = 5)
  ggsave(paste0(out_path, "_rf_cv_histogram_by_site.png"), gg_rf_histograms_by_site, width = 8, height = 5, dpi = 600)
  
  ### Procedure-Wise Performance
  if(procedure == "All") {
    
  gg_rf_roc_by_procedure <- 
    gg_input |>
      group_by(Procedure) |>
      roc_curve(truth = target, .pred_TRUE, event_level = "second") |> 
      ggplot(aes(x = 1 - specificity, y = sensitivity, group = Procedure, color = Procedure, linetype = fct_rev(Procedure))) +
      geom_path() + 
      scico::scale_color_scico_d(palette = "lipari", begin = 0.3, end = 0.7) + 
      ggtitle("ROC Curve: Machine Learning Model", subtitle = "Procedure-Wise Performance") + 
      theme(legend.position = c(0.8, 0.4))
  ggsave(paste0(out_path, "_rf_roc_curves_by_procedure.svg"), gg_rf_roc_by_procedure, width = 5, height = 5)
  ggsave(paste0(out_path, "_rf_roc_curves_by_procedure.png"), gg_rf_roc_by_procedure, width = 5, height = 5, dpi = 600)
  
  gg_rf_histograms_by_procedure <- 
    gg_input |>
      ggplot(aes(.pred_TRUE, fill = target)) +
      geom_histogram(breaks = seq(0, 1, by = 0.05), linewidth = 1, color = "black", closed = "left") + 
      stat_bin(aes(y = after_stat(count), label = ifelse(after_stat(count) < 4, "", after_stat(count)), color = target), position = position_stack(vjust = 0.5), breaks = seq(0, 1, by = 0.05), fontface = "bold", geom="text", show.legend = F) + 
      facet_wrap(~Procedure) + 
      scale_fill_manual(values = c("gray80", "#534C7A")) +
      scale_color_manual(values = c("black", "white")) +
      scale_x_continuous(name = "Predicted Probability", breaks = seq(0, 1, by = 0.25)) +
      ggtitle("Predicted Probability Distributions: Machine Learning Model", subtitle = "Predicting Significant Culture Growth from Gram Stain Results and Clinical Data") + 
      theme(axis.line.y.left = element_blank(), axis.title.y.left = element_blank(), legend.position = c("inside"), 
            legend.position.inside = c(0.75, 0.75), legend.title = element_blank())
  ggsave(paste0(out_path, "_rf_cv_histogram_by_procedure.svg"), gg_rf_histograms_by_procedure, width = 8, height = 5)
  ggsave(paste0(out_path, "_rf_cv_histogram_by_procedure.png"), gg_rf_histograms_by_procedure, width = 8, height = 5, dpi = 600)
  
  gg_npv_all_input <- 
    gg_input |>
      arrange(Site, Procedure, target, `Organism Abundance`, .pred_TRUE) |> 
      mutate(target = ifelse(target == "TRUE", 1, 0)) |>
      mutate(negatives = row_number(), total = n(), npv = 1 - cumsum(target) / total) |> 
      group_by(Site, Procedure) |>
      mutate(negatives_snp = row_number(), total_snp = n(), npv_snp = 1 - cumsum(target) / total_snp) |> 
      group_by(Site) |>
      mutate(negatives_s = row_number(), total_s = n(), npv_s = 1 - cumsum(target) / total_s) |>
      group_by(Procedure) |>
      mutate(negatives_p = row_number(), total_p = n(), npv_p = 1 - cumsum(target) / total_p) |> 
      ungroup() |> 
      mutate(Abundance = factor(case_when(`Organism Abundance` == "None" ~ "None", `Organism Abundance` %in% c("Rare", "Few") ~ "Rare/Few", `Organism Abundance` %in% c("Moderate", "Abundant") ~ "Moderate/Abundant", TRUE ~ "All"), levels = c("None", "Rare/Few", "Moderate/Abundant")))
 
  gg_all_npv_100 <- gg_npv_all_input |> filter(npv == 1) |> slice_tail(n = 1) |> transmute(negatives, prop = negatives/total, total, npv, Site = "ALL", Procedure = "ALL", Threshold = "1")
  gg_all_npv_99 <- gg_npv_all_input |> filter(npv > 0.985) |> slice_tail(n = 1) |> transmute(negatives, prop = negatives/total, total, npv, Site = "ALL", Procedure = "ALL", Threshold = "0.99")
  gg_all_npv_95 <- gg_npv_all_input |> filter(npv > 0.95) |> slice_tail(n = 1) |> transmute(negatives, prop = negatives/total, total, npv, Site = "ALL", Procedure = "ALL", Threshold = "0.95")
  
  gg_all_npv_100_p <- gg_npv_all_input |> group_by(Procedure) |> filter(npv_p == 1) |> slice_tail(n = 1) |> transmute(negatives = negatives_p, prop = negatives_p/total_p, total = total_p, npv = npv_p, Site = "ALL", Procedure, Threshold = "1")
  gg_all_npv_99_p <- gg_npv_all_input |> group_by(Procedure) |> filter(npv_p > 0.985) |> slice_tail(n = 1) |> transmute(negatives = negatives_p, prop = negatives_p/total_p, total = total_p, npv = npv_p, Site = "ALL", Procedure, Threshold = "0.99")
  gg_all_npv_95_p <- gg_npv_all_input |> group_by(Procedure) |> filter(npv_p > 0.95) |> slice_tail(n = 1) |> transmute(negatives = negatives_p, prop = negatives_p/total_p, total = total_p, npv = npv_p, Site = "ALL", Procedure, Threshold = "0.95")
  
  gg_all_npv_100_s <- gg_npv_all_input |> group_by(Site) |> filter(npv_s == 1) |> slice_tail(n = 1) |> transmute(negatives = negatives_s, prop = negatives_s/total_s, total = total_s, npv = npv_s, Site, Procedure = "ALL", Threshold = "1")
  gg_all_npv_99_s <- gg_npv_all_input |> group_by(Site) |> filter(npv_s > 0.985) |> slice_tail(n = 1) |> transmute(negatives = negatives_s, prop = negatives_s/total_s, total = total_s, npv = npv_s, Procedure = "ALL", Threshold = "0.99")
  gg_all_npv_95_s <- gg_npv_all_input |> group_by(Site) |> filter(npv_s > 0.95) |> slice_tail(n = 1) |> transmute(negatives = negatives_s, prop = negatives_s/total_s, total = total_s, npv = npv_s, Procedure = "ALL", Threshold = "0.95")
  
  gg_all_npv_100_snp <- gg_npv_all_input |> group_by(Site, Procedure) |> filter(npv_snp == 1) |> slice_tail(n = 1) |> transmute(negatives = negatives_snp, prop = negatives_snp/total_snp, total = total_snp, npv = npv_snp, Site, Procedure, Threshold = "1")
  gg_all_npv_99_snp <- gg_npv_all_input |> group_by(Site, Procedure) |> filter(npv_snp > 0.985) |> slice_tail(n = 1) |> transmute(negatives = negatives_snp, prop = negatives_snp/total_snp, total = total_snp, npv = npv_snp, Site, Procedure, Threshold = "0.99")
  gg_all_npv_95_snp <- gg_npv_all_input |> group_by(Site, Procedure) |> filter(npv_snp > 0.95) |> slice_tail(n = 1) |> transmute(negatives = negatives_snp, prop = negatives_snp/total_snp, total = total_snp, npv = npv_snp, Site, Procedure, Threshold = "0.95")
  
  sites_map = c("ALL", "S001", "S002", "S003", "S004")
  procedures_map = c("ALL", "ETA", "BAL")
  threshold_map = c("1", "0.99", "0.95")
  
  gg_npv_combined_input <- 
    bind_rows(gg_all_npv_100, gg_all_npv_99, gg_all_npv_95, 
              gg_all_npv_100_p, gg_all_npv_99_p, gg_all_npv_95_p, 
              gg_all_npv_100_s, gg_all_npv_99_s, gg_all_npv_95_s, 
              gg_all_npv_100_snp, gg_all_npv_99_snp, gg_all_npv_95_snp) |> 
    mutate(Site = factor(Site, levels = sites_map), 
           Procedure = factor(Procedure, levels = procedures_map), 
           Threshold = factor(Threshold, levels = threshold_map),
           false_negatives = round(negatives - npv * negatives),
           true_negatives = negatives - false_negatives,
           positives = total - negatives,
           fn_prop = false_negatives / total,
           tn_prop = true_negatives / total,
           percent = paste0(round(prop * 100), "%"),
           label = glue::glue("
                              {total} Cultures<br>
                              **{negatives} Canceled**<br>
                              ***{false_negatives} Significant***
                              ")) |> 
    pivot_longer(c("fn_prop", "tn_prop"), names_to = "Type", values_to = "Count") |>
    mutate(label = ifelse(Type == "tn_prop", label, ""),
           Type = fct_rev(factor(Type, levels = c("fn_prop", "tn_prop"))))
  
  gg_npv_by_site_and_procedure <-
    gg_npv_combined_input |>
      filter(Threshold == 0.99 & Site != "ALL" & Procedure == "ALL") |> 
      mutate(Site = str_replace_all(Site, "S00", "Site ")) |> 
      ggplot(aes(y = Count, x = Threshold, fill = fct_rev(Type))) + 
        geom_col() + 
        ggtext::geom_richtext(aes(label = label, x = Threshold), size = 6, hjust = 0, nudge_x = -0.48, vjust = 0, nudge_y = 0.01, fill = NA, label.color = NA) +
        facet_wrap(~Site, nrow = 1, strip.position = "bottom") + 
        scale_fill_manual(values = c("darkred", "grey80")) +
        scale_y_continuous(name = "Proportion of Cultures Canceled", limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        scale_x_discrete(name = "") +
        coord_cartesian(clip = "off") +
        ggtitle("Culture Cancelation Rates By Site", subtitle = "Machine Learning Approach: Set to NPV of 99%") + 
        theme(legend.position = "none", axis.title.x.bottom = element_blank(), plot.subtitle = element_text(margin = margin(0, 0, 8, 0)), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank())
  ggsave(paste0(out_path, "_cultures_canceled_by_npv_threshold.svg"), gg_npv_by_site_and_procedure, width = 9, height = 5)
  ggsave(paste0(out_path, "_cultures_canceled_by_npv_threshold.png"), gg_npv_by_site_and_procedure, dpi = 600, width = 9, height = 5)
  
  metrics_input <- 
    gg_npv_all_input |> 
    mutate(`Significant Growth` = factor(target == 1),
           `ML @ 99% NPV` = factor(.pred_TRUE > 0.04),
           `No Organisms Seen` = factor(Abundance != "None"))
  
  metrics_input |> npv(truth = `Significant Growth`, estimate = `ML @ 99% NPV`, event_level = "second")
  
  gg_diverging_bar_inputs <- 
    gg_npv_all_input |> 
      mutate(`Significant Growth` = target == 1,
             `ML @ 99% NPV` = .pred_TRUE > 0.04,
             `No Organisms Seen` = Abundance != "None") |>
      pivot_longer(c("ML @ 99% NPV", "No Organisms Seen"), names_to = "Type", values_to = "Value") |> 
      mutate(correct = case_when(!Value & !`Significant Growth` ~ "True Negatives",
                                 Value & !`Significant Growth` ~ "False Positives",
                                 !Value & `Significant Growth` ~ "False Negatives",
                                 Value & `Significant Growth` ~ "True Positives"),
             position = case_when(Type == "No Organisms Seen" ~ 1.15,
                                  Type == "ML @ 99% NPV" ~ 1.4)) |>
      count(Type, `Significant Growth`, Value, correct, position) |> 
      group_by(Type) |>
      mutate(total = sum(n)) |>
      group_by(Type, Value) |>
      mutate(pred_total = sum(n), pred_prop = round(pred_total/total * 100), subset_prop = round(n/pred_total * 100)) |>
      mutate(label = case_when(correct == "False Negatives" ~ glue("***{n} significant***"), 
                               correct == "True Negatives" ~ glue("*{total}* cultures<br>*{pred_total} canceled*")))

  gg_diverging_bar_all <- 
    ggplot() + 
      geom_col(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
               aes(x = fct_rev(Type), y = -n, fill = correct)) + 
      geom_col(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
               aes(x = fct_rev(Type), y = n, fill = correct)) + 
      ggtext::geom_richtext(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
                            aes(y = fct_rev(Type), x = 1, label = label), 
                            size = 8, hjust = 0, vjust = 2, fill = NA, label.color = NA) +
      ggtext::geom_richtext(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
                            aes(y = fct_rev(Type), x = 1, label = label), 
                            size = 8, hjust = 0, vjust = 0, fill = NA, label.color = NA) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      theme(legend.position = "none", axis.line.y.left = element_blank())
  
  gg_diverging_bar_inputs_by_procedure <- 
    gg_npv_all_input |> 
    mutate(`Significant Growth` = target == 1,
           `ML @ 99% NPV` = .pred_TRUE > 0.04,
           `No Organisms Seen` = Abundance != "None") |>
    pivot_longer(c("ML @ 99% NPV", "No Organisms Seen"), names_to = "Type", values_to = "Value") |> 
    mutate(correct = case_when(!Value & !`Significant Growth` ~ "True Negatives",
                               Value & !`Significant Growth` ~ "False Positives",
                               !Value & `Significant Growth` ~ "False Negatives",
                               Value & `Significant Growth` ~ "True Positives"),
           position = case_when(Type == "No Organisms Seen" ~ 1.15,
                                Type == "ML @ 99% NPV" ~ 1.4)) |>
    count(Type, `Significant Growth`, Value, Procedure, correct, position) |> 
    group_by(Type, Procedure) |>
    mutate(total = sum(n)) |>
    group_by(Type, Value, Procedure) |>
    mutate(pred_total = sum(n), pred_prop = round(pred_total/total * 100), subset_prop = round(n/pred_total * 100)) |>
    mutate(label = case_when(correct == "False Negatives" ~ glue("***{n} significant***"), 
                             correct == "True Negatives" ~ glue("*{total}* cultures<br>*{pred_total} canceled*")))
  
  gg_diverging_bar_all <- 
    ggplot() + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
             aes(y = fct_rev(Type), x = -n, fill = correct)) + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
             aes(y = fct_rev(Type), x = n, fill = correct)) + 
    ggtext::geom_richtext(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
                          aes(y = fct_rev(Type), x = 1, label = label), 
                          size = 8, hjust = 0, vjust = 2, fill = NA, label.color = NA) +
    ggtext::geom_richtext(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
                          aes(y = fct_rev(Type), x = 1, label = label), 
                          size = 8, hjust = 0, vjust = 0, fill = NA, label.color = NA) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    theme(legend.position = "none", axis.line.y.left = element_blank())
  
  ggplot() + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
             aes(x = position, y = -n, fill = correct), width = 0.095) + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
             aes(x = position, y = n, fill = correct), width = 0.095) + 
    scale_x_continuous(name = "", limits = c(1.1, 1.45), breaks = c(1.15, 1.3, 1.4), labels = c("No Organisms Seen", "ML @ 95% NPV", "ML @ 99% NPV")) +
    facet_wrap(~ Site, ncol = 1) + 
    coord_flip() + 
    theme(legend.position = "none", axis.line.y.left = element_blank())
  
  ggplot() + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "False Negatives"), 
             aes(x = position, y = -n, fill = correct), width = 0.095) + 
    geom_col(data = gg_diverging_bar_inputs |> filter(correct == "True Negatives"), 
             aes(x = position, y = n, fill = correct), width = 0.095) + 
    scale_x_continuous(name = "", limits = c(1.1, 1.45), breaks = c(1.15, 1.3, 1.4), labels = c("No Organisms Seen", "ML @ 95% NPV", "ML @ 99% NPV")) +
    facet_wrap(~ Procedure, ncol = 1) + 
    coord_flip() + 
    theme(legend.position = "none", axis.line.y.left = element_blank())
  
  }
  
  
  
  ### Get SHAP values from Final Fit
  library(shapviz)
  library(fastshap)
  library(glue)
  pfun <- function(object, newdata) { predict(object, new_data = newdata, type = "prob")[,2][[1]] }
  
  correct <- augment(rf_fit, train) |> filter(target == .pred_class)
  
  shap <- fastshap::explain(rf_fit, X = correct |> select(-target, -matches("pred")), pred_wrapper = pfun, shap_only = F)
  gg_shap <- 
    sv_importance(shapviz(shap), max_display = 20, kind = "beeswarm") +
    ggtitle("Feature Importance by SHAP Values", subtitle = "Correct Predictions from Training Set") +
    theme(legend.position = c(0.9, 0.4), axis.text.y.left = element_text(size = 12, face = "bold", color = "black"))
  ggsave(paste0(out_path, "_rf_shap_values.svg"), gg_shap, width = 8, height = 4)
  ggsave(paste0(out_path, "_rf_shap_values.png"), gg_shap, width = 8, height = 4, dpi = 600)
  
  list(cv = rf_cv_set, fit = rf_fit, shap = shap, train = train)

}

makeOrganismAbundanceNPVfigure <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  library(glue)
  
  input |> 
    ungroup() |>
    transmute(abundance = as.numeric(`Organism Abundance`), any_sig = factor(`Any Significant Growth`)) |> 
    arrange(abundance, any_sig) |>
    threshold_perf("abundance", truth = "any_sig", 
                   thresholds = seq(0.5, 5.5, by = 1), 
                   metrics = metric_set(npv, ppv, sensitivity, specificity), 
                   event_level = "second")
  
}

makeRachelFigures <- function(gram_input = tar_read(gram_culture_joined_by_subject) |> ungroup() |> select(-`Unique Organisms Grown`, -`Predominant Morphology`, -`No Organism Seen`), culture_input_by_organism = tar_read(culture_results_joined_by_organism), culture_input_by_subject = tar_read(culture_results_joined_by_subject)){
  
  ##### Gram Stain Results Tables #####
  gt_gram_stain_general_by_site <- 
    gram_input |> 
      rename("Predominant Morphology" = `Predominant Morphology (w/ Other)`) |>
      select(-SubjectID, -`Epithelial Cells`, -Age, -Sex, -c(`Gram-Positive Cocci`:`Other Morphology`)) |> 
      tbl_summary(by = Site, sort = list(`Predominant Morphology` = "frequency"), missing = "no") |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_site_per_specimen_table.html"))
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_site_per_specimen_table.pdf"))
  gt::gtsave(gt_gram_stain_general_by_site |> as_gt() |> tab_header(title = "Gram Stain Results", subtitle = "By Site, All Sources, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_site_per_specimen_table.docx"))
  
  gt_gram_stain_general_by_procedure <- 
    gram_input |> 
      rename("Predominant Morphology" = `Predominant Morphology (w/ Other)`) |>
      select(-SubjectID, -`Epithelial Cells`, -Age, -Sex, -c(`Gram-Positive Cocci`:`Other Morphology`)) |> 
      tbl_summary(by = Procedure, sort = list(`Predominant Morphology` = "frequency"), missing = "no") |> 
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      bold_p() |>
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_procedure_per_specimen_table.html"))
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_procedure_per_specimen_table.pdf"))
  gt::gtsave(gt_gram_stain_general_by_procedure |> as_gt() |> tab_header(title = "Gram Stain Metadata", subtitle = "By Collection Procedure, All Sites, One Row Per Subject"), 
             filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_by_procedure_per_specimen_table.docx"))
  
  gt_gram_stain_combo <- 
    tbl_merge(list(gt_gram_stain_general_by_site, gt_gram_stain_general_by_procedure), tab_spanner = c("By Site", "By Procedure")) |> as_gt() |> 
    tab_header(title = "Gram Stain Results", subtitle = "By Site and Collection Procedure")
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_combo_table.html"))
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_combo_table.png"), zoom = 6)
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_combo_table.pdf"))
  gt::gtsave(gt_gram_stain_combo, filename = paste0(figure_path, "rachel_asm_general_gram_stain_overview_combo_table.docx"))
  

  input
  
  ##### Make Contaminated Cultures Table #####
  gt_contaminated <- 
    culture_input_by_subject |> 
      tbl_summary(by = Contaminated, include = c("Site", "Procedure", "Unique Organisms: Significant")) |>
      italicize_levels() |>
      bold_labels() |> 
      add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
      add_overall() |>  
      modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")
  gt::gtsave(gt_contaminated |> as_gt() |> tab_header(title = "Contaminated Cultures", subtitle = "By Site and Collection Procedure"), 
             filename = paste0(figure_path, "/contaminated_cultures_by_site_and_procedure.html"))
  gt::gtsave(gt_contaminated |> as_gt() |> tab_header(title = "Contaminated Cultures", subtitle = "By Site and Collection Procedure"), 
             filename = paste0(figure_path, "/contaminated_cultures_by_site_and_procedure.docx"))
  
  ##### Of negative Gram stains, what grew? #####
  gram_input |> filter(`Organism Abundance` == "None" & `Any Significant Growth`) |> select(SubjectID) |> left_join(culture_input_by_organism) |> drop_na(Organism) |> filter(`Growth Significance` == "Significant")
  
  ##### Of Abundant Gram Positive Cocci, what grew? #####
  gram_input |> filter(`Organism Abundance` == "Abundant" & `Predominant Morphology (w/ Other)`== "Gram-positive cocci") |> select(SubjectID, `Any Concordant Growth`) |> count(`Any Concordant Growth`)
  
  ##### Signficant Growth by Predominant Morphology #####
  gg_significant_growth_by_predominant_morph <- 
    gram_input |> 
      mutate(
        `Predominant Morphology (w/ Other)` = ifelse(is.na(`Predominant Morphology (w/ Other)`), "None", as.character(`Predominant Morphology (w/ Other)`)), 
        Morphology = factor(`Predominant Morphology (w/ Other)`, 
                                 levels = c("Gram-positive cocci", "Yeast", "Mixed microorganism", "Gram-negative bacilli", "Other", "None"))) |>
      ggplot(aes(fill = fct_rev(Morphology), y = `Any Significant Growth`)) +
        geom_bar(position = "fill") +
        scico::scale_fill_scico_d(palette = "acton", begin = 1, end = 0) + 
        guides(fill = guide_legend(nrow = 1, reverse = T, direction = "horizontal")) +
        labs(title = "Significant Growth by Predominant Morphology", x = "Proportion of Specimens", y = "Significant Growth") +
        theme(legend.text.position = "top", legend.title = element_blank(), legend.position = c(0.75, 1.02))
    ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_morphology.png"), gg_significant_growth_by_predominant_morph, width = 12, height = 6, dpi = 600)
    ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_morphology_white.png"), gg_significant_growth_by_predominant_morph, width = 12, height = 6, dpi = 600, bg = "white")
    ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_morphology.svg"), gg_significant_growth_by_predominant_morph, width = 12, height = 6)
    
  ##### Significant Growth Count Figures #####
  gg_significant_growth_by_abundance <-
    gram_input |> 
      count(`Any Significant Growth`, `Organism Abundance`) |>
      mutate(total = sum(n), prop = round(n/sum(n), digits = 2) * 100, .by = `Organism Abundance`) |>
      filter(`Any Significant Growth`) |> 
      ggplot(aes(x = prop, y = `Organism Abundance`, alpha = `Organism Abundance`)) +
        geom_col(fill = "black") + 
        ggtext::geom_richtext(aes(label = glue("***{prop}%***<br>
                                               {n} of {total}"), x = prop), hjust = 0, fill = NA, label.color = NA, size = 4, color = "black") +
        scale_x_continuous(name = "Percent of Specimens with Significant Culture Growth", limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) + 
        scale_alpha_manual(values = seq(0.5, 1, by = 0.1)) +
        annotate("text", x = 95, y = 1, label = glue("p < 0.001, by Pearson's Chi-squared"), hjust = 1, vjust = 0, size = 4) +
        ggtitle("Significant Culture Growth by Organism Abundance", subtitle = "All Sites, All Procedures") +
        coord_cartesian(clip = "off") +
        theme(legend.position = "none", axis.title.y.left = element_blank(), axis.text.y.left = element_text(face = "bold", color = "black"))
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_all_procedures.png"), gg_significant_growth_by_abundance, width = 12, height = 5, dpi = 600)
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_all_procedures_white.png"), gg_significant_growth_by_abundance, width = 12, height = 5, dpi = 600, bg = "white")
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_all_procedures.svg"), gg_significant_growth_by_abundance, width = 12, height = 5)
  
  gg_significant_growth_by_abundance_per_procedure <-   
    gram_input |> 
      count(`Any Significant Growth`, `Organism Abundance`, Procedure) |>
      mutate(total = sum(n), prop = round(n/sum(n), digits = 2) * 100, .by = c(`Organism Abundance`, Procedure)) |>
      filter(`Any Significant Growth`) |> 
      ggplot(aes(x = prop, y = `Organism Abundance`, alpha = `Organism Abundance`, fill = Procedure)) +
        geom_col() + 
        ggtext::geom_richtext(aes(label = glue("***{prop}%***<br>{n}/{total}"), x = prop), hjust = 1, fill = NA, label.color = NA, size = 3, color = "black", alpha = 1) +
        scale_x_continuous(name = "Percent of Specimens with Significant Culture Growth", limits = c(0, 100), expand = expansion(mult = c(0, 0.01))) + 
        scale_alpha_manual(values = seq(0.5, 1, by = 0.1)) +
        scico::scale_fill_scico_d(palette = "lipari", begin = 0.25, end = 0.75) +
        facet_wrap(~Procedure, nrow = 2) + 
        annotate("text", x = 95, y = 1, label = glue("p < 0.001, by Pearson's Chi-squared"), hjust = 1, vjust = 0, size = 4) +
        ggtitle("Significant Culture Growth by Organism Abundance", subtitle = "All Sites, By Procedure") +
        coord_cartesian(clip = "off") +
        theme(legend.position = "none", axis.title.y.left = element_blank(), axis.text.y.left = element_text(face = "bold", color = "black"))
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_by_procedures.png"), gg_significant_growth_by_abundance_per_procedure, width = 12, height = 5, dpi = 600)
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_by_procedures_white.png"), gg_significant_growth_by_abundance_per_procedure, width = 12, height = 5, dpi = 600, bg = "white")
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_all_sites_by_procedures.svg"), gg_significant_growth_by_abundance_per_procedure, width = 12, height = 5)
  
  gg_significant_growth_by_abundance_per_site <-   
    gram_input |> 
      count(`Any Significant Growth`, `Organism Abundance`, Site) |>
      mutate(total = sum(n), prop = round(n/sum(n), digits = 2) * 100, .by = c(`Organism Abundance`, Site)) |>
      filter(`Any Significant Growth`) |> 
      ggplot(aes(x = prop, y = `Organism Abundance`, alpha = `Organism Abundance`, fill = Site)) +
      geom_col() + 
      ggtext::geom_richtext(aes(label = glue("***{prop}%***<br>{n} of {total}"), x = prop), hjust = 1, fill = NA, label.color = NA, size = 3, color = "black", alpha = 1) +
      scale_x_continuous(name = "Percent of Specimens with Significant Culture Growth", limits = c(0, 100), expand = expansion(mult = c(0.02, 0.02))) + 
      scale_alpha_manual(values = seq(0.5, 1, by = 0.1)) +
      scico::scale_fill_scico_d(palette = "lipari", begin = 0.1, end = 0.9) +
      facet_wrap(~Site, nrow = 2) + 
      annotate("text", x = 95, y = 1, label = glue("p < 0.001"), hjust = 1, vjust = 0, size = 4) +
      ggtitle("Significant Culture Growth by Organism Abundance", subtitle = "By Site, All Procedures") +
      coord_cartesian(clip = "off") +
      theme(legend.position = "none", axis.title.y.left = element_blank(), axis.text.y.left = element_text(face = "bold", color = "black"))
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_by_sites_all_procedures.png"), gg_significant_growth_by_abundance_per_site, width = 12, height = 5, dpi = 600)
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_by_sites_all_procedures_white.png"), gg_significant_growth_by_abundance_per_site, width = 12, height = 5, dpi = 600, bg = "white")
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_by_sites_all_procedures.svg"), gg_significant_growth_by_abundance_per_site, width = 12, height = 5)
  
  gg_significant_growth_by_morphology <-   
    gram_input |> 
      count(`Any Significant Growth`, `Organism Abundance`, `Predominant Morphology (w/ Other)`) |>
      mutate(total = sum(n), prop = round(n/sum(n), digits = 2) * 100, .by = c(`Organism Abundance`, `Predominant Morphology (w/ Other)`)) |>
      filter(`Any Significant Growth` & `Predominant Morphology (w/ Other)` %in% c("Gram-positive cocci", "Gram-negative bacilli", "Mixed microorganism", "Other")) |> 
      mutate(`Predominant Morphology (w/ Other)` = factor(`Predominant Morphology (w/ Other)`, levels = c("Gram-positive cocci", "Mixed microorganism", "Gram-negative bacilli", "Other"))) |> 
      ggplot(aes(x = prop, y = `Organism Abundance`, alpha = `Organism Abundance`, fill = `Predominant Morphology (w/ Other)`)) +
      geom_col() + 
      ggtext::geom_richtext(aes(label = glue("***{prop}%***<br>{n} of {total}"), x = prop), hjust = 1, fill = NA, label.color = NA, size = 3, color = "black", alpha = 1) +
      scale_x_continuous(name = "Percent of Specimens with Significant Culture Growth", limits = c(0, 100), expand = expansion(mult = c(0.02, 0.02))) + 
      scale_alpha_manual(values = seq(0.5, 1, by = 0.1)) +
      scico::scale_fill_scico_d(palette = "acton", begin = 0.1, end = 0.9) +
      facet_wrap(~`Predominant Morphology (w/ Other)`, nrow = 2) + 
      annotate("text", x = 95, y = 1, label = glue("p < 0.001"), hjust = 1, vjust = 0, size = 4) +
      ggtitle("Significant Culture Growth by Organism Abundance", subtitle = "All Sites, All Procedures, By Gram Stain Morphology") +
      coord_cartesian(clip = "off") +
      theme(legend.position = "none", axis.title.y.left = element_blank(), axis.text.y.left = element_text(face = "bold", color = "black"))
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_per_morphology.png"), gg_significant_growth_by_morphology, width = 12, height = 5, dpi = 600)
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_per_morphology_white.png"), gg_significant_growth_by_morphology, width = 12, height = 5, dpi = 600, bg = "white")
  ggsave(paste0(figure_path, "rachel_asm_significant_growth_by_gram_abundance_per_morphology.svg"), gg_significant_growth_by_morphology, width = 12, height = 5)
  
  ##### Organism Abundance for Cultures with Staph and Pseudomonas #####
  gt_staph_pseudo_significance_by_abundance <- 
    culture_input_by_organism |> 
      left_join(gram_input |> select(SubjectID, `Organism Abundance`)) |> 
      filter(Organism %in% c("Staphylococcus aureus", "Pseudomonas aeruginosa")) |>
      mutate(`Growth Significance` = fct_drop(`Growth Significance`)) |>
      tbl_strata(strata = "Organism", 
                 .tbl_fun = ~.x |> 
                   tbl_summary(by = "Growth Significance", include = c("Organism Abundance")) |> 
                   italicize_levels() |>
                   bold_labels() |> 
                   add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
                   bold_p() |>
                   modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*"))
  gt::gtsave(gt_staph_pseudo_significance_by_abundance |> as_gt() |> tab_header(title = "Significant Culture Growth By Gram Stain Abundance", subtitle = "Top Two Organisms"), 
             filename = paste0(figure_path, "/rachel_asm_staph_pseudo_significance_by_abundance.html"))
  gt::gtsave(gt_staph_pseudo_significance_by_abundance |> as_gt() |> tab_header(title = "Significant Culture Growth By Gram Stain Abundance", subtitle = "Top Two Organisms"), 
             filename = paste0(figure_path, "/rachel_asm_staph_pseudo_significance_by_abundance.docx"))
  
  
  gt_staph_pseudo_significance_by_abundance_regression <- 
    culture_input_by_organism |> 
      left_join(gram_input |> select(SubjectID, `Organism Abundance`)) |> 
      filter(Organism %in% c("Staphylococcus aureus", "Pseudomonas aeruginosa")) |>
      mutate(`Growth Significance` = fct_drop(`Growth Significance`),
             `Organism Abundance` = as.numeric(`Organism Abundance`)) |>
      select(Organism, `Organism Abundance`, `Growth Significance`) |>
      tbl_strata(strata = "Organism", 
                 .tbl_fun = ~.x |> 
                   tbl_uvregression(y = `Growth Significance`, method = glm, exponentiate = T, method.args = list(family = binomial)) |> 
                   italicize_levels() |>
                   bold_labels() |> 
                   bold_p())
  gt::gtsave(gt_staph_pseudo_significance_by_abundance_regression |> as_gt() |> tab_header(title = "Significant Culture Growth By Gram Stain Abundance", subtitle = "Top Two Organisms"), 
             filename = paste0(figure_path, "/rachel_asm_staph_pseudo_significance_by_abundance_regression.html"))
  gt::gtsave(gt_staph_pseudo_significance_by_abundance_regression |> as_gt() |> tab_header(title = "Significant Culture Growth By Gram Stain Abundance", subtitle = "Top Two Organisms"), 
             filename = paste0(figure_path, "/rachel_asm_staph_pseudo_significance_by_abundance_regression.docx"))
  
  ##### Culture Results from Negative Gram Stains ##### 
  gram_input |> 
    filter(`Organism Abundance` == "None") |> 
    left_join(culture_input_by_organism) |>
    filter(`Growth Significance` == "Significant") |>
    count(Organism, sort = T)
  
  ##### Organism Growth Rate By Procedure #####
  culture_input_by_subject |>
    tbl_summary(by = "Procedure", include = matches("Any|Significant")) |> 
    bold_labels() |> 
    italicize_levels() |> 
    add_p(test = everything() ~ "t.test", test.args = all_tests("t.test") ~ list(var.equal = T)) |> 
    bold_p()

}

makeCareyAnnFigures <- function(input = tar_read(gram_culture_joined_by_subject)){
  
  culture_raw <- read_csv("../../../data/culture_results_input.csv")
  
  culture_nurfi <- 
    culture_raw |> 
      select(SubjectID, LRT_REF_GRWTH_YN, LRT_REF_NF_AMT) |> 
      distinct()
  
  input_with_NURFI <- left_join(input, culture_nurfi)
  
  any_organism_on_gram <- 
    input_with_NURFI |> 
      ungroup() |> 
      mutate(growth_class = 
               case_when(
                 `Any Significant Growth` ~ "Pathogen",
                 LRT_REF_GRWTH_YN == "No" ~ "Negative",
                 LRT_REF_GRWTH_YN == "Yes" & !`Any Significant Growth` &
                   ((LRT_REF_NF_AMT %in% c(NA, "<1,000 CFU/mL", "1,000-10,000 CFU/mL", "10,000-100,000 CFU/mL") & Procedure == "ETA") | 
                   (LRT_REF_NF_AMT %in% c(NA, "<1,000 CFU/mL", "1,000-10,000 CFU/mL") & Procedure == "BAL")) ~ "Insignificant",
                 LRT_REF_GRWTH_YN == "Yes" & !`Any Significant Growth` &
                   ((LRT_REF_NF_AMT %in% c(">100,000 CFU/mL") & Procedure == "ETA") | 
                   (LRT_REF_NF_AMT %in% c("10,000-100,000 CFU/mL", ">100,000 CFU/mL") & Procedure == "BAL")) ~ "NURFI"))
  
  any_organism_on_gram |> 
    filter(!`No Organism Seen`) |> 
    select(Procedure, growth_class) |>
    mutate(`Culture Result` = fct_relevel(growth_class, "Negative", "Insignificant", "NURFI", "Pathogen")) |>
    drop_na() |> 
    tbl_summary(by = "Procedure", include = "Culture Result") |> 
      bold_labels() |>
      italicize_levels() |>
      add_p() |>
      bold_p() |> 
      as_gt() |> 
      tab_header(title = "Culture Summaries", subtitle = "Gram stains with any organism present")
  
  library(ggalluvial)
  any_organism_on_gram |> 
    transmute(Procedure, 
              `Organisms on Gram stain` = factor(ifelse(`No Organism Seen`, "No Organisms", "Organisms Present")), 
              `Culture Result` = factor(ifelse(growth_class == "Pathogen", "Pathogen Growth", "No/Insignificant Growth"), levels = c("No/Insignificant Growth", "Pathogen Growth")),
              fill_color = (`Culture Result` == "Pathogen Growth")) |> 
    drop_na() |> 
    count(Procedure, `Organisms on Gram stain`, `Culture Result`, fill_color, sort = T) |> 
      ggplot(aes(y = n, axis1 = `Organisms on Gram stain`, axis2 = `Culture Result`)) + 
        geom_alluvium(aes(fill = fill_color, color = fill_color), width = 2/3, alpha = 1) +
        geom_stratum(aes(fill = fill_color, color = fill_color), width = 2/3, fill = "grey75", color = "black") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
        facet_wrap(~Procedure, nrow = 2) + 
        scale_fill_manual(values = c("grey90", "#191841")) + 
        scale_color_manual(values = c("grey80", "#191841")) + 
        theme(legend.position = "none", axis.line.y.left = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x.bottom = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave("~/Downloads/signficant_growth_by_gram_sankey.svg", width = 8, height = 4)
}


