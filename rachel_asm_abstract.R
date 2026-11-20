##### Config ######
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(ggsankey))
extrafont::loadfonts(quiet = TRUE)
sysfonts::font_add_google("Comfortaa", "Comfortaa")
sysfonts::font_add_google("Nunito sans", "Nunito sans")
showtext::showtext_auto()

theme_ns <- theme(text = element_text(family = "Comfortaa"),
                  title = element_text(size = 14, face = "bold", margin = margin(0, 0, 8, 0)),
                  plot.subtitle = element_text(size = 12, face = "plain", hjust = 0, margin = margin(0, 0, 8, 0)),
                  axis.title = element_text(size = 12, face = "bold", margin = margin(4,4,4,4)),
                  axis.text = element_text(family = "Nunito sans"),
                  axis.title.x.bottom = element_text(face = "bold", margin = margin(4,0,0,0)),
                  axis.title.y.left = element_text(face = "bold", margin = margin(0,4,0,0)),
                  legend.title = element_text(face = "bold.italic", size = 12),
                  axis.line = element_line(),
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(), 
                  panel.background = element_blank(),
                  strip.text = element_text(size = 10, face = "bold.italic"),
                  strip.background = element_blank())
theme_set(theme_ns)

pattern_colors <- c(
  "#191841",
  "#5A589E",
  "#99CCFF",
  "#00B28E",
  "#FF5B5A",
  "#FCB952",
  "#58595B"
)

scale_color_pattern <- scale_color_manual(values = pattern_colors)
scale_fill_pattern <- scale_fill_manual(values = pattern_colors)

geoms_color <- c("point", "line", "text", "tile", "ribbon", "segment", "path", "step", "rect", "smooth", "errorbar", "crossbar", "linerange", "pointrange", "boxplot", "violin", "hline", "vline", "abline")
geoms_fill <- c("bar", "tile", "ribbon", "area", "rect", "polygon", "density_2d", "crossbar", "boxplot", "violin", "hex")

# Update geom defaults for color and fill
for (geom in geoms_color) {
  update_geom_defaults(geom, list(colour = pattern_colors[1]))
}

for (geom in geoms_fill) {
  update_geom_defaults(geom, list(fill = pattern_colors[1]))
}

theme_pattern <- function(data, ...){
  data %>% 
    tab_options(
      heading.background.color = pattern_colors[1],
      column_labels.background.color = "#DADCEE",
      footnotes.background.color = "#F7F7FB",
      ...
    ) |> 
    tab_style(
      style = cell_text(
        font = "Nunito sans"
      ), 
      locations = cells_body()
    ) |>
    tab_style(
      style = cell_text(
        font = "Nunito sans"
      ), 
      locations = cells_footnotes()
    ) |>
    tab_style(
      style = cell_text(
        font = "Comfortaa"
      ), 
      locations = cells_column_labels()
    ) |>
    tab_style(
      style = cell_text(
        font = "Comfortaa",
        weight = "bold"
      ), 
      locations = cells_title(groups = "title")
    ) |> 
    tab_style(
      style = cell_text(
        font = "Comfortaa",
        style = "italic"
      ), 
      locations = cells_title(groups = "subtitle")
    )
  
}

###### Read in Data ######
data_path = "~/git/Pattern/data/MedrioExport202312011541/"
figure_path = "../../../figures/multicenter_study/"

enrollment_data <- readxl::read_xlsx(paste0(data_path, "LRT Enrollment_Collection_Site.xlsx"), na = c("NA", "<NA>", "N/A", "n/a"))
collection_data <- readxl::read_xlsx(paste0(data_path, "LRT_Enrollment_Collection and culture set up specified.xlsx"), na = c("NA", "<NA>", "N/A", "n/a"))
culture_data <- 
  readxl::read_xlsx(paste0(data_path, "Reference_LRT_Culture_Results.xlsx"), na = c("NA", "<NA>", "N/A", "n/a")) |> 
  mutate(LRT_REF_RPT_BIN = factor(LRT_REF_RPT_BIN, levels = c("<1,000", "1,000-10,000", "10,000-100,000", ">100,000"), ordered = T))
repeat_culture_data <- readxl::read_xlsx(paste0(data_path, "Reference_LRT_Culture_Results_Repeat.xlsx") , na = c("NA", "<NA>", "N/A", "n/a"))
gram_stain_data <- 
  readxl::read_xlsx(paste0(data_path, "Reference_LRT_Gram Stain.xlsx"), na = c("NA", "<NA>", "N/A", "n/a")) |> 
  mutate(
    across(matches("PLS_AMT|ORG_AMT"), ~factor(., levels = c("None", "Rare", "Few", "Moderate", "Abundant"), ordered = T)),
    across(matches("EPI|CELL|RBC"), ~factor(., levels = c("None", "Few", "Abundant"), ordered = T)),
    LRT_REF_GS_MORPH = str_replace_all(pattern = " in clusters| in clusters | in pairs and chains | in pairs and chains| in pairs| in pairs | and pseudohyphae", replacement = "", LRT_REF_GS_MORPH),
    LRT_REF_GS_MORPH = fct_infreq(factor(LRT_REF_GS_MORPH)))
repeat_gram_stain_data <- readxl::read_xlsx(paste0(data_path, "Reference_LRT_Gram Stain_Repeat.xlsx"), na = c("NA", "<NA>", "N/A", "n/a"))

group_gram_morph_map <- read_delim(paste0(data_path, "../group_gram_morph_map.txt")) |> mutate(across(everything(), ~str_trim(.))) |> distinct() |> 
  mutate(Morphology = str_to_sentence(Morphology),
         gram = case_when(
           grepl("negative", `Gram stain`) ~ "Negative",
           grepl("positive", `Gram stain`) ~ "Positive",
           grepl("Acid", `Gram stain`, ignore.case = T) ~ "Acid-Fast",
           T ~ "Other")) |>
  distinct()

genus_group_map <- read_delim(paste0(data_path, "../genus_species_attribute_map.tsv")) |> mutate(across(everything(), ~str_trim(.))) |> distinct()

organism_gram_group_map <- left_join(genus_group_map, group_gram_morph_map)

###### Analysis ######

gram_stain_general_input <- 
  gram_stain_data |> 
  rowwise() |>
  mutate(
    TypeID = str_split(SubjectID, "-", simplify = T)[1],
    SiteID = str_split(SubjectID, "-", simplify = T)[2],
    SpecimenID = str_split(SubjectID, "-", simplify = T)[3],
    ReplicateID = str_split(SubjectID, "-", simplify = T)[4],
  ) |> ungroup() |> 
  group_by(SubjectID) |> 
  fill(everything(), .direction = "downup") |> 
  distinct()

culture_gram_join <- 
  full_join(
    gram_stain_results_input |> select(SubjectID, LRT_REF_GS_MORPH, LRT_REF_GS_SPECTYPE, matches("AMT"), Gram, Shape),
    culture_data_input |> select(SubjectID, LRT_REF_RPT_BIN, MALDI_ID, MALDI_ID_COLLAPSED, NUM_UNIQUE_ORGANISMS, max_abundance, culture_tag = `Gram stain`, culture_morph = Morphology, culture_gram = gram)) |>
  distinct() |> 
  mutate(gram_equals_culture = (LRT_REF_GS_MORPH == culture_tag) | (LRT_REF_GS_MORPH == "Mixed microorganism" & NUM_UNIQUE_ORGANISMS > 1))

gram_matches_culture <- 
  culture_gram_join |> 
  group_by(SubjectID, LRT_REF_GS_MORPH) |> 
  summarise(
    any_gram = any(!is.na(LRT_REF_GS_MORPH), na.rm = T), 
    any_culture = any(!is.na(culture_tag), na.rm = T), 
    gram_equals_culture,
    any_gram_equals_culture = any(gram_equals_culture, na.rm = T)) |>
  ungroup() |> 
  distinct()

gram_stain_concordance_by_subject <- 
  left_join(distinct(gram_stain_general_input |> select(SubjectID, LRT_REF_GS_SPECTYPE)), gram_matches_culture) |>
  group_by(SubjectID, LRT_REF_GS_SPECTYPE) |>
  summarise(no_gram = all(is.na(LRT_REF_GS_MORPH)),
            any_gram = !all(is.na(LRT_REF_GS_MORPH)), 
            any_gram_equals_culture = any(any_gram_equals_culture, na.rm = T)) |> 
  mutate(any_gram_equals_culture = ifelse(!any_gram, "No Organism Detected", as.character(any_gram_equals_culture)), 
         any_gram_equals_culture = ifelse(any_gram_equals_culture == "TRUE", "Concordant Culture Growth", any_gram_equals_culture), 
         any_gram_equals_culture = ifelse(any_gram_equals_culture == "FALSE", "No Concordant Culture Growth", any_gram_equals_culture)) |> ungroup()

gt_gram_stain_concordance_by_subject <- 
  gram_stain_concordance_by_subject |> 
  tbl_summary(by = LRT_REF_GS_SPECTYPE, include = c(no_gram, any_gram, any_gram_equals_culture), missing_text = "No Organism Detected") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
  bold_p() |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Gram Stain Concordance with Culture", subtitle = "All Sites, By Procedure Type")

gt::gtsave(gt_gram_stain_concordance_by_subject, filename = paste0(figure_path, "gram_stain_concordance_by_subject_by_procedure_table.html"))
gt::gtsave(gt_gram_stain_concordance_by_subject, filename = paste0(figure_path, "gram_stain_concordance_by_subject_by_procedure_table.docx"))







gram_results_by_subject <- 
  gram_stain_data |>
  group_by(SubjectID) |> 
  summarise(
    max_epi = max(LRT_REF_SQ_EPI_AMT, na.rm = T),
    max_poly = max(LRT_REF_PLS_AMT, na.rm = T),
    max_rbc = max(LRT_REF_RBC_AMT, na.rm = T),
    gram_NONE = all(is.na(LRT_REF_GS_MORPH), na.rm = T),
    gram_GPC = any(LRT_REF_GS_MORPH == "Gram-positive cocci", na.rm = T),
    gram_GNB = any(LRT_REF_GS_MORPH == "Gram-negative bacilli", na.rm = T),
    gram_MIXED = any(LRT_REF_GS_MORPH == "Mixed microorganism", na.rm = T),
    gram_GPB = any(LRT_REF_GS_MORPH == "Gram-positive bacilli", na.rm = T),
    gram_OTHER = ) |> 
  ungroup() |> 
  left_join(distinct(gram_stain_data |> distinct(SubjectID, LRT_REF_GS_SPECTYPE) |> drop_na()))

gt_gram_by_subject <- 
  tbl_summary(gram_results_by_subject |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, label = c("gram_GPC" ~ "Gram-Positive Cocci", "gram_GNB" ~ "Gram-Negative Bacilli", "gram_MIXED" ~ "Mixed Microorganism", "gram_GPB" ~ "Gram-Positive Bacilli", "max_poly" ~ "Maximum Neutrophils", "max_epi" ~ "Maximum Squamous Epithelial Cells", "max_rbc" ~ "Maximum Red Blood Cells"), statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing = "ifany", missing_text = "None") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p() |>  
  bold_p() |> 
  add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]")) |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
  modify_footnote(everything() ~ NA, abbreviation = T) |> 
  modify_footnote(update = everything() ~ NA) |>
  modify_footnote(update = everything() ~ "N, [95% CI]") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Gram Stain Results", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_gram_by_subject, filename = paste0(figure_path, "gram_stain_results_by_subject_table.html"))
gt::gtsave(gt_gram_by_subject, filename = paste0(figure_path, "gram_stain_results_by_subject_table.docx"))










culture_results_by_subject <- 
  culture_data_input |> 
    left_join(gram_results_by_subject |> select(SubjectID, LRT_REF_GS_SPECTYPE)) |> 
    group_by(SubjectID, LRT_REF_GS_SPECTYPE) |> 
    summarise(
      num_org = max(NUM_UNIQUE_ORGANISMS),
      max_normal_flora = max(LRT_REF_NF_AMT, na.rm = T),
      any_contam = any(REF_LRT_CONTAM_YN == "Yes", na.rm = T),
      any_growth = any(LRT_REF_GRWTH_YN == "Yes", na.rm = T),
      pseudomonas = any(MALDI_ID_GROUPED == "Pseudomonas aeruginosa", na.rm = T),
      staph = any(MALDI_ID_GROUPED == "Staphylococcus aureus", na.rm = T),
      yeast = any(MALDI_ID_GROUPED == "Yeast", na.rm = T),
      klebsiella = any(MALDI_ID_GROUPED == "Klebsiella pneumoniae", na.rm = T),
      steno = any(MALDI_ID_GROUPED == "Stenotrophomonas maltophilia", na.rm = T),
      serratia = any(MALDI_ID_GROUPED == "Serratia marcescens", na.rm = T),
      ecoli = any(MALDI_ID_GROUPED == "Escherichia coli", na.rm = T),
      acinetobacter = any(MALDI_ID_GROUPED == "Acinetobacter baumannii-calcoaceticus complex", na.rm = T),
      hflu = any(MALDI_ID_GROUPED == "Haemophilus influenzae", na.rm = T),
      culture_GPC = any(`Gram stain` == "Gram-positive cocci", na.rm = T),
      culture_GNB = any(`Gram stain` == "Gram-negative bacilli", na.rm = T),
      culture_MIXED = any(grepl("positive", `Gram stain`), na.rm = T) & any(grepl("negative", `Gram stain`), na.rm = T),
      culture_YEAST = any(`Gram stain` == "Yeast", na.rm = T),
      culture_GPB = any(`Gram stain` == "Gram-positive bacilli", na.rm = T)) |> 
    distinct() |> 
    ungroup()

gt_culture_by_subject <- 
  tbl_summary(culture_results_by_subject |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, 
              label = c("num_org" ~ "Number of Unique Organisms", "max_normal_flora" ~ "Abundance of Normal Flora", "any_contam" ~ "Contaminated?", "any_growth" ~ "Any Growth?", "culture_GPC" ~ "Gram-Positive Cocci", "culture_GNB" ~ "Gram-Negative Bacilli", "culture_MIXED" ~ "Mixed Microorganism", "culture_GPB" ~ "Gram-Positive Bacilli", "culture_YEAST" ~ "Yeast"), 
              statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing_text = "N/A") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
  bold_p() |> 
  add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]", all_continuous() ~ "[{conf.low} - {conf.high}]")) |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
  modify_footnote(everything() ~ NA, abbreviation = T) |> 
  modify_footnote(update = everything() ~ NA) |>
  modify_footnote(update = everything() ~ "N, (Percentage). [95% Confidence Interval]") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Culture Results", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_culture_by_subject, filename = paste0(figure_path, "culture_results_by_subject_table.html"))
gt::gtsave(gt_culture_by_subject, filename = paste0(figure_path, "culture_results_by_subject_table.docx"))

gt_culture_by_subject