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
    LRT_REF_GS_MORPH = fct_infreq(factor(LRT_REF_GS_MORPH))) |> 
    filter(!SubjectID %in% c("PR-S003-033-R0", "PR-S003-034-R0"))
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

gram_stain_results_input <- 
  gram_stain_data |> 
  select(SubjectID, matches("LRT"), VarGroup1Row) |> 
  group_by(SubjectID) |> 
  fill(everything(), .direction = "downup") |> 
  rowwise() |>
  mutate(
    TypeID = str_split(SubjectID, "-", simplify = T)[1],
    SiteID = str_split(SubjectID, "-", simplify = T)[2],
    SpecimenID = str_split(SubjectID, "-", simplify = T)[3],
    ReplicateID = str_split(SubjectID, "-", simplify = T)[4],
  ) |> ungroup() |> 
  distinct() |> 
  mutate(
    organism_detected = ifelse(LRT_REF_GS_ORG_AMT == "None", F, T), 
    Kingdom = 
      case_when(
        grepl("Gram", LRT_REF_GS_MORPH, ignore.case = T) ~ "Bacteria", 
        grepl("yeast|fung", LRT_REF_GS_MORPH, ignore.case = T) ~ "Fungus",
        grepl("mixed", LRT_REF_GS_MORPH, ignore.case = T) ~ "Mixed",
        T ~ "No Organism"),
    Gram = 
      case_when(
        grepl("pos", LRT_REF_GS_MORPH, ignore.case = T) ~ "Positive", 
        grepl("neg", LRT_REF_GS_MORPH, ignore.case = T) ~ "Negative",
        grepl("var", LRT_REF_GS_MORPH, ignore.case = T) ~ "Variable",
        T ~ NA), 
    Shape = 
      case_when(
        grepl(" cocci", LRT_REF_GS_MORPH, ignore.case = T) ~ "Cocci",
        grepl(" bacilli", LRT_REF_GS_MORPH, ignore.case = T) ~ "Bacilli",
        grepl("coccobacil|diplo", LRT_REF_GS_MORPH, ignore.case = T) ~ "Other"
      ), 
    Grouping = 
      case_when(
        grepl("clusters", LRT_REF_GS_MORPH, ignore.case = T) ~ "Clusters",
        grepl("chain", LRT_REF_GS_MORPH, ignore.case = T) ~ "Pairs/Chains",
      )
  )


culture_data_input <- 
  culture_data |> 
  select(SubjectID, matches("LRT"), VarGroup1Row) |> 
  select(-matches("HAEM|OPTO|BILE|PIC|_WS")) |>
  group_by(SubjectID) |> 
  fill(everything(), .direction = "downup") |> 
  rowwise() |>
  mutate(
    TypeID = str_split(SubjectID, "-", simplify = T)[1],
    SiteID = str_split(SubjectID, "-", simplify = T)[2],
    SpecimenID = str_split(SubjectID, "-", simplify = T)[3],
    ReplicateID = str_split(SubjectID, "-", simplify = T)[4],
  ) |> 
  mutate(
    MALDI_ID = ifelse(grepl("Other", LRT_REF_MALDI_ID, ignore.case = T), str_to_sentence(LRT_REF_ID_OTHER), LRT_REF_MALDI_ID)
  ) |> 
  group_by(SubjectID) |> 
  mutate(NUM_UNIQUE_ORGANISMS = length(do::unique_no.NA(MALDI_ID))) |> 
  ungroup() |>
  group_by(SubjectID, LRT_REF_MALDI_ID) |> 
  mutate(max_abundance = max(LRT_REF_RPT_BIN)) |> 
  left_join(organism_gram_group_map, by = join_by(MALDI_ID == `Reference testing identification `)) |> 
  distinct(SubjectID, MALDI_ID, max_abundance, .keep_all = T) |>
  ungroup() |> 
  mutate(MALDI_ID_GROUPED = factor(ifelse(is.na(`Group or alternative name`), MALDI_ID, `Group or alternative name`)), 
         MALDI_ID_COLLAPSED = fct_infreq(fct_lump_min(MALDI_ID_GROUPED, min = 5)))

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
  group_by(SubjectID) |>
  mutate(num_unique_tags = length(do::unique_no.NA(culture_tag))) |>
  ungroup() |>
  mutate(gram_equals_culture = (LRT_REF_GS_MORPH == culture_tag) | (LRT_REF_GS_MORPH == "Mixed microorganism" & num_unique_tags > 1))

gram_matches_culture <- 
  culture_gram_join |> 
  group_by(SubjectID, LRT_REF_GS_MORPH) |> 
  summarise(
    any_gram = any(!is.na(LRT_REF_GS_MORPH), na.rm = T), 
    any_culture = any(!is.na(culture_tag), na.rm = T), 
    culture_tag,
    gram_equals_culture,
    any_gram_equals_culture = any(gram_equals_culture, na.rm = T)) |>
  ungroup() |> 
  distinct()

gram_stain_concordance_by_subject <- 
  left_join(culture_data_input, gram_stain_data |> select(SubjectID, LRT_REF_GS_SPECTYPE) |> drop_na() |> distinct()) |>
  left_join(gram_matches_culture) |>
  group_by(SubjectID, LRT_REF_GS_SPECTYPE) |>
  summarise(no_gram = all(is.na(LRT_REF_GS_MORPH)),
            any_gram = !all(is.na(LRT_REF_GS_MORPH)), 
            any_culture = any(!is.na(culture_tag), na.rm = T),
            any_gram_equals_culture = any(any_gram_equals_culture, na.rm = T)) |> 
  mutate(any_gram_equals_culture = ifelse(!any_gram & !any_culture, "No organsims on Gram stain, without culture growth", as.character(any_gram_equals_culture)),
         any_gram_equals_culture = ifelse(!any_gram & any_culture, "No organsims on Gram stain, with culture growth", as.character(any_gram_equals_culture)),
         any_gram_equals_culture = ifelse(any_gram_equals_culture == "TRUE", "Organsims seen on Gram stain, with concordant culture growth", any_gram_equals_culture), 
         any_gram_equals_culture = ifelse(any_gram_equals_culture == "FALSE", "Organsims seen on Gram stain, without concordant culture growth", any_gram_equals_culture)) |> 
  mutate(any_gram_equals_culture = factor(any_gram_equals_culture, levels = c("No organsims on Gram stain, without culture growth", "Organsims seen on Gram stain, with concordant culture growth", "No organsims on Gram stain, with culture growth", "Organsims seen on Gram stain, without concordant culture growth"))) |>
  ungroup()

gt_gram_stain_concordance_by_subject <- 
  gram_stain_concordance_by_subject |> 
  tbl_summary(by = LRT_REF_GS_SPECTYPE, include = c(no_gram, any_gram, any_gram_equals_culture), missing_text = "No Organism Detected") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p() |>  
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
    gram_YEAST = any(LRT_REF_GS_MORPH == "Yeast", na.rm = T),
    gram_GPB = any(LRT_REF_GS_MORPH == "Gram-positive bacilli", na.rm = T),
    gram_OTHER = any(LRT_REF_GS_MORPH %in% c("Gram-negative coccobacilli", "Gram-negative diplococci", "Gram-variable bacilli", "Fungal elements", "Gram-positive coccobacilli"), na.rm = T)) |> 
  ungroup() |> 
  left_join(distinct(gram_stain_data |> distinct(SubjectID, LRT_REF_GS_SPECTYPE) |> drop_na()))

sums <- rowSums(gram_results_by_subject |> select(matches("gram_")), na.rm = T)

gram_results_by_subject_with_sums <- 
  gram_results_by_subject |> 
    mutate(sum = sums)

gram_by_subject_long_one_gram <- 
  gram_results_by_subject_with_sums |> 
    filter(sum == 1) |> 
    pivot_longer(matches("gram_"), names_to = "gram", values_to = "value") |> 
    filter(value) |> 
    distinct()

gram_by_subject_long_mixed <- 
  gram_results_by_subject_with_sums |> 
    filter(sum > 1) |> 
    pivot_longer(matches("gram_"), names_to = "gram", values_to = "value") |> 
    filter(value) |> 
    mutate(gram = "gram_MIXED") |> 
    distinct()

gram_by_subject_long <- bind_rows(gram_by_subject_long_one_gram, gram_by_subject_long_mixed) |> select(-sum, -value)
  
gt_gram_by_subject_long <- 
  tbl_summary(gram_by_subject_long |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, label = c("gram" ~ "Gram Finding", "max_poly" ~ "Maximum Neutrophils", "max_epi" ~ "Maximum Squamous Epithelial Cells", "max_rbc" ~ "Maximum Red Blood Cells"), statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing = "ifany", missing_text = "None") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
  bold_p() |> 
  add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]")) |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
#  modify_footnote(everything() ~ NA, abbreviation = T) |> 
#  modify_footnote(update = everything() ~ NA) |>
#  modify_footnote(update = everything() ~ "N, [95% CI]") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Gram Stain Results", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_gram_by_subject_long, filename = paste0(figure_path, "gram_stain_results_by_subject_long_table.html"))
gt::gtsave(gt_gram_by_subject_long, filename = paste0(figure_path, "gram_stain_results_by_subject_long_table.docx"))


gram_by_subject_wide <- bind_rows(gram_by_subject_long_one_gram, gram_by_subject_long_mixed) |> pivot_wider(names_from = gram, values_from = value) |> select(-sum) |> mutate(across(matches("gram"), ~ ifelse(is.na(.x), F, .x)))

gt_gram_by_subject_wide <- 
  tbl_summary(gram_by_subject_wide |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, label = c("max_poly" ~ "Maximum Neutrophils", "max_epi" ~ "Maximum Squamous Epithelial Cells", "max_rbc" ~ "Maximum Red Blood Cells"), statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing = "ifany", missing_text = "None") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>
  add_q(method = "BH") |> 
  bold_p() |> 
  add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]")) |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
  #  modify_footnote(everything() ~ NA, abbreviation = T) |> 
  #  modify_footnote(update = everything() ~ NA) |>
  #  modify_footnote(update = everything() ~ "N, [95% CI]") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Gram Stain Results", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_gram_by_subject_long, filename = paste0(figure_path, "gram_stain_results_by_subject_long_table.html"))
gt::gtsave(gt_gram_by_subject_long, filename = paste0(figure_path, "gram_stain_results_by_subject_long_table.docx"))







culture_results_by_subject <- 
  culture_data_input |> 
    left_join(gram_stain_data |> select(SubjectID, LRT_REF_GS_SPECTYPE) |> drop_na()) |> 
    group_by(SubjectID, LRT_REF_GS_SPECTYPE) |>
    mutate(NUM_UNIQUE_ORGANISMS = factor(ifelse(NUM_UNIQUE_ORGANISMS >= 3, ">=3", NUM_UNIQUE_ORGANISMS), levels = c("0", "1", "2", ">=3"), ordered = T)) |>
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

culture_data_input |>
  left_join(gram_stain_data |> select(SubjectID, LRT_REF_GS_SPECTYPE) |> drop_na()) |>
  group_by(SubjectID, LRT_REF_GS_SPECTYPE) |>
    summarise(any_ORG_REPORTED = any(!is.na(MALDI_ID), na.rm = T)) |> 
    ungroup() |> 
    count(LRT_REF_GS_SPECTYPE, any_ORG_REPORTED) |>
    group_by(LRT_REF_GS_SPECTYPE) |>
    mutate(prop = n / sum(n) * 100)

gt_culture_by_subject <- 
  tbl_summary(culture_results_by_subject |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, 
              label = c("num_org" ~ "Number of Unique Organisms", "any_contam" ~ "Contaminated?", "culture_GPC" ~ "Gram-Positive Cocci", "culture_GNB" ~ "Gram-Negative Bacilli", "culture_MIXED" ~ "Mixed Microorganism", "culture_GPB" ~ "Gram-Positive Bacilli", "culture_YEAST" ~ "Yeast"), 
              statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing_text = "N/A") |>
  italicize_levels() |>
  bold_labels() |> 
  add_p(test = everything() ~ "fisher.test", test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)) |>  
  bold_p() |> 
  add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]", all_continuous() ~ "[{conf.low} - {conf.high}]")) |> 
  modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
#  modify_footnote(everything() ~ NA, abbreviation = T) |> 
#  modify_footnote(update = everything() ~ NA) |>
#  modify_footnote(update = everything() ~ "N, (Percentage). [95% Confidence Interval]") |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Culture Results", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_culture_by_subject, filename = paste0(figure_path, "culture_results_by_subject_table.html"))
gt::gtsave(gt_culture_by_subject, filename = paste0(figure_path, "culture_results_by_subject_table.docx"))

gt_culture_by_subject



threshold_growth <-
  culture_data_input |> 
    left_join(gram_stain_data |> select(SubjectID, LRT_REF_GS_SPECTYPE) |> distinct() |> drop_na()) |> 
    mutate(threshold_growth_bin = 
             case_when(
               LRT_REF_GS_SPECTYPE == "BAL" & max_abundance %in% c("10,000-100,000", ">100,000") ~ "Significant",
               LRT_REF_GS_SPECTYPE == "ETA" & max_abundance %in% c(">100,000") ~ "Significant",
              LRT_REF_GRWTH_YN == "No" ~ "No Organisms Identified",
              LRT_REF_GRWTH_YN == "Yes" ~ "Insignificant"),
           threshold_growth_bin = factor(threshold_growth_bin, levels = c("No Organisms Identified", "Insignificant", "Significant"), ordered = T))
  
threshold_growth_by_subject <- 
  threshold_growth |> 
    group_by(SubjectID, LRT_REF_GS_SPECTYPE) |> 
    summarise(max_threshold_growth = max(threshold_growth_bin, na.rm = T))

gt_threshold_growth_by_subject <- 
  tbl_summary(threshold_growth_by_subject |> select(-SubjectID), by = LRT_REF_GS_SPECTYPE, 
              include = c("max_threshold_growth"), label = c("max_threshold_growth" ~ "Max Threshold Growth"), 
              statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing_text = "N/A") |>
    italicize_levels() |>
    bold_labels() |> 
    add_p() |>  
    bold_p() |> 
    add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]", all_continuous() ~ "[{conf.low} - {conf.high}]")) |> 
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*") |>
#    modify_footnote(everything() ~ NA, abbreviation = T) |> 
#    modify_footnote(update = everything() ~ NA) |>
#    modify_footnote(update = everything() ~ "N, (Percentage). [95% Confidence Interval]") |>
    as_gt() |> 
    theme_pattern() |>
    tab_header(title = "Threshold Growth", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_threshold_growth_by_subject, filename = paste0(figure_path, "threshold_results_by_subject_table.html"))
gt::gtsave(gt_threshold_growth_by_subject, filename = paste0(figure_path, "threshold_results_by_subject_table.docx"))



threshold_growth_joined_with_gram <- 
  left_join(threshold_growth_by_subject, gram_results_by_subject) |> 
    mutate(abundant_rbc = ifelse(max_rbc == "Abundant", "Abundant RBCs", "Few or No RBCs"))

gt_threshold_growth_by_procedure_and_rbc <-
  tbl_strata(
    data = threshold_growth_joined_with_gram |> ungroup(),
    strata = c(LRT_REF_GS_SPECTYPE),
    .tbl_fun = 
      ~ tbl_summary(
          .x, 
              by = abundant_rbc, 
              include = c("max_threshold_growth"), 
              label = c("max_threshold_growth" ~ "Max Threshold Growth"), 
              statistic = list(all_categorical() ~ "{n}, ({p}%)"), missing_text = "N/A") |>
    italicize_levels() |>
    bold_labels() |> 
    add_p() |>  
    bold_p() |> 
    add_ci(pattern = "{stat}<br>{ci}", statistic = list(all_categorical() ~ "[{conf.low}% - {conf.high}%]", all_continuous() ~ "[{conf.low} - {conf.high}]")) |> 
    modify_header(all_stat_cols() ~ "**{level}**<br>*N = {n}*")) |>
  as_gt() |> 
  theme_pattern() |>
  tab_header(title = "Threshold Growth by RBC Abundance", subtitle = "Per Subject, By Procedure Type")

gt::gtsave(gt_threshold_growth_by_procedure_and_rbc, filename = paste0(figure_path, "threshold_results_by_subject_and_rbc_abundance_table.html"))
gt::gtsave(gt_threshold_growth_by_procedure_and_rbc, filename = paste0(figure_path, "threshold_results_by_subject_and_rbc_abundance__table.docx"))

gt_threshold_growth_by_procedure_and_rbc



counts_of_orgs_from_no_gram_subjects <- 
  culture_data_input |> 
    filter(SubjectID %in% c(
      gram_stain_concordance_by_subject |> filter(any_gram_equals_culture == "No organsims on Gram stain, with culture growth") |> pluck("SubjectID"))) |>
      count(MALDI_ID_GROUPED, sort = T) |> 
      mutate(prop_of_subjects = n/length(gram_stain_concordance_by_subject |> filter(any_gram_equals_culture == "No organsims on Gram stain, with culture growth") |> pluck("SubjectID")))

significant_growth_pairs <- 
  threshold_growth |> 
    left_join(gram_results_by_subject |> distinct(SubjectID, LRT_REF_GS_SPECTYPE)) |> 
    group_by(SubjectID) |>
    mutate(num_significant_growth = n()) |>
    filter(num_significant_growth > 1)


org_counts <- significant_growth_pairs |> count(MALDI_ID_GROUPED, sort = T)
expected_props <- org_counts |> mutate(prop = n / length(unique(significant_growth_pairs$SpecimenID)))

organism_cooccurrence <- crossprod(table(significant_growth_pairs[c("SubjectID", "MALDI_ID_GROUPED")]))
## remove diagonal
diag(organism_cooccurrence) <- 0
## find maximum value in matrix
max_value <- max(organism_cooccurrence)
## find the row and column of the maximum value
max_value_row <- which(organism_cooccurrence == max_value, arr.ind = TRUE)[1]
max_value_col <- which(organism_cooccurrence == max_value, arr.ind = TRUE)[2]
## find the organism names
rownames(organism_cooccurrence)[max_value_row]
rownames(organism_cooccurrence)[max_value_col]

org_counts_BAL <- significant_growth_pairs |> filter(LRT_REF_GS_SPECTYPE == "BAL") |> count(MALDI_ID_GROUPED, sort = T)

organism_cooccurrence_BAL <- crossprod(table(org_counts_BAL[c("SubjectID", "MALDI_ID_GROUPED")]))
## remove diagonal
diag(organism_cooccurrence_BAL) <- 0
## find maximum value in matrix
max_value_BAL <- max(organism_cooccurrence_BAL)
## find the row and column of the maximum value
max_value_row_BAL <- which(organism_cooccurrence_BAL == max_value, arr.ind = TRUE)[1]
max_value_col_BAL <- which(organism_cooccurrence_BAL == max_value, arr.ind = TRUE)[2]
## find the organism names
rownames(organism_cooccurrence_BAL)[max_value_row]
rownames(organism_cooccurrence_BAL)[max_value_col]

org_counts_ETA <- significant_growth_pairs |> filter(LRT_REF_GS_SPECTYPE == "ETA") |> count(MALDI_ID_GROUPED, sort = T)

organism_cooccurrence_ETA <- crossprod(table(org_counts_ETA[c("SubjectID", "MALDI_ID_GROUPED")]))
## remove diagonal
diag(organism_cooccurrence_ETA) <- 0
## find maximum value in matrix
max_value_ETA <- max(organism_cooccurrence_ETA)
## find the row and column of the maximum value
max_value_row_ETA <- which(organism_cooccurrence_ETA == max_value, arr.ind = TRUE)[1]
max_value_col_ETA <- which(organism_cooccurrence_ETA == max_value, arr.ind = TRUE)[2]
## find the organism names
rownames(organism_cooccurrence_ETA)[max_value_row]
rownames(organism_cooccurrence_ETA)[max_value_col]


