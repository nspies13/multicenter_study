set.seed(12345)

# Configure local, non-versioned inputs and generated outputs. Environment
# variables make the workflow portable when this directory is used as a
# standalone repository.
data_path <- normalizePath(
  Sys.getenv("MULTICENTER_STUDY_DATA_DIR", unset = "../../../data"),
  mustWork = FALSE,
  winslash = "/"
)
figure_path <- normalizePath(
  Sys.getenv("MULTICENTER_STUDY_FIGURE_DIR", unset = "figures/multicenter_study"),
  mustWork = FALSE,
  winslash = "/"
)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)
data_path <- paste0(data_path, "/")
figure_path <- paste0(figure_path, "/")

theme_ns <- theme(text = element_text(family = "Helvetica"),
                  plot.title = element_text(size = 16, face = "bold.italic", hjust = 0, margin = margin(0, 0, 4, 0)),
                  plot.subtitle = element_text(size = 14, face = "plain", hjust = 0, margin = margin(0, 0, 2, 0)),
                  axis.title = element_text(size = 12, face = "bold", margin = margin(4,4,4,4)),
                  axis.title.x.bottom = element_text(face = "bold", margin = margin(4,0,0,0)),
                  axis.title.y.left = element_text(face = "bold", margin = margin(0,4,0,0)),
                  legend.title = element_text(face = "bold.italic", size = 12),
                  axis.line = element_line(),
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(), 
                  panel.background = element_blank(),
                  plot.background = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_blank(),
                  strip.text = element_text(size = 12, face = "bold.italic"),
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
theme_poster <- function(data, ...){
  
  data %>% 
    tab_options(
      heading.background.color = "#a51417ff",
      column_labels.background.color = "gray30",
      footnotes.background.color = "#F7F7FB",
      ...
    ) |> 
    tab_style(
      style = cell_text(
        font = "Helvetica"
      ), 
      locations = cells_body()
    ) |>
    tab_style(
      style = cell_text(
        font = "Helvetica"
      ), 
      locations = cells_footnotes()
    ) |>
    tab_style(
      style = cell_text(
        font = "Helvetica"
      ), 
      locations = cells_column_labels()
    ) |>
    tab_style(
      style = cell_text(
        font = "Helvetica",
        weight = "bold",
        size = "xx-large"
      ), 
      locations = cells_title(groups = "title")
    ) |> 
    tab_style(
      style = cell_text(
        font = "Helvetica",
        style = "italic",
        size = "x-large"
      ), 
      locations = cells_title(groups = "subtitle")
    )
  
}
