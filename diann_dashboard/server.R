# Load required libraries
library(shiny)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(plotly)
library(arrow)  # For reading parquet files

# Define PTM target information
TARGET_PTM <- list(
  Citrullination = list(name = "Citrullination",
                        path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib.parquet",
                        symbol = "UniMod:7", site = "R", mass_shift = 0.984016),
  Hypusine = list(name = "Hypusine",
                  path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_hypusine_report-lib.parquet",
                  symbol = "Hypusine", site = "K", mass_shift = 114.042927),
  Deoxyhypusine = list(name = "Deoxyhypusine",
                       path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_deoxyhypusine_report-lib.parquet",
                       symbol = "Deoxyhypusine", site = "K", mass_shift = 98.031300)
)

# Define Server
server <- function(input, output, session) {
  # Reactive values to store processed data
  processed_data <- reactiveValues(data = NULL)
  
  # Function to load and process data
  process_data <- function(selected_ptm) {
    ptm_info <- TARGET_PTM[[selected_ptm]]
    df <- read_parquet(ptm_info$path)
    metadata <- read_csv("../data/metadata.csv")
    
    dat <- df %>%
      left_join(metadata, by = "ID") %>%
      mutate(
        has_target_PTM = str_detect(Modified.Sequence, ptm_info$symbol),
        has_target_site = str_detect(Stripped.Sequence, ptm_info$site),
        num_PTM = str_count(Modified.Sequence, fixed(ptm_info$symbol))
      )
    return(dat)
  }
  
  observe({
    req(input$ptm_dropdown)
    data <- process_data(input$ptm_dropdown)
    processed_data$data <- data
  })
  
  total_intensity_data <- reactiveValues(data=NULL)
  
  observe({
    req(processed_data$data)
    total_intensity_data$data <-  processed_data$data %>%
      group_by(ID, assay, has_target_PTM) %>%
      summarise(total_intensity = sum(Precursor.Quantity, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(total_intensity = log10(total_intensity)) %>%
      left_join(metadata, by = "ID")
    
  })

  
  # Update Cancer Type Dropdown based on processed data
  observe({
    req(input$file_dropdown)
    data <- process_data(input$file_dropdown)
    processed_data$data <- data
    
    cancer_types <- unique(data$`Cancer Type`[!is.na(data$`Cancer Type`)])
    
    updateSelectInput(session, "cancer_type_dropdown",
                      choices = cancer_types, selected = cancer_types[1])
  })
  
  # Generate Histogram Plot
  output$histogram_plot <- renderPlotly({
    req(total_intensity_data$data, input$cancer_type_dropdown)
    total_intensity <- total_intensity_data$data %>%
      filter(has_target_PTM, `Cancer Type` %in% input$cancer_type_dropdown, assay %in% c("FT", "IgB"))
    
    unique_assays <- unique(total_intensity$assay)
    
    p <- ggplot(total_intensity, aes(x = total_intensity, fill = `Cancer Type`)) +
      geom_histogram(bins = input$bins_slider, alpha = 0.5, position = "identity") +
      facet_wrap(~assay) +
      labs(title = paste("Distribution of Total Intensity with", input$file_dropdown),
           x = "Total Intensity (log10)", y = "Frequency") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  
  output$matplotlib_plot <- renderPlotly({
    req(processed_data())
    
    total_intensity <- processed_data() %>%
      filter(!is.na(group))
    
    color_mapping <- list(
      "Case, Modified peptide" = "green",
      "Case, Unmodified peptide" = "orange",
      "Control, Modified peptide" = "red",
      "Control, Unmodified peptide" = "yellow"
    )
    
    total_intensity_labeled <- total_intensity_labeled %>%
      mutate(ptm_label = ifelse(has_target_PTM, "Modified peptide", "Unmodified peptide"),
             combined = paste(group, ptm_label, sep = ", "))
    
    
    p <- ggplot(total_intensity_labeled, aes(x = total_intensity, fill = combined)) +
      geom_histogram(bins = 100, alpha = 0.5, position = "identity") +
      labs(x = "Total intensity (log10)", y = "Frequency",
           title = paste("Distribution of total intensity of peptide with", TARGET_PTM[1],
                         "on site", TARGET_PTM[3])) +
      theme_minimal()
    
    ggplotly(p)
  })
  
}
