# Load required libraries
library(shiny)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(plotly)
library(arrow)

# Define PTM target information
TARGET_PTM <- list(
  Citrullination_test = list(name = "Citrullination_test",
                             path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib_sample.parquet",
                             symbol = "UniMod:7", site = "R", mass_shift = 0.984016),
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

metadata <- read_csv("../data/metadata.csv")

source('../R/helper_functions.R')

# ----------------- Helper Functions -----------------


# Function to process data given a ptm_choice
load_data <- function(ptm_choice) {
  ptm_info <- TARGET_PTM[[ptm_choice]]
  df <- read_parquet(ptm_info$path)
  df %>%
    mutate(
      has_target_PTM = str_detect(Modified.Sequence, ptm_info$symbol),
      has_target_site = str_detect(Stripped.Sequence, ptm_info$site),
      num_PTM = str_count(Modified.Sequence, fixed(ptm_info$symbol))
    )
}

# ----------------- Shiny Server -----------------

server <- function(input, output, session) {
  
  # Use reactiveVal for processed_data to allow manual updates
  raw_data <- reactiveVal(NULL)
  
  # Manually update processed_data when input$file_dropdown changes
  observeEvent(input$file_dropdown, {
    req(input$file_dropdown)
    new_data <- load_data(input$file_dropdown)
    raw_data(new_data)
    
    # Update Cancer Type Dropdown based on new processed data
    cancer_types <- unique(metadata$`Cancer Type`[!is.na(metadata$`Cancer Type`)])
    updateSelectInput(session, "cancer_type_dropdown",
                      choices = cancer_types, selected = cancer_types[1])
  })
  
  processed_data <- reactiveVal(NULL)
  
  observe({
    # q-value filter should only be applied to sequences with detected PTM
    # sequence without detected PTM will have PTM.Q.Value=0
    PTM_filtered <- raw_data() %>%
      filter(`PTM.Q.Value` <= input$PTM_qvalue_slider
             & `PTM.Site.Confidence` >= input$PTM_site_confidence_slider
             & has_target_PTM)
    normal_data <- raw_data() %>%
      filter(!has_target_PTM)
    combined <- rbind(PTM_filtered, normal_data)
    
    processed_data(combined)
  })
  
  # Total intensity data (log10 transformed)
  total_intensity_data <- reactive({
    req(processed_data())
    processed_data() %>%
      group_by(plate, ID, assay, evotip, well, has_target_PTM) %>%
      summarise(
        total_intensity = sum(Precursor.Quantity, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(log_total_intensity = log10(total_intensity)) %>%
      inner_join(metadata, by = "ID", suffix = c("", ".meta")) %>% 
      select(-ends_with(".meta")) %>% 
      filter(!is.na(`Cancer Type`))
  })
  
  # Protein-level aggregation
  protein_level_aggregation <- reactive({
    req(processed_data())
    processed_data() %>%
      filter(has_target_site) %>%
      group_by(plate, ID, assay, evotip, well, Protein.Group) %>%
      summarise(
        has_target_PTM = any(has_target_PTM), .groups='drop') %>%
      inner_join(metadata, by = "ID")
      
  })
  
  # Peptide-level aggregation
  peptide_level_aggregation <- reactive({
    req(processed_data())
    processed_data() %>%
      filter(has_target_site) %>%
      group_by(plate, ID, assay, evotip, well, Stripped.Sequence) %>%
      summarise(has_target_PTM = any(has_target_PTM),
                total_intensity = sum(Precursor.Quantity), .groups='drop') %>%
      inner_join(metadata, by = "ID")
  })
  
  
  
  # ----------------- Plot Outputs -----------------
  
  # Plot 1: Histogram of total intensity (combined case/control label)
  output$matplotlib_plot <- renderPlotly({
    req(total_intensity_data(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    
    plot_data <- total_intensity_data() %>%
      filter(!is.na(group)) %>%
      mutate(ptm_label = ifelse(has_target_PTM, "Modified peptide", "Unmodified peptide"),
             combined = paste(group, ptm_label, sep = ", "))
    
    plot_title <- paste("Distribution of total intensity of peptide with",
                        ptm_info$name, "on site", ptm_info$site)
    
    create_histogram_plot(plot_data,
                          x_var = "log_total_intensity",
                          fill_var = "combined",
                          bins = 100,
                          title = plot_title,
                          xlab = "Total intensity (log10)",
                          ylab = "Frequency")
  })
  
  # Plot 2: Protein-level case–control plot
  output$prot_prop_plot_case_control <- renderPlotly({
    req(protein_level_aggregation(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    plot_case_control(protein_level_aggregation(), ptm_info, type = "protein", bins = 50)
  })
  
  # Plot 3: Peptide-level case–control plot
  output$pep_prop_plot_case_control <- renderPlotly({
    req(peptide_level_aggregation(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    plot_case_control(peptide_level_aggregation(), ptm_info, type = "peptide", bins = 50)
  })
  
  # Plot 4: Histogram of total intensity filtered by Cancer Type
  output$histogram_plot <- renderPlotly({
    req(total_intensity_data(), input$cancer_type_dropdown, input$file_dropdown)
    
    plot_data <- total_intensity_data() %>%
      filter(has_target_PTM,
             `Cancer Type` %in% input$cancer_type_dropdown,
             assay %in% c("FT", "IgB"))
    
    plot_title <- paste("Distribution of Total Intensity with", input$file_dropdown)
    create_histogram_plot(plot_data,
                          x_var = "log_total_intensity",
                          fill_var = "`Cancer Type`",
                          bins = input$bins_slider,
                          title = plot_title,
                          xlab = "Total Intensity (log10)",
                          ylab = "Frequency",
                          facet_formula = "~ assay")
  })
  
  
  # Plot5: Peptide-level histogram by Cancer Type
  output$pep_prop_plot_cancer_type <- renderPlotly({
    req(peptide_level_aggregation(), input$cancer_type_dropdown, input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    
    plot_data <- peptide_level_aggregation() %>%
      filter(`Cancer Type` %in% input$cancer_type_dropdown,
             assay %in% c("FT", "IgB")) %>%
      group_by(ID, assay, group, `Cancer Type`) %>%
      summarise(mean_has_target_PTM = mean(has_target_PTM, na.rm = TRUE)) %>%
      ungroup()
    
    plot_title <- paste("Distribution of Proportion of Modified Peptide with",
                        ptm_info$name, "on site", ptm_info$site)
    create_histogram_plot(plot_data,
                          x_var = "mean_has_target_PTM",
                          fill_var = "`Cancer Type`",
                          bins = 50,
                          title = plot_title,
                          xlab = "Proportion of Modified Peptide",
                          ylab = "Frequency",
                          facet_formula = "~ assay")
  })
  
  # Plot 6: Protein-level histogram by Cancer Type
  output$prot_prop_plot_cancer_type <- renderPlotly({
    req(protein_level_aggregation(), input$cancer_type_dropdown, input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    
    plot_data <- protein_level_aggregation() %>%
      filter(`Cancer Type` %in% input$cancer_type_dropdown,
             assay %in% c("FT", "IgB")) %>%
      group_by(ID, assay, group, `Cancer Type`) %>%
      summarise(mean_has_target_PTM = mean(has_target_PTM, na.rm = TRUE)) %>%
      ungroup()
    
    plot_title <- paste("Distribution of Proportion of Modified Protein with",
                        ptm_info$name, "on site", ptm_info$site)
    create_histogram_plot(plot_data,
                          x_var = "mean_has_target_PTM",
                          fill_var = "`Cancer Type`",
                          bins = input$bins_slider,
                          title = plot_title,
                          xlab = "Proportion Modified",
                          ylab = "Frequency",
                          facet_formula = "~ assay")
  })
  
  
  # Plot 7: batch effect plot
  output$total_plate_batch_effect_plot <- renderPlotly({
    req(total_intensity_data(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    
    p <- total_intensity_data() %>% group_by(ID, plate) %>% 
      summarise(log_total_intensity = log10(sum(total_intensity)), .groups = 'drop') %>%
      ggplot(aes(x = factor(plate), y = log_total_intensity)) +
      geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +  # Violin plot for distribution
      geom_boxplot(width = 0.1, outlier.shape = NA) +  # Boxplot inside violin
      theme_minimal() +  
      labs(x = "Batch ID", y = "log10(Total intensity)", title = "Total intensity (all peptide + all assay) of samples grouped by plate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
    
    # Convert to interactive plotly plot, enabling hover only for outliers
    ggplotly(p, tooltip = "text")
  })
  
  # Plot 8: batch effect plot
  output$plate_batch_effect_plot <- renderPlotly({
    req(total_intensity_data(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    
    p <- total_intensity_data() %>% 
      ggplot(aes(x = factor(plate), y = log_total_intensity)) +
      geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +  # Violin plot for distribution
      geom_boxplot(width = 0.1, outlier.shape = NA) +  # Boxplot inside violin
      facet_wrap(~assay + has_target_PTM,
                 scales = "free_y",
                 labeller = labeller(has_target_PTM = c("TRUE" = "Modified peptides", "FALSE" = "Unmodified peptides"))) +
      theme_minimal() +  
      labs(x = "Batch ID", y = "log10(Total intensity)", title = "Total intensity of samples grouped by plate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed
    
    # Convert to interactive plotly plot, enabling hover only for outliers
    ggplotly(p, tooltip = "text")
  })
  
  output$violin_plot_conditions <- renderPlotly({
    req(total_intensity_data())  # Ensure data is available
    
    # Create Summation across Condition B (grouping by A)
    intensity_by_modification <- total_intensity_data() %>%
      group_by(ID, has_target_PTM, plate) %>%  # Condition A
      summarise(total_intensity_sum = sum(total_intensity, na.rm = TRUE), .groups = "drop")
    
    # Create Summation across Condition A (grouping by B)
    intensity_by_assay <- total_intensity_data() %>%
      group_by(ID, assay, plate) %>%  # Condition B (Modified vs. Unmodified)
      summarise(total_intensity_sum = sum(total_intensity, na.rm = TRUE), .groups = "drop")
    
    # Violin Plot for Condition A
    p1 <- ggplot(intensity_by_modification, aes(x = as.factor(plate), y = log10(total_intensity_sum))) +
      geom_violin(trim=FALSE, alpha = 0.5, fill='skyblue') +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      labs(title = "Total intensity (IgB + FT) of samples grouped by plate",
           x = "Plate", y = "Total Intensity") +
      facet_wrap(~has_target_PTM, 
                 scales = "free_y",
                 labeller = labeller(has_target_PTM = c("TRUE" = "Modified peptides", "FALSE" = "Unmodified peptides"))) +
      theme_minimal()
    
    # Violin Plot for Condition B
    p2 <- ggplot(intensity_by_assay, aes(x = as.factor(plate), y = log10(total_intensity_sum))) +
      geom_violin(trim=FALSE, alpha = 0.5, fill = "skyblue") +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      labs(title = "Total intensity (Modified + Unmodified) of samples grouped by plate",
           x = "Plate", y = "Total Intensity") +
      facet_wrap(~assay, scales = "free_y") +
      theme_minimal()
    
    # Convert ggplots to Plotly and display both
    subplot(ggplotly(p1), ggplotly(p2), nrows = 2, shareX = TRUE) 
  })
  
  
}
