# Load required libraries
library(shiny)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(plotly)
library(arrow)
library(data.table)

metadata <- read_csv("../data/metadata.csv")

# ----------------- Helper Functions -----------------

# Function to determine missing_fill value based on UI input
get_missing_fill_value <- function(missing_handling) {
  if (missing_handling == "fill with 0") {
    return(0)
  } else {
    return(NA)
  }
}

CAMPAGIN <- "RAINBOW"

# Function to process data given a ptm_choice
load_data <- function(ptm_choice) {
  ptm_info <- TARGET_PTM[[ptm_choice]]
  path <- glue::glue("../data/{CAMPAGIN}/UNIPROT_human_revi_2024_12_19_ProteinAG_{ptm_info$symbol}_report-lib.predicted.speclib.parquet")
  df <- read_parquet(path)
  df %>%
    mutate(
      has_target_PTM = str_detect(Modified.Sequence, paste0('\\(', ptm_info$symbol, '\\)', '|', '\\(', ptm_info$unimod, '\\)')),
      has_target_site = str_detect(Stripped.Sequence, ptm_info$site)
      # num_PTM = str_count(Modified.Sequence, fixed(ptm_info$symbol))
    )
}

# Function to create download handler
create_download_handler <- function(data_func, prefix, file_dropdown, normalization_dropdown) {
  downloadHandler(
    filename = function() {
      paste(prefix, TARGET_PTM[[file_dropdown]]$name, '_',
            normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      selected_table <- data_func()
      write.csv(selected_table, file, row.names = FALSE)
    }
  )
}

# Function to render data table
render_data_table <- function(data_func) {
  DT::renderDataTable({
    req(input$shared_cancer_type_dropdown)

    selected_table <- data_func()

    DT::datatable(selected_table,
                  options = list(scrollX = TRUE))
  })
}


# ----------------- Shiny Server -----------------

server <- function(input, output, session) {
  source('../R/helper_functions.R', local = TRUE)
  source('../R/PTM_definition.R', local = TRUE)
  source('../R/data_object.R', local = TRUE)
  source('../R/normalization.R', local = TRUE)
  source('../R/test.R', local = TRUE)
  source('../R/visualization.R', local = TRUE)

  # Use reactiveVal for processed_data to allow manual updates
  raw_data <- reactiveVal(NULL)

  # Manually update processed_data when input$file_dropdown changes
  observeEvent(input$file_dropdown, {
    req(input$file_dropdown)
    new_data <- load_data(input$file_dropdown)
    raw_data(new_data)
  })

  # Update Cancer Type Dropdown based on new processed data
  observe({
    cancer_types <- unique(metadata$`Cancer Type`[!is.na(metadata$`Cancer Type`)])
    updateSelectInput(session, "shared_cancer_type_dropdown", choices = cancer_types,
                      selected = cancer_types[1])

    # Update statistical analysis dropdowns
    updateSelectInput(session, "group1_cancer_types", choices = cancer_types,
                      selected = cancer_types[1])
    updateSelectInput(session, "group2_cancer_types", choices = cancer_types,
                      selected = cancer_types[2])
  })

  # Make group selections mutually exclusive
  observeEvent(input$group1_cancer_types, {
    req(input$group1_cancer_types)

    # Get all available cancer types
    all_cancer_types <- unique(metadata$`Cancer Type`[!is.na(metadata$`Cancer Type`)])

    if (input$auto_group2) {
      # Auto-select all cancer types not in group 1
      remaining_types <- setdiff(all_cancer_types, input$group1_cancer_types)
      updateSelectInput(session, "group2_cancer_types",
                       choices = all_cancer_types,
                       selected = remaining_types)
    } else {
      # Update group 2 choices to exclude group 1 selections
      available_for_group2 <- setdiff(all_cancer_types, input$group1_cancer_types)
      current_group2 <- intersect(input$group2_cancer_types, available_for_group2)

      updateSelectInput(session, "group2_cancer_types",
                       choices = available_for_group2,
                       selected = current_group2)
    }
  })

  observeEvent(input$group2_cancer_types, {
    req(input$group2_cancer_types)

    # Only update if auto_group2 is not checked to avoid circular updates
    if (!input$auto_group2) {
      # Get all available cancer types
      all_cancer_types <- unique(metadata$`Cancer Type`[!is.na(metadata$`Cancer Type`)])

      # Update group 1 choices to exclude group 2 selections
      available_for_group1 <- setdiff(all_cancer_types, input$group2_cancer_types)
      current_group1 <- intersect(input$group1_cancer_types, available_for_group1)

      updateSelectInput(session, "group1_cancer_types",
                       choices = available_for_group1,
                       selected = current_group1)
    }
  })

  # Handle auto-selection checkbox changes
  observeEvent(input$auto_group2, {
    all_cancer_types <- unique(metadata$`Cancer Type`[!is.na(metadata$`Cancer Type`)])

    if (input$auto_group2) {
      # Enable auto-selection mode
      if (!is.null(input$group1_cancer_types) && length(input$group1_cancer_types) > 0) {
        remaining_types <- setdiff(all_cancer_types, input$group1_cancer_types)
        updateSelectInput(session, "group2_cancer_types",
                         choices = all_cancer_types,
                         selected = remaining_types)
      }

      # Reset group 1 choices to all types
      updateSelectInput(session, "group1_cancer_types", choices = all_cancer_types)
    } else {
      # Disable auto-selection mode, restore manual selection
      if (!is.null(input$group1_cancer_types) && length(input$group1_cancer_types) > 0) {
        available_for_group2 <- setdiff(all_cancer_types, input$group1_cancer_types)
        current_group2 <- intersect(input$group2_cancer_types, available_for_group2)

        updateSelectInput(session, "group2_cancer_types",
                         choices = available_for_group2,
                         selected = current_group2)
      }
    }
  })

  # Processed data
  processed_data <- reactive({
    req(raw_data())
    PTM_filtered <- raw_data() %>%
      filter(`PTM.Q.Value` <= input$PTM_qvalue_slider,
             `PTM.Site.Confidence` >= input$PTM_site_confidence_slider,
             has_target_PTM)
    normal_data <- raw_data() %>%
      filter(!has_target_PTM)
    # Combine PTM and normal data
    # the unmodified peptides have to be processed separately because they have a PTM related confidence as 0
    rbind(PTM_filtered, normal_data)
  })

  # Total intensity data
  total_intensity_data <- reactive({
    req(processed_data())

    processed_data() %>%
      group_by(Run, has_target_PTM) %>%
      summarise(
        total_intensity = sum(Precursor.Quantity, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(log_total_intensity = log10(total_intensity)) %>%
      inner_join(generated_metadata(), by = 'Run', suffix = c('', '.meta')) %>%
      select(-contains('.meta'))
  })

  # Protein-level aggregation
  protein_level_aggregation <- reactive({
    req(processed_data())
    # processed_data() %>%
    #   filter(has_target_site) %>%
    #   group_by(Run, Protein.Group) %>%
    #   summarise(
    #     has_target_PTM = any(has_target_PTM), .groups='drop') %>%
    #   inner_join(generated_metadata(), by = 'Run', suffix = c('', '.meta')) %>%
    #   select(-contains('.meta'))

    dt <- as.data.table(processed_data())  # Convert to data.table

    dt <- dt[, .(
      Genes = Genes[1],
      Protein.Ids = Protein.Ids[1],
      has_target_PTM = has_target_PTM[1],
      total_intensity = sum(Precursor.Quantity, na.rm = TRUE)
    ), by = .(Run, Protein.Group)]

    merge(dt, generated_metadata(), by = "Run", suffixes = c("", ".meta"))


  })

  # Peptide-level aggregation
  # peptide level aggregate peptides with the same stripped sequence
  # no matter the difference in modified sequence (sites, number of PTMs)
  # for example, for peptide "PEPTIDE", its stripped sequence is "PEPTIDE". it have two modified peptides: "PEPTIDE(P)" and "PE(P)PTIDE". It will result in one peptide level entry with two modified peptides.
  peptide_level_aggregation <- reactive({
    req(processed_data())
    # processed_data() %>%
    #   filter(has_target_site) %>%
    #   group_by(Run, Stripped.Sequence) %>%
    #   summarise(has_target_PTM = any(has_target_PTM),
    #             total_intensity = sum(Precursor.Quantity), .groups='drop') %>%
    #   inner_join(generated_metadata(), by = 'Run', suffix = c('', '.meta')) %>%
    #   select(-contains('.meta'))

    dt <- as.data.table(processed_data())  # Convert to data.table

    dt <- dt[, .(
      Genes = Genes[1],
      Protein.Group = Protein.Group[1],
      Protein.Ids = Protein.Ids[1],
      has_target_PTM = has_target_PTM[1],
      total_intensity = sum(Precursor.Quantity, na.rm = TRUE)
    ), by = .(Run, Stripped.Sequence)]

    merge(dt, generated_metadata(), by = "Run", suffixes = c("", ".meta"))

  })


  # modified-peptide level
  modified_peptide_level <- reactive({

    dt <- as.data.table(processed_data())  # Convert to data.table

    dt[, .(
        Genes = Genes[1],
        Protein.Group = Protein.Group[1],
        Protein.Ids = Protein.Ids[1],
        has_target_PTM = has_target_PTM[1],
        Stripped.Sequence = Stripped.Sequence[1],
        total_intensity = sum(Precursor.Quantity, na.rm = TRUE)
      ), by = .(Run, Modified.Sequence)]

      })


  modified_pep_wide_format <- reactive({
    req(modified_peptide_level())
    PTM_name <- TARGET_PTM[[input$file_dropdown]]$name
    modified_peptide_level() %>%
      rename(!!PTM_name:= has_target_PTM) %>%
      # mutate(total_intensity = log10(total_intensity)) %>%
      tidyr::pivot_wider(id_cols=c(Protein.Group, Protein.Ids, Genes, Modified.Sequence, Stripped.Sequence, {{PTM_name}}),
                         names_from = Run,
                         values_from = total_intensity, values_fill = NA)
  })

  wide_format <- reactive({
    req(peptide_level_aggregation())

    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    PTM_name = ptm_info$name
    peptide_level_aggregation() %>%
      rename(!!PTM_name:= has_target_PTM) %>%
      # mutate(total_intensity = log10(total_intensity)) %>%
      tidyr::pivot_wider(id_cols=c(Protein.Group, Protein.Ids, Genes, Stripped.Sequence, {{PTM_name}}),
                  names_from = Run,
                  values_from = total_intensity, values_fill = NA)
  })

  protein_wide_format <- reactive({
    req(protein_level_aggregation())

    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    PTM_name = ptm_info$name
    protein_level_aggregation() %>%
      rename(!!PTM_name:= has_target_PTM) %>%
      # mutate(total_intensity = log10(total_intensity)) %>%
      tidyr::pivot_wider(id_cols=c(Protein.Group, Protein.Ids, Genes, {{PTM_name}}),
                  names_from = Run,
                  values_from = total_intensity, values_fill = NA)
  })

  annotated_data <- reactive({
    dat <- wide_format() %>% select(-c(Protein.Group, Protein.Ids, Genes, TARGET_PTM[[input$file_dropdown]]$name))
    dat_matrix <- as.matrix(dat[-1])
    # create row metadata
    row_metadata <- wide_format()[, c('Stripped.Sequence', TARGET_PTM[[input$file_dropdown]]$name, 'Protein.Group', 'Protein.Ids', 'Genes')]

    meta <- generated_metadata()
    # filter dat_matrix only keep those in metadata
    dat_matrix <- dat_matrix[, colnames(dat_matrix) %in% meta$Run]
    # Step 2: Reorder metadata based on column order in dat_matrix
    ordered_meta <- meta[match(colnames(dat_matrix), meta$Run), ]
    AnnotatedData$new(dat_matrix,
                                   col_metadata=ordered_meta,
                                   row_metadata=row_metadata)
  })

  protein_annotated_data <- reactive({
    dat <- protein_wide_format() %>% select(-c(Protein.Group, Protein.Ids, Genes, TARGET_PTM[[input$file_dropdown]]$name))
    dat_matrix <- as.matrix(dat)
    # create row metadata
    row_metadata <- protein_wide_format()[, c('Protein.Group', 'Protein.Ids', 'Genes', TARGET_PTM[[input$file_dropdown]]$name)]

    meta <- generated_metadata()
    # filter dat_matrix only keep those in metadata
    dat_matrix <- dat_matrix[, colnames(dat_matrix) %in% meta$Run]
    # Step 2: Reorder metadata based on column order in dat_matrix
    ordered_meta <- meta[match(colnames(dat_matrix), meta$Run), ]
    AnnotatedData$new(dat_matrix,
                                       col_metadata=ordered_meta,
                                       row_metadata=row_metadata)
  })

  modified_pep_annotated_data <- reactive({
    dat <- modified_pep_wide_format() %>% select(-c(`Modified.Sequence`, `Stripped.Sequence`, Protein.Group, Protein.Ids, Genes, TARGET_PTM[[input$file_dropdown]]$name))
    dat_matrix <- as.matrix(dat)
    # create row metadata
    row_metadata <- modified_pep_wide_format()[, c('Modified.Sequence', 'Stripped.Sequence', TARGET_PTM[[input$file_dropdown]]$name, 'Protein.Group', 'Protein.Ids', 'Genes')]

    meta <- generated_metadata()
    # filter dat_matrix only keep those in metadata
    dat_matrix <- dat_matrix[, colnames(dat_matrix) %in% meta$Run]
    # Step 2: Reorder metadata based on column order in dat_matrix
    ordered_meta <- meta[match(colnames(dat_matrix), meta$Run), ]
    AnnotatedData$new(dat_matrix,
                                       col_metadata=ordered_meta,
                                       row_metadata=row_metadata)
  })

  generated_metadata <- reactive({
    req(processed_data())

    processed_data()  %>%
    select(Run, IPAS, ID, plate, assay, evotip, well) %>%
    unique() %>%
    inner_join(metadata, by='ID') %>%
    filter(!is.na(`Cancer Type`))
  })


  # ----------------- Plot Outputs -----------------

  # Plot 1: Histogram of total intensity (combined case/control label)
  output$intensity_plot <- renderPlotly({
    req(total_intensity_data(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]

    plot_data <- total_intensity_data() %>%
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

  # Plot 3: Peptide-level case–control plot
  output$pep_prop_plot_case_control <- renderPlotly({
    req(peptide_level_aggregation(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    plot_case_control(peptide_level_aggregation(), generated_metadata(), ptm_info, type = "peptide", bins = 50)
  })

  # Plot 2: Protein-level case–control plot
  output$prot_prop_plot_case_control <- renderPlotly({
    req(protein_level_aggregation(), input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]
    plot_case_control(protein_level_aggregation(), generated_metadata(), ptm_info, type = "protein", bins = 50)
  })



  # Plot 4: Histogram of total intensity filtered by Cancer Type
  output$intensity_plot_cancer_type <- renderPlotly({
    req(total_intensity_data(), input$shared_cancer_type_dropdown, input$file_dropdown)

    plot_data <- total_intensity_data() %>%
      filter(has_target_PTM,
             `Cancer Type` %in% input$shared_cancer_type_dropdown,
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
    req(peptide_level_aggregation(), input$shared_cancer_type_dropdown, input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]

    plot_data <- peptide_level_aggregation() %>%
      filter(`Cancer Type` %in% input$shared_cancer_type_dropdown,
             assay %in% c("FT", "IgB")) %>%
      group_by(Run, `Cancer Type`, assay) %>%
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
    req(protein_level_aggregation(), input$shared_cancer_type_dropdown, input$file_dropdown)
    ptm_info <- TARGET_PTM[[input$file_dropdown]]

    plot_data <- protein_level_aggregation() %>%
      filter(`Cancer Type` %in% input$shared_cancer_type_dropdown,
             assay %in% c("FT", "IgB")) %>%
      group_by(Run, `Cancer Type`, assay) %>%
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


  #-------------------------modified peptide level------------------------------
  output$peptide_table <- DT::renderDataTable({
    req(input$shared_cancer_type_dropdown, input$missing_handling)

    # select column
    selected_table <- annotated_data()$create_filtered_view(
      assay = input$shared_assay_dropdown,
      cancer_type=input$shared_cancer_type_dropdown,
      normalization=input$shared_normalization_dropdown,
      missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()

    DT::datatable(selected_table,
                  options = list(
                    scrollX = TRUE
                  ))

  })

  output$modified_peptide_table <- DT::renderDataTable({
    req(input$shared_cancer_type_dropdown, input$missing_handling)

    # select column

    selected_table <- modified_pep_annotated_data()$create_filtered_view(
      assay = input$shared_assay_dropdown,
      cancer_type=input$shared_cancer_type_dropdown,
      normalization=input$shared_normalization_dropdown,
      missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()

    DT::datatable(selected_table,
                  options = list(
                    scrollX = TRUE
                  ))

  })

  output$protein_table <- DT::renderDataTable({
    req(input$shared_cancer_type_dropdown, input$missing_handling)

    # select column
    selected_table <- protein_annotated_data()$create_filtered_view(
      assay = input$shared_assay_dropdown,
      cancer_type=input$shared_cancer_type_dropdown,
      normalization=input$shared_normalization_dropdown,
      missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()

    DT::datatable(selected_table,
                  options = list(
                    scrollX = TRUE
                  ))

  })

  output$download_modified_peptide <- downloadHandler(
    filename = function() {
      paste("modified_peptide_level_", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Assuming peptide_data() is the reactive that contains your data for the table.
      # You can also use the data directly if it isn't reactive.
      req(input$missing_handling)
      selected_table <- modified_pep_annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      write.csv(selected_table, file, row.names = FALSE)
    }
  )

  output$download_modified_peptide_long <- downloadHandler(
    filename = function() {
      paste("modified_peptide_long_level_", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Assuming peptide_data() is the reactive that contains your data for the table.
      # You can also use the data directly if it isn't reactive.
      req(input$missing_handling)
      selected_table <- modified_pep_annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      # convert to long format
      selected_table <- selected_table %>%
        # this is bad hardcoded todo: fix it
        tidyr::pivot_longer(cols = 7:ncol(.),
                     names_to = "Run",
                     values_to = "intensity") %>%
        filter(!is.na(intensity))

      write.csv(selected_table, file, row.names = FALSE)
    }
  )


  output$download_peptide <- downloadHandler(
    filename = function() {
      paste("peptide_level_", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Assuming peptide_data() is the reactive that contains your data for the table.
      # You can also use the data directly if it isn't reactive.
      req(input$missing_handling)
      selected_table <- annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      write.csv(selected_table, file, row.names = FALSE)
    }
  )

  output$download_peptide_long <- downloadHandler(
    filename = function() {
      paste("peptide_level_long", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Assuming peptide_data() is the reactive that contains your data for the table.
      # You can also use the data directly if it isn't reactive.
      req(input$missing_handling)
      selected_table <- annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      # convert it into long format
      browser()
      selected_table <- selected_table %>%
        tidyr::pivot_longer(cols = 6:ncol(.),
                     names_to = "run",
                     values_to = "intensity") %>%
        filter(!is.na(intensity))

      write.csv(selected_table, file, row.names = FALSE)
    }
  )

  output$download_protein <- downloadHandler(
    filename = function() {
      paste("protein_level_", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$missing_handling)
      selected_table <- protein_annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      write.csv(selected_table, file, row.names = FALSE)
    }
  )

  output$download_protein_long <- downloadHandler(
    filename = function() {
      paste("protein_level_long_", TARGET_PTM[[input$file_dropdown]]$name, '_',
            input$shared_normalization_dropdown, '_',
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$missing_handling)
      selected_table <- protein_annotated_data()$create_filtered_view(
        assay = input$shared_assay_dropdown,
        cancer_type=input$shared_cancer_type_dropdown,
        normalization=input$shared_normalization_dropdown,
        missing_fill = get_missing_fill_value(input$missing_handling))$get_annotated_data()
      # convert it into long format
      selected_table <- selected_table %>%
        tidyr::pivot_longer(cols = 5:ncol(.),
                     names_to = "run",
                     values_to = "intensity") %>%
        filter(!is.na(intensity))

      write.csv(selected_table, file, row.names = FALSE)
    }
  )

  output$download_metadata <- downloadHandler(
    filename = function() {
      paste("metadata_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Assuming peptide_data() is the reactive that contains your data for the table.
      # You can also use the data directly if it isn't reactive.
      write.csv(generated_metadata(), file, row.names = FALSE)
    }
  )

  # Statistical analysis reactive values
  ttest_results_modified_peptide <- reactiveVal(NULL)
  ttest_results_peptide <- reactiveVal(NULL)
  ttest_results_protein <- reactiveVal(NULL)

  # Run modified peptide analysis when button is clicked
  observeEvent(input$run_modified_peptide_test, {
    req(input$group1_cancer_types, input$group2_cancer_types,
        input$stats_assay_dropdown, input$shared_normalization_dropdown,
        input$stats_test_type, input$missing_handling)

    # Validate that groups are different
    if (length(intersect(input$group1_cancer_types, input$group2_cancer_types)) > 0) {
      showNotification("Error: Groups cannot have overlapping cancer types!", type = "error")
      return()
    }

    if (length(input$group1_cancer_types) == 0 || length(input$group2_cancer_types) == 0) {
      showNotification("Error: Both groups must have at least one cancer type selected!", type = "error")
      return()
    }

    # Show progress
    showNotification("Running modified peptide analysis...", type = "message", duration = 2)

    tryCatch({
      # Run analysis for modified peptide level
      modified_peptide_results <- perform_statistical_analysis(
        modified_pep_annotated_data(),
        input$group1_cancer_types,
        input$group2_cancer_types,
        input$stats_assay_dropdown,
        input$shared_normalization_dropdown,
        input$stats_test_type,
        input$pvalue_threshold,
        missing_fill = get_missing_fill_value(input$missing_handling)
      )
      ttest_results_modified_peptide(modified_peptide_results)

      showNotification("Modified peptide analysis completed successfully!", type = "success")

    }, error = function(e) {
      showNotification(paste("Error in modified peptide analysis:", e$message), type = "error")
    })
  })

  # Run peptide analysis when button is clicked
  observeEvent(input$run_peptide_test, {
    req(input$group1_cancer_types, input$group2_cancer_types,
        input$stats_assay_dropdown, input$shared_normalization_dropdown,
        input$stats_test_type, input$missing_handling)

    # Validate that groups are different
    if (length(intersect(input$group1_cancer_types, input$group2_cancer_types)) > 0) {
      showNotification("Error: Groups cannot have overlapping cancer types!", type = "error")
      return()
    }

    if (length(input$group1_cancer_types) == 0 || length(input$group2_cancer_types) == 0) {
      showNotification("Error: Both groups must have at least one cancer type selected!", type = "error")
      return()
    }

    # Show progress
    showNotification("Running peptide analysis...", type = "message", duration = 2)

    tryCatch({
      # Run analysis for peptide level
      peptide_results <- perform_statistical_analysis(
        annotated_data(),
        input$group1_cancer_types,
        input$group2_cancer_types,
        input$stats_assay_dropdown,
        input$shared_normalization_dropdown,
        input$stats_test_type,
        input$pvalue_threshold,
        missing_fill = get_missing_fill_value(input$missing_handling)
      )
      ttest_results_peptide(peptide_results)

      showNotification("Peptide analysis completed successfully!", type = "success")

    }, error = function(e) {
      showNotification(paste("Error in peptide analysis:", e$message), type = "error")
    })
  })

  # Run protein analysis when button is clicked
  observeEvent(input$run_protein_test, {
    req(input$group1_cancer_types, input$group2_cancer_types,
        input$stats_assay_dropdown, input$shared_normalization_dropdown,
        input$stats_test_type, input$missing_handling)

    # Validate that groups are different
    if (length(intersect(input$group1_cancer_types, input$group2_cancer_types)) > 0) {
      showNotification("Error: Groups cannot have overlapping cancer types!", type = "error")
      return()
    }

    if (length(input$group1_cancer_types) == 0 || length(input$group2_cancer_types) == 0) {
      showNotification("Error: Both groups must have at least one cancer type selected!", type = "error")
      return()
    }

    # Show progress
    showNotification("Running protein analysis...", type = "message", duration = 2)

    tryCatch({
      # Run analysis for protein level
      protein_results <- perform_statistical_analysis(
        protein_annotated_data(),
        input$group1_cancer_types,
        input$group2_cancer_types,
        input$stats_assay_dropdown,
        input$shared_normalization_dropdown,
        input$stats_test_type,
        input$pvalue_threshold,
        missing_fill = get_missing_fill_value(input$missing_handling)
      )
      ttest_results_protein(protein_results)

      showNotification("Protein analysis completed successfully!", type = "success")

    }, error = function(e) {
      showNotification(paste("Error in protein analysis:", e$message), type = "error")
    })
  })

  # Statistical Analysis Table Outputs
  output$ttest_modified_peptide_results <- DT::renderDataTable({
    req(ttest_results_modified_peptide())

    results <- ttest_results_modified_peptide()

    # Remove the significant column
    results$significant <- NULL

    # Apply filters
    if (!is.null(input$auc_filter)) {
      results <- results[!is.na(results$auc) &
                        abs(results$auc - 0.5) >= input$auc_filter[1] |
                        abs(results$auc - 0.5) >= input$auc_filter[2], ]
    }

    if (!is.null(input$log2fc_filter)) {
      results <- results[!is.na(results$log2_fc) &
                        results$log2_fc >= input$log2fc_filter[1] &
                        results$log2_fc <= input$log2fc_filter[2], ]
    }

    if (!is.null(input$show_significant_only) && input$show_significant_only) {
      results <- results[!is.na(results$p_value) &
                        results$p_value < input$pvalue_threshold, ]
    }

    # Format numeric columns
    results$log2_fc <- round(results$log2_fc, 3)
    results$p_value <- round(results$p_value, 6)
    results$mean_group1 <- round(results$mean_group1, 2)
    results$mean_group2 <- round(results$mean_group2, 2)
    results$auc <- round(results$auc, 3)

    dt <- DT::datatable(results,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 25,
                    order = list(list(which(colnames(results) == "p_value") - 1, 'asc'))
                  ),
                  rownames = FALSE)

    # Make p-values green if below threshold
    dt %>%
      DT::formatStyle("p_value",
                      backgroundColor = DT::styleInterval(input$pvalue_threshold, c("lightgreen", "white")))
  })

  output$ttest_peptide_results <- DT::renderDataTable({
    req(ttest_results_peptide())

    results <- ttest_results_peptide()

    # Remove the significant column
    results$significant <- NULL

    # Apply filters
    if (!is.null(input$auc_filter)) {
      results <- results[!is.na(results$auc) &
                        abs(results$auc - 0.5) >= input$auc_filter[1] |
                        abs(results$auc - 0.5) >= input$auc_filter[2], ]
    }

    if (!is.null(input$log2fc_filter)) {
      results <- results[!is.na(results$log2_fc) &
                        results$log2_fc >= input$log2fc_filter[1] &
                        results$log2_fc <= input$log2fc_filter[2], ]
    }

    if (!is.null(input$show_significant_only) && input$show_significant_only) {
      results <- results[!is.na(results$p_value) &
                        results$p_value < input$pvalue_threshold, ]
    }

    # Format numeric columns
    results$log2_fc<- round(results$log2_fc, 3)
    results$p_value <- round(results$p_value, 6)
    results$mean_group1 <- round(results$mean_group1, 2)
    results$mean_group2 <- round(results$mean_group2, 2)
    results$auc <- round(results$auc, 3)

    dt <- DT::datatable(results,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 25,
                    order = list(list(which(colnames(results) == "p_value") - 1, 'asc'))
                  ),
                  rownames = FALSE)

    # Make p-values green if below threshold
    dt %>%
      DT::formatStyle("p_value",
                      backgroundColor = DT::styleInterval(input$pvalue_threshold, c("lightgreen", "white")))
  })

  output$ttest_protein_results <- DT::renderDataTable({
    req(ttest_results_protein())

    results <- ttest_results_protein()

    # Remove the significant column
    results$significant <- NULL

    # Apply filters
    if (!is.null(input$auc_filter)) {
      results <- results[!is.na(results$auc) &
                        abs(results$auc - 0.5) >= input$auc_filter[1] |
                        abs(results$auc - 0.5) >= input$auc_filter[2], ]
    }

    if (!is.null(input$log2fc_filter)) {
      results <- results[!is.na(results$log2_fc) &
                        results$log2_fc >= input$log2fc_filter[1] &
                        results$log2_fc <= input$log2fc_filter[2], ]
    }

    if (!is.null(input$show_significant_only) && input$show_significant_only) {
      results <- results[!is.na(results$p_value) &
                        results$p_value < input$pvalue_threshold, ]
    }

    # Format numeric columns
    results$log2_fc <- round(results$log2_fc, 3)
    results$p_value <- round(results$p_value, 6)
    results$mean_group1 <- round(results$mean_group1, 2)
    results$mean_group2 <- round(results$mean_group2, 2)
    results$auc <- round(results$auc, 3)

    dt <- DT::datatable(results,
                  options = list(
                    scrollX = TRUE,
                    pageLength = 25,
                    order = list(list(which(colnames(results) == "p_value") - 1, 'asc'))
                  ),
                  rownames = FALSE)

    # Make p-values green if below threshold
    dt %>%
      DT::formatStyle("p_value",
                      backgroundColor = DT::styleInterval(input$pvalue_threshold, c("lightgreen", "white")))
  })

  # Statistical Analysis Download Handlers
  output$download_ttest_modified_peptide <- downloadHandler(
    filename = function() {
      group1_str <- paste(input$group1_cancer_types, collapse = "_")
      group2_str <- paste(input$group2_cancer_types, collapse = "_")
      paste("ttest_modified_peptide_", group1_str, "_vs_", group2_str, "_",
            input$stats_assay_dropdown, "_", input$stats_normalization_dropdown, "_",
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(ttest_results_modified_peptide())
      write.csv(ttest_results_modified_peptide(), file, row.names = FALSE)
    }
  )

  output$download_ttest_peptide <- downloadHandler(
    filename = function() {
      group1_str <- paste(input$group1_cancer_types, collapse = "_")
      group2_str <- paste(input$group2_cancer_types, collapse = "_")
      paste("ttest_peptide_", group1_str, "_vs_", group2_str, "_",
            input$stats_assay_dropdown, "_", input$stats_normalization_dropdown, "_",
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(ttest_results_peptide())
      write.csv(ttest_results_peptide(), file, row.names = FALSE)
    }
  )

  output$download_ttest_protein <- downloadHandler(
    filename = function() {
      group1_str <- paste(input$group1_cancer_types, collapse = "_")
      group2_str <- paste(input$group2_cancer_types, collapse = "_")
      paste("ttest_protein_", group1_str, "_vs_", group2_str, "_",
            input$stats_assay_dropdown, "_", input$stats_normalization_dropdown, "_",
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(ttest_results_protein())
      write.csv(ttest_results_protein(), file, row.names = FALSE)
    }
  )

  # Volcano Plot Outputs
  output$volcano_modified_peptide <- renderPlotly({
    req(ttest_results_modified_peptide())
    
    results <- ttest_results_modified_peptide()
    
    # Calculate -log10(p_value) and handle p_value = 0 cases
    results$neg_log10_pval <- ifelse(results$p_value == 0, 
                                     max(-log10(results$p_value[results$p_value > 0]), na.rm = TRUE) + 1,
                                     -log10(results$p_value))
    
    # Create significance categories for coloring
    results$significance <- "Not Significant"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc > 0] <- "Upregulated"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc < 0] <- "Downregulated"
    
    # Create the volcano plot
    p <- ggplot(results, aes(x = log2_fc, y = neg_log10_pval, color = significance)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("Not Significant" = "gray", 
                                   "Upregulated" = "red", 
                                   "Downregulated" = "blue")) +
      geom_hline(yintercept = -log10(input$pvalue_threshold), linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(x = "Log2 Fold Change", 
           y = "-Log10(P-value)",
           title = "Volcano Plot - Modified Peptide Level",
           color = "Significance") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$volcano_peptide <- renderPlotly({
    req(ttest_results_peptide())
    
    results <- ttest_results_peptide()
    
    # Calculate -log10(p_value) and handle p_value = 0 cases
    results$neg_log10_pval <- ifelse(results$p_value == 0, 
                                     max(-log10(results$p_value[results$p_value > 0]), na.rm = TRUE) + 1,
                                     -log10(results$p_value))
    
    # Create significance categories for coloring
    results$significance <- "Not Significant"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc > 0] <- "Upregulated"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc < 0] <- "Downregulated"
    
    # Create the volcano plot
    p <- ggplot(results, aes(x = log2_fc, y = neg_log10_pval, color = significance)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("Not Significant" = "gray", 
                                   "Upregulated" = "red", 
                                   "Downregulated" = "blue")) +
      geom_hline(yintercept = -log10(input$pvalue_threshold), linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(x = "Log2 Fold Change", 
           y = "-Log10(P-value)",
           title = "Volcano Plot - Peptide Level",
           color = "Significance") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$volcano_protein <- renderPlotly({
    req(ttest_results_protein())
    
    results <- ttest_results_protein()
    
    # Calculate -log10(p_value) and handle p_value = 0 cases
    results$neg_log10_pval <- ifelse(results$p_value == 0, 
                                     max(-log10(results$p_value[results$p_value > 0]), na.rm = TRUE) + 1,
                                     -log10(results$p_value))
    
    # Create significance categories for coloring
    results$significance <- "Not Significant"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc > 0] <- "Upregulated"
    results$significance[!is.na(results$p_value) & results$p_value < input$pvalue_threshold & 
                        !is.na(results$log2_fc) & results$log2_fc < 0] <- "Downregulated"
    
    # Create the volcano plot
    p <- ggplot(results, aes(x = log2_fc, y = neg_log10_pval, color = significance)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("Not Significant" = "gray", 
                                   "Upregulated" = "red", 
                                   "Downregulated" = "blue")) +
      geom_hline(yintercept = -log10(input$pvalue_threshold), linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(x = "Log2 Fold Change", 
           y = "-Log10(P-value)",
           title = "Volcano Plot - Protein Level",
           color = "Significance") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = c("x", "y"))
  })
}
