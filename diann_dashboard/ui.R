library(shiny)
library(plotly)

source('../R/PTM_definition.R')

# First page: PTM Analysis tab
ptmAnalysisTab <- tabPanel(
  title = "Overview",
  fluidPage(
    titlePanel("PTM Analysis in Shiny"),

    # Sidebar (Filters) at the Top
    fluidRow(
      column(4,
             selectInput("file_dropdown", "Select PTM File",
                         choices = names(TARGET_PTM),
                         selected = "Citrullination_test")),
      column(4,
             sliderInput("PTM_qvalue_slider", "PTM Q-value threshold", min=0, max=0.1,
                         value=0.95, step = 0.01)),
      column(4,
             sliderInput("PTM_site_confidence_slider", "Location confidence", min=0, max=1,
                         value=0, step = 0.01)),

    ),

    # Main Panel with Plots
    fluidRow(
      column(12, plotlyOutput("intensity_plot"))
    ),

    fluidRow(
      column(12, plotlyOutput('pep_prop_plot_case_control'))
    ),

    fluidRow(
      column(12, plotlyOutput('prot_prop_plot_case_control'))
    ),

    fluidRow(
      column(8,
             selectInput("cancer_type_dropdown", "Select Cancer Types",
                         choices = NULL, selected = c('Breast', 'MERIT Control'), multiple = TRUE)
      ),
      column(4,
             sliderInput("bins_slider", "Number of Bins", min = 10, max = 100, value = 30, step = 1)
      )
    ),

    fluidRow(
      column(12, plotlyOutput("intensity_plot_cancer_type"))
    ),

    fluidRow(
      column(12, plotlyOutput('pep_prop_plot_cancer_type'))
    ),

    fluidRow(
      column(12, plotlyOutput('prot_prop_plot_cancer_type'))
    )
  )
)

# Additional pages defined as tabPanel
batchEffect <- tabPanel(
  title = "Batch effect",
  fluidPage(
    plotlyOutput("total_plate_batch_effect_plot", height = "400px"),
    plotlyOutput("plate_batch_effect_plot", height = "800px"),
    plotlyOutput("violin_plot_conditions", height = "800px"),
  )
)

modifiedPeptideLevel <- tabPanel(
  title = "Modified peptide level",
  fluidPage(
    fluidRow(
      # dt table
      selectInput("cancer_type_dropdown_output", "Select Cancer Types",
                  choices = NULL, selected = c('Breast', 'MERIT Control'), multiple = TRUE),
      selectInput("assay_dropdown_output", "Assay",
                  choices = c('FT', 'IgB'), selected = c('FT','IgB'), multiple = TRUE),
      selectInput("normalization_dropdown_output", "normaliztion",
                  choices = c('none', 'median', 'plate_median'), selected = 'none'),
    ),

    # warning
    tags$div(
      class = "alert alert-warning",
      "Warning: normalization are performed separately for different assays."
    ),
    DT::dataTableOutput("modified_peptide_table"),
    downloadButton("download_modified_peptide", "Download Data"),
    downloadButton("download_modified_peptide_long", "Download Long Format Data"),
    downloadButton("download_metadata", "Download Metadata")
  )
)

peptideLevel <- tabPanel(
  title = "Peptide level",
  fluidPage(
    # warning
    tags$div(
      class = "alert alert-warning",
      "Warning: normalization are performed separately for different assays."
    ),
    DT::dataTableOutput("peptide_table"),
    downloadButton("download_peptide", "Download Data"),
    downloadButton("download_peptide_long", "Download Long Format Data")
  )
)

proteinLevel <- tabPanel(
  title = "Protein level",
  fluidPage(
    fluidRow(
      # dt table
      selectInput("cancer_type_dropdown_protein", "Select Cancer Types",
                  choices = NULL, selected = c('Breast', 'MERIT Control'), multiple = TRUE),
      selectInput("assay_dropdown_protein", "Assay",
                  choices = c('FT', 'IgB'), selected = c('FT','IgB'), multiple = TRUE),
      selectInput("normalization_dropdown_protein", "normaliztion",
                  choices = c('none', 'median', 'plate_median'), selected = 'none'),
    ),

    # warning
    tags$div(
      class = "alert alert-warning",
      "Warning: normalization are performed separately for different assays."
    ),
    DT::dataTableOutput("protein_table"),
    downloadButton("download_protein", "Download Data"),
    downloadButton("download_protein_long", "Download Long Format Data")
  )
)

statisticalAnalysis <- tabPanel(
  title = "Statistical Analysis",
  fluidPage(
    titlePanel("T-test Analysis Between Cancer Types"),

    fluidRow(
      column(4,
             selectInput("group1_cancer_types", "Group 1 Cancer Types",
                         choices = NULL, multiple = TRUE,
                         selected = NULL),
             tags$small("Select cancer types for Group 1")
      ),
      column(4,
             selectInput("group2_cancer_types", "Group 2 Cancer Types",
                         choices = NULL, multiple = TRUE,
                         selected = NULL),
             tags$small("Select cancer types for Group 2"),
             br(),
             checkboxInput("auto_group2", "Auto-select Group 2",
                          value = FALSE),
             tags$small("When checked, Group 2 will be all cancer types not in Group 1")
      ),
      column(4,
             selectInput("stats_test_type", "Statistical Test Type",
                         choices = c("T-test" = "t_test", "Wilcoxon" = "wilcoxon"),
                         selected = "t_test"),
             selectInput("stats_assay_dropdown", "Assay",
                         choices = c('FT', 'IgB'), selected = 'FT'),
             selectInput("stats_normalization_dropdown", "Normalization",
                         choices = c('none', 'median', 'plate_median'), selected = 'none')
      )
    ),

    fluidRow(
      column(12,
             tags$div(
               class = "alert alert-info",
               "T-test will be performed for each peptide/modified peptide between the two selected groups.
                Log2 fold change is calculated as log2(mean_group1 / mean_group2)."
             )
      )
    ),

    fluidRow(
      column(6,
             actionButton("run_ttest", "Run T-test Analysis",
                         class = "btn-primary btn-lg"),
             br(), br()
      ),
      column(6,
             numericInput("pvalue_threshold", "P-value Threshold",
                         value = 0.05, min = 0.001, max = 1, step = 0.001)
      )
    ),

    fluidRow(
      column(12,
             tabsetPanel(
               tabPanel("Modified Peptide Level Results",
                       DT::dataTableOutput("ttest_modified_peptide_results"),
                       br(),
                       downloadButton("download_ttest_modified_peptide", "Download Results")
               ),
               tabPanel("Peptide Level Results",
                       DT::dataTableOutput("ttest_peptide_results"),
                       br(),
                       downloadButton("download_ttest_peptide", "Download Results")
               ),
               tabPanel("Protein Level Results",
                       DT::dataTableOutput("ttest_protein_results"),
                       br(),
                       downloadButton("download_ttest_protein", "Download Results")
               )
             )
      )
    )
  )
)

# Main UI using navbarPage to combine multiple pages
ui <- navbarPage(
  "DIANN PTM dashboard",
  ptmAnalysisTab,
  batchEffect,
  modifiedPeptideLevel,
  peptideLevel,
  proteinLevel,
  statisticalAnalysis
)
