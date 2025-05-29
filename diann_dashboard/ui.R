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
    h2("This is another page of your dashboard"),
    p("Add more content here.")
  )
)

# Main UI using navbarPage to combine multiple pages
ui <- navbarPage(
  "DIANN PTM dashboard",
  ptmAnalysisTab,
  batchEffect,
  modifiedPeptideLevel,
  peptideLevel,
  proteinLevel
)
