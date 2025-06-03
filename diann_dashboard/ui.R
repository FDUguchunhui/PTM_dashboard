library(shiny)
library(plotly)
library(bslib)

source('../R/PTM_definition.R')

# Custom CSS for styling and collapsible sidebar
custom_css <- "
<style>
  /* Main styling */
  .main-header {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 20px;
    margin-bottom: 20px;
    border-radius: 10px;
    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
  }

  .main-header h1 {
    margin: 0;
    font-weight: 300;
    font-size: 2.5em;
  }

  /* Sidebar styling */
  .custom-sidebar {
    background: linear-gradient(180deg, #f8f9fa 0%, #e9ecef 100%);
    border-radius: 10px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    padding: 20px;
    margin-right: 15px;
    transition: all 0.3s ease;
    min-height: 100vh;
  }

  .sidebar-collapsed {
    margin-left: -280px;
  }

  .sidebar-toggle {
    position: fixed;
    top: 20px;
    left: 20px;
    z-index: 1000;
    background: #667eea;
    color: white;
    border: none;
    border-radius: 50%;
    width: 50px;
    height: 50px;
    font-size: 18px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
    transition: all 0.3s ease;
  }

  .sidebar-toggle:hover {
    background: #5a67d8;
    transform: scale(1.1);
  }

  .main-content {
    transition: all 0.3s ease;
    padding-left: 20px;
  }

  .main-content.expanded {
    margin-left: -280px;
  }

  /* Card styling for plots */
  .plot-card {
    background: white;
    border-radius: 10px;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    padding: 20px;
    margin-bottom: 20px;
    transition: transform 0.2s ease;
  }

  .plot-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
  }

  /* Tab styling */
  .nav-tabs .nav-link {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 8px 8px 0 0;
    margin-right: 5px;
    transition: all 0.3s ease;
  }

  .nav-tabs .nav-link.active {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    border-color: #667eea;
  }

  .nav-tabs .nav-link:hover {
    background: rgba(102, 126, 234, 0.1);
  }

  /* Form controls styling */
  .form-control, .form-select {
    border-radius: 8px;
    border: 2px solid #e9ecef;
    transition: border-color 0.3s ease;
  }

  .form-control:focus, .form-select:focus {
    border-color: #667eea;
    box-shadow: 0 0 0 0.2rem rgba(102, 126, 234, 0.25);
  }

  /* Button styling */
  .btn-primary {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    border: none;
    border-radius: 8px;
    padding: 10px 20px;
    font-weight: 500;
    transition: all 0.3s ease;
  }

  .btn-primary:hover {
    transform: translateY(-1px);
    box-shadow: 0 4px 8px rgba(102, 126, 234, 0.3);
  }

  /* Alert styling */
  .alert {
    border-radius: 8px;
    border: none;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
  }

  .alert-warning {
    background: linear-gradient(135deg, #ffeaa7 0%, #fdcb6e 100%);
    color: #2d3436;
  }

  .alert-info {
    background: linear-gradient(135deg, #74b9ff 0%, #0984e3 100%);
    color: white;
  }

  /* Responsive design */
  @media (max-width: 768px) {
    .custom-sidebar {
      margin-left: -100%;
    }

    .main-content.expanded {
      margin-left: 0;
    }

    .sidebar-toggle {
      display: block;
    }
  }
</style>
"

# JavaScript for sidebar toggle functionality
sidebar_js <- "
<script>
$(document).ready(function() {
  var sidebarCollapsed = false;

  $('#sidebar-toggle').click(function() {
    if (sidebarCollapsed) {
      $('.custom-sidebar').removeClass('sidebar-collapsed');
      $('.main-content').removeClass('expanded');
      $(this).html('<i class=\"fas fa-bars\"></i>');
    } else {
      $('.custom-sidebar').addClass('sidebar-collapsed');
      $('.main-content').addClass('expanded');
      $(this).html('<i class=\"fas fa-times\"></i>');
    }
    sidebarCollapsed = !sidebarCollapsed;
  });
});
</script>
"

# First page: PTM Analysis tab
ptmAnalysisTab <- tabPanel(
  title = "Overview",
  div(class = "plot-card",
    titlePanel("PTM Analysis Dashboard"),
  ),

  # Main Panel with Plots
  div(class = "plot-card",
    plotlyOutput("intensity_plot")
  ),

  div(class = "plot-card",
    plotlyOutput('pep_prop_plot_case_control')
  ),

  div(class = "plot-card",
    plotlyOutput('prot_prop_plot_case_control')
  ),

  div(class = "plot-card",
    fluidRow(
      column(4,
             sliderInput("bins_slider", "Number of Bins", min = 10, max = 100, value = 30, step = 1))
    ),
    plotlyOutput("intensity_plot_cancer_type")
  ),

  div(class = "plot-card",
    plotlyOutput('pep_prop_plot_cancer_type')
  ),

  div(class = "plot-card",
    plotlyOutput('prot_prop_plot_cancer_type')
  )
)

# Additional pages defined as tabPanel
batchEffect <- tabPanel(
  title = "Batch Effect",
  div(class = "plot-card",
    plotlyOutput("total_plate_batch_effect_plot", height = "400px")
  ),
  div(class = "plot-card",
    plotlyOutput("plate_batch_effect_plot", height = "800px")
  ),
  div(class = "plot-card",
    plotlyOutput("violin_plot_conditions", height = "800px")
  )
)

modifiedPeptideLevel <- tabPanel(
  title = "Modified Peptide Level",
  div(class = "plot-card",
    tags$div(
      class = "alert alert-warning",
      tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 8px;"),
      "Warning: Normalization is performed separately for different assays."
    ),
    DT::dataTableOutput("modified_peptide_table"),
    br(),
    div(class = "text-center",
      downloadButton("download_modified_peptide", "Download Data", class = "btn-primary"),
      downloadButton("download_modified_peptide_long", "Download Long Format", class = "btn-primary"),
      downloadButton("download_metadata", "Download Metadata", class = "btn-primary")
    )
  )
)

peptideLevel <- tabPanel(
  title = "Peptide Level",
  div(class = "plot-card",
    tags$div(
      class = "alert alert-warning",
      tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 8px;"),
      "Warning: Normalization is performed separately for different assays."
    ),
    DT::dataTableOutput("peptide_table"),
    br(),
    div(class = "text-center",
      downloadButton("download_peptide", "Download Data", class = "btn-primary"),
      downloadButton("download_peptide_long", "Download Long Format", class = "btn-primary")
    )
  )
)

proteinLevel <- tabPanel(
  title = "Protein Level",
  div(class = "plot-card",
    tags$div(
      class = "alert alert-warning",
      tags$i(class = "fas fa-exclamation-triangle", style = "margin-right: 8px;"),
      "Warning: Normalization is performed separately for different assays."
    ),
    DT::dataTableOutput("protein_table"),
    br(),
    div(class = "text-center",
      downloadButton("download_protein", "Download Data", class = "btn-primary"),
      downloadButton("download_protein_long", "Download Long Format", class = "btn-primary")
    )
  )
)

statisticalAnalysis <- tabPanel(
  title = "Statistical Analysis",
  div(class = "plot-card",
    div(class = "main-header",
      h3("T-test Analysis Between Cancer Types", style = "margin: 0; color: white;")
    ),

    fluidRow(
      column(4,
             selectInput("group1_cancer_types", "Group 1 Cancer Types",
                         choices = NULL, multiple = TRUE,
                         selected = 'Breast'),
             tags$small("Select cancer types for Group 1", class = "text-muted")
      ),
      column(4,
             selectInput("group2_cancer_types", "Group 2 Cancer Types",
                         choices = NULL, multiple = TRUE,
                         selected = 'Gastric'),
             tags$small("Select cancer types for Group 2", class = "text-muted"),
             br(), br(),
             checkboxInput("auto_group2", "Auto-select Group 2",
                          value = FALSE),
             tags$small("When checked, Group 2 will be all cancer types not in Group 1", class = "text-muted")
      ),
      column(4,
             selectInput("stats_test_type", "Statistical Test Type",
                         choices = c("T-test" = "t_test", "Wilcoxon" = "wilcoxon"),
                         selected = "t_test"),
             selectInput("stats_assay_dropdown", "Assay",
                         choices = c('FT', 'IgB'), selected = 'FT')
      )
    ),

    tags$div(
      class = "alert alert-info",
      tags$i(class = "fas fa-info-circle", style = "margin-right: 8px;"),
      "T-test will be performed for each peptide/modified peptide between the two selected groups.
       Log2 fold change is calculated as log2(mean_group1 / mean_group2)."
    ),

    fluidRow(
      column(6,
             actionButton("run_test", "Run Analysis",
                         class = "btn-primary btn-lg"),
             br(), br()
      )

    ),

    # Filtering controls
    fluidRow(
      column(3,
             sliderInput("auc_filter", "Abs AUC from 0.5",
                        min = 0, max = 0.5, value = 0, step = 0.01),
             tags$small("Filter results by Area Under Curve (AUC) values", class = "text-muted")
      ),
      column(3,
             sliderInput("log2fc_filter", "Log2 Fold Change Filter",
                        min = -10, max = 10, value = c(-10, 10), step = 0.1),
             tags$small("Filter results by log2 fold change values", class = "text-muted")
      ),

      column(3,
             sliderInput("pvalue_threshold", "P-value Threshold",
                          value = 0.05, min = 0.001, max = 0.1, step = 0.001),
             tags$small("Filter results by p-value", class = "text-muted")
      ),

      column(3,
             br(),
             checkboxInput("show_significant_only", "Show significant results only",
                          value = FALSE),
             tags$small("Show only results below p-value threshold", class = "text-muted")
      )
    ),

    tabsetPanel(
      tabPanel("Modified Peptide Level Results",
              DT::dataTableOutput("ttest_modified_peptide_results"),
              br(),
              downloadButton("download_ttest_modified_peptide", "Download Results", class = "btn-primary")
      ),
      tabPanel("Peptide Level Results",
              DT::dataTableOutput("ttest_peptide_results"),
              br(),
              downloadButton("download_ttest_peptide", "Download Results", class = "btn-primary")
      ),
      tabPanel("Protein Level Results",
              DT::dataTableOutput("ttest_protein_results"),
              br(),
              downloadButton("download_ttest_protein", "Download Results", class = "btn-primary")
      )
    )
  )
)

# Main UI with collapsible sidebar
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bg = "#ffffff",
    fg = "#2c3e50",
    primary = "#667eea",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter")
  ),

  # Include custom CSS and JavaScript
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"),
    HTML(custom_css),
    HTML(sidebar_js)
  ),

  # Sidebar toggle button
  tags$button(
    id = "sidebar-toggle",
    class = "sidebar-toggle",
    HTML('<i class="fas fa-bars"></i>')
  ),

  # Main header
  div(class = "main-header",
    h1("DIANN PTM Dashboard", style = "margin: 0;"),
    p("Advanced Post-Translational Modification Analysis Platform", style = "margin: 5px 0 0 0; opacity: 0.9;")
  ),

  fluidRow(
    # Sidebar
    column(3,
      div(class = "custom-sidebar",
        h4("Global Settings", style = "color: #667eea; margin-bottom: 20px;"),
        tags$hr(style = "border-color: #667eea;"),

        # PTM File Selection and Filters
        div(style = "margin-bottom: 15px;",
          selectInput("file_dropdown", "Select PTM File",
                      choices = names(TARGET_PTM),
                      selected = "Citrullination_test")
        ),

        div(style = "margin-bottom: 15px;",
          sliderInput("PTM_qvalue_slider", "PTM Q-value threshold",
                      min = 0, max = 0.5, value = 0.1, step = 0.01)
        ),

        div(style = "margin-bottom: 15px;",
          sliderInput("PTM_site_confidence_slider", "Location confidence",
                      min = 0, max = 1, value = 0, step = 0.01)
        ),

        tags$hr(style = "border-color: #e9ecef;"),

        # Shared dropdowns that will be used across tabs
        div(style = "margin-bottom: 15px;",
          selectInput("shared_cancer_type_dropdown", "Select Cancer Types",
                      choices = NULL,
                      selected = c('Breast', 'MERIT Control'),
                      multiple = TRUE)
        ),

        div(style = "margin-bottom: 15px;",
          selectInput("shared_assay_dropdown", "Assay",
                      choices = c('FT', 'IgB'),
                      selected = c('FT','IgB'),
                      multiple = TRUE)
        ),

        div(style = "margin-bottom: 15px;",
          selectInput("shared_normalization_dropdown", "Normalization",
                      choices = c('none', 'median', 'plate_median'),
                      selected = 'none')
        ),

        tags$hr(style = "border-color: #e9ecef;"),
        tags$small("These settings apply to Modified Peptide Level, Peptide Level, and Protein Level tabs.",
                  class = "text-muted")
      )
    ),

    # Main content
    column(9,
      div(class = "main-content",
        navbarPage(
          "",
          ptmAnalysisTab,
          batchEffect,
          modifiedPeptideLevel,
          peptideLevel,
          proteinLevel,
          statisticalAnalysis
        )
      )
    )
  )
)
