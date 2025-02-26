library(shiny)
library(plotly)

# First page: PTM Analysis tab
ptmAnalysisTab <- tabPanel(
  title = "Overview",
  fluidPage(
    titlePanel("PTM Analysis in Shiny"),
    
    # Sidebar (Filters) at the Top
    fluidRow(
      column(4,
             selectInput("file_dropdown", "Select PTM File",
                         choices = c("Citrullination_test", "Citrullination", "Hypusine", "Deoxyhypusine"),
                         selected = "Citrullination_test")
      )
    ),
    
    fluidRow(
      column(12, 
             sliderInput("PTM_qvalue_slider", "PTM Q-value threshold", min=0.9, max=1, 
                         value=0.95, step = 0.01))
    ),
    
    # Main Panel with Plots
    fluidRow(
      column(12, plotlyOutput("matplotlib_plot"))
    ),
    
    fluidRow(
      column(12, plotlyOutput('pep_prop_plot_case_control'))
    ),
    
    fluidRow(
      column(12, plotlyOutput('prot_prop_plot_case_control'))
    ),
    
    fluidRow(
      column(8,
             selectInput("cancer_type_dropdown", "Select Cancer Types", choices = NULL, multiple = TRUE)
      ),
      column(4,
             sliderInput("bins_slider", "Number of Bins", min = 10, max = 100, value = 30, step = 1)
      )
    ),
    
    fluidRow(
      column(12, plotlyOutput("histogram_plot"))
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
    h2("This is another page of your dashboard"),
    p("Add more content here.")
  )
)

peptideLevel <- tabPanel(
  title = "Peptide level",
  fluidPage(
    h2("This is another page of your dashboard"),
    p("Add more content here.")
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
  peptideLevel,
  proteinLevel
)
