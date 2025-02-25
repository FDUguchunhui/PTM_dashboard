# Load required libraries
library(shiny)
library(plotly)

# Define UI
ui <- fluidPage(
  titlePanel("PTM Analysis in Shiny"),
  
  # Sidebar (Filters) at the Top
  fluidRow(
    column(4,
           selectInput("file_dropdown", "Select PTM File",
                       choices = c("Citrullination", "Hypusine", "Deoxyhypusine"),
                       selected = "Citrullination")
    ),
  ),
  
  # Main Panel with Plots
  
  fluidRow(
    column(12, plotlyOutput("matplotlib_plot")),
  ),
  
  fluidRow(
    column(12, plotlyOutput('prop_plot_case_control')),
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
  )
)
