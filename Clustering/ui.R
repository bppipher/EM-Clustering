#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#



library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)

# Create a bootstrap fluid layout
ui <- navbarPage("Dynamic EM Clustering",
                 tabPanel("Clustering",
                          fluidPage(
                            
                            #App Title ----
                            titlePanel("Dynamic Expectation-Maximization Clustering plot"),
                            
                            # Main Row for displaying plot ----
                            fluidRow(
                              
                              plotOutput(
                                "clusterPlot", "100%", "600px"
                              )
                            ),
                            # Row layout with input and output ----
                            fluidRow(
                              column(2,
                                     actionButton("clear",
                                                  "Clear Points",
                                                  icon = icon("fas fa-eraser"))
                              ),
                              column(2,
                                     checkboxInput("RUN","Check to run EM", value = FALSE)
                              )
                            ),
                            fluidRow(
                              column(3,
                                     wellPanel("Number of Points: ", verbatimTextOutput("numPoints", placeholder = TRUE))
                              ),
                              column(2,
                                     wellPanel("Number of Clusters?",numericInput("noclusters", label = NULL, value = 2, min = 2, max = 10, step = 1))
                              ),
                              tableOutput("factorVScluster")
                            )
                          )
                 ),
                 tabPanel("Data",
                          fluidPage(
                            
                            # App title ----
                            titlePanel("Input Data"),
                            
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                
                                # Input: Method to add data ----
                                radioButtons("upload", "Which method to input data?",
                                             choices = c(Preset = "preset",
                                                         Upload = "upload",
                                                         Generate = "generate"),
                                             selected = NULL),
                                
                                # Condition: Chose to upload
                                conditionalPanel(condition = "input.upload == 'upload'",
                                                 # Horizontal line ----
                                                 tags$hr(),
                                                 # Input: Select a file ----
                                                 fileInput("file1", "Choose CSV File",
                                                           multiple = FALSE,
                                                           accept = c("text/csv",
                                                                      "text/comma-separated-values,text/plain",
                                                                      ".csv")
                                                 ),
                                                 # Horizotnal line ----
                                                 tags$hr(),
                                                 # Input: Checkbox if file has header
                                                 checkboxInput("header", "Header", TRUE),
                                                 # Input: Select separator
                                                 radioButtons("sep", "Separator",
                                                              choices = c(Comma = ",",
                                                                          Semicolon = ";",
                                                                          Tab = "\t"),
                                                              selected = ","),
                                                 # Input: Select quote
                                                 radioButtons("quote", "Quote",
                                                              choices = c(None = "",
                                                                          "Double Quote" = '"',
                                                                          "Single Quote" = "'"),
                                                              selected = '"')
                                ),
                                
                                # Condition: Chose to use preset data
                                conditionalPanel(condition = "input.upload == 'preset'",
                                                 # Horizontal line
                                                 tags$hr(),
                                                 # Input: Choose dataset ----
                                                 selectInput("dataset", "Or choose a dataset:",
                                                             choices = c("None",
                                                                         "faithful",
                                                                         "iris (full)",
                                                                         "iris (sepal only)")
                                                 )
                                                 
                                ),
                                # Input: Select Columns ----
                                conditionalPanel(condition = "input.upload !='Generate' && output.availabledata == 'YES'",
                                                 numericInput('xcol', 'Enter desired x column:', NA),
                                                 numericInput('ycol', 'Enter desired y column:', NA),
                                                 numericInput('factorcol', 'Enter Factor column (if present):', NA)
                                ),
                                
                                # Condition: Chose to generate data
                                conditionalPanel(condition = "input.upload == 'generate'",
                                                 # Horizontal line
                                                 tags$hr(),
                                                 numericInput("numclusts",
                                                              "Input No. of Clusters",
                                                              value = 2,
                                                              min = 2,
                                                              max = 10),
                                                 numericInput("pointstogen",
                                                              "Input No. of Points to generate",
                                                              value = 100,
                                                              min = 1,
                                                              max = 10000)
                                ),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                conditionalPanel(condition = "output.availabledata == 'YES'",
                                                 # Input: Choose to push data to cluster plot
                                                 radioButtons("disp", "Display",
                                                              choices = c(Head = "head",
                                                                          All = "all"),
                                                              selected = "head"),
                                                 # Horizontal line ---
                                                 tags$hr(),
                                                 # Input: Choose to push data to cluster plot
                                                 checkboxInput("jit", "Jitter data?", value = FALSE),
                                                 actionButton("pushtoclust", "Export to Cluster Plot?",
                                                              icon = icon("fas fa-angle-double-right"))
                                )
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(
                                
                                conditionalPanel(condition = "input.upload == 'generate'", 
                                                 # Input: Creating matrix arguments for clusters
                                                 uiOutput("clustcovs")
                                ),
                                
                                # Output: Data file ----
                                splitLayout(cellWidths = c("70%", "30%"),
                                  tableOutput("contents"),
                                  tableOutput("factors")
                                )
                              )
                              
                            )
                          )
                 ),
                 tabPanel("How-to",
                          fluidPage(
                            strong("Instructions:"),
                            p("First head to `Data` and either upload a csv or use a preset.
                              Be sure to properly set the csv parameters before uploading as if
                              you attempt to load a csv with different formatting than specified
                              the app may fail.  (eg, uploading data without headers when 'Header' is checked.)
                              It is also worth noting that generating Mixed Bivariate Normal data is not currently
                              working, mostly due to limitations to Shiny's \"Input\" interfaces, although the functions to do Mixed MVN data generation are present in the code.",
                              style = "font-family: 'times'; font-si16pt"),
                            p("Ensure that all data is numerical, or that you only have a single factor included.
                              If you have a factor be sure to mark it as such or the applet will behave unexpectedly.
                              Only one factor is supported in the app.  If you do not mark a factor, or have non-numerical data,
                              the app may fail.  Uploading a factor is not necessary, but if you do so the app will compare
                              the cluster levels to the factor levels using a frequency table.  
                              The app uses all available data, and works off of a Mixed Multivariate Normal Expectation-Maximization algorithm utilizing Bayes rule, and is not clustering strictly in two dimensions despite the plot only supporting two dimensions.",
                              style = "font-family: 'times'; font-si16pt"),
                            p("Once data has been loaded, and you've been sure to mark the factor if present,
                              you can click the `Export to Cluster Plot?` button to send it to the Clustering tab.
                              An optional setting is to jitter the data, useful for datasets like `iris` which are rounded.  Jittering data can help the convergence
                              and stability of the app.",
                              p("From here you can head to the Clustering tab, select the number of clusters
                                you're interested in fitting towards and then select the checkbox to run.  The generated plot will color code the points to
                                the cluster each point had the highest responsibility in at convergence.  If a factor was selected a frequency table will display.
                                ",
                                style = "font-family: 'times'; font-si16pt"),
                              p("A good demo would be to upload only two columns of `iris` manually, with factors, and attempt to fit 3 clusters.
                                Note the EM algorithms accuracy (or lack thereof) in correctly classifying the flower species into clusters with only a two-dimensional viewpoint.
                                You can then use the entire `iris` data set, with the factor, and compare how the EM algorithm then has notable success in the higher-dimension setting.  The algorithm makes a final class assignment by rounding, and uses the computed Mean and Variance-Covariance matrix of the data as an initialization.",
                                style = "font-family: 'times'; font-si16pt"),
                              p("When changing datasets it is recommended to uncheck `Check to run EM` before pushing new dataset to the Cluster plot.",
                                style = "font-family: 'times'; font-si16pt"),
                              style = "font-family: 'times'; font-si16pt"),
                            p("Shiny has a known bug where you can override NumericInput's minimum and maximum, which can also cause the app to fail.",
                              style = "font-family: 'times'; font-si16pt"),
                            p("The app was originally written using GGplot and Plot_ly but they introduced noticable lag to the app once uploaded to shinyapps.io",
                              style = "font-family: 'times'; font-si16pt")
                            
                          )
                 )
)
return(ui)

