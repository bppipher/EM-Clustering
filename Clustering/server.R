
library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)

source("clusterscripts.R")
source("EMscripts.R")

server <- function(input, output, session) {
  
  #########################################
  ######### Clustering            #########
  #########################################
  
  # Create a spot where we can store additional
  # reactive values for this session
  val <- reactiveValues(x = NULL, y = NULL,
                        xlimits = c(-10,10), ylimits = c(-10,10),
                        xlabels = "X", ylabels = "Y")
  
  # Count the number of points
  output$numPoints <- renderText({
    length(val$x)
  })
  # Reset the points on button click
  observeEvent(input$clear,{
    if (input$clear > 0){
      isolate({
        updata <- NULL
        compdata <- NULL
        factordata <- NULL
        val$x <- NULL
        val$y <- NULL
        val$xlimits <- c(-10,10)
        val$ylimits <- val$xlimits
        val$xlabels <- "X"
        val$ylabels <- "Y"
      })
    }
  })
  
  # Generate the plot of the clustered points
  # Generic R plot version
  output$clusterPlot <- renderPlot({
    
    tryCatch({
      # Try to cluster
      if (length(val$x) <= 1){
        stop("We can't cluster less than 2 points")
      }
      if(input$RUN){
        nclust = input$noclusters
        plotdata = data.matrix(compdata())
        EM = ExpectationMaximization(plotdata, nclust)
        x = as.data.frame(plotdata)
        Cluster = as.factor(EM$Cluster)
        levels(Cluster) = paste("Cluster", 1:nclust)
        x = cbind(x, Cluster)
        if(!is.na(input$factorcol)){
          Factor = updata()[,input$factorcol]
          plot(x = val$x,
               y = val$y,
               col = (seq(1:nclust)+1)[x$Cluster],
               pch = (seq(1:nlevels(Factor))+16)[updata()[,input$factorcol]],
               xlim = val$xlimits,
               ylim = val$ylimits,
               xlab = val$xlabels,
               ylab = val$ylabels)
          output$factorVScluster <- renderTable({
            if(!is.na(input$factorcol) && input$RUN){
              factor = input$factorcol
              table(Factor = updata()[,input$factorcol], Cluster = x$Cluster)
              
            }
          })
        } else {
          plot(x = val$x,
               y = val$y,
               col = (seq(1:nclust)+1)[x$Cluster],
               xlim = val$xlimits,
               ylim = val$ylimits,
               xlab = val$xlabels,
               ylab = val$ylabels)
        }
      } else {
        plot(val$x, val$y,
             xlab = val$xlabels, ylab = val$ylabels,
             xlim = val$xlimits, ylim = val$ylimits)
      }
    }, error = function(warn){
      # Otherwise just plot the points and instructions
      plot(val$x, val$y,
           xlim = val$xlimits, ylim = val$ylimits,
           xlab = val$xlabels, ylab = val$ylabels)
      text(0, 0, "Unable to create clusters.")
    }
    )
  })
  
  #########################################
  #########  File Uploading       #########
  #########################################
  avail.columns <- reactive(NULL)
  updata <- reactive({
    df = NULL
    if(input$upload == "upload"){
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      req(input$file1)
      
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df <- read.csv(input$file1$datapath,
                         header = input$header,
                         sep = input$sep,
                         quote = input$quote)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
    }
    if(input$upload == "preset"){
      df <- switch(input$dataset,
                   "None" = NULL,
                   "faithful" = faithful,
                   "iris (full)" = iris,
                   "iris (sepal only)" = iris[,-c(3,4)])
    }
    if(input$upload == "generate"){
      df <- NULL
      # df <- rmvmixnorm(n = input$pointstogen,
      #                  mus = ,
      #                  covs = ,
      #                  probs = ,
      #                  nclust = input$numclusts)
    }
    # Update select input immediately after clicking on the action button.
    if(!is.null(df)){
      updateNumericInput(session, 'xcol', value = 1, min = 1, max = ncol(df))
      updateNumericInput(session, 'ycol', value = if(ncol(df)>1){2}else{1}, min = 1, max = ncol(df))
      updateNumericInput(session, 'factorcol', value = NULL, min = 1, max = ncol(df))
    }
    
    return(df)
    
  })
  compdata <- reactive({
    retdata = NULL
    if(!is.na(input$factorcol)){
      factor = -input$factorcol
    } else {
      factor = TRUE
    }
    if(!is.null(updata())){
      if(input$jit){
        retdata = as.data.frame(lapply(updata()[,factor],jitter))
      } else {
        retdata = updata()[,factor]
      }
    }
    return(retdata)
  })
  
  output$contents <- renderTable({
    if(!is.na(input$factorcol)){
      factor = -input$factorcol
    } else {
      factor = TRUE
    }
    if(input$disp == "head") {
      return(head(updata())[,factor])
    }
    else {
      return(updata()[,factor])
    }
  })
  factordata = reactive({NULL})
  output$factors <- renderTable({
    if(!is.na(input$factorcol)){
      factordata = as.data.frame(updata()[,input$factorcol])
      colnames(factordata) = "Factors"
      if(input$disp == "head"){
        return(head(factordata))
      } else {
        return(factordata)
      }
    }
  })
  
  output$availabledata <- reactive({
    if(!is.null(updata())){
      "YES"
    }
  })
  outputOptions(output, "availabledata", suspendWhenHidden = FALSE)
  
  output$clustcovs <- renderUI({
    clustcovsnames <- c("X Variance","Y Variance","XY Covariance", "X Mean", "Y Mean", "Cluster Probability")
    lapply(seq(input$numclusts), function(j){
      sidebarPanel(#column(width=3,
        #lapply(seq(6),function(i){
        #numericInput(inputId = paste0("clustergenEntry_",i,"_Cluster_",j),label = clustcovsnames[i],value = i)
        #})
        numericInput(inputId = paste0("clustergenEntry_","XVar","_Cluster_",j),
                     label = "X Variance", value = 1),
        numericInput(inputId = paste0("clustergenEntry_","YVar","_Cluster_",j),
                     label = "Y Variance", value = 1),
        numericInput(inputId = paste0("clustergenEntry_","XYCov","_Cluster_",j),
                     label = "XY Covariance", value = 0),
        numericInput(inputId = paste0("clustergenEntry_","XMean","_Cluster_",j),
                     label = "X Mean", value = 0),
        numericInput(inputId = paste0("clustergenEntry_","YMean","_Cluster_",j),
                     label = "Y Mean", value = 0),
        numericInput(inputId = paste0("clustergenEntry_","Prob","_Cluster_",j),
                     label = "Cluster Probability", value = 1/input$numclusts)
      )
    })
  })
  EMData = reactive({NULL})
  observeEvent(input$pushtoclust,{
    if(input$pushtoclust > 0 && !is.null(compdata())){
      isolate({
        newdata <- data.matrix(compdata())
        val$x <- newdata[ , input$xcol]
        val$y <- newdata[ , input$ycol]
        val$xlimits <- c(min(val$x)*1.04, max(val$x)*1.04)
        val$ylimits <- c(min(val$y)*1.04, max(val$y)*1.04)
        if(!is.null(colnames(newdata)[input$xcol])){
          val$xlabels <- colnames(newdata)[input$xcol]
        } else{
          val$xlabels <- "X"
        }
        if(!is.null(colnames(newdata)[input$ycol])){
          val$ylabels <- colnames(newdata)[input$ycol]
        } else {
          val$ylabels <- "Y"
        }
      })
    }
  }
  )
}

return(server)

