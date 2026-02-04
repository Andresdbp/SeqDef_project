#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(viridis)
library(ape)
library(phytools)

# Define UI 
ui <- fluidPage(
    # Application title
    titlePanel("Visualizing Sequencing Distribution"),
    # Sidebar
    sidebarLayout(
        sidebarPanel(
          fileInput("file", label = h3("File input"),
                    accept=c('csv', 'comma-separated-values','.csv')),
          fluidRow(
            column(12, uiOutput("taxlabel"))
          ),
          fluidRow(
            column(6, uiOutput("taxstart")),
            
            column(6, uiOutput("taxend")),
          ),
          uiOutput("select"),
          uiOutput("levels"),
          fluidRow(
            column(12, uiOutput("datlabel"))
          ),
          fluidRow(
            column(6, uiOutput("datstart")),
            column(6, uiOutput("datend")),
          ),
          uiOutput("trait"),
          uiOutput("handleNA"),
          br(),
          uiOutput("plot"),
          hr(),
          verbatimTextOutput("messageNA")
        ),
        # Plot
        mainPanel(
          tags$head(tags$style("#treePlot{height:95vh !important;width:80vh !important;}")),
          tableOutput("message"),
          plotOutput("treePlot")
        )
    )
)

# Define Server logic
server <- function(input, output) {
  # Functions
  input_file <- reactive({
    if (is.null(input$file)) {
      return("")
    }
    return(read.csv(file = input$file$datapath, as.is=F))
  })
  taxstart <- reactive({
    if (is.na(input$taxstart)) {
      return("")
    }
    return(input$taxstart)
  })
  taxend <- reactive({
    if (is.null(input$taxend)) {
      return(NA)
    }
    return(input$taxend)
  })
  datstart <- reactive({
    if (is.na(input$datstart)) {
      return("")
    }
    return(input$datstart)
  })
  datend <- reactive({
    if (is.na(input$datend)) {
      return("")
    }
    return(input$datend)
  })
  tax.name <- reactive({
    return(input$levels)
    })
  levels <- reactive({
    dat <- input_file()
    choices <- c()
    colstart <- last_selection()
    colend <- taxend()
    counter <- 1
    for (i in colstart:colend) {
      if (length(unique(dat[,i])) > 1) {
        choices[counter] <- colnames(dat)[i]
        counter <- counter + 1
      }
    }
    return(choices)
    })
  traits <- reactive({
    dat <- input_file()
    choices <- c()
    colstart <- datstart()
    colend <- datend()
    counter <- 1
    for (i in colstart:colend) {
      if (length(unique(dat[,i])) > 1) {
        choices[counter] <- colnames(dat)[i]
        counter <- counter + 1
      }
    }
    return(choices)
  })
  dat.name <- reactive({ input$trait})
  handleNA <- reactive({ input$handleNA })
  parse <- function(i) {
    dat <- input_file()
    # Get selected values from previous selectInputs
    selected_values <- sapply(seq(input$taxstart, i - 1), function(j) {
      input[[paste0("select_", j)]]
    })
    # Check if selected values are not NULL or empty
    if (!all(sapply(selected_values, function(val) is.null(val) || val == "All"))) {
      # Filter rows based on selected values
      filtered_data <- dat
      for (j in seq(input$taxstart, i - 1)) {
        if (!is.null(input[[paste0("select_", j)]]) && input[[paste0("select_", j)]] != "All") {
          filtered_data <- filtered_data[filtered_data[, j] == input[[paste0("select_", j)]], ]
        }
      }
      return(as.character(sort(unique(filtered_data[, i]))))
    } else {
      # Return all unique values if selected values are NULL or empty
      return(as.character(sort(unique(dat[, i]))))
    }
  }
  filtered_data <- reactive({
    dat <- input_file()
    # Filter data based on selected values from selectInputs
    for (i in seq(input$taxstart, input$taxend)) {
      if (input[[paste0("select_", i)]] != "All") {
        dat <- dat[dat[, i] == input[[paste0("select_", i)]], ]
      }
    }
    return(dat)
  })
  last_selection <- reactive({
    last_non_all <- 1
    for (i in seq(input$taxstart, (input$taxend))) {
      if (input[[paste0("select_", i)]] != "All") {
        last_non_all <- i + 1
      }
    }
    if (last_non_all == 1) {
      last_non_all <- 2
    }
    if (last_non_all >= input$taxend) {
      last_non_all <- input$taxend
    }
    last_non_all <- pmin(last_non_all, input$taxend)
    return(last_non_all)
  })
  
  # Logic behind UI elements
  observeEvent(
    eventExpr = input[["file"]],
    handlerExpr = {
      req(input_file())
      dat <- input_file()
      output$taxlabel <- renderUI({
        h3("Taxonomy Columns")
      })
      output$taxstart <- renderUI({
        numericInput("taxstart", label = "First Column", value = 1)
      })
      output$taxend <- renderUI({
        numericInput("taxend", label = "Last Column", value = taxend())
      })
    }
  )
  observeEvent(
    eventExpr = input[["taxend"]],
    handlerExpr = {
      req(input$taxend)
      dat <- input_file()
      output$select = renderUI({
        lapply(seq(input$taxstart, input$taxend), function(i) {
          selected_value <- if (!is.null(input[[paste0("select_", i)]])) input[[paste0("select_", i)]] else "All"
          selectInput(
            inputId = paste0("select_", i),
            label = colnames(dat)[i],
            choices = c("All", parse(i)),
            selected = selected_value
          )
        })
      })
    }
  )
  
  observeEvent(
    eventExpr = input[["select_1"]],
    handlerExpr = {
      req(levels())
      output$levels <- renderUI({
        selectInput("levels", label = "Select level to plot:"
                    , choices = levels())
      })
    }
  )
  observeEvent(
    eventExpr = input[["levels"]],
    handlerExpr = {
      req(levels())
      output$datlabel <- renderUI({
        h3("Data Columns")
      })
      output$datstart <- renderUI({
        numericInput("datstart", label = "First Column", value = (taxend() + 1))
      })
      output$datend <- renderUI({
        numericInput("datend", label = "Last Column", value = (taxend() + 1))
      })
    }
  )
  observeEvent(
    eventExpr = input[["datend"]],
    handlerExpr = {
      req(datend())
      output$trait <- renderUI({
        selectInput("trait", label = h3("Select Trait")
                    , choices = traits())
      })
    }
  )
  observeEvent(
    eventExpr = input[["trait"]],
    handlerExpr = {
      req(traits())
      output$handleNA <- renderUI({
        radioButtons("handleNA", label = h3("How to handle missing data:"),
                     choices = list("Remove NA" = 1, "Turn NA into 0" = 2), 
                     selected = 1)
      })
    }
  )
  observeEvent(
    eventExpr = input[["trait"]],
    handlerExpr = {
      req(traits())
      output$plot <- renderUI({
        actionButton(
          inputId = "plot",
          label = "Generate Plot"
        )
      })
    }
  )
  observeEvent(
    eventExpr = input[["plot"]],
    handlerExpr = {
      dat <- filtered_data()
      tax <- tax.name()
      trait <- dat.name()
      tax.index <- which(colnames(dat) == tax)
      trait.index <- which(colnames(dat) == trait)
      handleNA <- handleNA()
      if (length(unique(dat[,tax.index])) > 1) {
        dataNA <- sum(is.na(dat[,trait.index]))
        if (handleNA == 1) { # Remove NA
          dat <- dat[!is.na(dat[,trait.index]),]
          actionNA <- "removed"
        }
        if (handleNA == 2) { # Make NA into 0
          dat[is.na(dat[,trait.index]),trait.index] <- 0
          actionNA <- "turned into 0"
        }
        dat[,tax.index] <- as.character(dat[,tax.index])
        counts <- as.data.frame(table(dat[,tax.index]))
        N <- counts$Freq
        names(N) <- counts$Var1
        set.seed(1)
        frm <- reformulate(paste(colnames(dat)[1:tax.index],collapse="/"))
        dat_unique <- do.call(rbind, lapply(unique(dat[,tax.index]), function(value) {
          dat[dat[,tax.index] == value, ][1, ]
        }))
        for (i in 1:length(dat_unique[,tax.index])) {
          dat_unique[,trait.index][i] <- mean(dat[,trait.index][which(dat[,tax.index] == dat_unique[i,tax.index])])
        }
        dat_unique[,tax.index] <- as.factor(dat_unique[,tax.index])
        tr1 <- as.phylo(frm, data = dat_unique, collapse=TRUE)
        tr1$edge.length <- rep(1, nrow(tr1$edge))
        tr1 <- force.ultrametric(tr1, method="extend")
        tip.label <- tr1$tip.label
        clade.label <- tr1$tip.label
        tr2 <- multi2di(tr1)
        tr2$edge.length[tr2$edge.length==0] <- 0.0001
        depth <- sapply(tr2$tip.label, function(x, y) {
          0.6 * y$edge.length[which(tr2$edge[, 2] == which(y$tip.label == x))]
        }, y = tr2)
        N.sort <- c()
        for (i in 1:length(tip.label)) {
          N.sort[i] <- N[names(N) == tip.label[i]]
        }
        N.sort <- N.sort*2
        names(N.sort) <- tip.label
        trans<-data.frame(tip.label,clade.label,N.sort,depth)
        tr2 <- phylo.toBackbone(tr2, trans)
        
        # Plot the tree
        output$treePlot <- renderPlot({
          layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))
          cols <- viridis(100, option="A")[round(1 + ((dat_unique[,trait.index] - min(dat_unique[,trait.index])) / (max(dat_unique[,trait.index]) - min(dat_unique[,trait.index]))) * 99)]
          if (max(dat_unique[,trait.index]) == min(dat_unique[,trait.index])) {
            cols <- viridis(100, option="A")[rep(x = 50, times = length(dat_unique[,trait.index]))]
          }
          par(mar = rep(0, 4))
          plot(tr2, col=cols)
          colfunc <- colorRampPalette(viridis(100, option="A"))
          legend_image <- as.raster(matrix(colfunc(100), ncol=1))
          par(mar = rep(0.2, 4))
          plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
          text(x=0.9, y = 0.8, labels = paste(colnames(dat_unique)[trait.index]), cex = 1)
          text(x=1.3, y = seq(0.75,0.25,l=5), labels = signif(seq(max(dat_unique[,trait.index]), min(dat_unique[,trait.index]), l=5), digits = 2))
          rasterImage(legend_image, 0.5, 0.75, 1, 0.25)
        })
        if (dataNA < 1) {
          output$messageNA <- renderText({})
        }
        if (dataNA == 1) {
          output$messageNA <- renderPrint({paste(dataNA, "data point was", actionNA)})
        }
        if (dataNA > 1) {
          output$messageNA <- renderPrint({paste(dataNA, "data points were", actionNA)})
        }
        output$message <- renderText({})
      }
      if (length(unique(dat[,tax.index])) == 1) {
        output$message <- renderText({paste0("Single data point selected: ", dat[,tax.index], ". ", colnames(dat)[trait.index], ": " , dat[,trait.index])})
        output$treePlot <- renderPlot({})
      }
    }
  )
}
  



# Run the application 
shinyApp(ui = ui, server = server)
