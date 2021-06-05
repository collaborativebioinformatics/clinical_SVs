library(shiny)
library(tippy)
library(shinythemes)
library(tidyverse)
library(tidygraph)

sv <- function(){

  ## head
  head <- tagList(
    tags$link(
      href = "assets/pushbar.css",
      rel="stylesheet",
      type="text/css"
    ),
    tags$link(
      href = "assets/custom.css",
      rel = "stylesheet",
      type = "text/css"
    ),
    tags$script(
      src = "assets/pushbar.js"
    ),
    tags$link(
      href = "assets/please-wait.css",
      rel = "stylesheet",
      type = "text/css"
    ),
    tags$script(
      src = "assets/please-wait.min.js"
    ),
    tags$script(
      src = "https://unpkg.com/micromodal/dist/micromodal.min.js"
    ),
    tags$link(
      rel = "stylesheet",
      href = "https://use.fontawesome.com/releases/v5.7.2/css/all.css",
      integrity = "sha384-fnmOCqbTlWIlj8LyTjo7mOUStjsKC4pOpQbqyi7RrhN7udi9RwhKkMHpvLbHG9Sr",
      crossorigin = "anonymous"
    ),
    tags$script(
      src = "assets/custom.js"
    ),
    tags$style(
      paste0(".pushbar{background-color:", "white", ";}")
    )
  )

   # particles_json <- jsonlite::fromJSON(
   #  system.file("assets/particlesjs-config.json")
   # )

  svtypes = c('DEL', 'DUP', 'INV', 'INS')
  svsize.max = 1e9

  ui <- shinyUI(
    tagList(tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),
            navbarPage(
              title = "SVCall",
              fluid = TRUE,
              #inverse = inverse,
              windowTitle = "SVCall",
              header = head,
              theme = shinythemes::shinytheme( "paper"),
              id = "tabs",
              tabPanel(
                "HOME",
                shinyjs::useShinyjs(),
                div(
                  class = "container",
                  style = "min-height:90vh;",
                  div(
                    style = "width: 100%; height: 300px; position: relative;z-index:-9;",
                    div(
                      id = "particles-target",
                      style = "position: absolute; top: 0; bottom: 0; right: 0; left: 0;"
                    ),
                    div(
                      style = "padding-top:60px;",
                      h1("\u0043\u006c\u0069\u006e\u0069\u0063\u0061\u006c\u0053\u0056", class = "center"),
                      h3("Clinically Reportable Structural Variant Calls", class = "center")
                    )
                  ),
                   # particlesjs::particles(particles_json, target_id = "particles-target", element_id = "particles"),
                  #img(src = 'ClinicalSVsLogo.png', height = '100px', width = '100px'),
                  tabsetPanel(
                    type = "tabs",
                    tabPanel(
                      "SEARCH",
                      fluidRow(
                        column(
                          1,
                          br(),
                          actionButton("opts", "", icon = icon("plus")),
                          tippy_this("opts", "More options")
                        ),
                        column(
                          9,
                          textInput("q", "", width = "100%", placeholder = "Enter your search query here."),
                          tippy_this("q", "Your search query")
                        ),
                        column(
                          2,
                          br(),
                          actionButton(
                            "submit",
                            "Search",
                            icon = icon("search"),
                            width = "100%",
                            class = "btn btn-primary"
                          )
                        )
                      ),
                      div(
                        id = "options",
                        style = "display:none;",
                        h3("Options")
                      )
                    ),
                    tabPanel(
                      "LOAD",
                      fileInput(
                        "file",
                        label = "Choose one or more previously downloaded VCF file (.vcf.gz)",
                        accept = c(".RData", ".rdata",".vcf.gz"),
                        placeholder = " No file selected",
                        multiple = TRUE,
                        width = "80%"
                      )
                    )
                  ),
                  br(),
                  br(),
                  div(
                    style = "position:fixed;bottom:0px;right:43%;",
                    p(
                      class = "center",
                      "Visit",
                      a(
                        "https://github.com/collaborativebioinformatics/clinical_SVs",
                        href = "https://github.com/collaborativebioinformatics/clinical_SVs",
                        target = "_blank"
                      ),
                      "for more information."
                    )
                  )
                )
              ),
              tabPanel(
                "Annotated",
                #networks_ui("networks")
                #textOutput("text1"),
                splitLayout(cellWidths = c("25%", "25%", "25%", "25%"),
                            numericInput('size.min', 'Minimum SV size (bp)', 0, 0),
                            numericInput('size.max', 'Maximum SV size (bp)', svsize.max, svsize.max),
                            downloadButton('downloadData', 'Download annotated VCF'),
                            downloadButton('downloadPlot', 'Download Plots')
                ),
                checkboxGroupInput('svtypes', "SV type", svtypes, svtypes, inline = TRUE)
              )
            )
    ))
  server <- function(input, output, session){
    output$text1 <- renderText({
      paste("Developing")
    })

    observeEvent(input$submit, {
      updateTabsetPanel(session = session, inputId = "tabs", selected = "Annotated")
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(data, con)
      }
    )

    output$downloadPlot <- downloadHandler(
      filename = function() { paste(input$dataset, '.png', sep='') },
      content = function(file) {
        device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = plotInput(), device = device)
      }
    )
  }

  shinyApp(ui, server)
}

####
