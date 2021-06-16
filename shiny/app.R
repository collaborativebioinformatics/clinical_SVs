library(shiny)
library(tippy)
library(shinythemes)
library(tidyverse)
library(tidygraph)
library(shinyjs)
library(particlesjs)
annot.rdata <- 'annotation_data.RData'
out.vcf <- 'clinical-sv-annotated.vcf'
out.csv <- 'clinical-sv-table.csv'
load(annot.rdata)

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

  particles_json <- jsonlite::fromJSON(
    "assets/particlesjs-config.json"
  )

  svtypes = c('DEL', 'DUP', 'INV', 'INS')
  svsize.max = 1e9

  ui <- shinyUI(
    tagList(
      tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),
            fluidPage(
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
                   particlesjs::particles(particles_json, target_id = "particles-target", element_id = "particles"),
                  #img(src = 'assets/ClinicalSVsLogo.png', height = '100px', width = '100px'),
                  tabsetPanel(
                    type = "tabs",
                    tabPanel(
                      "Annotate",
                      fluidRow(
                        # column(
                        #   1,
                        #   br(),
                        #   actionButton("opts", "", icon = icon("plus")),
                        #   tippy_this("opts", "More options")
                        # ),
                        # column(
                        #   9,
                        #   textInput("q", "", width = "100%", placeholder = "Enter your search query here."),
                        #   tippy_this("q", "Your search query")
                        # ),
                        column(
                             9,
                          fileInput(
                            "file",
                            label = "Upload the VCF file (.vcf)",
                            accept = c(".vcf","text/csv","text/comma-separated-values, text/plain"),
                            placeholder = " No file selected",
                            multiple = TRUE,
                            width = "100%"
                          ),
                        # column(
                        #   2,
                        #   br(),
                          actionButton(
                            "submit",
                            "Anotate",
                            icon = icon("search"),
                            width = "100%",
                            class = "btn btn-primary"
                          # )
                          )
                        )
                      )
                      ,
                      div(
                        id = "options",
                        style = "display:none;",
                        h3("Options")
                      )
                    )
                  ),
                  br(),
                  #textOutput("text1"),
                  #textOutput("text2"),
                  splitLayout(
                              downloadButton('downloadvcf', 'Download annotated VCF'),
                              downloadButton('downloadcsv', 'Download annotated CSV'),
                              downloadButton('downloadPlot', 'Download Plots')
                  ),
                  a(id = "toggleAdvanced", "Show/hide advanced info"),
                  shinyjs::hidden(div(id = "advanced",
                  checkboxGroupInput('svtypes', "SV type", svtypes, svtypes, inline = TRUE)
                  )),
                  fluidRow(column(7,dataTableOutput('newvcf'))
                  ),
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
              )
              # tabPanel(
              #   "Annotated",
                #networks_ui("networks")

            )
    ))
  server <- function(input, output, session){

    observe({
      shinyjs::toggleState("submit", !is.null(input$file) && input$file != "")
    })

    output$text1 <- renderText({
      paste0(input$file)
    })
    output$text2 <- renderText({
      readSVvcf(paste0(input$file$datapath), out.fmt='vcf', keep.ids=TRUE)
    })


    output$newvcf <- renderDataTable({
      req(input$file$datapath)


      ## read VCF
      suppressWarnings(suppressMessages(library(sveval)))
      suppressWarnings(suppressMessages(library(GenomicRanges)))
      vcf.o = readSVvcf(input$file$datapath, out.fmt='vcf', keep.ids=TRUE)

      ## make sure chromosomes are in the form 'chrX'
      if(all(!grepl('chr', seqlevels(vcf.o)))){
        seqlevels(vcf.o) = paste0('chr', seqlevels(vcf.o))
      }

      ## annotate gene overlapped by SVs
      source('annotate_genes.R')
      vcf.o = annotate_genes(vcf.o, genc)

      ## annotate frequency
      source('annotate_frequency.R')
      vcf.o = annotate_frequency(vcf.o, gnomad)

      ## annotate known clinical SVs
      source('annotate_known_clinical_SVs.R')
      vcf.o = annotate_known_clinical_SVs(vcf.o, clinsv)

      ## clinical ranks, to order the SVs and select top 5 for example
      source('annotate_clinical_score.R')
      vcf.o = annotate_clinical_score(vcf.o,genc)

      ## write annotated VCF


      ## write tables
      svs <- tibble(gene=info(vcf.o)$GENE,
                   variant_id=names(vcf.o),
                   chr=as.character(seqnames(vcf.o)),
                   start=start(vcf.o),
                   end=end(vcf.o),
                   size=abs(unlist(lapply(info(vcf.o)$SVLEN, '[', 1))),
                   frequency=info(vcf.o)$AF,
                   svtype=info(vcf.o)$SVTYPE,
                   clinsv=info(vcf.o)$CLINSV,
                   clinrk=info(vcf.o)$CLINRK) %>%
        arrange(clinrk)
      return(svs)
      })


    observeEvent(input$submit, {
      updateTabsetPanel(session = session, inputId = "tabs", selected = "Annotated")
    })

    output$downloadvcf <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.vcf', sep='')
      },
      content = function(con) {
       writeVcf(vcf.o, con)
       }
    )

    output$downloadcsv <- downloadHandler(
      filename = function() {
        paste('data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(svs,con)
      }
    )

    shinyjs::onclick("toggleAdvanced", shinyjs::toggle(id = "advanced", anim = TRUE))

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
