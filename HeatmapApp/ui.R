require(shiny)
require(shinythemes)
require(visNetwork)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(colourpicker)
library(formattable)

##

#fluidPage(theme = shinytheme("cosmo"),
fluidPage(theme = shinytheme("journal"),
          
          # Application title
          #titlePanel("LION/web | Lipid Ontology heatmap for lipidomics"),
          br(),
          
          sidebarLayout(
            
            
            sidebarPanel( width = 4,
                          
                          tabsetPanel(id = "method", type = "tabs", selected = "Data upload",
                          
                                      tabPanel(id = "Settings", title = "", icon = icon("cog", lib = "font-awesome"),
                                               br(),
                                               br(),
                                               
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "RemoveRedundantTerms",
                                                   value = TRUE,
                                                   label = "Exclude parent terms with same associations as child",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = 'When a LION-term contains the same lipids as its parent (for instance, glycerophosphocholines and diacylglycerophosphocholines), the parent term (here glycerophosphocholines) will be excluded in analysis. '
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                      materialSwitch(
                                                        inputId = "SmartMatching",
                                                        value = TRUE,
                                                        label = "Use 'smartmatching' (beta)",
                                                        status = 'primary',
                                                        right = TRUE
                                                      ),
                                                      content = "If possible, unmatched lipids are associated with related LION-terms: `TG(O-16:0/18:2/22:6)`	is associated with `alkyldiacylglycerols`, `C16:0`, `C18:2`, and `C22:6`"
                                               ), 
                                               br(),
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                      materialSwitch(
                                                        inputId = "FAprediction",
                                                        value = FALSE,
                                                        label = "Predict fatty acid assocations (beta)",
                                                        status = 'primary',
                                                        right = TRUE
                                                      ),
                                                      content = "For sum-formatted phospholipids, e.g. PC(34:1), predict the most likely fatty acid assocations, based on reported fatty acid compositions. For example, PC(34:1) will be associated to C16:0 and C18:1. Only use for datasets of mammalian origin"
                                               ), 
                                               br(),
                                          
                                               
                                               popify(placement = "right", title = "Info", options=list(container="body"),
                                                 materialSwitch(
                                                   inputId = "LIONselection",
                                                   label = "Preselect LION-terms for analysis",
                                                   status = 'primary',
                                                   right = TRUE
                                                 ),
                                                 content = "By default, all LION-terms will be used in the analysis. With this option, you can pre-select LION-terms to use in the enrichment analysis."
                                               ), 
                                               em(""),
                                               shinyTree("tree",checkbox = TRUE,search = TRUE),
                                               br()
                                      ), 
                                      
                                      tabPanel("Data upload", 
                                               br(),
                                               br(),
                                               "Choose CSV File:",
                                               fluidRow(
                                                 column(
                                                   offset = 0,
                                                   #style='margin-left:2%;' ,
                                                   width = 10,
                                                   popify(
                                                     placement = "bottom",
                                                     title = "File-input info",
                                                     fileInput(
                                                       "file1",
                                                       label = NULL,
                                                       multiple = FALSE,
                                                       accept = c("text/csv",
                                                                  "text/comma-separated-values,text/plain",
                                                                  ".csv")
                                                     ),
                                                     content = 'Format your dataset as comma seperated value files (.csv), with the first column reserved for metabolites and the other columns for numeric data (containing decimal points). Use double column headers; with row 1 containing condition identifiers and row 2 containing sample identifiers. Submit at least duplicates per condition. Dataset should be normalized before submission. Download a dataset below for an example.'
                                                     
                                                   )
                                                 ),
                                                 column(
                                                   offset = 0,
                                                   width = 1,
                                                   style = "margin-top: 5px;",
                                                   align = "center",
                                                   popify(
                                                     placement = "right",
                                                     title = "Lipid nomenclature",
                                                     options = list(container = "body"),
                                                     el = icon(name = "question", lib = "font-awesome", "fa-2x"),
                                                     content = 'Format lipids in LIPIDMAPS notation style: a class-prefix followed by (summed) fatty acid(s) surrounded by parentheses. Examples are: PC(32:1); PE(18:1/16:0); SM(d18:1/18:0); TAG(54:2); etc. Check www.lipidmaps.org for more examples. LION will try to reformat alternative notation styles into LIPIDMAPS format.'
                                                   )
                                                 )
                                               ),
                                               
                                               downloadLink("examplePre1",
                                                            "example set 1 (organelle fractions) [1]"),
                                               
                                               br(),
                                               downloadLink("examplePre2",
                                                            "example set 2 (CHO-k1 incubated with several FFAs) [2]"),
                                               br(),
                                               downloadLink("examplePre3", "example set 3 (CHO-k1 incubated with AA) [2]"),
                                               br(),
                                               
                                               em(tags$small("[1] adapted from Andreyev AY et al, 2010, [2] from Molenaar MR et al, 2019")),br(),
                                               #em(tags$small("[2] from Molenaar MR et al, 2019")),
                                               br(),
                                               # br(),
                                               br(),
                                               uiOutput("SubmitUI"),
                                               br(),
                                               uiOutput("selectTresholdsUI")
                                  
                                               
                                               
                                               
                                               
                                      )

                        
                          
                          

            )),
            
            mainPanel(
              tabsetPanel(id = "tabs", type = "tabs", 
                          tabPanel("General information", 
                                   h3("LION/web | LION-PCA heatmap module"),
                                   br(),
                                   img(src="LIONicon heatmap.png", align = "left", height = '200px'), #, width = '550px'),
                                   br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                                   p('The Lipid Ontology (LION) enrichment analysis web application (LION/web) is a novel bioinformatics tool for lipidomics that',
                                    'enables users to search for enriched LION-terms in lipidomic subsets. LION-terms contain detailed lipid classification',
                                    'by LIPIDMAPS, biophysical data, lipid functions and organelle associations.'),
                                   p("Any comments, questions or suggestions? Do you want to add lipid annotations? Please use our contact form."),
                                   br(),
                                   p("Please cite:"),
                                   a('LION/web: a web-based ontology enrichment tool for lipidomic data analysis.',
                                     href="https://doi.org/10.1093/gigascience/giz061", target="_blank"),br(),
                                   tags$b("Gigascience. 2019 Jun 1;8(6). pii: giz061. doi: 10.1093/gigascience/giz061."),
                                   br(),
                                   em("Martijn R. Molenaar, Aike Jeucken, Tsjerk A. Wassenaar, Chris H. A. van de Lest, Jos F. Brouwers, J. Bernd Helms."), 
                                   br(),
                                   br(),
                                   br(),a('Lipidomic profiling of rat hepatic stellate cells during activation reveals a two-stage process accompanied by increased levels of lysosomal lipids.',
                                     href="https://doi.org/10.1016/j.jbc.2023.103042", target="_blank"),br(),
                                   tags$b("J Biol Chem. 2023 Feb 18;299(4):103042. doi: 10.1016/j.jbc.2023.103042."),
                                   br(),
                                   em("Martijn R Molenaar, Maya W Haaker, A Bas Vaandrager, Martin Houweling, J Bernd Helms."), 
                                   br(),
                                   br(),
                                   br(),
                                   p("Division Cell Biology, Metabolism & Cancer, Department Biomolecular Health Sciences, Universiteit Utrecht, The Netherlands"),
                                   'v. 2023.04.15',
                                   
                                   br(),
                                   br()
                                   
                                   ),
                          tabPanel("LION input", 
                                   br(),
                                   #em(   
                                   htmlOutput("mapping_percentage"),    #textOutput("mapping_percentage")  ),
                                   br(),
                                   formattableOutput("value"),        #tableOutput("value"),
                                   uiOutput("downloadInputCNTRL"),
                                   
                                   br()
                                   
                                   
                          ),
                        
                          tabPanel("LION heatmap", 
                                   br(),
                                   plotOutput("ontology.graph", height = "auto"), #, height = 650),
                                   
                                   fluidRow(
                                     br(),
                                     br(),
                                     br(),
                                     bsCollapse(id = "figureOptions", 
                                                bsCollapsePanel("Figure export", 
                                                                column(offset = 0,width = 5,
                                                                splitLayout(
                                                                       numericInput(inputId = "figureHeight", label = "height (inch)", min = 1, max = 49,value = 9.06),
                                                                       numericInput(inputId = "figureWidth", label = "width (inch)", min = 1, max = 49, value = 14.3)
                                                                )),
                                                                br(),
                                                                column(offset = 0,width = 9,
                                                                       downloadButton("downloadPlotSVG", "Download heatmap as SVG"), 
                                                                       downloadButton("downloadPlotPNG", "Download heatmap as PNG"),
                                                                       downloadButton("downloadTable", "Download LION-term information"))
                                                                )
                                     )
                                    
                                   ),
                                   
                                   br()
                                   )
                          
                          
                                   
                                   
              )
            )
          )
)
