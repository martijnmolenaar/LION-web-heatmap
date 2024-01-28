options(stringsAsFactors = FALSE,shiny.sanitize.errors = TRUE)

source(file  = "data/20191010 LION_tree_structure.R")

#### libraries

require(shiny)
require(data.table)
require(ggplot2)
require(ggthemes)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(httr)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(colortools)
library(colourpicker)
library(cowplot)
library(formattable)
library(sortable)

## loading lipid ontology data
require(RSQLite)
require(topOnto)
require('topOnto.LION.db')
topOnto::initONT('LION')


associationFile  <-  "data/20190704 LION_association.txt"

LIONterms_rules <- read.csv(file = 'data/20191008 LIONterms_rules.csv', header = T)
LIONterms_rules$RULE1[LIONterms_rules$RULE1 == ""] <- "-"
LIONterms_FAs <- read.csv(file = 'data/20191008 LIONterms_FAs.csv', header = T)

FA_composition_table <- read.csv(file = 'data/FA_composition_table.csv')

## define functions

source('data/20191008 LIONweb functions.R')


## read associations
lipidID2TERM <- readMappings(file = associationFile)   ## topOnto function, but removes spaces


# Define server logic for random distribution application
function(input, output, session) {
  
  
  ### hide tabs at start-up
  hideTab(inputId = "tabs", target = "LION input")
  hideTab(inputId = "tabs", target = "LION heatmap")
  
  showNotification(ui = "",
                   action =  p("By using this app you agree with the", a('Terms of Usage.',
                                                                         href="https://martijnmolenaar.github.io/lipidontology.com/faq.html#basics", 
                                                                         target="_blank")),
                   duration = 30, type = "default")
  
    ## pre-processing with CSVs:
  
  input_data <- reactive({
    
    req(input$file1)
    hideTab(inputId = "tabs", target = "LION input")
    hideTab(inputId = "tabs", target = "LION heatmap")
    
    updateTabsetPanel(session, "tabs",
                      selected = "General information"
    )
    
    file_location <- input$file1$datapath
    
    df <- try(read.csv(file_location,
                       header = FALSE,
                       sep = "," #, quote = ""
    ))
    
    
    
    
    ## error handling
    
    errors <- NULL
    if (!(is.data.frame(df))) {
      df <- data.frame(errors = df[1], column = 0)
      errors <- c(errors, "ERROR: File type not supported")
    } else {
      if (dim(df)[2]  == 1) {
        errors <- c(errors, "ERROR: No commas found to seperate columns")
      }
      if (dim(df)[2]  > 1 & dim(df)[2]  < 3) {
        errors <- c(errors, "ERROR: Unkown error")
      }
      if (sum(is.na(df[, -1]))   > 0) {
        errors <- c(errors, "ERROR: There are missing values")
      }
      if (sum( df[-c(1,2),-1] == "" ) > 0) {
        errors <- c(errors, "ERROR: Dataset contains empty cells")
      }
      if (all(table(as.character(df[1, -1])) < 2)) {
        errors <- c(errors, "ERROR: Some or all conditions are n < 2")
      }
      if (any(duplicated(df[2, -1]))) {
        errors <-
          c(errors, "ERROR: One or more samples names are not unique")
      }
      if (sum(apply(df[-1,-c(1:2)],1,function(row){all(row == 0)})) > 0) {
        errors <-
          c(errors, "ERROR: Dataset contains rows with only zeros")
      }
    }
    errors <- paste(errors, collapse = "<br>")
    
    ## end error handling
    
   
    #if(!(is.null(errors))){
    if(errors==""){
      input_data <- list(df = df,
         matrix = sapply(df[-c(1,2),-1], as.numeric),
         conditions = unique(as.character(df[1,-1])),
         samples = as.character(df[2,-1]),
         meta = df[1:2,-1],
         IDs = df[-c(1,2),1],
         errors = errors)
      
      input_data$species_scaled <- t(scale(t(input_data$matrix)))
      input_data$species_scaled[  is.nan(input_data$species_scaled)] <- 0   ## data with zero SD all to 0
      
      rownames(input_data$species_scaled) <- convertLipidNames(input_data$IDs)
      colnames(input_data$species_scaled) <- input_data$samples
      
      input_data$samples_scaled <- scale(input_data$matrix)
      rownames(input_data$samples_scaled) <- convertLipidNames(input_data$IDs)
      colnames(input_data$samples_scaled) <- input_data$samples
      
      
    } else {
      input_data <-list(df = NULL,
           matrix = NULL,
           conditions = NULL,
           samples = NULL,
           meta = NULL,
           IDs = NULL,
           errors = errors,
           samples_scaled = NULL,
           species_scaled = NULL)
    }
    

    return(input_data)
    
  })
  
  output$SubmitUI <- renderUI({
    
    input_data <- input_data()
    
    if (input_data$errors == "") {
    
      fluidRow(
        column(
          offset = 0,
          width = 12,
          
          popify(
            placement = "bottom",
            title = "Info",
            checkboxInput(
              inputId = "normalization",
              label = "normalize signals as percentage",
              value = FALSE
            ),
            content = 'If checked, lipid signals are expressed as percentage of the total signal per sample. This is useful when input data is not normalized.'
          ),
          
          popify(
            placement = "top",
            title = "Info",
            actionButton("submitB", "  Submit Data", 
                         icon = icon("lightbulb-o", lib = "font-awesome")),  
            content = 'After clicking "submit", LION/web matches the input lipids to the LION-database. Only matched identifiers will be used in the  analysis. Subsequently, LION/web selects the most fluctuating LION-terms by assessing enrichment over a set number of principle components (similar approach as: Wagner F, 2015. PLoS ONE). The selected LION-terms will be presented in a heatmap of the provided dataset.')
        
          
   
        )
      )
    } else {     ### if there are errors...
      HTML(paste('<font color="red">',
                 input_data$errors,
                 "</font>",sep=""))
    }
   
  })
  
  output$PCA_screeplot <- renderPlot({
    PCA_scree_plot <- dataBlock()$PCA_scree_plot
    
    PCA_scree_plot
   
  })
  
  
  output$selectTresholdsUI <- renderUI({
    
    dataBlock <- dataBlock()
    
    isolate(input_data <- input_data())
  
    
    wellPanel(
      tabsetPanel(id = "heatmap_units", selected = "LION-terms", type = "pills",
        tabPanel("LION-terms", 
                 br(),
                 popify(
                   placement = "right",
                   title = "Number of principal components",
                   options = list(container = "body"),
                   el = plotOutput("PCA_screeplot", height = 250),
                   content = 'This plot shows the number of derived LION-terms (left) and the cumulative variance explained (right) as a function of the number of principal components. LION-PCA performs LION-term enrichment based on the loadings of a given number of principal components, which can be set below.'
                 ),
                 
                 br(),
                 br(),
                 numericInput("PCAn", label = "number of principal components", value = min(3,dataBlock$max_pcas),
                              min = 1, max = dataBlock$max_pcas)),
        tabPanel("lipid species", br())
      ),
      
      br(),
      
      bsCollapse(id = "heatmapColorOptions",
                 bsCollapsePanel(title = "hierarchical clustering options", 
                                 # column(
                                 #   offset = 0,
                                 #   width = 5,
                                   materialSwitch(
                                     inputId = "ClusterSamples",
                                     label = "cluster samples",
                                     status = 'primary',
                                     value = FALSE,
                                     right = TRUE
                                   ),
                                   materialSwitch(
                                     inputId = "ClusterTerms",
                                     label = "cluster terms",
                                     status = 'primary',
                                     value = TRUE,
                                     right = TRUE
                                   ),
                                   #uiOutput("clusterTermsUI"),
                                   wellPanel(fluidRow(
                                     column(width = 12,
                                            materialSwitch(
                                              inputId = "CutClusters",
                                              label = "cut dendrogram",
                                              status = 'primary',
                                              value = TRUE,
                                              right = TRUE
                                            ),
                                            numericInput(
                                              "clusternr",
                                              label = "number of clusters", width = '90%',
                                              value = 5,
                                              min = 1,
                                              max = 20
                                            )
                                     )
                                   ))
                                   
                                 #)
                                 
                 ),
                 bsCollapsePanel(title = "legend colors conditions", 
                                 column(
                                   offset = 0,
                                   width = 5,
                                   lapply(1:length(input_data$conditions), function(i) {
                                     conditions_i <- input_data$conditions[i]
                                     colourInput(
                                       inputId = paste("color", i, sep = ""),
                                       label = conditions_i, #showColour = TRUE, 
                                       value = colorRampPalette(c("lightgray", "darkblue","black"))(length(input_data$conditions))[i]  #"red"
                                     )
                                   })
                                 )
                                 
                 ),
                 bsCollapsePanel(title = "re-order conditions",
                                 rank_list(
                                     text = "Re-order conditions",
                                     labels = as.list(input_data$conditions),
                                     input_id = "ranked_conditions"
                                   )
                                 ),
                 bsCollapsePanel(title = "heatmap colors",
                                 column(
                                   offset = 0,
                                   width = 5,
                                   lapply(c("down-regulated","up-regulated"), function(i) {
                                     conditions_i <- input_data$conditions[i]
                                     colourInput(
                                       inputId = i,
                                       label = i, #showColour = TRUE, 
                                       
                                       value = ifelse(i == "down-regulated", "#DBB51A", "red") #"red"
                                     )
                                   })
                                 )
                                 )
                 ),
      
      br(),
      popify(
          placement = "top",
          title = "Info",
          actionButton("showPlot", "  Show Heatmap", 
                       icon = icon("chart-line", lib = "font-awesome")),  
          content = 'Generate heatmap'),
      br(),
      br()
      
    )
   
  })
  
  

  ## examples
  
  output$examplePre1 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set - adapted from Andreyev AY et al 2010","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set - adapted from Andreyev AY et al 2010.csv", file)
    },
    contentType = "text/csv"
    )
  output$examplePre2 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set FA incorporation","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set FA incorporation.csv", file)
    },
    contentType = "text/csv"
  )
  output$examplePre3 <- downloadHandler(
    filename <- function() {
      paste("lipidomics set membrane fluidity after AA incubation","csv",sep=".")
    },
    
    content <- function(file) {
      file.copy("data/lipidomics set membrane fluidity after AA incubation.csv", file)
    },
    contentType = "text/csv"
  )
  
  
  observeEvent(input$submitB, {
    ## first unhide tabs
    showTab(inputId = "tabs", target = "LION input")
    # showTab(inputId = "tabs", target = "LION heatmap")
    
    ## goto LION input
    updateTabsetPanel(session, "tabs",
                      selected = "LION input")
  })
  
  observeEvent(input$showPlot, {
    
    showTab(inputId = "tabs", target = "LION heatmap")
    visHeatmapBlock()
    
    ## goto LION input
    updateTabsetPanel(session, "tabs",
                      selected = "LION heatmap")
  })
  
  ### heatmap analysis
  dataBlock <- eventReactive(input$submitB, {  
                
    errorhandling <- NULL
    
    input_data <- isolate(input_data())
    
    if(isolate(input$normalization)){  ## normalization option switched on
      input_data$matrix[] <- 
        apply(input_data$matrix,2,function(i){
          i / sum(i) * 100
        })
      
      input_data$samples_scaled[] <- scale(input_data$matrix)
      input_data$species_scaled[] <- t(scale(t(input_data$matrix)))
      input_data$species_scaled[  is.nan(input_data$species_scaled)] <- 0
    }
    
    ###  table 1
    
    withProgress(message = 'progress:', value = 0, {
      incProgress(.1, detail = paste("mapping input data"))
      
      lipidExistance <- data.frame(input = paste(sprintf("[#%04d]", 1:length(input_data$IDs)),input_data$IDs))
      lipidExistance$'simplified input' <-  convertLipidNames(gsub("^\\[#\\d+\\] ","",lipidExistance$input))
      
      
      ##### new 20191106
      lipidExistance_list <-
        lapply(1:dim(lipidExistance)[1], function(row_i) {
          
          lipid_i <- lipidExistance[row_i, 2]   ## lipid_i is lipid of this iteration
          
          LION_ID <-
            unlist(lipidID2TERM[lipid_i == names(lipidID2TERM)])
          
          if (is.null(LION_ID)) {
            ### if input is LION:xxxx
            LION_ID <-
              unlist(lipidID2TERM[lipid_i == unlist(lipidID2TERM)])
            LION_ID <-
              LION_ID[!grepl("^SLM:|^LM..\\d+", names(LION_ID))]   ## remove SwissLipids/LIPIDMAPS IDs
          }
          
          if (!is.null(LION_ID)) {
            output_df <-
              data.frame(
                input = lipidExistance[row_i, 1],
                'simplified input' = lipid_i,
                name = names(LION_ID),
                LION = LION_ID,
                match = "direct"
              )
          } else {
            ### no direct matching is found
            if (isolate(input$SmartMatching)) {
              ### 'smartmatching' is on
              lipid_index <-
                list(
                  generalized = simplifyLipidnames(lipid_i),
                  headgroup = getLipidHeadgroup(lipid_i),
                  linkage = getLipidLinkage(lipid_i),
                  FAs = getFAs(lipid_i)
                )
              
              generalized_LION_ID <-
                unlist(lipidID2TERM[lipid_index$generalized == names(lipidID2TERM)])
              if (!is.null(generalized_LION_ID)) {
                terms <-
                  data.frame(name = names(lipidID2TERM)[lipid_index$generalized == names(lipidID2TERM)],
                             LION = generalized_LION_ID)
              } else {
                ## no generalized lipid found
                LIONterms_rules_i <-
                  LIONterms_rules[lipid_index$headgroup == LIONterms_rules$RULE1 ,]
                
                if (all(LIONterms_rules_i$RULE2 == "")) {
                  terms <- LIONterms_rules_i[, 1:2]
                } else {
                  terms <-
                    LIONterms_rules_i[lipid_index$linkage == LIONterms_rules_i$RULE2,][, 1:2]
                }
              }
              
              
              terms <-
                rbind(terms, LIONterms_FAs[LIONterms_FAs$name %in% lipid_index$FAs,])
              if (dim(terms)[1] > 0) {
                output_df <-
                  data.frame(input = lipidExistance[row_i, 1],
                             'simplified input' = lipid_i,
                             terms,
                             match = "smart matching")
                
              } else {
                ## no match by smart matching
                output_df <-
                  data.frame(
                    input = lipidExistance[row_i, 1],
                    'simplified input' = lipid_i,
                    name = "not found",
                    LION = "not found",
                    match = ""
                  )
              }
              
              
            } else {
              ### 'smartmatching' is off, and no match found
              output_df <-
                data.frame(
                  input = lipidExistance[row_i, 1],
                  'simplified input' = lipid_i,
                  name = "not found",
                  LION = "not found",
                  match = ""
                )
            }
            
            
          }
          if(isolate(input$FAprediction)){    ## predict FAs if nescerrary
            
            lipid_index <-
              list(
                generalized = simplifyLipidnames(lipid_i),
                headgroup = getLipidHeadgroup(lipid_i),
                FAs = getFAs(lipid_i)
              )
            
            if(lipid_index$generalized == lipid_i &   length(lipid_index$FAs) == 1){   ## is FA prediction applicable?
              
              predicted_FAs <- predict_FAs(headgroup = lipid_index$headgroup, 
                                           summedFA = gsub("C","",lipid_index$FAs), 
                                           composition_table = FA_composition_table )
              
              if(dim(LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,])[1]>0){    ## any result?
                output_df <-
                  rbind(output_df,
                        data.frame(input = lipidExistance[row_i, 1],
                                   'simplified input' = lipid_i,
                                   LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,],
                                   match = "FA-prediction"))
              } 
              
              
            }
            
            
          }
          
          return(output_df)
        })
      
      
      matching_statistics <- data.frame(total = length(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                        matched = sum(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                        percent = mean(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})    ) * 100)
      
      
      lipidExistance <- do.call("rbind",lipidExistance_list)
      
      colnames(lipidExistance) <- c('input','simplified input','LION name','LION ID', "match")
      
    
      color_df <- data.frame(match = c("direct","smart matching", "FA-prediction",""),
                             color = c("#4f4d4d","#9c3c2d","#c4942b","#bdbbbf"))
      
      lipidExistance_toview <-
        do.call("rbind", lapply(lipidExistance_list, function(lipidExistance_i) {
          
          color_pattern <- color_df$color[match(lipidExistance_i$match, color_df$match)]
          
          data.frame(
            input = unique(lipidExistance_i$input),
            name = paste("<font color='",color_pattern,"'>",lipidExistance_i$name,"</font>" , sep ="", collapse = "<br>"),
            LION = paste("<a href='",
                         'https://bioportal.bioontology.org/ontologies/LION/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F',
                         gsub(":","_",lipidExistance_i$LION),
                         "' style='color: ", color_pattern,
                         "' target='_blank'>",
                         lipidExistance_i$LION,
                         "</a>",
                         sep = "", collapse = "<br>")
            
            
          )
          
         
        }))
      
      colnames(lipidExistance_toview) <- c('input','LION name','LION ID')
      lipidExistance_toview[[3]] <- gsub(" href=.+Fnot found'","",lipidExistance_toview[[3]])
      
      ##### end new 20191006
      
      lipidExistance$'simplified input' <- NULL        ## remove this column, not interesting for user
      
      
      ## use matched lipids as assocation table 
      lipidExistance_feasable <- lipidExistance[lipidExistance$`LION ID` != "not found",]
      
      matched_lipidID2TERM <- 
        sapply(unique(lipidExistance_feasable$input), function(input){
          lipidExistance_feasable$`LION ID`[lipidExistance_feasable$input == input]
        }, simplify = FALSE)
      
   
    })
    
    
    if (!any(dim(lipidExistance) == c(1, 1))) {
      ### perform when there is input (dim 1 by 1 >> no input)
      
      withProgress(message = 'progress:', value = 0, {
        # Generate an HTML table view of the ontology data
        
        incProgress(.5, detail = paste("submitting data"))
        
        #### heatmap insert
        
        
        
        pca_data <-
          prcomp(input_data$samples_scaled,
                 center = TRUE,
                 scale = TRUE)
        PCA_df <- as.data.frame(pca_data$rotation)
        PCA_df$samples <- rownames(PCA_df)
        pca_data$explained_variance <-
          pca_data$sdev / sum(pca_data$sdev)
        
        lipidIDrank <-
          rank(pca_data$x[, 1])  ### mock ranking, need it for ONTdata construction
        
        #### convert input to LION-friendly formatting
        
        names(lipidIDrank) <- unique(lipidExistance$input)
        
        mySel <- function(allScore) {
          return(rep(TRUE, length(lipidIDrank)))
        }
        
        ONTdata <- new(
          ### making ontology object
          "topONTdata",
          ontology = "LION",
          allGenes = lipidIDrank,
          annot = annFUN.gene2GO,
          gene2GO = matched_lipidID2TERM,
          geneSelectionFun = mySel
        )
        
        LUT <-
          data.frame(ID = names(ONTdata@termName),
                     ### look-up table for LION-IDs and descriptions
                     Discription = ONTdata@termName)
        
        #
        max_pcas <- min(c(dim(pca_data$rotation)[2],10))
        
        incProgress(.7, detail = paste("performing PCA"))

        LION_term_by_PCA <-
          ## doing analysis of PC1, PC2, PCn..., with ks2 (2-tailed KS) algorithm to assess both directions
          
          
          sapply(as.data.frame(pca_data$x[,1:max_pcas]),function(pc_values) {
            ONTdata@allScores <- rank(pc_values)
            
            resultFis <-
              runTest(ONTdata,
                      algorithm = "classic",
                      statistic = "ks2")
            
            
            resultFis@score <- abs(resultFis@score)   ### ks2 is now signed
            
            to_display <-
              GenTable(ONTdata,
                       'p-value' = resultFis,
                       topNodes = 2000)
            
            to_display <-
              to_display[grep("LION", to_display$TERM.ID), ]
            
            
            
            to_display$Term <-
              LUT$Discription[match(to_display$TERM.ID, LUT$ID)]        ### otherwise, some names are abbrev.
            to_display <-
              to_display[to_display$Annotated > 2, ]                         ### remove terms with 1 or 2 lipids
            
            to_display$`p-value` <-
              gsub("< 1e", "< 1.0e", to_display$`p-value`)           ### '< 1e-30' cannot be understood
            to_display$`p-value` <-
              gsub("< ", "", to_display$`p-value`)                   ### '< 1e-30' cannot be understood
            
            to_display$'FDR q-value' <-
              ### correction for multiple comparison
              format(p.adjust(as.numeric(to_display$`p-value`), "fdr"), digits = 3)
            
            colnames(to_display) <-
              ### for readibility
              c(
                "Term ID",
                "Discription",
                "Annotated",
                "Significant",
                "Expected",
                "p-value",
                "FDR q-value"
              )
            to_display <-
              to_display[, c(1, 2, 3, 6, 7)]
            
            
            ## select significant terms or otherwise the first 3 terms
            if (length(to_display$`Term ID`[as.numeric(to_display$`p-value`) < 0.05]) == 0) {
              list(
                LION_IDs = to_display$`Term ID`[1:3],
                ONTdata = ONTdata,
                table = to_display
              )
            } else {
              list(
                LION_IDs = to_display$`Term ID`[as.numeric(to_display$`p-value`) < 0.05],
                ONTdata = ONTdata,
                table = to_display
              )
            }
            
            
          }, simplify = FALSE)
        
        if (isolate(input$LIONselection)) {  ### switch for LION-selection
          TermsOfInterest <-  get_selected(isolate(input$tree), format = "names")
          TermsOfInterest <- unique(c(unlist(TermsOfInterest),
                                      unlist(sapply(TermsOfInterest, function(element) {
                                        attr(element, "ancestry")
                                      })))) ## vector with LION-names
          
          TermsOfInterest <- LUT$ID[match(TermsOfInterest, LUT$Discription)]
          
          nr_of_LION_terms_by_PCA <-
            sapply(1:length(LION_term_by_PCA), function(i){
              length(unique(unlist(sapply(1:i, function(j){  
                LION_term_by_PCA[[j]]$LION_IDs[LION_term_by_PCA[[j]]$LION_IDs %in%   TermsOfInterest]
                
              }))))
            })
        } else {
          nr_of_LION_terms_by_PCA <-
            sapply(1:length(LION_term_by_PCA), function(i){
              length(unique(unlist(sapply(1:i, function(j){  
                LION_term_by_PCA[[j]]$LION_IDs
                
              }))))
            })
        }
        
        
        incProgress(.8, detail = paste("calculating PCA statistics"))
        
        
        LION_pca_df <-
           data.frame(pca = factor(paste("PC",1:max_pcas,sep=""), levels = paste("PC",1:max_pcas,sep="")),
                      cumvar = cumsum(pca_data$explained_variance)[1:max_pcas]*100,
                   nr_of_LIONs = nr_of_LION_terms_by_PCA)
        
        PCA_scree_plot <-
          plot_grid(
            ggplot(data = LION_pca_df,
                   aes(x = pca, y = nr_of_LIONs)) +
              labs(title = "significant LION-terms", x = "", y = "# LION-terms") +
              geom_bar(stat = "identity", fill = "#C93838", width = .5) +
              theme_linedraw() +
              theme(
                axis.title.y = element_text(face = "bold", size = 9),
                axis.title.x = element_text(face = "bold", size = 9), axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(face = "bold", size = 11, hjust = .5)
              ),
            ggplot(data = LION_pca_df,
                   aes(x = pca, y = cumvar)) +
              labs(title = "PCA scree plot", x = "", y = "cumulative variance (%)") +
              geom_bar(stat = "identity", width = .5) +
              theme_linedraw() +
              theme(
                axis.title.y = element_text(face = "bold", size = 9),
                axis.title.x = element_text(face = "bold", size = 9),axis.text.x = element_text(angle = 45, hjust = 1),
                plot.title = element_text(face = "bold", size = 11, hjust = .5)
              )
          )
        
        
         incProgress(.9, detail = paste("output to browser"))
        
        list(
          ONTdata = ONTdata,
          lipidExistance = lipidExistance,
          lipidExistance_toview = lipidExistance_toview,
          matching_statistics = matching_statistics,
          nr_of_LION_terms_by_PCA = nr_of_LION_terms_by_PCA,
          PCA_scree_plot = PCA_scree_plot,
          LION_pca_df = LION_pca_df,
          max_pcas = max_pcas,
          LION_term_by_PCA = LION_term_by_PCA
          
          
        )
        
      })
    } else {
      ### empty if there is no input
      list(
        ONTdata = NULL,
        lipidExistance = NULL,
        lipidExistance_toview = NULL,
        matching_statistics = NULL,
        nr_of_LION_terms_by_PCA = NULL,
        PCA_scree_plot = NULL,
        LION_pca_df = NULL,
        max_pcas = NULL,
        LION_term_by_PCA = NULL
      
      )
      
    }  ### if perform when there is input (dim 1 by 1 >> no input)
    
    
  })
  
  
  visHeatmapBlock <- reactive({    
    input$showPlot
    
    isolate(dataBlock <- dataBlock())
    LION_term_by_PCA <- dataBlock$LION_term_by_PCA
    ONTdata <- dataBlock$ONTdata
    lipidExistance <- dataBlock$lipidExistance
    
    isolate(input_data <- input_data())
    
    if(isolate(input$normalization)){  ## normalization option switched on
      input_data$matrix[] <- 
        apply(input_data$matrix,2,function(i){
          i / sum(i) * 100
        })
      
      input_data$samples_scaled[] <- scale(input_data$matrix)
      input_data$species_scaled[] <- t(scale(t(input_data$matrix)))
      input_data$species_scaled[is.nan(input_data$species_scaled)] <- 0
    }
    
    #browser()
    
    
    if(isolate(input$heatmap_units != "lipid species")){
  
    withProgress(message = 'progress:', value = 0, {
      # Generate an HTML table view of the ontology data
      
      incProgress(.2, detail = paste("select PCA-derived terms"))
    
    ### grab LION-terms that are enriched in both PCA directions
    LION_term_by_PCA <-
      unique(unlist(sapply(LION_term_by_PCA[1:isolate(input$PCAn)], function(n) {
        n[['LION_IDs']][!is.na(n[['LION_IDs']])]
      }, simplify = FALSE)))
    
    
    LUT <-
      data.frame(ID = names(ONTdata@termName),
                 ### look-up table for LION-IDs and descriptions
                 Discription = ONTdata@termName)
    

    incProgress(.6, detail = paste("remove redundant terms"))
    ####  limiting by LION-term selection
    if (isolate(input$LIONselection)) {
      ### switch for LION-selection
      TermsOfInterest <-
        get_selected(isolate(input$tree), format = "names")
      TermsOfInterest <- unique(c(unlist(TermsOfInterest),
                                  unlist(sapply(TermsOfInterest, function(element) {
                                    attr(element, "ancestry")
                                  })))) ## vector with LION-names
      ## convert to LION IDs and limit
      LION_term_by_PCA <-
        LION_term_by_PCA[LION_term_by_PCA %in% LUT$ID[match(TermsOfInterest, LUT$Discription)]]
    }
    
    ## limiting redundant/similar parent terms
    
    if (isolate(input$RemoveRedundantTerms)) {
      ### switch for LION-selection
      ONT_DAG <- ONTdata@graph
      ONT_DAG_lev <- buildLevels(ONT_DAG)
      DAG.env <- ONT_DAG_lev$nodes2level
      
      DAGlevel <-   sapply(LION_term_by_PCA, function(id) {
        DAG.env[[id]]
      })
      
      lipidsInTerms <- genesInTerm(ONTdata)
      
      TermContent <-
        sapply(LION_term_by_PCA, function(id) {
          lipidsInTerms[[id]]
        })
      
      test_similarity <-
        sapply(names(TermContent), function(term_i) {
          sapply(names(TermContent), function(term_j) {
            if (term_i != term_j) {
              if (length(TermContent[[term_i]]) == length(TermContent[[term_j]])) {
                same <- all(TermContent[[term_i]] %in% TermContent[[term_j]])
              } else {
                same <- FALSE
              }
              
              
              outputList <- list(
                term_a = term_i,
                term_b = term_j,
                isSame = same,
                term_a_level = DAGlevel[names(DAGlevel) == term_i],
                term_b_level = DAGlevel[names(DAGlevel) == term_j]
              )
              
              
              outputList$remove <- ""
              
              if ((outputList$term_a_level - outputList$term_b_level) == 1) {
                outputList$remove <- outputList$term_b
              }
              if ((outputList$term_a_level - outputList$term_b_level) == -1) {
                outputList$remove <- outputList$term_a
              }
              
              outputList
              
            } else {
              list(
                term_a = term_i,
                term_b = term_j,
                remove = "",
                isSame = FALSE
              )
            }
          }, simplify = FALSE)
        })
      

      test_similarity <-
        test_similarity[unlist(lapply(test_similarity, function(n) {
          n[['isSame']]
        }))]
      
      test_similarity <-
        unique(unlist(lapply(test_similarity, function(n) {
          n[['remove']]
        })))
      
      LION_term_by_PCA <-
        LION_term_by_PCA[!LION_term_by_PCA %in% test_similarity]
      
    }
    
    
    
    ###
    lipidsInTerms <- genesInTerm(ONTdata)
    
    incProgress(.8, detail = paste("generate matrix for heatmap"))
    
    #browser()
    LION_term_by_PCA_list <-
      ### mean signal per sample per LION-term
      sapply(LION_term_by_PCA, function(term) {
        colMeans(input_data$species_scaled[unique(lipidExistance$input) %in% lipidsInTerms[[term]], ])
      }, simplify = FALSE)
    
    if(length(LION_term_by_PCA_list)>1){
   
    
    ### LION term names
    
    names(LION_term_by_PCA_list) <-
      paste(LUT$Discription[match(names(LION_term_by_PCA_list), LUT$ID)],
            " (",
            names(LION_term_by_PCA_list),
            ")",
            sep = "")
    ### contruct matrix
    LION_term_by_PCA_matrix <-
      do.call("rbind", LION_term_by_PCA_list)
    
    ## reorder by condition re-ordering
    LION_term_by_PCA_matrix <- do.call("cbind",sapply(isolate(input$ranked_conditions), function(condition_i){
      LION_term_by_PCA_matrix[,input_data$meta[1,] == condition_i]
    }, simplify = FALSE))
    
    
    ### annotation matrix for pheatmap
    annot <-
      data.frame(samples = factor(unlist(input_data$meta[1,][match(colnames(LION_term_by_PCA_matrix), input_data$meta[2,])])))
    
    
    rownames(annot) <- colnames(LION_term_by_PCA_matrix)
    
    color_scheme <-  list(samples =
                            sapply(1:length(input_data$conditions), function(i) {
                              isolate(input[[paste0("color", i)]])
                            }))
    names(color_scheme$samples) <- input_data$conditions
    color_scheme$samples <- color_scheme$samples[    match(isolate(input$ranked_conditions), names(color_scheme$samples))]
    
    ### scale again
    LION_term_by_PCA_matrix <-
      t(scale(t(LION_term_by_PCA_matrix)))
    
    
      LION_term_by_PCA_matrix[is.nan(LION_term_by_PCA_matrix)] <- 0
    
    ScaleColors <-
      #colorRampPalette(c("#DBB51A", "lightgray", "red"))(100)[seq(from = -max(abs(c(
      colorRampPalette(c(isolate(input$`down-regulated`), "lightgray", isolate(input$`up-regulated`)))(100)[seq(from = -max(abs(c(
        min(LION_term_by_PCA_matrix),
        max(LION_term_by_PCA_matrix)
      ))),
      to = max(abs(c(
        min(LION_term_by_PCA_matrix),
        max(LION_term_by_PCA_matrix)
      ))),
      length.out = 100) >= min(LION_term_by_PCA_matrix) &
        seq(from = -max(abs(c(
          min(LION_term_by_PCA_matrix),
          max(LION_term_by_PCA_matrix)
        ))),
        to = max(abs(c(
          min(LION_term_by_PCA_matrix),
          max(LION_term_by_PCA_matrix)
        ))),
        length.out = 100) <= max(LION_term_by_PCA_matrix)]
    
    incProgress(.9, detail = paste("generate heatmap"))
    
      heatmap_plot <-
        pheatmap(
          mat = LION_term_by_PCA_matrix,
          color = ScaleColors ,
          ## palet for data
          cluster_rows = isolate(input$ClusterTerms),
          cluster_cols = isolate(input$ClusterSamples),
          cutree_rows = min(dim(LION_term_by_PCA_matrix)[1], ifelse(isolate(input$CutClusters),isolate(input$clusternr),NA)),
          #cutree_rows = ifelse(isolate(input$CutClusters),isolate(input$clusternr),NA),
          annotation_colors =  color_scheme,
          annotation_col = annot,
          fontsize = 8,
          gaps_col = (which(!duplicated(annot)) - 1)[-1]
        )                          ## separate samples
 

  
    
    to_ggplot_display <- as.ggplot(heatmap_plot)
    #### end heatmap insert
    
    ### report
    
    lipidsInTerms_dataset <-
      lipidsInTerms[match(gsub("\\)|\\(", "", regmatches(
        names(LION_term_by_PCA_list),
        regexpr("\\(LION:.+\\)$", names(LION_term_by_PCA_list))
      )),  names(lipidsInTerms))]
    
    names(lipidsInTerms_dataset) <- names(LION_term_by_PCA_list)
    
    lipidsInTerms_dataset <-
      data.frame(
        term = names(lipidsInTerms_dataset),
        nr = sapply(lipidsInTerms_dataset, function(term) {
          length(term)
        }),
        lipids = sapply(lipidsInTerms_dataset, function(term) {
          paste(gsub("^\\[.+\\] ", "", term), collapse = "; ")
        })
      )
    
    names(lipidsInTerms_dataset) <- NULL
    
    } else {
      to_ggplot_display <- ggplot()+
        labs(title = "No LION-term meets the criteria")+
        theme_minimal()
      lipidsInTerms_dataset <- NULL
    }
    ###
    })
      
    } else {   ## when input$heatmap_units == "lipid species"
      ### contruct matrix
      LION_term_by_PCA_matrix <- input_data$species_scaled
      
      ## reorder by condition re-ordering
      LION_term_by_PCA_matrix <- do.call("cbind",sapply(isolate(input$ranked_conditions), function(condition_i){
        LION_term_by_PCA_matrix[,input_data$meta[1,] == condition_i]
      }, simplify = FALSE))
      
      
      ### annotation matrix for pheatmap
      annot <-
        data.frame(samples = factor(unlist(input_data$meta[1,][match(colnames(LION_term_by_PCA_matrix), input_data$meta[2,])])))
      
      
      rownames(annot) <- colnames(LION_term_by_PCA_matrix)
      
      color_scheme <-  list(samples =
                              sapply(1:length(input_data$conditions), function(i) {
                                isolate(input[[paste0("color", i)]])
                              }))
      names(color_scheme$samples) <- input_data$conditions
      color_scheme$samples <- color_scheme$samples[    match(isolate(input$ranked_conditions), names(color_scheme$samples))]
      
      ScaleColors <-
        colorRampPalette(c(isolate(input$`down-regulated`), "lightgray", isolate(input$`up-regulated`)))(100)[seq(from = -max(abs(c(
          min(LION_term_by_PCA_matrix),
          max(LION_term_by_PCA_matrix)
        ))),
        to = max(abs(c(
          min(LION_term_by_PCA_matrix),
          max(LION_term_by_PCA_matrix)
        ))),
        length.out = 100) >= min(LION_term_by_PCA_matrix) &
          seq(from = -max(abs(c(
            min(LION_term_by_PCA_matrix),
            max(LION_term_by_PCA_matrix)
          ))),
          to = max(abs(c(
            min(LION_term_by_PCA_matrix),
            max(LION_term_by_PCA_matrix)
          ))),
          length.out = 100) <= max(LION_term_by_PCA_matrix)]
      
      heatmap_plot <-
        pheatmap(
          mat = LION_term_by_PCA_matrix,
          color = ScaleColors ,
          cluster_rows = isolate(input$ClusterTerms),
          cluster_cols = isolate(input$ClusterSamples),
          cutree_rows = min(dim(LION_term_by_PCA_matrix)[1], ifelse(isolate(input$CutClusters),isolate(input$clusternr),NA)),
          annotation_colors =  color_scheme,
          annotation_col = annot,
          fontsize = 8,
          gaps_col = (which(!duplicated(annot)) - 1)[-1]
        )                          ## separate samples
      
      
      to_ggplot_display <- as.ggplot(heatmap_plot)
      #### end heatmap insert
      
      ### report
      lipidsInTerms_dataset <- NULL
      
      
    }
   
    list(to_ggplot_display = to_ggplot_display,
         lipidsInTerms_dataset = lipidsInTerms_dataset)
    
  })
      
  ### mapping percentage
  output$mapping_percentage <- renderText({          
    cat("before mapping %..")
    if (is.null(dataBlock()$lipidExistance)){  
      output_text <- "Unexpected input"  
    } else {
      
      matching_statistics <- dataBlock()$matching_statistics
      
      color_df <- data.frame(match = c("direct matching","smart matching", "fatty acid-association prediction",""),
                             color = c("#4f4d4d","#9c3c2d","#c4942b","#bdbbbf"))
      
      #paste("<br><br><font color='",color_df$color,"'>",color_df$match,"</font>" , sep ="", collapse = "<br>")
     
      output_text <- paste("<i>",
        matching_statistics$matched, " out of ",  matching_statistics$total,  " ", "(", round(matching_statistics$percent, digits = 2), "%", ")",
        " identifiers are matched to LION. </i><br><br>",   
        paste("<font color='",color_df$color,"'>",color_df$match,"</font>" , sep ="", collapse = "<br>"),
        sep = "")
    }
    cat("mapping % done\n")
    output_text
  })   
      
      
      
      ###
  output$downloadInputCNTRL <- renderUI({
       if(is.null(dataBlock()$lipidExistance)){} else {
         downloadButton("downloadInput", "Download input table")
       }
    })
  
 
 
  output$value <- renderFormattable({
    table <- dataBlock()$lipidExistance_toview   #lipidExistance
    formattable(table,align = c("l","l","l"))
  })     #, sanitize.text.function=identity)
      
      output$downloadInput <- downloadHandler(
        filename = function() {
          paste("LION-input-job",isolate(input$submitB), ".csv", sep="")
          
        },
        content = function(file) {
          write.csv(dataBlock()$lipidExistance, 
                    file, row.names = FALSE,quote = TRUE)
        }
      )
      
      
      output$downloadTable <- downloadHandler(
        filename = function() {
          paste("LION-analysis-job",isolate(input$submitB), ".csv", sep="")
          
        },
        content = function(file) {
          write.csv(visHeatmapBlock()$lipidsInTerms_dataset, 
                    file, row.names = FALSE,col.names = TRUE, quote = TRUE)
        }
      )
      
    
      output$downloadPlotPNG <- downloadHandler(
        filename = function() { 
          paste("LION-heatmap-plot-job",isolate(input$submitB),  '.png', sep='') },
        content = function(file) {
          ggsave(file,visHeatmapBlock()$to_ggplot_display,   device = 'png', dpi=300, width = input$figureWidth, height =input$figureHeight, units = "in")
        }
      )
      output$downloadPlotSVG <- downloadHandler(
        filename = function() { 
          paste("LION-heatmap-plot-job",isolate(input$submitB),  '.svg', sep='') },
        content = function(file) {
          ggsave(file,visHeatmapBlock()$to_ggplot_display,   device = svg,  width = input$figureWidth, height =input$figureHeight, units = "in")
        }
      )
      
      height_plot <- reactive({
          input$showPlot

        if(isolate(input$heatmap_units == "lipid species")){
          h <- isolate(max(c(50,length(isolate(input_data()$IDs))), na.rm = TRUE) * 9)
          
        } else {
          h <- isolate(max(c(50,dataBlock()$nr_of_LION_terms_by_PCA[input$PCAn]), na.rm = TRUE) * 9)
        }
         
        return(h)
        
      })
     
        
      output$ontology.graph <- renderPlot({
        
        visHeatmapBlock <- visHeatmapBlock()
        visHeatmapBlock$to_ggplot_display
      }, height = height_plot)
      
      
      ## email
      observe({
        if(is.null(input$send) || input$send==0) return(NULL)
        from <- isolate(input$from)
        to <- "xxxx" 
        subject <- isolate(input$subject)
        msg <- isolate(input$message)
        
        sendEmail(subject = subject, from = from, mail_message = msg)
        showModal(modalDialog(
          footer = modalButton("Ok"),
          size = "m",
          easyClose = TRUE,
          title = paste("Confirmation:",subject), 
          p("Thank you for contacting us. Your message was sent successfully.")
          
        ))
      })
      
      
      
      output$tree <- renderTree({        ### LION-tree
        LIONstructure
      })
      
      
}
