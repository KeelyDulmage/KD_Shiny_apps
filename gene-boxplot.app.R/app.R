# Shiny web app to visualize UMI gene expression data as boxplots

library(shiny)
library(pheatmap)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

filtdf <- readRDS()

# Define UI for application that draws boxplots for genes
ui <- fluidPage(
    titlePanel("UMI gene expression data"),

    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "sel_genes",label = h4("Gene names"),
                        value = "Cdh1, Epcam, Ptprc"),
            actionButton("go", "Go"),
            helpText("Separate gene IDs by comma"),
            radioButtons("radio", h4("Y Axis Scale"),
                         choices = list("MI Counts" = "MI_count", "Log 2 MI Counts" = "Log2_MI_Count",
                                        "Log 10 MI Counts" = "Log10_MI_Count"), selected = "Log2_MI_Count"),
            htmlOutput("missing_genes")),  

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("g_boxplots", inline = F)
                )
            )
    )
# Define server logic required to draw boxplots
server <- function(input, output) {
        # gene selection here based on sel_genes
        # and also calculate what's missing
    
    get_data <- eventReactive(input$go,{
    sel_geness <- unique(as.character(unlist(strsplit(input$sel_genes, ", "))))
    sel_geness <- tolower(sel_geness)
    sel_genes_2 <- gsub("(^[[:alpha:]])", "\\U\\1", sel_geness, perl=TRUE)
    sel_genes_3 <- paste0("^",c("Cell_Type", sel_genes_2),"$", collapse = "|")
    
    sel_table <- filtdf[, grep(sel_genes_3, colnames(filtdf)), drop = F] # grab genes of interest from the table
    
    missing <- setdiff(sel_genes_2, colnames(sel_table))
    
    # make data long for plotting.  Note, this assumes that genes are in all columns after the first one
    
    sel_table.long <- sel_table %>% 
        tidyr::gather(Gene, MI_count, colnames(sel_table)[2]:colnames(sel_table)[ncol(sel_table)]) %>%
        dplyr::mutate(Log2_MI_Count = log2(MI_count + 1), Log10_MI_Count = log10(MI_count + 1))
    
    obj <- list(data = sel_table.long, miss_list = missing, ngenes = length(sel_genes_2))
    
    obj
    
    })
    
    ngene_height <- function(){
        obj <- get_data()
        return(250*(2 + obj$ngenes) %/% 3)
    }

    
    output$g_boxplots <- renderPlot({
        plot.obj <- get_data()
        # ggplot object here
        ggplot(plot.obj$data, aes_string(x = "Cell_Type", y = paste(input$radio))) +
            geom_boxplot(aes(fill = Cell_Type), col = "gray40", alpha = 0.4, outlier.shape = NA) +
            geom_jitter(width = .1, size = 0.5) +
            facet_wrap(~Gene, ncol = 3, scales = "free") +
            theme(axis.text.x  = element_text(angle = 90), text = element_text(size = 14, face = "bold"))
    }, height = ngene_height)
    
    output$missing_genes <- renderUI({
        text.obj <- get_data()
        if(length(text.obj$miss_list) >= 1){
        HTML(paste0("<b>", "The following genes were not found: ", paste(text.obj$miss_list, collapse = ", "), 
              '<br/>', "Check that you are using A) Mouse ID, and B) Gene name rather than protein.",'<br/>', 
              "If gene is still missing it is likely to have been filtered due to lack of expression.","<b>")
        )
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
