#' Author: Anne Muldoon
#' Title: app.R for final_project
#' Purpose: Create an interactive Rshiny app for explorinig data & visulization 
#'          of RNAseq Data
if (!require("colourpicker")) install.packages("colourpicker")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("DT")) install.packages("DT")
if (!require("sortable")) install.packages("sortable")
library(shiny)
library(dplyr)
library(bslib)
library(ggplot2)
library(pheatmap)
library(colourpicker) # you might need to install this
library(sortable)
options(shiny.maxRequestSize = 15 * 1024^2)
#dev.off()
ui <- fluidPage(
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
                 sidebarPanel(
                   fileInput("metadata", "Upload metadata file",
                             accept=c('.tsv')
                             ),
                   actionButton("load_metadata", "Load Metadata"),
                   radioButtons("feature",
                   label = "Choose which part of the experiment to show 
                   metadata summary for",
                   choices = c('Cell Differentiation','Explanted Cardiac Myocytes',
                               'Heart Ventricle','Heart Apex'))),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("summary"),
                            tableOutput("filtered")),
                   tabPanel("Table",
                            tableOutput("sortable")
                            ),
                   tabPanel("Plots",
                            plotOutput("categorical",height = "400px", width = "600px")
                            )
                   )
               )
      ),
      tabPanel("Counts",
               #This component allows the user to choose different gene 
               # filtering thresholds and assess their effects using diagnostic
               # plots of the counts matrix.
               sidebarPanel(
                 fileInput("norm_counts", "Upload normalized counts file",
                           accept=c('text/csv', 'text/comma-separated- values,text/plain', 
                                    '.csv')),
                 actionButton("load_data", "Load Data"),
                 uiOutput("variance"),
                 uiOutput("non_zero")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table",
                            tableOutput("counts_table")
                   ),
                   tabPanel("Scatter Plot",
                            plotOutput("counts_scatter_1",height = "400px", width = "600px"),
                            plotOutput("counts_scatter_2",height = "400px", width = "600px")
                   ),
                   tabPanel("Clustered Heatmap",
                            plotOutput("counts_heatmap",height = "400px", width = "600px")
                   ),
                   tabPanel("Scatter Plot of PCA",
                            sidebarPanel(
                              uiOutput("components_x"),
                              uiOutput("components_y")
                            ),
                            mainPanel(
                              plotOutput("counts_pca",height = "400px", width = "600px")
                            )
                   )
                 )
               )
      ),
      tabPanel("Differential Expression",
               sidebarPanel(
                 fileInput("deseq", "Upload deseq2 results file",
                           accept=c('.csv')
                 ),
                 actionButton("load_results", "Load Results")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            DT::dataTableOutput("deseq_summary"),
                   ),
                   tabPanel("Volcano Plot",
                            sidebarPanel(
                              uiOutput("color1"),
                              uiOutput("color2"),
                              uiOutput("logFC"),
                              uiOutput("padj")
                            ),
                            plotOutput("volcano",height = "400px", width = "600px")
                   ),
                 )
               )
      ),
      tabPanel("FGSEA",
               sidebarPanel(
                  fileInput("fgsea_res", "Upload fgsea data file",
                            accept=c('.csv')
                 ),
                  actionButton("load_david", "Load fgsea file")
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Top Pathways",
                             sidebarPanel(
                               uiOutput("top_pathways")
                             ),
                             plotOutput("DAVID",height = "400px", width = "600px")
                    ),
                    tabPanel("Table",
                             sidebarPanel(
                                uiOutput("fgsea_padj"),
                              radioButtons("fgsea_filter",
                                           label = "Choose to show positive, negative, or all NES pathways",
                                           choices = c('Positive','Negative','All')),
                                actionButton("apply_fgsea_filter", "Apply Filter"),
                                downloadButton("download_button", label = "Download", class = NULL)
                             ),
                             mainPanel(
                               DT::dataTableOutput("fgsea_table")
                             )
                    ),
                    tabPanel("Scatter Plot",
                             sidebarPanel(
                               uiOutput("scatter_fgsea_padj"),
                             ),
                             mainPanel(
                               plotOutput("fgsea_scatter",height = "400px", width = "600px")
                             )
                             #plot("fgsea_scatter")
                    )
                 )
                 )
        )
      )
    )
)
server <- function(input, output, session){ 
  rv <- reactiveValues(metadata_df = NULL, counts_df = NULL, 
                       deseq_df=NULL, fgsea_df=NULL)
  
  observeEvent(input$load_metadata, {
    req(input$metadata)
    tsv<-read.delim(req(input$metadata$datapath, header = TRUE,sep='\t'))
    rv$metadata_df<-as.data.frame(tsv)
  })
  observeEvent(input$load_data, {
    req(input$norm_counts)
    csv<-read.csv(req(input$norm_counts$datapath, header = TRUE))
    rv$counts_df<-as.data.frame(csv)
  })
  observeEvent(input$load_results, {
    req(input$deseq)
    csv<-read.csv(req(input$deseq$datapath, header = TRUE))
    rv$deseq_df<-as.data.frame(csv)
  })
  observeEvent(input$load_david, {
    req(input$fgsea_res)
    csv<-read.csv(req(input$fgsea_res$datapath, header = TRUE))
    rv$fgsea_df<-as.data.frame(csv)
    
    output$fgsea_padj<-renderUI({
      sliderInput(inputId = "padj_fgsea", min = 0, max = 1,label = 
                    "Select the p_adjusted vals to filter:", 
                  value = 0.5, step = 0.1)
    })
  })
  
  observeEvent(input$apply_fgsea_filter, {
    req(input$padj_fgsea, rv$fgsea_df,input$fgsea_filter)

    filtered_data <- rv$fgsea_df %>% 
      filter(padj <= input$padj_fgsea)
    if(input$fgsea_filter=='Positive'){
      filtered_data<-filtered_data%>%filter(NES>=0)
    } else if (input$fgsea_filter=='Negative'){
      filtered_data<-filtered_data%>%filter(NES<0)
    } else{
      filtered_data<-filtered_data
    }
    output$fgsea_table <- DT::renderDataTable({
      datatable(filtered_data, options = list(pageLength = 10, searchable = TRUE))
      })
  })
  filtered_data <- reactive({
    req(rv$fgsea_df)
    req(input$fgsea_scatter_padj)
    filtered<-rv$fgsea_df %>%
      mutate(log_padj = -log10(padj)) %>%
      mutate(color = ifelse(padj <= input$fgsea_scatter_padj, "highlight", "grey"))
  })

  output$sortable<-renderTable({
    req(rv$metadata_df)
    rv$metadata_df
  })
  filter_summary <- reactive({
    req(rv$metadata_df, input$feature)
    feature <- input$feature
    df<-rv$metadata_df
    if(feature=='Cell Differentiation'){
      df<-df%>%filter(startsWith(sample_name,'ESC')|
                      startsWith(sample_name,'MES')|
                      startsWith(sample_name,'CP')|
                      startsWith(sample_name,'CM'))
    } else if(feature=='Explanted Cardiac Myocytes'){
      df<-df%>%filter(startsWith(sample_name,'ex'))
    } else if(feature=='Heart Ventricle'){
      df<-df%>%filter(startsWith(sample_name,'vP')|
                      startsWith(sample_name,'vA')|
                      startsWith(sample_name,'iP'))
    }else{
      df<-df%>%filter(startsWith(sample_name,'vD')|
                      startsWith(sample_name,'iD'))
    }
    return(df)
  })
  output$filtered<-renderTable(filter_summary())
  build_summary<-reactive({
    df<-filter_summary()
    new_df<-data.frame('Number_of_Samples'=nrow(df),
                       Timepoints=n_distinct(df$timepoint))
    #%>%rename('Number of Samples'='Number_of_Samples')
      #rename('Number of Samples'=Num)
    if(any(startsWith(df$sample_name,'v')==TRUE)){
      new_df<-new_df%>%mutate('Replicates'='2: heart, 3: isolated cardiomyocytes')
    }else{
        new_df<-new_df%>%mutate('Replicates'=max(df$replicate))
    }
    return(new_df)
  })
  output$summary<-renderTable({
    req(build_summary())
    build_summary()
    })
  categorical_plot<-reactive({
    req(rv$metadata_df)
    df<-rv$metadata_df%>%
      filter(!is.na(sample_name), !is.na(description), !is.na(timepoint))%>%
      mutate(Experiment=case_when(startsWith(description,'emb')
                                  |startsWith(description,'mes')
                                  |startsWith(description,'car') ~ "Cell Differentiation",
                                  startsWith(description,'exp') ~ "Explanted Cardiac Myocytes",
                                  endsWith(description,'ventrical') ~ "Heart Ventricle",
                                  endsWith(description,'surgery') ~ "Heart Apex Surgery"))
    df_cell<-df%>%filter(startsWith(sample_name,'ESC')|
                           startsWith(sample_name,'MES')|
                           startsWith(sample_name,'CP')|
                           startsWith(sample_name,'CM'))%>%dplyr::select(timepoint,Experiment)
    df_myo<-df%>%filter(startsWith(sample_name,'ex'))%>%dplyr::select(timepoint,Experiment)
    df_vent<-df%>%filter(startsWith(sample_name,'vP')|
                           startsWith(sample_name,'vA')|
                           startsWith(sample_name,'iP'))%>%dplyr::select(timepoint,Experiment)
    df_apex<-df%>%filter(startsWith(sample_name,'vD')|
                           startsWith(sample_name,'iD'))%>%dplyr::select(timepoint,Experiment)
    total_df<-bind_rows(df_cell,df_myo,df_vent,df_apex)%>%
      group_by(timepoint, Experiment) %>%
      summarise(num_samples = n(), .groups = "drop")%>%arrange(timepoint)
    total_df$timepoint <- factor(total_df$timepoint, levels = c("day0", "day1",
                                                                "day2","day3",
                                                                "day4","day7",
                                                                "day10","week8"))
    if(nrow(total_df) > 0) {
      plt<-ggplot(total_df, aes(x = timepoint, y = num_samples, fill = Experiment)) +
        geom_bar(stat = "identity", position = "stack") +  # Use 'dodge' to separate the bars
        labs(
          title = "Number of Samples per Timepoint",
          x = "Timepoint",
          y = "Number of Samples",
          fill = "Experiment"
        ) +
        theme_minimal()
    } else {
      plt <- ggplot() + labs(title = "No valid data available")
    }
    return(plt)
  })
  output$categorical<-renderPlot(categorical_plot())
  
  output$variance<-renderUI({
    sliderInput(inputId = "var", min = 0, max = 100,label = 
                  "Include genes with at least percentile variance of:", 
                value = 50, step = 1)
  })
  #non-zero: max value should be total count of samples (reactive?)
  output$non_zero<-renderUI({
    sliderInput(inputId = "no_zero", min = 0, max = 36,label = 
                  "Include genes with number of samples being non-zero:", 
                value = 18, step = 1)
  })
  output$components_x<-renderUI({
    sliderInput(inputId = "pc_x", min = 1, max = 20,label = 
                  "Choose which principal component to plot on the x-axis:", 
                value = 1, step = 1)
  })
  output$components_y<-renderUI({
    sliderInput(inputId = "pc_y", min = 1, max = 20,label = 
                  "Choose which principal component to plot on the y-axis:", 
                value = 1, step = 1)
  })
  output$color1<-renderUI({
    colourInput("color_1", label="Base point color", "#22577A")
  })
  #get input from colorInput, will feed this to volcano plot
  output$color2<-renderUI({
    colourInput("color_2", label="Highlight point color", "#FFCF56")
  })
  output$logFC<-renderUI({
    sliderInput(inputId = "log_FC", min = 0, max = 10,label = 
                  "Color signifance based on log2FoldChange value (up/downregulation):", 
                value = 5, step = 1)
    ##has to be positive and. negative in plot
  })
  output$padj<-renderUI({
    sliderInput(inputId = "p_adj", min = -50, max = -1,label = 
                  "Select the magnitude of the p adjusted coloring:", 
                value = -25, step = 1)
  })
  output$top_pathways<-renderUI({
    sliderInput(inputId = "paths", min = 1, max = 10,label = 
                  "Top Number of Pathways to plot:", 
                value = 5, step = 1)
  })
  
  output$download_button<-downloadHandler(
    filename = function() {
      paste('filtered_data-', Sys.Date(), '.csv', sep='')
      },
    content = function(con) {
      req(rv$fgsea_df)
      req(input$padj_fgsea)
      req(input$fgsea_filter)
      df <- rv$fgsea_df %>% filter(padj < input$padj_fgsea)
      if(input$fgsea_filter=='Positive'){
        df<-df%>%filter(NES>=0)
      } else if (input$fgsea_filter=='Negative'){
        df<-df%>%filter(NES<0)
      } else{
        df<-df
      }
      write.csv(df, con)
    }
  )
  output$scatter_fgsea_padj<-renderUI({
    sliderInput(inputId = "fgsea_scatter_padj", min = 0, max = 1,label = 
                  "Select the value of the p adjusted coloring:", 
                value = 0.5, step = 0.1)
  })
  calculate_filters <- reactive({#,var,non_zero) {
    #number of samples
    # total number of genes
    # number and % of genes passing current filter
    # number and % of genes not passing current filter
    req(rv$counts_df)
    req(input$var)
    req(input$no_zero)
    
    variance<-isolate(input$var)
    non_zero<-isolate(input$no_zero)
    
    old_df<-rv$counts_df%>%select(Gene)
    new_df<-rv$counts_df%>%
      select(-Gene,-GeneID,-Coordinates)%>%
      mutate(non_zero_counts=rowSums(. != 0),
             variances=apply(.,MARGIN=1,var))
    percentile <- quantile(new_df$variances, probs = variance / 100,na.rm=TRUE)
    new_df$pass_filter <- new_df$non_zero_counts >= non_zero & new_df$variances >= percentile
    total_df<-bind_cols(old_df,new_df)
    total_df<-as.data.frame(total_df)
    return(total_df)
  })
  make_count_table<-function(df){
    req(rv$counts_df)
    new_df<-data.frame('Number_of_Genes'=nrow(df))%>%
      mutate('Number of Samples'=as.integer(ncol(df)-4))%>%
      mutate('Number of Filtered Genes'=sum(df$pass_filter==TRUE))%>%
      mutate('Percent of Filtered Genes'=((sum(df$pass_filter==TRUE))/nrow(df))*100)%>%
      mutate('Number of Genes not passing filter'=sum(df$pass_filter==FALSE))%>%
      mutate('Percent of Genes not passing filter'=((sum(df$pass_filter==FALSE))/nrow(df))*100)
    #%>%rename('Number of Genes'='Number_of_Genes')
    return(new_df)
  }
  draw_scatterplot_1<-function(df){
    plot_df<-df%>%select(-Gene,-non_zero_counts,-variances,-pass_filter)
    df<-df%>%mutate(median_count=apply(plot_df,MARGIN=1,median))
    plt<-ggplot(df, aes(x = median_count, y = non_zero_counts, color = pass_filter)) +
      geom_point() +  # Use 'dodge' to separate the bars
      labs(
        title = "Median Count vs Non-zero Counts",
        x = "Median Count",
        y = "Non-Zero Counts",
        fill = "pass_filter"
      ) +
      theme_minimal()

    return(plt)
  }
  draw_scatterplot_2<-function(df){
    plot_df<-df%>%select(-Gene,-non_zero_counts,-variances,-pass_filter)
    df<-df%>%mutate(median_count=apply(plot_df,MARGIN=1,median))
    plt<-ggplot(df, aes(x = median_count, y = log(variances), color = pass_filter)) +
      geom_point() +  # Use 'dodge' to separate the bars
      labs(
        title = "Median Count vs log(Variance)",
        x = "Median Count",
        y = "log(Variance)",
        fill = "pass_filter"
      ) +
      theme_minimal()
    return(plt)
  }
  make_counts_heatmap<-function(df){
    #looks better than transformed counts
    #if i have time: split it into 4based on. experiment
    print(head(df))
    df<-as.data.frame(df)%>%filter(pass_filter==TRUE)%>%
      dplyr::select(-Gene,-non_zero_counts,-variances,-pass_filter)
    print(df)
    plt<-pheatmap(df,
                  scale = "row",
                  color = colorRampPalette(c("blue", "white", "red"))(100),
                  main = "Gene Expression Heatmap",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE,
                  show_colnames = TRUE)
    # plt<-pheatmap(df,
    #             scale = "row",  # Normalize the data by rows (genes)
    #             col = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette from blue to red
    #             main = "Gene Expression Heatmap",  # Title
    #             xlab = "Samples",  # X-axis label
    #             ylab = "Genes",  # Y-axis label
    #             margins = c(5, 10),  # Adjust margins for labels
    #             cexCol = 0.8,  # Control the size of the column labels
    #             cexRow = 0.8,  # Control the size of the row labels
    #             dendrogram = "both"  # Cluster both rows and columns
    # )
    return(plt)
  }
  make_counts_PCA<-reactive({
    req(input$pc_x)
    req(input$pc_y)
    xaxis<-isolate(input$pc_x)
    yaxis<-isolate(input$pc_y)
    df<-calculate_filters()%>%
      filter(pass_filter==TRUE)%>%
      select(-Gene,-non_zero_counts,-variances,-pass_filter)
    pca_results <- prcomp(t(df), center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca_results$x)
    colnames(pca_df) <- paste0("PC", 1:ncol(pca_df))
    pca_df<-pca_df%>%mutate(Experiment=rownames(pca_df))
    
    ggplot(pca_df, aes(x = .data[[paste0("PC", xaxis)]], y = .data[[paste0("PC", yaxis)]],color=Experiment)) +
      geom_point() +  # Color points based on PC1 values
      labs(title = "PCA of Filtered Genes",
           x = paste("Principal Component ", xaxis),
           y = paste("Principal Component ",yaxis)) +
      theme_minimal()
  })

  output$counts_table<-renderTable(make_count_table(calculate_filters()))
  output$counts_scatter_1<-renderPlot(draw_scatterplot_1(calculate_filters()))
  output$counts_scatter_2<-renderPlot(draw_scatterplot_2(calculate_filters()))
  output$counts_heatmap<-renderPlot(make_counts_heatmap(calculate_filters()))
  output$counts_pca<-renderPlot(make_counts_PCA())
  
  output$deseq_summary<-DT::renderDataTable({
    req(rv$deseq_df)
    df<-rv$deseq_df%>%relocate(GeneID,.before=baseMean)
    datatable(df,
              options = list(
                pageLength = 5,          # Number of rows to show per page
                lengthChange = FALSE,    # Disable the option to change page length
                autoWidth = TRUE,        # Automatically adjust column widths
                columnDefs = list(
                  list(targets = 0:2, searchable = TRUE)  # Make columns searchable
                )
              ))
  })
  volcano_plot <-function(df, color1, color2,LFC,p) {
    truth<-df$padj<10**p & abs(df$log2FoldChange)>LFC#& (dataf$log2FoldChange<-LFC|dataf$log2FoldChange>LFC)
    
    plt<-ggplot(df,aes(x=log2FoldChange,y=-log10(padj),color=truth)) +
      geom_point()+scale_colour_manual(values = c(color1,color2))+
      labs(color = "Gene is Differentially Expressed")+
      ggtitle("Volcano Plot of Differentially Expressed Genes")
    return(plt)
  }
  #maybe add a table showing top ~10 DEGs
  #also maybe add other diff express
  #also maybe add venn diagram
  output$volcano<-renderPlot({
    req(input$color_1)
    req(input$color_2)
    req(rv$deseq_df)
    df<-rv$deseq_df%>%
      filter(!is.na(log2FoldChange), !is.na(padj))
    volcano_plot(df,input$color_1, input$color_2,input$log_FC,input$p_adj)
  })
  plot_pathways<-function(fgsea_results, num_paths){
    fgsea_results<-fgsea_results%>%
      arrange(NES)%>%
      select(pathway,padj,NES)
    plt_high<- top_n(fgsea_results, num_paths, NES)
    plt_low<-top_n(fgsea_results, -num_paths, NES)
    plt_res<-plt_high%>%bind_rows(plt_low)%>%mutate(Color = ifelse(NES <0, "red","blue"))
    plt_res$pathway <- factor(plt_res$pathway, levels = plt_res$pathway)
    #%>%slice(-10:-1)
    plt<-plt_res%>%ggplot(aes(y=pathway,x=NES,fill=Color)) +
      geom_bar(stat = "identity")+
      theme(axis.text=element_text(size=10),axis.title=element_text(size=10))+
      scale_fill_identity(guide = FALSE)
    return(plt)
  }
  output$DAVID<-renderPlot({
    req(input$paths)
    req(rv$fgsea_df)
    df<-rv$fgsea_df%>%filter(!is.na(pathway),!is.na(padj),!is.na(padj))
    plot_pathways(df,input$paths)
  })

  output$fgsea_table <- DT::renderDataTable({
    req(rv$fgsea_df)
    datatable(rv$fgsea_df, options = list(pageLength = 10, searchable = TRUE))
  })
  output$fgsea_scatter<-renderPlot({
    df<-filtered_data()
    
    ggplot(df, aes(x = NES, y = log_padj, color = color)) +
      geom_point(size = 3) +
      scale_color_manual(values = c("highlight" = "blue", "grey" = "grey")) +
      labs(x = "NES (Normalized Enrichment Score)", y = "-log10(Adjusted P-value)", 
           title = "Scatter Plot of NES vs -log10(padj)") +
      theme_minimal() +
      theme(legend.position = "none")  
  })

}

shinyApp(ui = ui, server = server)





