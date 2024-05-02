# 1. Declare functions
{
  write_gephi_gexf <- function(GO_object,label_indexes,output_file) {
    
    df_nodes_all=GO_object$result
    df_nodes_all$label=""
    set_base=which(!is.na(df_nodes_all$community))
    df_nodes_all$label[set_base[label_indexes]]=df_nodes_all$description[set_base[label_indexes]]
    df_nodes_all=df_nodes_all[,c(1:2,ncol(df_nodes_all),3:(ncol(df_nodes_all)-1))]
    GO_object_out=GO_object
    GO_object_out$result=df_nodes_all
    
    df_nodes=GO_object$result
    df_nodes=df_nodes[which(!is.na(df_nodes$community)),]
    df_nodes$label=""
    df_nodes$label[label_indexes]=df_nodes$description[label_indexes]
    df_nodes=df_nodes[,c(1:2,16,3:15)]
    df_edges=GO_object$network
    
    nodes <- data.frame(ID=df_nodes$GO_ID,Label=df_nodes$label)
    edges=df_edges[,1:2]
    edgesWeight=as.numeric(df_edges$Weight)
    
    #df_edges$edge_ID=c(1:nrow(df_edges))
    #edgesId=data.frame(df_edges$edge_ID)
    #edgesAtt=data.frame(df_edges[,4,3])
    df_nodes$significance <- -log10(df_nodes$fdr)
    nodesAtt=df_nodes
    
    # Escribir el archivo .gexf
    write.gexf(nodes=nodes,edges=edges,edgesWeight=edgesWeight,nodesAtt=nodesAtt,output=paste0(output_file,".gexf"))
    return(GO_object_out)
  }
  output_preparer <- function(data_input) {
    data_input$a="A"
    data_input$b="B"
    
    for(i in 1:nrow(data_input)) {
      data_input$a[i] <- paste(data_input$Levels[[i]], collapse = ", ")
      data_input$b[i] <- paste(data_input$Members[[i]][1][[1]], collapse=", ")
    }
    data_input$Levels=data_input$a
    data_input$Members=data_input$b
    data_input=data_input[,1:(ncol(data_input)-2)]
    return(data_input)
  }
  build_network <- function(tab, threshold) {
    links <- data.frame(node_source = character(0), node_target = character(0), P = numeric(0), stringsAsFactors = FALSE)
    
    for (i in 1:(nrow(tab))) {
      for (j in (i):nrow(tab)) {
        elements1 <- unlist(strsplit(tab$Members[i], ","))
        elements2 <- unlist(strsplit(tab$Members[j], ","))
        
        common_elements <- intersect(elements1, elements2)
        P <- length(common_elements) / min(length(elements1), length(elements2))
        
        if (P > threshold) {
          links <- rbind(links, c(tab[i, "Ontology"], tab[j, "Ontology"], P))
        }
      }
    }
    
    colnames(links) <- c("GO_source", "GO_target", "Weight")
    return(links)
  }
  GO=function(query,background_matrix,metadata_matrix,bg=NA,min_GO_level=0,max_GO_level=15,test_mode="enrichment",pvalue_correction_method="BH",min_term_size = 1,max_term_size = Inf,th_hit=0.1,th_hit_OR=2,threshold_links =0.5,merge_duplicates=TRUE){
    
    ## Quality check 1: Are there cols, or raws, duplicated in the background matrix? (In the rownames, dupes are not allowed). Are there genes in the query.
    
    unicos=length(unique(colnames(background_matrix)))
    todos=ncol(background_matrix)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated terms in the background matrix:",todos," submitted, ",unicos," are unique"))
      background_matrix=background_matrix[,which(!duplicated(colnames(background_matrix)))]
    }
    
    unicos=length(unique(rownames(background_matrix)))
    todos=nrow(background_matrix)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated genes in the background matrix:",todos," submitted, ",unicos," are unique"))
      background_matrix=background_matrix[which(!duplicated(rownames(background_matrix))),]
    }
    
    unicos=length(unique(query))
    todos=length(query)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated genes in query:",todos," submitted, ",unicos," are unique"))
      query=query[which(!duplicated(query))]
    }
    
    ## Si había background:
    ## 1. Completo QC: bg no contiene dupes, y todo el query vive dentro del bg.
    ## 2. Subset the matrix with annotations row-wise: only genes within background should be in. AND expand it to include genes in BG not annotated.
    
    ## Si NO había background:
    ## 1. Define the background como el set de genes con alguna anotación EN CUALQUIER NIVEL. (Idea: si un gen tiene una anotación fuera de los niveles target, aunque no la tenga en ninguno de los términos de los niveles target, podría haberla tenido, y por tanto puede y debe entrar en el bg)
    ## 2. Subset the query to live in that background.
    
    ## Finally in both cases:
    ## 3. Subset the annotation matrix to contain the terms that will  be tested, and only those. Declare also query_vector.
    ## 4. locate terms containing exactly the same set of genes in the defined universe through a uniqueness_index attribute for terms
    
    
    if(length(bg)>1)
    {
      ## 1. completamos QC: dupes in BG??
      unicos=length(unique(bg))
      todos=length(bg)
      if(unicos!=todos)
      {
        print(paste("Warning: repeated genes in bg:",todos," submitted, ",unicos," are unique"))
        bg=bg[which(!duplicated(bg))]
      }
      ## 1. Completamos QC: Is query within bg?
      within=length(which(query %in% bg))
      todos=length(query)
      if(within!=todos)
      {
        print(paste("Warning: from ",todos," submitted in the query, only ",within," are in the bg, and therefore will be tested"))
        query=query[which(query %in% bg)]
      }
      #2. Subset the matrix with annotations row-wise:only genes within background should be in...
      background_matrix=background_matrix[which(rownames(background_matrix) %in% bg),]
      extra_genes=bg[which(!bg %in% rownames(background_matrix))]
      #2...AND expand it to include genes in BG not annotated.
      appendix <- data.frame(matrix(FALSE, nrow = length(extra_genes), ncol = ncol(background_matrix)))
      colnames(appendix) <- colnames(background_matrix)
      rownames(appendix) <- extra_genes
      background_matrix=rbind(background_matrix,appendix)
      
      # Check that 2 worked well:
      background_matrix=background_matrix[order(rownames(background_matrix)),]
      bg=bg[order(bg)]
      #print("Vio bg")
    }else{
      # 1. define bg from annotated genes in all levels.
      genes_presence_vector=apply(background_matrix,1,sum)
      annotated_genes=rownames(background_matrix)[which(genes_presence_vector>0)]
      background_matrix=background_matrix[annotated_genes,]
      background_matrix=background_matrix[order(rownames(background_matrix)),]
      bg=rownames(background_matrix)
      ## 2. Subset the query to live in that background.
      query=query[which(query %in% bg)]
      #print("NO VIO bg")
    }
    
    ##3. Subset the annotation matrix to contain the terms that will  be tested, and only those.
    
    #query_vector=data.frame(background_matrix[,1])
    #query_vector[,1]=FALSE
    #colnames(query_vector)="query"
    #rownames(query_vector)=rownames(background_matrix)
    #query_vector[which(rownames(query_vector) %in% query),1]=TRUE
    
    
    ## filter terms: 1st according to levels:
    GO_levels <- min_GO_level:max_GO_level
    mode <- switch(test_mode, "enrichment" = "greater", "depletion"="less", "two.sided"="two.sided")
    terms_within_levels <- c()
    for (i in 1:nrow(metadata_matrix)) {
      if (sum(metadata_matrix$Levels[i][[1]] %in% GO_levels)>0) {
        terms_within_levels <- c(terms_within_levels, i)
      }
    }
    
    metadata_matrix <- metadata_matrix[terms_within_levels,]
    background_matrix=background_matrix[,terms_within_levels]
    ## filter terms: 2nd according to term size:
    
    number_genes_per_term=apply(background_matrix,2,sum)
    metadata_matrix$Size_in_universe=number_genes_per_term
    set_passing=which(number_genes_per_term>=min_term_size & number_genes_per_term<=max_term_size)
    
    metadata_matrix <- metadata_matrix[set_passing,]
    background_matrix=background_matrix[,set_passing]
    
    #4. Identify terms that have exactly the same annotations using a uniqueness_index column
    
    colapsos=apply(background_matrix,2,function(x){paste(x,collapse="")})
    set_first_instances=which(!duplicated(colapsos))
    metadata_matrix$uniqueness_index=0
    metadata_matrix$uniqueness_index[set_first_instances]=c(1:length(set_first_instances))
    
    dupe_instances=unique(colapsos[which(duplicated(colapsos))])
    
    for(i in 1:length(dupe_instances))
    {
      set=which(colapsos==dupe_instances[i])
      value=max(metadata_matrix$uniqueness_index[set])
      metadata_matrix$uniqueness_index[set]=value
    }
    
    ## Now, we load number of genes in query and percentage of the total term that that represents
    query_matrix <- background_matrix[which(rownames(background_matrix) %in% query),]
    
    number_genes_vector <- c()
    percentage_genes_vector <- c()
    
    for (i in 1:nrow(metadata_matrix)) {
      
      number_of_genes <- sum(query_matrix[,metadata_matrix[i,2]])
      percentage <- (number_of_genes/metadata_matrix$Size_in_universe[i]) * 100
      
      number_genes_vector <- c(number_genes_vector, number_of_genes)
      percentage_genes_vector <- c(percentage_genes_vector, percentage)
    }
    
    metadata_matrix$Number_of_genes <- number_genes_vector
    metadata_matrix$Percentage <- percentage_genes_vector
    
    ## Omitiré filtrado por num genes en la intersección, percentage of genes in the interseccion y tb por terminos con tamaño minimo.
    pval_vector <- c()
    OR_vector <- c()
    members_list <- list()
    
    for (i in  1:nrow(metadata_matrix)){
      
      in_both <- sum(query_matrix[,metadata_matrix$Category[i]])
      in_term_out_query <- metadata_matrix$Size_in_universe[i] - in_both
      out_term_in_query <- length(query) - in_both
      out_both <- nrow(background_matrix) - in_both - in_term_out_query - out_term_in_query
      
      contingency_table <- data.frame(
        c(in_both, in_term_out_query),
        c(out_term_in_query,out_both)
      )
      
      test=fisher.test(contingency_table,alternative=mode)
      pval_vector[i] <- test$p.value
      OR_vector[i] <- test$estimate
      
      members <- rownames(query_matrix)[which(query_matrix[,metadata_matrix$Category[i]])]
      members_list[[i]] <- list(members)
    }
    metadata_matrix$pvalue <- pval_vector
    metadata_matrix$Odds_Ratio <- OR_vector
    metadata_matrix$Members <- members_list
    
    #And now we apply the pvalue correction
    print("Nominal p values distribution summary:")
    print(summary(metadata_matrix$pvalue))
    #metadata_matrix <- metadata_matrix[which(is.na(metadata_matrix$pvalue) == FALSE & is.infinite(metadata_matrix$pvalue) == FALSE),]
    
    
    if(pvalue_correction_method=="qval"){
      set=which(!duplicated(metadata_matrix$uniqueness_index))
      metadata_matrix$P.adj=2
      qs=qvalue(metadata_matrix$pvalue[set])
      pi0=qs$pi0
      metadata_matrix$P.adj[set]=qs$qvalues
    }else{
      set=which(!duplicated(metadata_matrix$uniqueness_index))
      metadata_matrix$P.adj=2
      metadata_matrix$P.adj[set]=p.adjust(metadata_matrix$pvalue[set], method=pvalue_correction_method)
      pi0=1
    }
    ## For the sets that are dupes, I copy the same FDR along the entire term: these are not many different terms, only one that is repeated.
    uniqueness_indexes_w_several_terms=unique(metadata_matrix$uniqueness_index[which(duplicated(metadata_matrix$uniqueness_index))])
    for(i in uniqueness_indexes_w_several_terms){
      set=which(metadata_matrix$uniqueness_index==i)
      metadata_matrix$P.adj[set]=min(metadata_matrix$P.adj[set])
    }
    
    ## Order things:
    background_matrix=background_matrix[,order(metadata_matrix$P.adj)]
    metadata_matrix=metadata_matrix[order(metadata_matrix$P.adj),]
    
    length(which(colnames(background_matrix)!=metadata_matrix$Category))
    rownames(metadata_matrix)=metadata_matrix$Ontology
    metadata_matrix <- output_preparer(metadata_matrix)
    metadata_matrix=metadata_matrix[,c(1:11,13,12)]
    
    #And finally we output the results: the background matrix containing only the genes and terms considered in the analyses, the metadata table with the results and pi0
    
    hits_table=metadata_matrix[which(metadata_matrix$P.adj<th_hit & metadata_matrix$Odds_Ratio>th_hit_OR),]
    if(merge_duplicates) hits_table=hits_table[which(!duplicated(hits_table$uniqueness_index)),]
    
    
    network <- build_network(hits_table, threshold=threshold_links)
    #print(network)
    graph <- graph.data.frame(network, directed = FALSE)
    community <- cluster_louvain(graph)
    
    #length(which(community$names==hits_table$Ontology))
    hits_table$community=community$membership
    hits_table=hits_table[order(hits_table$community,hits_table$P.adj),]
    hits_table$index=c(1:nrow(hits_table))
    hits_table=hits_table[,c(1:12,14:15,13)]
    
    ## This is selected manually: these indexes are actually not good for this case
    
    results_output=merge(metadata_matrix,hits_table,by=colnames(metadata_matrix),all=TRUE)
    results_output=results_output[,c(1:12,14:15,13)]
    
    colnames(results_output)=c("GO_ID","description","levels","size_in_annotation","ontology_group","size_in_universe","uniqueness_index","number_of_genes_in_query","percentage_of_genes_in_query","p", "OR","fdr","community","node_numeric_id","genes_in_query")
    results_output=results_output[,c(1,2,3,5,14,13,7,4,6,8,9,11,10,12,15)]
    results_output=results_output[order(results_output$community,results_output$fdr),]
    
    factor_original <- factor(results_output$uniqueness_index,levels=unique(results_output$uniqueness_index))
    niveles <- levels(factor_original)
    niveles_numeros <- seq_along(niveles)
    factor_numeros <- factor(factor_original, levels = niveles, labels = niveles_numeros)
    
    results_output$uniqueness_index=factor_numeros
    results_output=results_output[order(results_output$uniqueness_index,results_output$fdr),]
    
    output <- list(background_matrix=background_matrix,result=results_output,pi0=pi0,network=network)
    return(output)
    
  }
  heatmaplotter <- function (namelist,data,feature_data,columns=NULL,name){
    if (is.null(columns)) {
      data_mod<-data[which(rownames(data) %in% namelist),]
    }else { data_mod<-data[which(rownames(data) %in% namelist),columns]
    }
    data_mod<-data_mod[match(namelist,rownames(data_mod)),]
    cols=colnames(data_mod)
    
    data_mod=t(apply(data_mod,1,scale))
    colnames(data_mod)=cols
    data_mod<-data_mod[match(namelist,rownames(data_mod)),]
    
    # meanvec<-apply(data_mod, 1, mean,na.rm=TRUE)
    #for(i in 1:nrow(data_mod)){
    #     data_mod[i,]<-data_mod[i,]-meanvec[i]
    #}
    data_mod<-as.matrix(data_mod)
    chunk=feature_data[which(feature_data$Locus.Tag %in% rownames(data_mod)),]
    chunk<-chunk[match(rownames(data_mod),chunk$Locus.Tag),]
    rownames(data_mod)=chunk$Gene_name
    
    data_expanded=melt(data_mod)
    colnames(data_expanded)=c("Gene","Sample","Z")
    
    
    
    
    finalplot<-ggplot(data_expanded, aes(x=Sample,y=Gene, fill= Z)) +
      geom_tile(color = "white",
                lwd = 0.5,
                linetype = 1) +
      scale_fill_gradient2(low="blue", high="red")+coord_fixed()+ggtitle(name)
    return(finalplot)
  }
  names_select=function(tab,attribute,values,features)
  {
    genes=tab$genes_in_query[which(tab[,attribute] %in% values)]
    genes <- paste(genes, collapse = ", ")
    genes <- unique(unlist(strsplit(genes, ", ")))
    print(length(genes))
    features=features[which(features$Locus.Tag %in% genes),]
    if(length(genes)!=nrow(features))
    {
      print("Se pierden alguno de los ",length(genes)," genes")
    }
    return(features)
  }
  
}

# 2. Load GO annotation objects
{
  #load("Inputs/GO_items/GO_in_R_for_loading")
  load("boolean_matrix_bp")
  load("BP_metadata_matrix")
  load("boolean_matrix_cc")
  load("CC_metadata_matrix")
  
  ## Merge ontologies
  
  dim(GO_gene_term_matrix_cc) #THis matrix contains information about terms related to cellular component and their genes
  dim(GO_gene_term_matrix_bp) #This matrix contains information about terms related to biological process and their genes
  dim(CC_metadata_matrix)     #This matrix contains metadata about the cellular component terms
  dim(BP_metadata_matrix)     #This matrix contains metadata about the biological process terms
  
  ## Los genes están en el mismo orden
  length(which(rownames(GO_gene_term_matrix_cc)==rownames(GO_gene_term_matrix_bp)))
  ## los términos en cada caso cuadran entre la matriz y la tabla de metadatos
  length(which(colnames(GO_gene_term_matrix_cc)==CC_metadata_matrix$Category))
  length(which(colnames(GO_gene_term_matrix_bp)==BP_metadata_matrix$Category))
  
  GO_gene_term_matrix_all=cbind(GO_gene_term_matrix_bp,GO_gene_term_matrix_cc)
  
  ## Add the Ontology_group column:
  BP_metadata_matrix$Ontology_group="Biological process"
  CC_metadata_matrix$Ontology_group="Cell compartment"
  ALL_metadata_matrix=rbind(BP_metadata_matrix,CC_metadata_matrix)   #This way we have merged both terms
}


#Example pipeline

{
  enrichment <- GO(query = query_genes,                           #Here we enter the query
                     background_matrix = GO_gene_term_matrix_all, #Here we load the matrix     
                     metadata_matrix = ALL_metadata_matrix,       #Here we load the metadata matrix
                     bg=background,                               #Here we load the background, THIS IS OPTIONAL, YOU CAN RUN IT WITHOUT BACKGROUND
                     min_GO_level = 4,                            #Here you set the minimum GO level you have to use
                     max_GO_level = 6,                            #Here you set the maximum GO level you have to use
                     test_mode = "enrichment",                    #Here you select the mode of the test
                     pvalue_correction_method = "qval",  #Here you set the pvalue correction method
                     min_term_size= 3,                   #Here you set the minimum size a term needs to have to be taken into account
                     max_term_size = Inf,                #Here you set the maximum allowed term size
                     th_hit=0.1,                         #Here you select the threshold for fdr
                     th_hit_OR=2,                        #Here you select the threshold for Odds Ratio
                     threshold_links = 0.5,              #And here you select  the threshold for connectivity, to create the networks
                     merge_duplicates=TRUE               #There are terms that have exactly the same genes. This option removes them and leaves only one
  )
  
  res_enrichment=enrichment$result
  enrichment=enrichment[which(!is.na(enrichment$community)),]

  #You can set the terms with labels this way: label_indexes_enrichment=c(1,5,22,29,34,39,44,49,54:59)
  
  #And this way you can write your gephi network
  
  
  
  
}