library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(biomaRt)
library(ggplot2)


read_host_genecounts <- function(reports_dir){
  print(list.files(reports_dir) )
  list.files(reports_dir) %>% 
    lapply(function(x){
      sample_name <- stringr::str_match(x, "^(.*)_\\d+_reads_per_gene.star.tab$")[,2]
      fname <- paste(reports_dir, "/", x, sep="")
      if (file.size(fname) == 0) return(NULL)
      cbind(sample_name=sample_name,
            read.csv(fname, stringsAsFactors=F, sep='\t', header=FALSE)[,c(1,2)],
            stringsAsFactors=F)
    }) %>%
    do.call(what=dplyr::bind_rows) %>%
    pivot_wider(names_from = sample_name, values_from = V2) %>% 
    dplyr::slice(-1:-4) %>%
    mutate(V1 = sapply( strsplit( V1, split="\\." ), "[", 1 )) %>%
    set_rownames(.$V1)
}


filter_protein_coding <- function(reports, hgnc=TRUE){
  genenames<- rownames(reports) #sapply( strsplit( rownames(reports), split="\\." ), "[", 1 )
  
  if(file.exists("/Users/katrina.kalantar/code/earli-study/reference/genemap.csv")){
    genemap <- read.csv("/Users/katrina.kalantar/code/earli-study/reference/genemap.csv")
  } else{
    genenames %>% 
      getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description"),
             filters = "ensembl_gene_id",
             values = .,
             mart) -> 
      genemap
    write.csv(genemap, "/Users/katrina.kalantar/code/earli-study/reference/genemap.csv")
  }
  
  
  #create an index in which you are rearranging ensmbl gene ids from genemap into same order as they are in genenames
  idx <- match(genenames, genemap$ensembl_gene_id )
  
  cbind(as.character(genemap$hgnc_symbol[ idx ]), 
        as.character(genemap$description[ idx ]), 
        as.character(genemap$gene_biotype[ idx ]), 
        as.character(genemap$ensembl_gene_id[ idx ])) %>% 
    as.data.frame() %>% 
    `colnames<-`(c("hgnc_symbol", "description", "gene_biotype", "ensembl_gene_id") ) %>%
    filter(gene_biotype == "protein_coding") %>%
    filter(!is.na(ensembl_gene_id)) -> pc
  
  reports %>% 
    filter(V1 %in% pc$ensembl_gene_id) -> reports_pc
  
  if(hgnc){
    print('hello')
    reports_pc %>% mutate(V1 = pc$hgnc_symbol[match(.$V1, pc$ensembl_gene_id)]) %>%
      as.data.frame() -> reports_pc
    rownames(reports_pc) <- make.unique(as.character(reports_pc$V1))
    reports_pc  %>% 
      subset(., select = -c(V1)) -> reports_pc
  }
  
  return(reports_pc)
}

apply_genecount_qc_filters <- function(genecounts, 
                         threshold_total_sum = 1000000, 
                         threshold_nonzero_count = 10000){
  genecounts %>% 
    dplyr::select(which(colSums(.) > threshold_total_sum)) %>%
    dplyr::select(which(colSums(. > 0) > threshold_nonzero_count))
}

plot_numeric_bar <- function(numeric_array){
  mm = melt(numeric_array)
  mm$names = rownames(mm)
  ggplot(aes(x=names, y = value), data = mm) +
    geom_bar(stat = 'identity') +
    coord_flip()
}



assign_classifier_groups<-function(metadata, 
                                   classifier_variable, 
                                   true_positive_label, 
                                   true_negative_label, 
                                   train_fract = .5, test_fract = .2, valid_fract = .3){
  
  variable = levels(metadata[,c(classifier_variable)])
  only_validation = variable[!(variable %in% c(true_positive_label, true_negative_label))]
  metadata$classifier <- 0
  
  lapply(variable, function(x){
    total <- table(metadata[,c(classifier_variable)])[x]
    if(!is.na(total)){
      labs <- c()
      if(x %in% c(true_positive_label, true_negative_label)){
        labs <- sample(c(rep("train",round(total*train_fract)),
                         rep("test",round(total*test_fract)),
                         rep("valid",total-(round(total*train_fract) + round(total*test_fract)))))        
      }else{
        labs <- sample(c(rep("valid",total)))     
      }

      metadata %>%
        filter((!!sym(classifier_variable)) == x) %>%
        mutate(classifier = labs) -> labeled_subset
      
      return(labeled_subset)
    }
  }) %>%
    bind_rows(., .id = "column_label") ->
    metadata
}

get_eset_object <- function(metadata, genecounts, metadata_match_var, metadata_var_of_interest, classes_of_interest){
  
  for(i in seq(1:length(metadata_var_of_interest))){
    mvoi <- metadata_var_of_interest[i]
    coi <- classes_of_interest[[i]]
    metadata %>%
      filter(!!sym(mvoi) %in% coi) ->
      metadata
  }
  rownames(metadata) <- metadata[,c(metadata_match_var)]
  
  genecounts %>% 
    dplyr::select(which(colnames(.) %in% metadata[,c(metadata_match_var)])) -> filtered_genecounts
  
  metadata <- metadata[colnames(filtered_genecounts),]
  
  eset <- ExpressionSet(assayData = as.matrix(filtered_genecounts),
                        phenoData = AnnotatedDataFrame(metadata))
}
