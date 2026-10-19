library(dplyr)
library(Seurat)
library(ggplot2)
library(xlsx)
library(gprofiler2) # Replaced biomaRt
library(Matrix)
library(progress)
future::plan("multicore", workers = 4)

##########################################USER INPUT##################################################################
file.list<- c("Neurons_zebrafish_ES_complete_081926.rds")
species.name.to.convert<-"hsapiens" 
##########################################################################################################################
###Common species inputs are the following
#Human-hsapiens
#Mouse -mmusculus
#Zebrafish - drerio
#Rat-rnorvegicus
#Drosophila -dmelanogaster
#C.elegans- celegans
#Chicken -ggallus
#Macaque -mmulatta
#Chimpanzee -ptroglodytes
#Pig -sscrofa

##########################################################################################################
#######################The converting function ##############################################################
count_matrix_converter <- function(data.input, species){
  genes <- row.names(data.input)
  pb <- progress_bar$new(
    format = " Converting matrix [:bar] :percent eta: :eta",
    total = length(genes), clear = FALSE, width= 60)
  
  # UPDATED: Using gprofiler2 instead of biomaRt
  ortho_raw <- gorth(
    query = genes,
    source_organism = "drerio",
    target_organism = species,
    numeric_ns = ""
  )
  
 
  ortho <- data.frame(
    ensembl_id_input = ortho_raw$input_ensg,                        
    zebrafish_gene = ortho_raw$input,        
    human_gene = ortho_raw$ortholog_name,
    ensemble_id_output=ortho_raw$ortholog_ensg
  )
  
  write.csv(ortho, "orthologue_conversion_file.csv")
  
  c=2
  mat1 <- Matrix(0, nrow = 1, 
                 ncol = length(colnames(data.input)), 
                 sparse = TRUE)
  c2=2
  mat2 <- Matrix(0, nrow = 1, 
                 ncol = length(colnames(data.input)), 
                 sparse = TRUE)
  
  for (i in 1:length(genes)){
    pb$tick()
    check1<- which(genes[i]== ortho[,2])
    if (length(check1)==1){
      check2<- nchar(ortho[,3][check1])
      check3<- is.na(ortho[,3][check1])
      if(check2>0 && check3==FALSE){
        mat1<-rbind(mat1, data.input[i,])
        row.names(mat1)[c]<-ortho[,3][check1]
        c=c+1
      }
    }else if (length(check1)>1){
      for (j in check1){
        check4<-nchar(ortho[,3][j])
        check5<-is.na(ortho[,3][j])
        if(check4>0 && check5==FALSE){
          mat2<-rbind(mat2,data.input[i,])
          row.names(mat2)[c2]<-ortho[,3][j]
          c2=c2+1
        }
      }
    }
  }
  
  mat1<-mat1[-1,]
  mat2<-mat2[-1,]
  comb.mat<- rbind(mat1, mat2)
  colnames(comb.mat)<- colnames(data.input)
  
  # Replaced duplicate purging loops with Matrix aggregation
  gene_names <- rownames(comb.mat)
  unique_genes <- unique(gene_names)
  
  pb_agg <- progress_bar$new(
    format = " Aggregating duplicate genes [:bar] :percent eta: :eta",
    total = 1, clear = FALSE, width= 60)
  
  agg_mat <- Matrix::fac2sparse(factor(gene_names, levels = unique_genes))
  comb.mat_clean <- agg_mat %*% comb.mat
  
  rownames(comb.mat_clean) <- unique_genes
  colnames(comb.mat_clean) <- colnames(data.input)
  
  pb_agg$tick()
  
  return(comb.mat_clean)
}




############################################################################################################
######################The loop for all the rds files########################################################

for ( filen in file.list ){
  status<- paste0("Converting ", filen, " to ",species.name.to.convert," orthologues....")
  print(status)
  

  filename<- filen
  zebra<- readRDS(file = filename)
  
  zebra.data<- zebra@assays$RNA@layers$counts 
  genenames<-row.names(zebra@assays$RNA@features@.Data)
  row.names(zebra.data)<-genenames
  colnames(zebra.data)<-row.names(zebra@assays$RNA@cells@.Data)
  
  # The biomaRt connection step (useMart) was completely removed from here
  
  species<- species.name.to.convert
  convert.data<-  count_matrix_converter(zebra.data,species)
  new_assay <- CreateAssay5Object(counts = convert.data)
  assay_name<-paste0("RNA_",species)
  zebra[[assay_name]]<- new_assay
  new.filename<-paste0("Converted_",species,"_",filename)
  saveRDS(zebra, file = new.filename)
}