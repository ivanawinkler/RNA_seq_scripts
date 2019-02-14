#####################################################################
#libraries
#####################################################################

library(jsonlite)
library("DESeq2")
library(biomaRt)
#-------------------------------------------------------------------
#organising metadata
#------------------------------------------------------------------

#defining the dataframe
metadata=read_json("metadata.cart.2017-04-13T11-14-43.896288.json",null="na",.na="na")

#substituting null to NA
for (i in seq(from=1,to=length(metadata))){
  if (is.null(metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_death)){metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_death=NA}
  if (is.null(metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up)){metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up=NA}
}

#defining columns of metadata table
metadata_table=as.data.frame(metadata[[1]]$associated_entities[[1]]$entity_submitter_id,na.rm=T)
metadata_table[1,2]=metadata[[1]]$cases[[1]]$samples[[1]]$sample_type_id
if (metadata_table[1,2]=="01"| metadata_table[1,2] == "02"){metadata_table[1,3]="tumor"} else {metadata_table[1,3]="control"}
metadata_table[1,4]=metadata[[1]]$cases[[1]]$diagnoses[[1]]$tumor_stage
metadata_table[1,5]=metadata[[1]]$cases[[1]]$diagnoses[[1]]$morphology
metadata_table[1,6]=metadata[[1]]$cases[[1]]$diagnoses[[1]]$days_to_death
metadata_table[1,7]=metadata[[1]]$cases[[1]]$samples[[1]]$initial_weight
metadata_table[1,8]=metadata[[1]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_normal_cells
metadata_table[1,9]=metadata[[1]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_stromal_cells
metadata_table[1,10]=metadata[[1]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_tumor_cells
metadata_table[1,11]=metadata[[1]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_necrosis
metadata_table[1,12]=metadata[[1]]$file_name
metadata_table[1,13]=metadata[[1]]$cases[[1]]$diagnoses[[1]]$vital_status
metadata_table[1,14]=metadata[[1]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up
metadata_table=data.frame(lapply(metadata_table, as.character), stringsAsFactors=FALSE)

#filling the columns
for (i in seq(from=2,to=length(metadata))){
  metadata_table[i,1]=metadata[[i]]$associated_entities[[1]]$entity_submitter_id
  metadata_table[i,2]=metadata[[i]]$cases[[1]]$samples[[1]]$sample_type_id
  if (metadata_table[i,2]=="01"| metadata_table[i,2] == "02"){metadata_table[i,3]="tumor"} else {metadata_table[i,3]="control"}
  metadata_table[i,4]=metadata[[i]]$cases[[1]]$diagnoses[[1]]$tumor_stage
  metadata_table[i,5]=metadata[[i]]$cases[[1]]$diagnoses[[1]]$morphology
  metadata_table[i,6]=metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_death
  metadata_table[i,7]=metadata[[i]]$cases[[1]]$samples[[1]]$initial_weight
  metadata_table[i,8]=metadata[[i]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_normal_cells
  metadata_table[i,9]=metadata[[i]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_stromal_cells
  metadata_table[i,10]=metadata[[i]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_tumor_cells
  metadata_table[i,11]=metadata[[i]]$cases[[1]]$samples[[1]]$portions[[1]]$slides[[1]]$percent_necrosis
  metadata_table[i,12]=metadata[[i]]$file_name
  metadata_table[i,13]=metadata[[i]]$cases[[1]]$diagnoses[[1]]$vital_status
  metadata_table[i,14]=metadata[[i]]$cases[[1]]$diagnoses[[1]]$days_to_last_follow_up
  
}
colnames(metadata_table)=c("barcode", "sample_type","tissue_status","tumor_stage","morphology","days_to death","initial_weight","percent_normal_cells","percent_stromal_cells",
                           "percent_tumor_cells","percent_necrosis","file_name","vital_status","days_to last_follow_up")
metadata_table[,12]=gsub(".gz","",metadata_table[,12])

#generating metadata table for Deseq2
colData=metadata_table[,c(1:3)]
colData=colData[order(colData$barcode),]
rownames(colData)=colData[,1]
colData=(colData[,c(2,3)])
colnames(colData)=c("tissue_type","condition")
colData$condition=factor(colData$condition,levels=c("control","tumor"))

#---------------------------------------------------------------------------------------------------------------------------------
#Function that reads in files from STAR in one dataframe and metadata-takes as arguments: mypath (path to folder containing 
#read count files); mypath_metadata (path to file containing metadata-dataframe-columns: sample names, conditions and type(SR or PE))
#----------------------------------------------------------------------------------------------------------------------------------


file_importeR=function(mypath){
  filenames=list.files(mypath, full.names=TRUE) #contains names of files in the folder
  mydataframefile=read.table(filenames[1],header=F,sep="\t",row.names = NULL) #reads in the first file
  mydataframefile=mydataframefile[c(1:(dim(mydataframefile)[1]-5)),]
  colnames(mydataframefile)=c("V1",filenames[1])
  for (x in seq(from=2, to=length(filenames))){
    mydataframe=read.table(filenames[x],header=F,sep="\t",row.names = NULL) #reads in other files
    mydataframe=mydataframe[1:(dim(mydataframe)[1]-5),]
    colnames(mydataframe)=c("V1",filenames[x])
    mydataframefile=merge(mydataframefile,mydataframe,by.x="V1",by.y="V1") #merges the other files
  }
  rownames(mydataframefile)=mydataframefile[,1]
  mydataframefile=mydataframefile[,2:length(mydataframefile)]
  assign("countData",mydataframefile,envir = .GlobalEnv) #returns the dataframe to global environment
}
  
file_importeR("./lihc_rnaseq")

#remove name of the folder- DON'T forget (you won't be able to match) 
colnames(countData)=gsub("./lihc_rnaseq/","",colnames(countData))

#changing names in countdata table to barcodes
names_list=setNames(as.list(metadata_table[,12]),metadata_table[,1])
for (i in seq(from=1, to=length(countData))){
  for (j in seq(from=1, to=length(names_list))){
  if (colnames(countData)[i]==names_list[[j]]){colnames(countData)[i]=names(names_list[j])}
}
}
countData=countData[,sort(names(countData))] 

#---------------------------------------------------------------------------------------------------------------------------------
#DEG function - takes as arguments: counData (dataframe containing gene counts with sample_names as columns and genes as rownames); 
#               colData (metadata-dataframe-columns: sample names, conditions and type(SR or PE)); condition_list (list of paired 
#               conditions which are compared); cook (logical: T (default)-cook distance filtering turned on)
#---------------------------------------------------------------------------------------------------------------------------------
multiple_deseq2_analyseR=function(countData,colData,cook=T){
  
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  dds$conditions=factor(dds$condition,levels = c("control","tumor"))
  dds=dds[rowSums(counts(dds))>1,]
  #differential expression analysis
  dds=DESeq(dds)
  res=results(dds,cooksCutoff = cook)
  return(res) #returns the list to global environment #returns the list to global environment
}

results_tcga=multiple_deseq2_analyseR(countData = countData, colData = colData)
results_tcga_df=as.data.frame(results_tcga)
#-----------------------------------------------------------------------------------------------------------------------------------
#Function that returns dataframe with ensembl_ids, gene description,entrez_id and official gene symbol- it takes as an argument filter_id
#(identifier used for mapping;default ensembl_id) and res- DEG dataframe
#-------------------------------------------------------------------------------------------------------------------------------------

ensembleR=function(res){
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
  filters=listFilters(ensembl)
  ensemnl_list=getBM(attributes=c("ensembl_gene_id","description","entrezgene","external_gene_name"),filters="ensembl_gene_id",values = gsub("\\..*","",as.vector(row.names(res))),mart=ensembl)
  ensemnl_list=ensemnl_list[!duplicated(ensemnl_list$ensembl_gene_id),]
  row.names(res)=gsub("\\..*","",row.names(res))
  full_list=merge(res,ensemnl_list,by.x="row.names",by.y="ensembl_gene_id")
  colnames(full_list)[1]="EnsemblID"
  assign("full_list",full_list,envir = .GlobalEnv) #returns the list to global environment #returns the list to global environment
}

ensembleR(results_tcga_df)
#--------------------------------------------------------------------------------------------------------------------------------------
#Plotting function-returns in form of pdf file cook distance plot, dispersio plot, heatmap of top 50 disregulated genes, heatmap of 
#sample to sample distances and PCA plot. Arguments: colData (df with metadata) and countData (count matrix)
#--------------------------------------------------------------------------------------------------------------------------------------

plotteR=function(colData,countData){
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  
  #filtering out genes with 0 counts
  dds=dds[rowSums(counts(dds))>1,]
  dds=DESeq(dds)
  
  #generatig normalised counts
  rld=varianceStabilizingTransformation(dds,blind=F)
  rld_df=assay(rld)
  #graph plotting
  pdf("RNA_seq_plots",onefile = T)
  boxplot(log10(assays(dds)[["cooks"]]),range=0,las=2,main="Cooks distance")
  plotDispEsts(dds,main="Dispersion plot")
  library("pheatmap")
  select=order(rowMeans(counts(dds,normalized=T)),decreasing = T)[1:20]
  pheatmap(assay(rld)[select,], cluster_rows=F, show_rownames=T,
           cluster_cols=T,main="Heatmap of top 50 disregulated genes")
  sampleDists <- dist(t(assay(rld)))
  #Heatmap of sample_to_sample_distance
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(colnames(assay(rld)))
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,main="Heatmap of sample_to_sample_distance")
  #PCA analysis
  #PCA plotting function
  library(rgl)
  plotPCA_f=function(x,y, nGroup,v) {
    n=ncol(x) 
    if(!(n %in% c(2,3))) { # check if 2d or 3d
      stop("x must have either 2 or 3 columns")
    }
    
    fit=hclust(dist(x), method="complete") # cluster
    groups= cutree(fit, k=nGroup)
    if(n == 3) { # 3d plot
      plot3d(x, col=v, type="s", size=1, axes=F)
      #text3d(x, text=rownames(model$x[,1:2]),adj=1.5)
      axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=TRUE)
      grid3d("x")
      grid3d("y")
      grid3d("z")
    } else { # 2d plot
      maxes= apply(abs(x), 2, max)
      rangeX= c(-maxes[1], maxes[1])
      rangeY= c(-maxes[2], maxes[2])
      plot(x, col=colData$condition, pch=1, xlab=c(colnames(x)[1],y[1]), ylab=c(colnames(x)[2],y[2]), xlim=rangeX, ylim=rangeY,main="PCA plot")
      lines(c(0,0), rangeX*2)
      lines(rangeY*2, c(0,0))
      #text(x, labels=rownames(model$x[,1:2]), pos= 1 )
    }
  }
  model= prcomp(t(assay(rld)), scale=TRUE)
  names(vector_3d)=rownames(colData)
  vector_3d=as.numeric(colData$condition)
  PoV <- (model$sdev^2/sum(model$sdev^2))*100
  PoV=as.character(PoV)
  PoV=paste(PoV,"%")
  plotPCA_f(model$x[,1:2], PoV, 2,vector_3d)
  plotPCA_f(model$x[,1:3],PoV, 4,vector_3d)
  
  dev.off()
}

#----------------------------------------------------------------------------------------------------------------------------------------
#Function that exports the result df as .csv file. It takes as arguments:result df, cutoff(numeric cutoff for padj filtering-all genes that
#have pvalue adjusted below the threshold are exported to the .csv file), file_name-(string, name of the exported file)
#--------------------------------------------------------------------------------------------------------------------------------------

exporteR=function(res_df,cutoff,file_name){
  padj_vector=colnames(res_df)[grepl("padj",names(res_df))]
  res_tmp=res_df[res_df[,padj_vector[1]]<cutoff,]
  res_output_f = res_tmp[!is.na(res_tmp[,padj_vector[1]]),]
  if (length(padj_vector)>1){
    for (i in seq(from=2,to=length(padj_vector))){
      res_tmp=res_df[res_df[,padj_vector[i]]<cutoff,]
      res_output = res_tmp[!is.na(res_tmp[,padj_vector[i]]),]
      res_output_f=rbind(res_output_f,res_output)
    }
  }
  res_output_f=res_output_f[!duplicated(res_output_f$Gene_id),]
  write.csv(res_output_f,file_name)
}

write.table(full_list,"TCGA_full_list",quote = F,row.names = F,sep="+")
write.table(rld_df,"rld_tcga",quote = F,sep=",")
write.table(metadata_table,"metadata_tcga",quote = F,sep=",")
