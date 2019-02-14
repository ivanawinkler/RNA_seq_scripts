#---------------------------------------------------------------------------------------------------------------------------------
#Function that reads in files from STAR in one dataframe and metadata-takes as arguments: mypath (path to folder containing 
#read count files); mypath_metadata (path to file containing metadata-dataframe-columns: sample names, conditions and type(SR or PE))
#----------------------------------------------------------------------------------------------------------------------------------
      
file_importeR=function(mypath,mypath_metadata){
  colData=read.csv(mypath_metadata,row.names = 1,)#loading countData
  filenames=list.files("./DESeq2_input", full.names=TRUE) #contains names of files in the folder
  mydataframefile=read.table(filenames[1],header=F,sep="\t",row.names = NULL) #reads in the first file
  mydataframefile=mydataframefile[5:dim(mydataframefile)[1],c(1,3)]
  colnames(mydataframefile)=c("V1",filenames[1])
  for (x in seq(from=2, to=length(filenames))){
    mydataframe=read.table(filenames[x],header=F,sep="\t",row.names = NULL) #reads in other files
    mydataframe=mydataframe[5:dim(mydataframe)[1], c(1,3)]
    colnames(mydataframe)=c("V1",filenames[x])
    mydataframefile=merge(mydataframefile,mydataframe,by.x="V1",by.y="V1") #merges the other files
  }
  colnames(mydataframefile)=sub(paste(mypath,"/",sep=""),"",colnames(mydataframefile))
  colnames(mydataframefile)=sub("_Reads.out.tab","",colnames(mydataframefile))
  rownames(mydataframefile)=mydataframefile[,1]
  countData=as.matrix(mydataframefile[2:17])
  countData=countData[, rownames(colData)]
  assign("countData",countData,envir = .GlobalEnv) #returns the dataframe to global environment
  assign("colData",colData,envir = .GlobalEnv) 
}

#---------------------------------------------------------------------------------------------------------------------------------
#DEG function - takes as arguments: counData (dataframe containing gene counts with sample_names as columns and genes as rownames); 
#               colData (metadata-dataframe-columns: sample names, conditions and type(SR or PE)); condition_list (list of paired 
#               conditions which are compared); cook (logical: T (default)-cook distance filtering turned on)
#---------------------------------------------------------------------------------------------------------------------------------
multiple_deseq2_analyseR=function(countData,colData,condition_list,cook=T){
  library("DESeq2")
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  
  #filtering out genes with 0 counts
  #dds=dds[rowSums(counts(dds))>1,]
  
  #differential expression analysis
  dds=DESeq(dds)
  res_list=list()
  for (i in seq(from=1, to=length(condition_list))){
    res=results(dds,cooksCutoff=T,contrast=c("condition",condition_list[[i]]))
    colnames(res)=sapply(colnames(res),function(x) paste(x,sapply(condition_list[i], "[", c(1)),sep="_"))
    res_list[[i]]=res
  }
  
  #merging all result tables in one
  res_temp=merge(as.data.frame(res_list[[1]]),as.data.frame(res_list[[2]]),by.x="row.names",by.y="row.names")
  colnames(res_temp)[1]="Gene_id"
  for (i in seq(from=3, to=length(condition_list))){
  res_temp=merge(res_temp,as.data.frame(res_list[[i]]),by.x="Gene_id",by.y="row.names")
  }
  assign("res_temp",res_temp,envir = .GlobalEnv) #returns the list to global environment #returns the list to global environment
}
#---------------------------------------------------------------------------------------------------------------------------------
#Function that returns dataframe with identifed outliers (pvalue=NA) only- it takes as arguments : res_temp (DEG dataframe with cook 
# filtering on) and res_temp_out (DEG dataframe with cook filtering off)
#----------------------------------------------------------------------------------------------------------------------------------
outlieR=function(res_temp,res_temp_out){
  res_temp_out=res_temp_out[is.na(res_temp$pvalue_nodul),]
  assign("res_outlier",res_temp_out,envir = .GlobalEnv) #returns the list to global environment
}
#-----------------------------------------------------------------------------------------------------------------------------------
#Function that returns dataframe with ensembl_ids, gene description,entrez_id and official gene symbol- it takes as an argument filter_id
#(identifier used for mapping;default ensembl_id) and res- DEG dataframe
#-------------------------------------------------------------------------------------------------------------------------------------
ensembleR=function(res){
  library(biomaRt)
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  filters=listFilters(ensembl)
  ensemnl_list=getBM(attributes=c("ensembl_gene_id","description","entrezgene","external_gene_name"),filters="ensembl_gene_id",values = gsub("\\..*","",as.vector(res[,1])),mart=ensembl)
  ensemnl_list=ensemnl_list[!duplicated(ensemnl_list$ensembl_gene_id),]
  res[,1]=gsub("\\..*","",res[,1])
  full_list=merge(res,ensemnl_list,by.x="Gene_id",by.y="ensembl_gene_id")
  assign("full_list",full_list,envir = .GlobalEnv) #returns the list to global environment #returns the list to global environment
}
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
  rld=rlog(dds,blind=F)

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
    text3d(x, text=rownames(model$x[,1:2]),adj=1.5)
    axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=TRUE)
    grid3d("x")
    grid3d("y")
    grid3d("z")
  } else { # 2d plot
    maxes= apply(abs(x), 2, max)
    rangeX= c(-maxes[1], maxes[1])
    rangeY= c(-maxes[2], maxes[2])
    plot(x, col=colData$condition, pch=19, xlab=c(colnames(x)[1],y[1]), ylab=c(colnames(x)[2],y[2]), xlim=rangeX, ylim=rangeY,main="PCA plot")
    lines(c(0,0), rangeX*2)
    lines(rangeY*2, c(0,0))
    text(x, labels=rownames(model$x[,1:2]), pos= 1 )
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



