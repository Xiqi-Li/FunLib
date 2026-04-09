# =============================================================================
# GENE EXPRESSION ANALYSIS FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   1. check_low_expression          - Flag low-expression samples
#   2. prepare_unsupervised_data     - Feature selection for clustering (MAD/CV/DQ/GUMBEL)
#   3. dge_limma                     - Differential expression via limma
#   4. dge_edgeR                     - Differential expression via edgeR
#   5. dge_DESeq                     - Differential expression via DESeq2
#   6. run_dge_gsea_pipeline         - Complete DGE+GSEA pipeline with plots
#   7. plot_volcano_with_annotations - Annotated volcano plot
#   8. unsupervised_analysis         - Multi-method unsupervised (UMAP/tSNE/PCA/MDS)
#   9. consensus_immunedeconvolute   - Multi-method immune deconvolution
#  10. geneset_activity              - ssGSEA/GSVA/zscore geneset scoring
#
# SOURCE PROVENANCE:
#   mRNA.R (primary):    woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02)
#   funcsInPembro.R:     Pembro/funcsInPembro.R (2025-07-02)
# =============================================================================


# SOURCE: funcsInPembro.R (2025-07-02)
# NOTE: This version returns a named list (sample_pct_low_expression + outlier_samples).
#       The mRNA.R version (2026-03-02) returns only sample_pct_low_expression vector
#       and adds na.rm=T to the sum() call. The Pembro version is used here as it is
#       more informative.
check_low_expression<-function(expressions,low_expression_cutoff=1,low_expression_pct_outlier_factor=1.5){
  tatal_low_expression<-sum(expressions<low_expression_cutoff)
  average_pct_low_expression<-tatal_low_expression/(nrow(expressions)*ncol(expressions))*100
  cat("On average, every sample have ",average_pct_low_expression,"% low expression genes.\n\n")
  sample_pct_low_expression<-colSums(expressions<low_expression_cutoff)/nrow(expressions)*100
  outlier_sample_pct_low_expression<-sample_pct_low_expression[sample_pct_low_expression>average_pct_low_expression*low_expression_pct_outlier_factor]
  if(length(outlier_sample_pct_low_expression)==0){
    cat("distribution of low expression gene for every sample is similar.")
  } else{
    cat("sample",paste(names(outlier_sample_pct_low_expression),collapse = ","),"looks having too many low expression genes: ",paste(outlier_sample_pct_low_expression,collapse=","),"% respectively.\n\n")
  }
  return(
    list(sample_pct_low_expression=sample_pct_low_expression,
         outlier_samples=outlier_sample_pct_low_expression))
}


# SOURCE: mRNA.R (2026-03-02)
#' prepare_unsupervised_data
#' extract qualified subset data with most variance for unsupervised analysis

#' @param expressions, dataframe, no default, gene expressions table,usually should optimally be clean and filtered by prepare_clean_RNA_sample_info_and_protein_expressions,must be provided
#' @param method, character, default="MAD", methods for calculation of variance, and etc, and later used by filtering, CV for coefficient of variance, DQ for quantile of both mean and variance (dual quantile), GUMBLE for gumbel distribution of mad, ATC for accumulated correlation with other genes

#' @param mad_top_n, numeric, default=-1, number of genes to be kept with top mad values, if it is -1 keep all genes
#' @param mad_top_quantile, numeric, default=0.75, percentage of genes to be kept with top mad values if mad_top_n is -1

#' @param cv_top_n, numeric, default=-1, number of genes to be kept with top cv values, if it is -1 keep all genes
#' @param cv_top_mean_quantile, numeric, default=0.5, percentage of genes to be kept with top mad values if cv_top_n is -1

#' @param dq_top_mean_quantile, numeric, default=0.5, genes to be kept with mean value quantile above defined dq_top_mean_quantile
#' @param dq_top_var_quantile, numeric, default=0.5,  genes to be kept with variance value quantile above defined dq_top_var_quantile

#' @param gumbel_p_cutoff, numeric, default=0.1, genes with p value of MAD gumbel distribution less than gumbel_p_cutoff will be kept

#' @param atc_top_n, numeric, default=-1, number of genes to be kept with top atc values, if it is -1 keep all genes
#' @param atc_cutoff, numeric, default=0.2, genes to be kept with atc values more than atc_cutoff

#' @param remove_outlier, logic, default FALSE, whether remove outlier value before calculation of variance, and etc

#' @return data.frame. data for unsupervised analysis.
#' @export
#'

prepare_unsupervised_data<-function(expressions,method=c("MAD","CV","DQ","GUMBEL","ATC","GUMBELATC"),mad_top_n=-1,mad_top_quantile=0.75,cv_top_n=-1,cv_top_mean_quantile=0.5,dq_top_mean_quantile=0.5,dq_top_var_quantile=0.5,gumbel_p_cutoff=0.1,atc_top_n=-1,atc_cutoff=0.2,remove_outlier=F){
  method=match.arg(method)
  #MAD
  if(method=='MAD'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=mad(dat,na.rm=T);return(mad)})
    if(mad_top_n==-1){
      results=expressions[mads>stats::quantile(mads,mad_top_quantile),]
    } else{
      results=expressions[order(mads,decreasing = T)[1:mad_top_n],]
    }
  }
  #GUMBEL
  if(method=='GUMBEL'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=stats::mad(dat,na.rm=T);return(mad)})
    scale=sqrt(var(mads,na.rm=T)*6/(pi^2))
    mu=mean(mads,na.rm=T)+scale*digamma(1)
    threshold=ordinal::qgumbel(1-gumbel_p_cutoff,mu,scale)
    results=expressions[mads>=threshold,]
  }
  #CV
  if(method=="CV"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    #'      rowMeans(expressions,na.rm=T)
    top_mean_filter=stats::quantile(means,cv_top_mean_quantile,na.rm=T)
    expressions_=expressions[means>=top_mean_filter,]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    cvs=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);cv=goeveg::cv(dat,na.rm=T);return(cv)})
    if(cv_top_n==-1){
      results=expressions_
    } else {
      results=expressions_[order(cvs,decreasing = T)[1:min(cv_top_n,nrow(expressions_))],]
    }
  }
  #DQ
  if(method=="DQ"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    top_mean_filter=stats::quantile(means,dq_top_mean_quantile,na.rm=T)
    expressions_=expressions[(means>=top_mean_filter & (!is.na(means))),]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    vars=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);var=stats::var(dat,na.rm=T);return(var)})
    top_var_filter=stats::quantile(vars,dq_top_var_quantile)
    results<-expressions_[(!is.na(vars)) & vars>=top_var_filter,]
  }
  if(method=="ATC"){
    atcs<-cola::ATC(expressions)
    if(atc_top_n!=-1){
      results<-expressions[order(atcs,decreasing = T)<=atc_top_n,]
    } else {
      results<-expressions[atcs>=atc_cutoff,]
    }
  }
  if(method=="GUMBELATC"){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=stats::mad(dat,na.rm=T);return(mad)})
    scale=sqrt(var(mads,na.rm=T)*6/(pi^2))
    mu=mean(mads,na.rm=T)+scale*digamma(1)
    threshold=ordinal::qgumbel(1-gumbel_p_cutoff,mu,scale)
    gumbel_results=expressions[mads>=threshold,]

    atcs<-cola::ATC(expressions)
    if(atc_top_n!=-1){
      atc_results<-expressions[order(atcs,decreasing = T)<=atc_top_n,]
    } else {
      atc_results<-expressions[atcs>=atc_cutoff,]
    }
    results<-expressions[intersect(rownames(gumbel_results),rownames(atc_results)),]
  }
  results<-results[rowSums(is.na(as.matrix(results)))<(0.25*ncol(results)),]
  return(results)
}


# SOURCE: funcsInPembro.R (2025-07-02)
# NOTE: This version adds robust_lm and weight parameters, and uses arrayWeights().
#       The mRNA.R version (2026-03-02) lacks these parameters; use this version for
#       compatibility with run_dge_gsea_pipeline().
#'
#' @description differential gene expression analysis using limma algorithm.
#'
#' @param expressions \code{data.frame()}. Gene expressions table.
#' @param is_rawcount \code{logical()}. Whether expressions is raw count or not. Default set to FALSE.
#' @param is_logged \code{logical()}. Whether provided expressions was log2 transformed or not. Default set to TRUE.
#' @param normalize \code{logical()}. Normalize provided expressions by library size. Default set to FALSE.
#' @param sample_frequency_threshold \code{numeric()}. Genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed. Default 0.5.
#' @param clinic_info \code{data.frame()}. Clinical information table.
#' @param ID_col \code{character()}. Column name for Sample ID (or Patient ID), should be consistent with column names of expressions.
#' @param group_col \code{character()}. Column name for group factor that differential gene expression will be conducted on.
#' @param covariate_col \code{character()}. column name for covariate factor that effects should be removed from model.
#' @param block_col \code{character()}. Column name for block factor if test is conducted on block model, such as paired test.
#' @param contrasts \code{vector(mode="character")} Specific contrasts if preferred, elements should be exactly same as group factor.
#' @param method \code{character(1)}. Method to conduct test. One of:
#' \itemize{
#'  \item "limma_trend" for non raw count expressions
#'  \item "limma_voom" for raw count expressions.
#' }
#' Default set as "limma_trend".
#' @param robust_lm \code{logical(1)}. Whether to use robust linear model. Default set to FALSE.
#' @param weight \code{logical(1)}. Whether to use weights in linear model. Default set to FALSE.
#'
#' @import DESeq2 statmod
#' @return \code{list()}, contains expressions, method, design, contrasts, test, and statistics of limma test.
#' @export
#'
#' @examples
#' \dontrun{
#' data(woodman)
#' clean_log2_protein_expressions=log2_expressions
#' results<-dge_limma(
#'   clean_log2_protein_expressions,
#'   clinical_info = clean_RNA_sample_info,
#'   ID_col = "Sample_ID",
#'   group_col = "APOLLO_TIMPOINT",
#'   contrasts = c("TP2-Baseline"),
#'   method ="limma_trend")
#' }
#'
#'
dge_limma<-function(expressions,is_rawcount=FALSE,is_logged=T,normalize=FALSE,sample_frequency_threshold=0.5,
                    clinic_info,
                    ID_col,
                    group_col,
                    covariate_col,
                    block_col,
                    contrasts,
                    method=c("limma_trend","limma_voom"),
                    robust_lm = FALSE,
                    weight = FALSE){
  require(limma)
  require(edgeR)
  require(DESeq2)
  stopifnot("ID column was not found in clinic_info"=ID_col %in% colnames(clinic_info))
  stopifnot("group column was not found in clinic_info"=group_col %in% colnames(clinic_info))
  clinic_info[[group_col]]<-factor(clinic_info[[group_col]])
  if(!missing(covariate_col)) {
    stopifnot("covariate column was not found in clinic_info" = covariate_col %in% colnames(clinic_info))
    clinic_info[[covariate_col]]<-factor(clinic_info[[covariate_col]])
    if(any(levels(clinic_info[[covariate_col]]) %in% levels(clinic_info[[group_col]]))){clinic_info[[covariate_col]]<-factor(paste(clinic_info[[covariate_col]],"_",sep=""))}
    formula<-paste(paste("~",covariate_col,sep="+"),group_col,sep="+")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c("Base",paste(levels(clinic_info[[covariate_col]])[-1],levels(clinic_info[[covariate_col]])[1],sep="-"),paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    if(!missing(contrasts)){
      contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
      colnames(contrast.matrix)<-contrasts
      rownames(contrast.matrix)<-colnames(design)
      for(contrast in contrasts){
        if(contrast %in% rownames(contrast.matrix)){
          contrast.matrix[contrast,contrast]<-1
        } else {
          contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
          contrast.matrix[contrast_,contrast]<-c(1,-1)
        }
      }

    }
  } else {
    formula<-paste("~0+",group_col,sep="")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-levels(clinic_info[[group_col]])
    contrast.matrix<-makeContrasts(contrasts=contrasts,levels = design)
  }
  if(!missing(block_col)) {
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    clinic_info[[block_col]]<-factor(clinic_info[[block_col]])
  }
  stopifnot("expression colnames do not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  method=match.arg(method)
  if(is_rawcount){
    dge <- DGEList(counts=expressions)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    if(method=="limma_trend"){
      logCPM <- cpm(dge, log=TRUE, prior.count=3)
      is_rawcount<-FALSE
      is_logged<-TRUE
      expressions<-logCPM
    } else if(method=="limma_voom"){
      v <- voom(dge, design, plot=TRUE, normalize="quantile")
      expressions<-v$E
      if(!missing(block_col)){
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        v <- voom(v, design, plot = TRUE, block = clinic_info[[block_col]], correlation = cor$consensus)
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        fit <- lmFit(v, design, block = clinic_info[[block_col]], correlation = cor$consensus)
      } else {
        fit <- lmFit(v, design)
      }
    }
  }
  if(!is_rawcount) {
    if(!is_logged){
      expressions<-log2(expressions+1)
    }
    if(normalize){
      library_size<-apply(expressions,2,sum)
      expressions<-t(t(expressions)/library_size)*mean(library_size)
    }
    expressions<-expressions[(rowSums(expressions==0))<(sample_frequency_threshold*ncol(expressions)),]
    cat("non raw_count data will be analyzed with limma\n")
    if(!missing(block_col)){
      corfit <- duplicateCorrelation(expressions,design,block=clinic_info[[block_col]])
      aw <- arrayWeights(expressions, design)
      fit <- lmFit(expressions,design,block=clinic_info[[block_col]],correlation=corfit$consensus,
                   method = if (robust_lm) "robust" else "ls", weights=if (weight) aw else NULL)
    } else {
      aw <- arrayWeights(expressions, design)
      fit <- lmFit(expressions,design,
                   method = if (robust_lm) "robust" else "ls", weights=if (weight) aw else NULL)
    }
  }
  fit<-contrasts.fit(fit,contrast.matrix)
  fit<-eBayes(fit,robust = T)
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    statistics=topTable(fit,coef=contrast,number=nrow(expressions))
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_statistics<-topTable(fit,number = nrow(expressions))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:length(contrasts))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method=method,design=design,contrast.matrix=contrast.matrix,fit=fit,statistics=contrast_statistics))
}


# SOURCE: mRNA.R (2026-03-02)
#' dge_edgeR
#' differential gene expression analysis for raw count expressions using edgeR algorithm

#' @param expressions, dataframe, no default, gene expressions table,must be provided

#' @param sample_frequency_threshold, numeric, default 0.5, genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed

#' @param clinic_info, dataframe, no default, clinic information table, must be provide
#' @param ID_col, character, no default, column name for Sample ID (or Patient ID), should be consistent with column names of expressions
#' @param group_col, character, no default, column name for group factor that differential gene expression will be conducted based on
#' @param covariate_col, character, no default, column name for covariate factor that effects should be removed from model
#' @param block_col, character, no default, column name for block factor if test is conducted on block model, such as paired test

#' @param contrasts, string vector, no default, specific contrasts if prefered, elements should be exactly same as group factor

#'
#' @return list. ontains expressions,methond, design,contrasts,test, and statistics of edgeR test
#' @export
#'

dge_edgeR<-function(expressions,sample_frequency_threshold=0.5,clinic_info,ID_col,group_col,covariate_col,block_col,contrasts){
  stopifnot("expression colnames wass not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  if(!missing(block_col)){
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    expressions_<-as.data.frame(t(expressions))
    expressions_[["Sample_ID"]]<-rownames(expressions_)
    expressions_[["block_col"]]<-clinic_info[[block_col]]
    sample_block_indices<-expressions_ %>% dplyr::group_by(block_col) %>% dplyr::group_indices()
    sample_id_df<-data.frame(Sample_ID=rownames(expressions_),block_col=clinic_info[[block_col]],sample_block_indices=sample_block_indices)
    unique_sample_id_df<-sample_id_df[(1:nrow(expressions_))[!duplicated(sample_id_df$sample_block_indices)],]
    expressions_<-expressions_ %>% dplyr::group_by(block_col) %>% dplyr::summarise_all(mean)
    expressions_<-as.data.frame(expressions_)
    rownames(expressions_)<-unique_sample_id_df$Sample_ID[match(expressions_$block_col,unique_sample_id_df$block_col)]
    expressions_<-floor(expressions_[,-1])
    expressions<-data.frame(t(expressions_))
    clinic_info<-clinic_info[match(unique_sample_id_df$Sample_ID,clinic_info[[ID_col]]),]
  }
  stopifnot("ID column was not found in clinic_info"=ID_col %in% colnames(clinic_info))
  stopifnot("group column was not found in clinic_info"=group_col %in% colnames(clinic_info))
  clinic_info[[group_col]]<-factor(clinic_info[[group_col]])
  if(!missing(covariate_col)) {
    stopifnot("covariate column was not found in clinic_info"=covariate_col %in% colnames(clinic_info))
    clinic_info[[covariate_col]]<-factor(clinic_info[[covariate_col]])
    if(any(levels(clinic_info[[covariate_col]]) %in% levels(clinic_info[[group_col]]))){clinic_info[[covariate_col]]<-factor(paste(clinic_info[[covariate_col]],"_",sep=""))}
    formula<-paste(paste("~",covariate_col,sep="+"),group_col,sep="+")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c("Base",paste(levels(clinic_info[[covariate_col]])[-1],levels(clinic_info[[covariate_col]])[1],sep="-"),paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    if(!missing(contrasts)){
      contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
      colnames(contrast.matrix)<-contrasts
      rownames(contrast.matrix)<-colnames(design)
      for(contrast in contrasts){
        if(contrast %in% rownames(contrast.matrix)){
          contrast.matrix[contrast,contrast]<-1
        } else {
          contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
          contrast.matrix[contrast_,contrast]<-c(1,-1)
        }
      }

    }
  } else {
    formula<-paste("~0+",group_col,sep="")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-levels(clinic_info[[group_col]])
    contrast.matrix<-limma::makeContrasts(contrasts=contrasts,levels = design)
  }
  dge <- edgeR::DGEList(counts=expressions)
  keep <- edgeR::filterByExpr(dge,design)
  keep<-keep & (rowSums(dge$counts<=1)<(sample_frequency_threshold*ncol(dge$counts)))
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge<-edgeR::estimateCommonDisp(dge,design,robust=T)
  fit<-edgeR::glmQLFit(dge,design)

  expressions<-edgeR::cpm(dge, log=TRUE, prior.count=3)
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    test<-edgeR::glmQLFTest(fit,contrast=contrast.matrix[,contrast])
    statistics=as.data.frame(edgeR::topTags(test,n=nrow(expression)))
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_test<-edgeR::glmQLFTest(fit,contrast=contrast.matrix)
    F_statistics<-as.data.frame(edgeR::topTags(F_test,n = nrow(expressions)))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:(length(contrasts)+1))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method="edgeR",design=design,contrast.matrix=contrast.matrix,fit=fit,statistics=contrast_statistics))
}


# SOURCE: mRNA.R (2026-03-02)
#' dge_DESeq
#' differential gene expression analysis for raw count expressions using DESeq2 algorithm

#' @param expressions, dataframe, no default, gene expressions table,must be provided

#' @param sample_frequency_threshold, numeric, default 0.5, genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed

#' @param clinic_info, dataframe, no default, clinic information table, must be provide
#' @param ID_col, character, no default, column name for Sample ID (or Patient ID), should be consistent with column names of expressions
#' @param group_col, character, no default, column name for group factor that differential gene expression will be conducted based on
#' @param covariate_col, character, no default, column name for covariate factor that effects should be removed from model
#' @param block_col, character, no default, column name for block factor if test is conducted on block model, such as paired test

#' @param contrasts, string vector, no default, specific contrasts if prefered, elements should be exactly same as group factor

#'
#' @return list. contains expressions,methond, design,contrasts,test, and statistics of DESeq2 test
#' @export
#'

dge_DESeq<-function(expressions,sample_frequency_threshold=0.5,clinic_info,ID_col,group_col,covariate_col,block_col,contrasts){
  stopifnot("expression colnames wass not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  if(!missing(block_col)){
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    expressions_<-as.data.frame(t(expressions))
    expressions_[["Sample_ID"]]<-rownames(expressions_)
    expressions_[["block_col"]]<-clinic_info[[block_col]]
    sample_block_indices<-expressions_ %>% dplyr::group_by(block_col) %>% dplyr::group_indices()
    sample_id_df<-data.frame(Sample_ID=rownames(expressions_),block_col=clinic_info[[block_col]],sample_block_indices=sample_block_indices)
    unique_sample_id_df<-sample_id_df[(1:nrow(expressions_))[!duplicated(sample_id_df$sample_block_indices)],]
    expressions_<-expressions_ %>% dplyr::group_by(block_col) %>% dplyr::summarise_all(mean)
    expressions_<-as.data.frame(expressions_)
    rownames(expressions_)<-unique_sample_id_df$Sample_ID[match(expressions_$block_col,unique_sample_id_df$block_col)]
    expressions_<-floor(expressions_[,-1])
    expressions<-data.frame(t(expressions_))
    clinic_info<-clinic_info[match(unique_sample_id_df$Sample_ID,clinic_info[[ID_col]]),]
  }
  stopifnot("ID column was not found in clinic_info"=ID_col %in% colnames(clinic_info))
  stopifnot("group column was not found in clinic_info"=group_col %in% colnames(clinic_info))
  clinic_info[[group_col]]<-factor(clinic_info[[group_col]])
  if(!missing(covariate_col)) {
    stopifnot("covariate column was not found in clinic_info"=covariate_col %in% colnames(clinic_info))
    clinic_info[[covariate_col]]<-factor(clinic_info[[covariate_col]])
    if(any(levels(clinic_info[[covariate_col]]) %in% levels(clinic_info[[group_col]]))){clinic_info[[covariate_col]]<-factor(paste(clinic_info[[covariate_col]],"_",sep=""))}
    formula<-paste(paste("~",covariate_col,sep="+"),group_col,sep="+")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c("Base",paste(levels(clinic_info[[covariate_col]])[-1],levels(clinic_info[[covariate_col]])[1],sep="-"),paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    if(!missing(contrasts)){
      contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
      colnames(contrast.matrix)<-contrasts
      rownames(contrast.matrix)<-colnames(design)
      for(contrast in contrasts){
        if(contrast %in% rownames(contrast.matrix)){
          contrast.matrix[contrast,contrast]<-1
        } else {
          contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
          contrast.matrix[contrast_,contrast]<-c(1,-1)
        }
      }

    }
  } else {
    formula<-paste("~",group_col,sep="")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c(levels(clinic_info[[group_col]])[1],paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
    colnames(contrast.matrix)<-contrasts
    rownames(contrast.matrix)<-colnames(design)
    for(contrast in contrasts){
      if(contrast %in% rownames(contrast.matrix)){
        contrast.matrix[contrast,contrast]<-1
      } else {
        contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
        contrast.matrix[contrast_,contrast]<-c(1,-1)
      }
    }
  }
  formula<-gsub("0\\+","",formula)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = expressions,
                                colData = clinic_info,
                                design = as.formula(formula))
  keep <- rowSums(BiocGenerics::counts(dds)) >= ncol(expressions)*2
  keep<-keep & (rowSums(BiocGenerics::counts(dds)<=1)<(sample_frequency_threshold*ncol(BiocGenerics::counts(dds))))
  dds <- dds[keep,]
  dds <- DESeq2::DESeq(dds)
  expressions<-tryCatch({SummarizedExperiment::assay(DESeq2::vst(dds,blind = F))},error=function(e){SummarizedExperiment::assay(DESeq2::vst(dds))})
  print(dim(expressions))
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    statistics=as.data.frame(DESeq2::results(dds,contrast=contrast.matrix[,contrast],cooksCutoff = F))[,-1]
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_statistics<-as.data.frame(DESeq2::results(dds,cooksCutoff = F))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:length(contrasts))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method="DESeq2",design=design,contrast.matrix=contrast.matrix,fit=dds,statistics=contrast_statistics))
}


# SOURCE: funcsInPembro.R (2025-07-02)
run_dge_gsea_pipeline <- function(
    sample_info,
    log2_expressions,
    phenotype,
    contrast,
    pathways,
    logFC_thresh = log2(1.5),
    show_boxplots = TRUE,
    show_volcano = TRUE,
    show_heatmap = TRUE,
    run_gsea = TRUE,
    make_datatables = TRUE,
    gsea_rank_metric = "t",
    gsea_rank_p_thresh = NULL,
    robust_lm=F,
    weight=F
) {
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(DT)
  library(fgsea)
  library(ComplexHeatmap)

  # Load reference data
  oncoKB <- read.csv(system.file("extdata/cancerGeneList_oncoKB_lastupdate05192023.tsv", package = "expr"), sep = "\t", header = TRUE, check.names = FALSE)
  rownames(oncoKB) <- oncoKB$`Hugo Symbol`
  oncoKB <- oncoKB %>%
    mutate(
      `Is Oncogene` = `Is Oncogene` == "Yes",
      `Is Tumor Suppressor Gene` = `Is Tumor Suppressor Gene` == "Yes",
      `OncoKB Annotated` = `OncoKB Annotated` == "Yes"
    )

  cancerGeneCensus <- read.csv(system.file("extdata/concensus_Census_allWed_Jan_18_20_33_56_2023.tsv", package = "expr"), sep = "\t", header = TRUE, check.names = FALSE)
  rownames(cancerGeneCensus) <- cancerGeneCensus$`Gene Symbol`
  cancerGeneCensus$Hallmark <- cancerGeneCensus$Hallmark != ""

  load(system.file("extdata/protein_coding_ensemble2symbol.RData", package = "expr"))

  # Setup group info
  group <- data.frame(as.factor(sample_info[, phenotype]))
  dimnames(group) <- list(rownames(sample_info), phenotype)
  split_groups <- setNames(group[[phenotype]], sample_info$Sample_ID)
  group <- data.frame(Sample_ID = sample_info$Sample_ID, group)
  group <- group[complete.cases(group),]

  # DGE
  comp_limma <- dge_limma(
    expressions = log2_expressions[, group$Sample_ID],
    is_rawcount = FALSE,
    is_logged = TRUE,
    normalize = FALSE,
    sample_frequency_threshold = 0.5,
    clinic_info = group,
    ID_col = "Sample_ID",
    group_col = phenotype,
    contrasts = contrast,
    method = "limma_trend",
    robust_lm = robust_lm,
  )

  stats <- comp_limma$statistics %>% arrange(get(sprintf("%s:P.Value", contrast)))
  stats <- cbind(stats, protein_coding_ensemble2symbol[match(rownames(stats), protein_coding_ensemble2symbol$Gene_Symbol), -1])

  tmp <- apply(oncoKB[stats$Gene_Symbol, c("OncoKB Annotated", "Is Oncogene", "Is Tumor Suppressor Gene")], 1, function(x) {
    x[1] <- !x[1]
    paste(c("not annotated", "oncogene", "TSG")[x], collapse = "; ")
  })
  tmp <- replace(tmp, tmp == "NA; NA; NA", NA)
  stats$oncoKB <- tmp

  tmp <- apply(cancerGeneCensus[stats$Gene_Symbol, c("Tier", "Hallmark")], 1, function(x) {
    paste(paste("Tier", x[1]), ifelse(x[2], "Hallmark", ""))
  })
  tmp <- replace(tmp, tmp == "TierNA NA", NA)
  stats$cancerGeneCensus <- trimws(tmp)

  sigRow <- stats[[sprintf("%s:P.Value", contrast)]] < 0.05 & abs(stats[[sprintf("%s:logFC", contrast)]]) > logFC_thresh

  dge_dt <- NULL
  if (make_datatables && sum(sigRow) > 0) {
    dge_dt <- stats[sigRow,] %>%
      DT::datatable(filter = "top", options = list(initComplete = JS("function(settings, json) { $(this.api().table().header()).css({'font-size': '12px'}); }"))) %>%
      DT::formatRound(columns = names(which(sapply(.$x$data, is.numeric))), digits = 2) %>%
      DT::formatStyle(columns = colnames(.$x$data), `font-size` = '12px')
  }


  # Optional Plots
  DGE_plots <- NULL
  if (show_boxplots && sum(sigRow) > 0) {
    DGE_plots <- ggplotly(
      plotBoxPlot(
        data_mx = log2_expressions,
        sampleAttr = group,
        MOI = rownames(stats[sigRow,])[1:min(12, sum(sigRow))],
        header4x = phenotype,
        header4ID = "Sample_ID",
        header4color = phenotype
      ) + theme(legend.position = "none") + labs(title = "Top differentially expressed genes")
    )
  }

  volcano_plot <- NULL
  if (show_volcano) {
    dff <- stats %>%
      mutate(diffexpressed = ifelse(
        get(sprintf("%s:P.Value", contrast)) < 0.05,
        ifelse(get(sprintf("%s:logFC", contrast)) > logFC_thresh, "Up",
               ifelse(get(sprintf("%s:logFC", contrast)) < -logFC_thresh, "Down", "No")),
        "No"
      )) %>%
      mutate(`-log10(P)` = -log10(get(sprintf("%s:P.Value", contrast))),
             label = Gene_Symbol)
    dff$label[!(dff$`-log10(P)` > -log10(0.05) & dff$diffexpressed != "No")] <- NA

    volcano_plot <- ggplot(dff, aes_string(x = sprintf("%s:logFC", contrast), y = "`-log10(P)`", col = "diffexpressed", label = "label")) +
      geom_point() +
      geom_text_repel() +
      theme_minimal() +
      scale_color_manual(name = "expression", values = c("dodgerblue3", "gray50", "firebrick3")) +
      geom_vline(xintercept = c(-logFC_thresh, logFC_thresh), linetype = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) +
      labs(title = sprintf("Differentially expressed genes (p<0.05, |log2FC|> %s)", logFC_thresh),
           x = sprintf("%s:logFC", contrast))
  }

  h.supervised <- NULL
  if (show_heatmap && sum(sigRow) > 1) {
    topAnno <- getHeatMapAnnotation(
      sampleAttr = sample_info,
      nGroupMax = 20,
      tracks = c("CohortCat", "cancerTypesI", "clinical_benefit", "gender", "organ", "TIL (0-3)", "ImmuneScore", "lowexpressgene_pct", "EIScat")
    )
    h.supervised <- Heatmap(t(scale(t(log2_expressions[rownames(stats[sigRow,]), ]))),
                            name = "expression",
                            cluster_columns = cluster_within_group(t(scale(t(log2_expressions[rownames(stats[sigRow,]), ]))), sample_info[[phenotype]]),
                            column_title = "DGE between baseline groups with/without clinical benefit",
                            show_row_names = FALSE,
                            show_column_names = FALSE,
                            top_annotation = topAnno,
                            column_split = 2)
  }


  # GSEA
  rank <- switch(
    gsea_rank_metric,
    "logFC" = {
      if (!is.null(gsea_rank_p_thresh)) {
        a <- stats %>%
          dplyr::filter(!is.na(.data[[sprintf("%s:logFC", contrast)]]) & .data[[sprintf("%s:P.Value", contrast)]]<gsea_rank_p_thresh)
        setNames(a[[sprintf("%s:logFC", contrast)]], rownames(a))
      } else {
        a=stats %>%
          dplyr::filter(!is.na(.data[[sprintf("%s:logFC", contrast)]]))
        setNames(a[[sprintf("%s:logFC", contrast)]], rownames(a))
      }
    },
    "SNR" = {
      # Compute Signal-to-Noise Ratio (mean difference divided by pooled SD)
      groups <- unlist(strsplit(contrast, "-"))
      group1 <- groups[1]
      group2 <- groups[2]
      group_labels <- factor(sample_info[[phenotype]], levels = c(group1, group2))
      names(group_labels)=sample_info$Sample_ID
      group1_samples <- names(group_labels[group_labels == group1])
      group2_samples <- names(group_labels[group_labels == group2])

      exprs_group1 <- log2_expressions[, group1_samples]
      exprs_group2 <- log2_expressions[, group2_samples]

      mean1 <- rowMeans(exprs_group1, na.rm = TRUE)
      mean2 <- rowMeans(exprs_group2, na.rm = TRUE)
      sd1 <- apply(exprs_group1, 1, sd, na.rm = TRUE)
      sd2 <- apply(exprs_group2, 1, sd, na.rm = TRUE)

      snr_rank <- (mean1 - mean2) / (sd1 + sd2)
      names(snr_rank) <- rownames(log2_expressions)
      snr_rank
      # snr[is.na(snr)] <- 0
      # sort(snr, decreasing = TRUE)
    },
    "signed_pval" = {
      logfc <- stats[[sprintf("%s:logFC", contrast)]]
      pval <- stats[[sprintf("%s:P.Value", contrast)]]
      signed_score <- -log10(pval) * sign(logfc)
      signed_score[is.na(signed_score)] <- 0
      setNames(signed_score, rownames(stats))
    },
    "t" = {
      stats[[sprintf("%s:t", contrast)]][is.na(stats[[sprintf("%s:t", contrast)]])] <- 0
      setNames(stats[[sprintf("%s:t", contrast)]], rownames(stats))
    },
    stop("Invalid gsea_rank_metric. Choose from 'logFC', 'SNR', 'signed_pval', or 't'.")
  )


  gseaRess <- gseaRess.sig <- gseaPlots <- list()
  if (run_gsea) {
    for (geneSet in names(pathways)) {
      pathway <- pathways[[geneSet]]
      # rank <- stats %>%
      #   arrange(get(sprintf("%s:logFC", contrast))) %>%
      #   filter(!!rlang::sym(sprintf("%s:P.Value", contrast)) < 0.05) %>%
      #   select(sprintf("%s:logFC", contrast)) %>%
      #   tibble::rownames_to_column("x") %>%
      #   tibble::deframe()

      gseaRes <- fgsea(
        pathways = pathway,
        stats = rank,
        minSize = min(sapply(pathway, length)),
        gseaParam = 1
      )
      gseaRess[[geneSet]] <- gseaRes

      if (make_datatables) {
        gseaRess.sig[[geneSet]] <- gseaRes %>%
          filter(pval < 0.05) %>%
          arrange(pval) %>%
          DT::datatable(filter = "top", options = list(initComplete = JS("function(settings, json) { $(this.api().table().header()).css({'font-size': '12px'}); }"))) %>%
          DT::formatRound(columns = names(which(sapply(.$x$data, is.numeric))), digits = 2) %>%
          DT::formatStyle(columns = colnames(.$x$data), `font-size` = '12px')
      }

      fgseaRes <- gseaRes %>% filter(pval < 0.05) %>% arrange(NES)
      gseaPlots[[geneSet]] <- plotGseaTable(pathway[fgseaRes$pathway], rank, fgseaRes, gseaParam = 1)
    }
  }

  return(list(
    DGE_stats = stats,
    DGE_datatable = dge_dt,
    Boxplot = DGE_plots,
    Volcano = volcano_plot,
    Heatmap = h.supervised,
    GSEA_results = gseaRess,
    GSEA_datatables = gseaRess.sig,
    GSEA_plots = gseaPlots
  ))
}


# SOURCE: funcsInPembro.R (2025-07-02)
#' Plot Annotated Volcano Plot with GSEA and Pathway Highlights
#'
#' This function creates an enhanced volcano plot from a differential gene expression table,
#' highlighting genes from selected pathways, optionally marking key genes, and annotating
#' GSEA results (NES and p-values). Genes that belong to multiple pathways will be represented
#' with multiple color-coded lineranges.
#'
#' @param dge_table A data frame containing differential expression results. Must include columns:
#'   \code{Gene_Symbol}, \code{Yes-No:logFC}, and \code{Yes-No:P.Value}.
#' @param gsea_df A data frame of GSEA results with columns \code{pathway}, \code{NES}, and \code{padj}.
#' @param pathways A named list of gene sets (e.g., from MSigDB), structured as \code{pathways$hallmark}.
#' @param pathway_labels A named character vector mapping display group names (e.g., "IFNG") to
#'   GSEA pathway names (e.g., "HALLMARK_INTERFERON_GAMMA_RESPONSE").
#' @param pathway_colors A named character vector assigning colors to each pathway label (e.g., \code{c("IFNG" = "#ea838e")}).
#' @param key_genes Optional character vector of gene symbols to highlight as key genes.
#' @param key_color Color for key genes if \code{"key"} is not already in \code{pathway_colors}. Default is \code{"#0afdfe"}.
#' @param rect_color_up Background rectangle color for upregulated significant genes. Default is \code{"#ea838e"}.
#' @param rect_color_down Background rectangle color for downregulated significant genes. Default is \code{"#827fd6"}.
#' @param step_divisor Numeric divisor used to control vertical spacing of GSEA annotation bars. Default is \code{6}.
#' @param linerange_linewidth Line width for \code{geom_linerange()} used to mark gene hits in each pathway. Default is \code{0.1}.
#' @param logFC_col Column name in \code{dge_table} for log fold change values. Default is \code{"Yes-No:logFC"}.
#' @param pval_col Column name in \code{dge_table} for p-values. Default is \code{"Yes-No:P.Value"}.
#' @param fontsize Font size for text annotations. Default is \code{3}.
#' @param stepscale Numeric scale factor for adjusting the vertical spacing of GSEA annotations. Default is \code{1}.
#'
#' @return A \code{ggplot} object showing the volcano plot with pathway and GSEA annotations.
#'
#' @examples
#' plot_volcano_with_annotations(
#'   dge_table = my_dge,
#'   gsea_df = my_gsea,
#'   pathways = pathways,
#'   pathway_labels = c("IFNG" = "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
#'   pathway_colors = c("IFNG" = "#ea838e"),
#'   key_genes = c("STAT1", "IRF1")
#' )
#'
#' @import ggplot2
#' @import ggrepel
#' @export


plot_volcano_with_annotations <- function(
    dge_table,
    logFC_col ="Yes-No:logFC",
    pval_col = "Yes-No:P.Value",
    logFCthresh = log2(1.5),
    gsea_df,
    pathways,
    pathway_labels,
    pathway_colors,
    key_genes = NULL,
    key_color = "#0afdfe",
    rect_color_up = "#ea838e",
    rect_color_down = "#827fd6",
    step_divisor = 6,
    linerange_linewidth = 0.1,
    fontsize = 3,
    stepscale=1
){
  diff <- dge_table
  colnames(diff)[colnames(diff) == logFC_col] <- "logFC"
  colnames(diff)[colnames(diff) == pval_col] <- "P"

  # Initialize gene_col
  diff$gene_col <- "all"

  # Assign pathway groups
  for (group in names(pathway_labels)) {
    path_name <- pathway_labels[[group]]
    genes_in_pathway <- pathways[[path_name]]
    diff$gene_col[diff$Gene_Symbol %in% genes_in_pathway] <- group
  }

  # Add key genes if provided
  if (!is.null(key_genes)) {
    diff$gene_col[diff$Gene_Symbol %in% key_genes] <- "key"
  }

  # Ensure key color is in color list
  if (!"key" %in% names(pathway_colors) && !is.null(key_genes)) {
    pathway_colors["key"] <- key_color
  }

  # Background and category-specific data
  data_bg <- diff[diff$gene_col == "all", ]
  data_geneset <- diff[diff$gene_col != "all", ]
  data_key <- if (!is.null(key_genes)) diff[diff$Gene_Symbol %in% key_genes, ] else NULL

  # Set up base plot with up/down highlight rects
  p <- ggplot(diff, aes(x = logFC, y = -log10(P))) +
    annotate("rect", fill = rect_color_down, alpha = 0.3,
             xmin = min(diff[diff$P < 0.05 & diff$logFC < 0, "logFC"]),
             xmax = if (is.null(logFCthresh)) max(diff[diff$P < 0.05 & diff$logFC < 0, "logFC"]) else -logFCthresh,
             ymin = -log10(0.05), ymax = max(-log10(diff$P))) +
    annotate("rect", fill = rect_color_up, alpha = 0.3,
             xmin = if (is.null(logFCthresh)) min(diff[diff$P < 0.05 & diff$logFC > 0, "logFC"]) else logFCthresh,
             xmax = max(diff[diff$P < 0.05 & diff$logFC > 0, "logFC"]),
             ymin = -log10(0.05), ymax = max(-log10(diff$P))) +
    geom_point(data = data_bg, shape = 21, color = "black", alpha = 0.2, stroke = 0.3) +
    geom_point(data = data_geneset, aes(fill = gene_col), shape = 21, alpha = 0.8, color = "black", stroke = 0.3) +
    scale_fill_manual(values = pathway_colors) +
    geom_vline(xintercept = 0, color = "#b2b2b2", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "#b2b2b2", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), color = "#ff7d82", linewidth = 0.6, linetype = "dotted") +
    theme_classic() +
    theme(legend.position = "none")

  # Add key gene labels
  if (!is.null(data_key) && nrow(data_key) > 0) {
    p <- p + geom_label_repel(
      data = data_key,
      aes(label = Gene_Symbol),
      size = 6, fill = "white", alpha = 0.8,
      box.padding = unit(0.35, "lines"),
      point.padding = 0.5,
      segment.colour = "#4c4b5e", segment.size = 0.5,
      min.segment.length = 0
    )
  }


  # Add lineranges and NES/padj for each pathway
  yScale <- range(-log10(diff$P))
  step <- diff(yScale) / step_divisor
  xScale <- range(diff$logFC)

  linerange_str <- ""
  subset_datas=  extra_layers=list()
  for (group in names(pathway_labels)) {
    path_name <- pathway_labels[[group]]
    genes <- pathways[[path_name]]
    subset_data <- diff[diff$Gene_Symbol %in% genes, ]
    subset_data$gene_col2=ifelse(subset_data$Gene_Symbol %in% genes, group, "all")
    subset_datas[[group]]=subset_data

    ymin_val <- -which(names(pathway_labels) == group) * step
    ymax_val <- ymin_val + 0.8 * step

    geom_line <- sprintf(
      "geom_linerange(data =subset_datas[['%s']], aes(x = logFC, ymin = %s, ymax = %s), color = '%s', size = 0.5, linewidth = %s)",
      group,
      ymin_val,
      ymax_val,
      pathway_colors[[group]],
      linerange_linewidth
    )

    # Append with "+" if not the first
    if (nchar(linerange_str) > 0) {
      linerange_str <- paste0(linerange_str, " +\n", geom_line)
    } else {
      linerange_str <- geom_line
    }


    # Label left
    extra_layers[[length(extra_layers) + 1]] <- annotate(
      "text",
      x = 1.1 * xScale[1],
      y = ymin_val + 0.4 * step,
      label = gsub("HALLMARK_", "", path_name),
      color = pathway_colors[[group]],
      hjust = 0,
      size = fontsize
    )

    # Label right (NES, pval)
    extra_layers[[length(extra_layers) + 1]] <- annotate(
      "text",
      x = 1.1 * xScale[2],
      y = ymin_val + 0.4 * step,
      label = sprintf(
        "NES: %.2f\nfdr: %s",
        unlist(gsea_df[gsea_df$pathway == path_name, "NES"]),
        formatC(unlist(gsea_df[gsea_df$pathway == path_name, "padj"]), format = "e", digits = 2)
      ),
      color = pathway_colors[[group]],
      lineheight = stepscale*step,
      hjust = 1,
      size = fontsize
    )

  }
  p_final <- Reduce(`+`, c(list(p), extra_layers))
  # Evaluate the complete geom_linerange block
  p <- eval(parse(text = paste0("p_final+\n",linerange_str)))
  return(p)
}


# SOURCE: mRNA.R (2026-03-02)
#' unsupervised_analysis
#'
#' unsupervised analysis of expression
#'
#' @param dat data.frame. data for unsupervised analysis.
#' @param labels character. labels for samples.
#' @param scaled character. scale direction.
#' @param run_umap logic. run umap method or not.
#' @param umap_params list. contain parameters for umap.
#' @param umap_config list. contain configuration for umap.
#' @param run_tsne logic. run umap method or not.
#' @param tsne_params list. contain parameters for tsne.
#' @param run_pca logic. run umap method or not.
#' @param pca_params list. contain parameters for pca
#' @param run_mds logic. run umap method or not.
#' @param mds_params list. contain parameters for mds.
#' @param run_heatmap logic. run umap method or not.
#' @param heatmap_params list. contain parameters for heatmap.
#' @param cols character. colors for labels.
#'
#' @return list. contain results for each method.
#' @export
#'

unsupervised_analysis<-function(dat,labels=NULL,scaled=c("row","column","none"),run_umap=T,umap_params=NULL,umap_config=umap::umap.defaults,run_tsne=T,tsne_params=NULL,run_pca=T,pca_params=NULL,run_mds=T,mds_params=NULL,run_heatmap=T,heatmap_params=NULL,cols){
  scaled=match.arg(scaled)
  if(scaled=='row'){
    dat<-na.omit(t(scale(t(dat))))
  }
  if(scaled=="column"){
    dat<-scale(dat)
  }
  if(missing(cols)){
    cols<-c("#006400","#00008b","#b03060","#ff4500","#ffd700","#7fff00","#00ffff","#ff00ff","#6495ed","#ffdab9")
    if(!is.null(labels)){
      if(length(unique(labels))>10){
        cols=rep("grey",length(unique(labels)))
        names(cols)<-unique(labels)
      } else {
        cols<-sample(cols,length(unique(labels)))
        names(cols)<-unique(labels)
      }
    }
  }
  if(run_umap){
    umap_res<-do.call(umap::umap,c(list(d=t(dat),config = umap_config),umap_params))
    layout<-data.frame(umap_res$layout)
    colnames(layout)<-c("L1","L2")
    umap_plot<-ggplot2::ggplot(data=layout,mapping = ggplot2::aes(x=L1,y=L2))+ggplot2::geom_point()+ggplot2::labs(title="umap_plot")
    if(!is.null(labels)){
      layout<-cbind(layout,Label=labels)
      umap_plot<-ggplot2::ggplot(data=layout,mapping = ggplot2::aes(x=L1,y=L2,col=Label))+ggplot2::geom_point()+ggplot2::labs(title="umap_plot")+ggplot2::scale_color_manual(values = cols,guide=F)
    }
  }
  if(run_tsne){
    dat_tsne<-unique(t(dat))
    dat_tsne<-as.matrix(dat_tsne)
    tsne_res<-do.call(Rtsne::Rtsne,c(list(X=dat_tsne, perplexity = min(floor((nrow(dat_tsne) - 1) / 3),100)),tsne_params))
    Y=data.frame(tsne_res$Y)
    colnames(Y)<-c("Y1","Y2")
    tsne_plot<-ggplot2::ggplot(data=Y,mapping = ggplot2::aes(x=Y1,y=Y2))+ggplot2::geom_point()+ggplot2::labs(title="tsne_plot")
    if(!is.null(labels)){
      Y<-cbind(Y,Label=labels)
      tsne_plot<-ggplot2::ggplot(data=Y,mapping = ggplot2::aes(x=Y1,y=Y2,col=Label))+ggplot2::geom_point()+ggplot2::labs(title="tsne_plot")+ggplot2::scale_color_manual(values = cols,guide=F)
    }
  }

  if(run_pca){
    pca_res<-do.call(PCAtools::pca,c(list(mat=t(dat)),pca_params))
    loadings<-data.frame(pca_res$loadings)
    colnames(loadings)<-paste("L",1:ncol(loadings),sep="")
    pca_plot<-ggplot2::ggplot(data=loadings,mapping=ggplot2::aes(x=L1,y=L2))+ggplot2::geom_point()+ggplot2::labs(title="pca_plot")
    if(!is.null(labels)){
      loadings<-cbind(loadings,Label=labels)
      pca_plot<-ggplot2::ggplot(data=loadings,mapping = ggplot2::aes(x=L1,y=L2,col=Label))+ggplot2::geom_point()+ggplot2::labs(title="pca_plot")+ggplot2::scale_color_manual(values = cols,guide=F)
    }
  }


  if(run_mds){
    distances<-dist(t(dat))
    fit_mds<-do.call(stats::cmdscale,c(list(d=distances,eig = T,k=2),mds_params))
    points<-data.frame(fit_mds$points)
    colnames(points)<-paste("P",1:ncol(points),sep="")
    mds_plot<-ggplot2::ggplot(data=points,mapping=ggplot2::aes(x=P1,y=P2))+ggplot2::geom_point()+ggplot2::labs(title="mds_plot")
    if(!is.null(labels)){
      points<-cbind(points,Label=labels)
      mds_plot<-ggplot2::ggplot(data=points,mapping = ggplot2::aes(x=P1,y=P2,col=Label))+ggplot2::geom_point()+ggplot2::labs(title="mds_plot")+ggplot2::scale_color_manual(values = cols,guide=F)
    }
  }

  if(run_heatmap){
    mads<-apply(dat,1,mad,na.rm=T)
    heatmap_dat<-dat[names(sort(mads,decreasing = T)),]
    heatmap_plot<-ggplotify::as.ggplot(do.call(ComplexHeatmap::Heatmap,c(list(matrix=heatmap_dat,show_row_names=F,show_column_names=F),heatmap_params)))+ggplot2::labs(title = "heatmap")
    if(!is.null(labels)){
      top_anno<-ComplexHeatmap::HeatmapAnnotation(df=data.frame(label=labels),col=list(label=cols))
      heatmap_plot<-ggplotify::as.ggplot(do.call(ComplexHeatmap::Heatmap,c(list(matrix=heatmap_dat,name="expressions",show_row_names=F,show_column_names=F,top_annotation=top_anno),heatmap_params)))+ggplot2::labs(title = "heatmap")
    }
  }
  gridExtra::grid.arrange(grobs=list(ggplot2::ggplotGrob(umap_plot),ggplot2::ggplotGrob(tsne_plot),ggplot2::ggplotGrob(pca_plot),ggplot2::ggplotGrob(mds_plot),ggplot2::ggplotGrob(heatmap_plot)),layout_matrix=rbind(c(1,2,3),
                                                                                                                                                              c(4,5,5)))
  return(list(umap=umap_res,tsne=tsne_res,pca=pca_res,mds=fit_mds))
}


# SOURCE: mRNA.R (2026-03-02)
consensus_immunedeconvolute<-function(expressions,methods=c("abis","bindea","cibersort","consensustme","danaher","davoli","dcq","deconseq","epic","ImmuCellAI","mcpcounter","quantiseq","timer","xcell"),celltype_mapping,bindea_reference,cibersort_reference,consensustme_indication,danaher_reference,davoli_reference,deconseq_reference,timer_indication,method_frequency_cutoff=2,background_noise=0.00001){
  results<-list()
  immune_abis<-NULL
  immune_bindea<-NULL
  immune_cibersort<-NULL
  immune_consensustme<-NULL
  immune_danaher<-NULL
  immune_davoli<-NULL
  immune_dcq<-NULL
  immune_deconseq<-NULL
  immune_epic<-NULL
  immune_immucellai<-NULL
  immune_mcpcounter<-NULL
  immune_quantiseq<-NULL
  immune_timer<-NULL
  immune_xcell<-NULL
  if("abis" %in% methods) {
    cat("running abis\n")
    immune_abis<-immunedeconv::deconvolute(expressions,method="abis",arrays=FALSE)
    immune_abis[,-1][immune_abis[,-1]<0]<-0
    immune_abis[,-1]<-t(t(immune_abis[,-1])/colSums(immune_abis[,-1]))
    assertthat::are_equal(colnames(expressions),colnames(immune_abis)[-1])
    results[["abis"]]<-immune_abis
  }

  if("bindea" %in% methods){
    cat("running bindea\n")
    bindea_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Bindea" & (!is.na(celltype_mapping$cell_type)),]
    bindea<-lapply(unique(bindea_reference[["Cell_Type"]]),function(ct){return(bindea_reference[["Symbol"]][bindea_reference[["Cell_Type"]]==ct])})
    names(bindea)<-unique(bindea_reference[["Cell_Type"]])
    immune_bindea<-corto::ssgsea(inmat=expressions,groups=bindea,minsize = 0)
    immune_bindea<-immune_bindea[bindea_mapping$method_cell_type,]
    rownames(immune_bindea)<-bindea_mapping$cell_type
    colnames(immune_bindea)<-colnames(expressions)
    immune_bindea[immune_bindea<0]<-0
    immune_bindea<-t(t(immune_bindea)/colSums(immune_bindea))
    immune_bindea<-data.frame(cell_type=rownames(immune_bindea),immune_bindea,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_bindea)[-1])
    results[["bindea"]]<-immune_bindea
  }

  if("cibersort" %in% methods){
    cat("running cibersort\n")
    cibersort_mapping<-celltype_mapping[celltype_mapping$method_dataset=="cibersort" & (!is.na(celltype_mapping$cell_type)),]
    immune_cibersort<-CIBERSORT(X=expressions,Y=cibersort_reference,perm=10)
    immune_cibersort<-immune_cibersort[cibersort_mapping$method_cell_type,]
    rownames(immune_cibersort)<-cibersort_mapping$cell_type
    immune_cibersort<-immune_cibersort[,colnames(expressions)]
    immune_cibersort<-data.frame(cell_type=rownames(immune_cibersort),immune_cibersort,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_cibersort)[-1])
    results[["cibersort"]]<-immune_cibersort
  }

  if("consensustme" %in% methods){
    cat("running consensustme\n")
    consensustme_mapping<-celltype_mapping[celltype_mapping$method_dataset=="consensus_tme",]
    immune_consensustme<-immunedeconv::deconvolute_consensus_tme(expressions,indications=rep(consensustme_indication,ncol(expressions)),method="ssgsea")
    consensustme_mapping<-consensustme_mapping[match(rownames(immune_consensustme),consensustme_mapping$method_cell_type),]
    rownames(immune_consensustme)<-consensustme_mapping$cell_type
    immune_consensustme[immune_consensustme<0]<-0
    immune_consensustme<-t(t(immune_consensustme)/colSums(immune_consensustme))
    immune_consensustme<-data.frame(cell_type=rownames(immune_consensustme),immune_consensustme,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_consensustme)[-1])
    results[["consensustme"]]<-immune_consensustme
  }

  if("danaher" %in% methods){
    cat("running danaher\n")
    danaher_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Danaher" & (!is.na(celltype_mapping$cell_type)),]
    danaher<-lapply(unique(danaher_reference[["Cell_Type"]]),function(ct){return(danaher_reference[["Symbol"]][danaher_reference[["Cell_Type"]]==ct])})
    names(danaher)<-unique(danaher_reference[["Cell_Type"]])
    immune_danaher<-corto::ssgsea(inmat=expressions,groups=danaher,minsize = 0)
    danaher_mapping<-danaher_mapping[danaher_mapping$method_cell_type %in% rownames(immune_danaher),]
    immune_danaher<-immune_danaher[danaher_mapping$method_cell_type,]
    rownames(immune_danaher)<-danaher_mapping$cell_type
    colnames(immune_danaher)<-colnames(expressions)
    immune_danaher[immune_danaher<0]<-0
    immune_danaher<-t(t(immune_danaher)/(colSums(immune_danaher,na.rm=T)+background_noise))
    immune_danaher<-data.frame(cell_type=rownames(immune_danaher),immune_danaher,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_danaher)[-1])
    results[["danaher"]]<-immune_danaher
  }

  if("davoli" %in% methods){
    cat("running davoli\n")
    davoli_mapping<-celltype_mapping[celltype_mapping$method_dataset=="Davoli" & (!is.na(celltype_mapping$cell_type)),]
    davoli<-lapply(unique(davoli_reference[["Cell_Type"]]),function(ct){return(davoli_reference[["Symbol"]][davoli_reference[["Cell_Type"]]==ct])})
    names(davoli)<-unique(davoli_reference[["Cell_Type"]])
    immune_davoli<-corto::ssgsea(inmat=expressions,groups=davoli,minsize = 0)
    davoli_mapping<-davoli_mapping[davoli_mapping$method_cell_type %in% rownames(immune_davoli),]
    immune_davoli<-immune_davoli[davoli_mapping$method_cell_type,]
    rownames(immune_davoli)<-davoli_mapping$cell_type
    colnames(immune_davoli)<-colnames(expressions)
    immune_davoli[immune_davoli<0]<-0
    immune_davoli<-t(t(immune_davoli)/(colSums(immune_davoli,na.rm=T)+background_noise))
    immune_davoli<-data.frame(cell_type=rownames(immune_davoli),immune_davoli,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_davoli)[-1])
    results[["davoli"]]<-immune_davoli
  }

  if("dcq" %in% methods){
    cat("running dcq\n")
    dcq_mapping=celltype_mapping[celltype_mapping$method_dataset=="dcq",]
    immune_dcq<-ADAPTS::estCellPercent.DCQ(refExpr=ADAPTS::LM22,geneExpr=expressions)
    rownames(immune_dcq)<-trimws(gsub("\\.+"," ",rownames(immune_dcq)))
    dcq_mapping<-dcq_mapping[na.omit(match(rownames(immune_dcq),dcq_mapping$method_cell_type)),]
    immune_dcq<-immune_dcq[dcq_mapping$method_cell_type,]
    rownames(immune_dcq)<-dcq_mapping$cell_type
    immune_dcq<-data.frame(cell_type=rownames(immune_dcq),immune_dcq/100,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_dcq)[-1])
    results[["dcq"]]<-immune_dcq
  }
  if("deconseq" %in% methods){
    cat("running deconseq\n")
    deconseq_mapping<-celltype_mapping[celltype_mapping$method_dataset=="deconseq" & (!is.na(celltype_mapping$cell_type)),]
    immune_deconseq<-as.data.frame(t(DeconRNASeq_(as.data.frame(expressions),deconseq_reference)$out.all))
    immune_deconseq<-immune_deconseq[deconseq_mapping$method_cell_type,]
    rownames(immune_deconseq)<-deconseq_mapping$cell_type
    colnames(immune_deconseq)<-colnames(expressions)
    immune_deconseq<-data.frame(cell_type=rownames(immune_deconseq),immune_deconseq,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_deconseq)[-1])
    results[["deconseq"]]<-immune_deconseq

  }
  if("epic" %in% methods){
    cat("running epic\n")
    immune_epic<-immunedeconv::deconvolute(expressions,method="epic",tumor=T,scale_mrna = T)
    assertthat::are_equal(colnames(expressions),colnames(immune_epic)[-1])
    results[["epic"]]<-immune_epic
  }
  try(if("ImmuCellAI" %in% methods){
    cat("running ImmuCellAI\n")
    immune_immucellai<-ImmuCellAI::ImmuCellAI_new(expressions,data_type = "rnaseq",group_tag = FALSE,response_tag = FALSE)
    immune_immucellai<-data.frame(cell_type=rownames(as.data.frame(t(immune_immucellai$Sample_abundance))),as.data.frame(t(immune_immucellai$Sample_abundance)))
    assertthat::are_equal(colnames(expressions),colnames(immune_immucellai)[-1])
    results[["immucellai"]]<-immune_immucellai
  },silent = T)

  if("mcpcounter" %in% methods){
    cat("running mcpcounter\n")
    mcp_counter_mapping<-celltype_mapping[celltype_mapping$method_dataset=="mcp_counter" & (!is.na(celltype_mapping$cell_type)),]
    immune_mcpcounter<-immunedeconv::deconvolute_mcp_counter(expressions)
    immune_mcpcounter<-t(t(immune_mcpcounter)/colSums(immune_mcpcounter))
    mcp_counter_mapping<-mcp_counter_mapping[match(rownames(immune_mcpcounter),mcp_counter_mapping$method_cell_type),]
    rownames(immune_mcpcounter)<-mcp_counter_mapping$cell_type
    immune_mcpcounter<-data.frame(cell_type=rownames(immune_mcpcounter),immune_mcpcounter,check.names = F)
    assertthat::are_equal(colnames(expressions),colnames(immune_mcpcounter)[-1])
    results[["mcpcounter"]]<-immune_mcpcounter
  }

  if("quantiseq" %in% methods){
    cat("running quantiseq\n")
    immune_quantiseq<-immunedeconv::deconvolute(expressions,method='quantiseq',tumor=T,arrays = F,scale_mrna = T)
    assertthat::are_equal(colnames(expressions),colnames(immune_quantiseq)[-1])
    results[["quantiseq"]]<-immune_quantiseq
  }

  if("timer" %in% methods){
    cat("running timer\n")
    immune_timer<-immunedeconv::deconvolute(expressions,method="timer",indications = rep(timer_indication,ncol(expressions)))
    assertthat::are_equal(colnames(expressions),colnames(immune_timer)[-1])
    results[["timer"]]<-immune_timer
  }

  if("xcell" %in% methods){
    cat("running xcell\n")
    immune_xcell<-immunedeconv::deconvolute(expressions,method="xcell",tumor=T,arrays=F)
    assertthat::are_equal(colnames(expressions),colnames(immune_xcell)[-1])
    results[["xcell"]]<-immune_xcell
  }

  immune_consensus<-rbind(immune_abis,immune_bindea,immune_cibersort,immune_consensustme,immune_danaher,immune_davoli,immune_dcq,immune_epic,immune_mcpcounter,immune_quantiseq,immune_timer,immune_xcell)
  method_frequency<-table(immune_consensus$cell_type)
  results[["method_frequency"]]<-method_frequency
  consensus_cell_type<-setdiff(names(method_frequency[method_frequency>=method_frequency_cutoff]),"uncharacterized cell")
  immune_consensus<-immune_consensus[immune_consensus$cell_type %in% consensus_cell_type,]
  immune_consensus[,-1]<-immune_consensus[,-1]+background_noise
  immune_consensus<-immune_consensus %>% dplyr::group_by(cell_type) %>% dplyr::summarise_all(function(x){exp(mean(log(x),na.rm=T))})
  immune_consensus<-immune_consensus[rowSums(immune_consensus[,-1]<background_noise)<(ncol(expressions)/3),]
  results[["consensus"]]<-immune_consensus
  return(results)
}


# SOURCE: mRNA.R (2026-03-02)
geneset_activity<-function(expressions,scale=c("row","column","none"),genesets,methods=c("ssgsea","gsva","zscore"),aggregate=T){
  scale=match.arg(scale)
  if(scale=="row"){
    expressions<-t(scale(t(expressions)))
  }
  if(scale=="column"){
    expressions=t(scale(expressions))
  }
  ssgsea_activity=NULL
  gsva_activity=NULL
  zscore_activity=NULL
  agg_activity=NULL
  if("ssgsea" %in% methods){
    ssgsea_activity=corto::ssgsea(expressions,genesets)
    colnames(ssgsea_activity)<-colnames(expressions)
  }
  if("gsva" %in% methods){
    gsva_activity=GSVA::gsva(expressions,genesets)
    colnames(gsva_activity)<-colnames(expressions)
  }
  if("zscore" %in% methods){
    zscore_activity<-expressions[FALSE,]
    for(gs in names(genesets)){
      gs_activity<-colMeans(expressions[intersect(rownames(expressions),genesets[[gs]]),],na.rm=T)
      zscore_activity<-rbind(zscore_activity,gs_activity)
    }
    rownames(zscore_activity)<-names(genesets)
    colnames(zscore_activity)<-colnames(expressions)
  }
  if(aggregate){
    agg_activity=data.frame(matrix(0,nrow=length(genesets),ncol = ncol(expressions)))
    rownames(agg_activity)<-names(genesets)
    colnames(agg_activity)<-colnames(expressions)
    for(method in methods){
      agg_activity<-agg_activity+list(ssgsea_activity=ssgsea_activity,gsva_activity=gsva_activity,zscore_activity=zscore_activity)[[paste(method,"activity",sep="_")]]
    }
    agg_activity<-agg_activity/length(methods)
  }
  return(list(ssgsea_activity=ssgsea_activity,gsva_activity=gsva_activity,zscore_activity=zscore_activity,agg_activity=agg_activity))
}
