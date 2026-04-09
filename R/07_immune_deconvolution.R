# =============================================================================
# IMMUNE DECONVOLUTION ALGORITHM IMPLEMENTATIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   DeconRNASeq_        - DeconRNASeq with bug fixes (bulk RNA immune deconvolution)
#   CoreAlg             - SVM core algorithm for CIBERSORT
#   doPerm              - Permutation test for CIBERSORT p-value
#   CIBERSORT           - CIBERSORT immune cell fraction estimation
#   immu_deconvolution  - High-level wrapper: deconvolute + test + plot
#
# NOTE: consensus_immunedeconvolute (multi-method consensus wrapper) is in
#       03_expression_analysis.R
#
# SOURCE PROVENANCE:
#   DeconRNASeq_/CoreAlg/doPerm/CIBERSORT:
#     woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02)
#   immu_deconvolution:
#     Pembro/funcsInPembro.R (2025-07-02)
# =============================================================================

#' DeconRNASeq_
#'
#' a replacement for DeconRNASeq, fix some bugs
#'
#' @param datasets data.frame. expression data.frame
#' @param signatures data.frame. reference signature data.frame
#' @param proportions data.frame. proportion matrix from different tissue/cell types.
#' @param checksig logic. whether the condition number of signature matrix should be checked, default = FALSE
#' @param known.prop logic. whether the proportions of cell types have been known in advanced for proof of concept, default = FALSE
#' @param use.scale logic. whether the data should be centered or scaled, default = TRUE
#' @param fig logic. whether to generate the scatter plots of the estimated cell fractions vs. the true proportions of cell types, default = TRUE
#'
#' @return list. a list containing estimated cell type fraction matrix for all the mixture samples and svd calculated PCA on the mixture samples to estimate the number of pure sources according to the cumulative R2
#' @export
#'

DeconRNASeq_<-function (datasets, signatures, proportions = NULL, checksig = FALSE,
                        known.prop = FALSE, use.scale = TRUE, fig = TRUE)
{
  if (is.null(datasets))
    stop(" Missing the mixture dataset, please provide a tab-delimited text file for mixture samples.")
  if (is.null(signatures))
    stop(" Missing the signature dataset, please provide a tab-delimited text file for pure tissue/cell types.")
  if (is.null(proportions) && known.prop)
    stop(" Missing the known proprotions, please provide a tab-delimited text file containing known fractions for pure tissue/cell types.")
  x.signature <- signatures
  x.data <- datasets
  if (is.data.frame(x.signature) == F)
    stop("signature datasets must be a dataframe")
  if (sum(is.na(x.signature)) > 0)
    stop("signature data cannot have NAs. please exclude or impute missing values.")
  if (is.data.frame(x.data) == F)
    stop("mixture datasets must be a dataframe")
  if (sum(is.na(x.data)) > 0)
    stop("mixture data cannot have NAs. please exclude or impute missing values.")
  numofg <- nrow(x.signature)
  Numofx <- ncol(x.signature)
  if (numofg < Numofx)
    stop("The number of genes is less than the number of cell types, which means less independent equations than unknowns.")
  x.data.temp <- pcaMethods::prep(x.data, scale = "none", center = TRUE)
  x.data.pca <- pcaMethods::pca(x.data.temp, method = "svd", center = FALSE, nPcs = Numofx)
  out.pca <- summary(x.data.pca)
  Var <- pcaMethods::R2cum(x.data.pca)
  numofmix <- order(Var > 0.99, decreasing = T)[1]
  if (fig) {
    if (Numofx != numofmix) {
      cat("\n Attention: the number of pure cell types =",
          Numofx, " defined in the signature matrix;\n")
      cat("\n PCA results indicate that the number of cell types in the mixtures =",
          numofmix, "\n")
    }
  }
  if (checksig && numofg >= 40) {
    step <- seq(20, numofg, by = 20)
    sig.cond <- sapply(step, function(x) kappa(scale(x.signature[1:x,
    ])))
    DeconRNASeq::condplot(step, sig.cond)
  }
  common.signature <- rownames(x.signature) %in% rownames(x.data)
  common.data <- rownames(x.data) %in% rownames(x.signature)
  x.data <- x.data[common.data, ]
  x.signature <- x.signature[common.signature, ]
  x.subdata <- x.data[rownames(x.signature), ]
  Numofx <- ncol(x.signature)
  if (use.scale) {
    AA <- scale(x.signature)
  }
  else {
    AA <- x.signature
  }
  EE <- rep(1, Numofx)
  FF <- 1
  GG <- diag(nrow = Numofx)
  HH <- rep(0, Numofx)
  out.all <- c()
  for (i in colnames(x.subdata)) {
    BB <- x.subdata[, i]
    if (use.scale) {
      BB <- scale(BB)
    }
    out <- limSolve::lsei(AA, BB, EE, FF, GG, HH)
    out.all <- rbind(out.all, out$X)
  }
  mean.rmse <- 0
  rmse <- c()
  if (known.prop) {
    x.proportions <- proportions
    x.proportions <- x.proportions[colnames(x.data), ]
    parray <- ggplot2::ggplot()
    length(parray) <- ncol(out.all)
    for (i in 1:ncol(out.all)) {
      A.error <- DeconRNASeq::rmse(x.proportions[, i], out.all[, i])
      rmse <- c(rmse, A.error)
      x <- out.all[, i]
      y <- x.proportions[, i]
      xlabel <- paste("estimated ", colnames(x.proportions)[i])
      ylabel <- paste("actual ", colnames(x.proportions)[i])
      main <- sprintf("Scatter plot of proportions,\n RMSE = %.3f",
                      A.error)
      parray[[i]] <- ggplot2::ggplot(data.frame(x, y), ggplot2::aes(x, y)) +
        ggplot2::geom_point(alpha = 0.3) + ggplot2::labs(title = main) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red",
                    size = 1) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)
    }
    if (fig) {
      DeconRNASeq::multiplot(plotlist = parray, cols = 2)
    }
    mean.rmse <- mean(rmse)
  }
  if (known.prop) {
    return(list(out.all = out.all, out.pca = out.pca, out.rmse = mean.rmse))
  }
  else {
    return(list(out.all = out.all, out.pca = out.pca))
  }
}


#' CoreAlg
#' Core algorithm
#' @param X matrix. cell-specific gene expression
#' @param y matrix. mixed expression per sample
#' @export
#' @details Source provenance: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02).
# SOURCE: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02)
CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' doPerm
#' do permutations
#' @param perm numeric. Number of permutations.
#' @param X matrix. cell-specific gene expression.
#' @param Y matrix. mixed expression per sample.
#' @export
#' @details Source provenance: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02).
# SOURCE: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02)
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' CIBERSORT
#' Main functions
#' @param X matrix. ig_matrix file path to gene expression from isolated cells
#' @param Y matrix. mixture_file heterogenous mixed expression
#' @param perm numeric. perm Number of permutations
#' @param QN logic. Perform quantile normalization or not (TRUE/FALSE)
#' @export
#' @details Source provenance: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02).
# SOURCE: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R (2026-03-02)
CIBERSORT <- function(X, Y, perm=0, QN=TRUE){

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

# SOURCE: Pembro/funcsInPembro.R (2025-07-02)
#' High-level wrapper: deconvolute + test + plot
#'
#' High-level wrapper: deconvolute + test + plot
#' @param project_dir Output plotting or file parameter used by the existing implementation.
#' @param gene_expressions Function argument documented from the legacy interface.
#' @param sample_info Input data frame or matrix.
#' @param ID_col Column name used by the existing implementation.
#' @param cf_algorithm Function argument documented from the legacy interface.
#' @param log2_transform_cf Function argument documented from the legacy interface.
#' @param background Function argument documented from the legacy interface.
#' @param mode Option controlling how the function runs.
#' @param test_method Function argument documented from the legacy interface.
#' @param paired Logical flag controlling optional behavior.
#' @param group_col Column name used by the existing implementation.
#' @param contrasts Function argument documented from the legacy interface.
#' @param pvalue_cutoff Numeric tuning parameter used by the existing implementation.
#' @param top_anno Function argument documented from the legacy interface.
#' @param ... Additional arguments passed through to downstream functions.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: Pembro/funcsInPembro.R (2025-07-02).
#'
#' @examples
#' \dontrun{
#' immu_deconvolution(...)
#' }
#' @export
immu_deconvolution<-function(project_dir,gene_expressions,sample_info,ID_col,cf_algorithm='xcell',log2_transform_cf=TRUE,background=0.00001,mode="supervised",test_method=c("ttest","limma"),paired=F,group_col,contrasts,pvalue_cutoff=0.05,top_anno=NULL,...){
  require(immunedeconv)
  test_method<-match.arg(test_method)
  if(!dir.exists(project_dir)){
    dir.create(project_dir)
  }
  stopifnot(all(colnames(gene_expressions)==sample_info[[ID_col]]))
  cellular_fraction<-deconvolute(gene_expressions,cf_algorithm,tumor = T)
  cellular_fraction<-as.data.frame(cellular_fraction)
  rownames(cellular_fraction)<-cellular_fraction$cell_type
  cellular_fraction<-cellular_fraction[,-1]
  write.csv(cellular_fraction,file.path(project_dir,"cellular_fraction.csv"))
  cf<-cellular_fraction
  if(log2_transform_cf){
    log2_cellular_fraction<-log2(cellular_fraction+background)
    cf<-log2_cellular_fraction
    write.csv(log2_cellular_fraction,file.path(project_dir,"log2_cellular_fratcion.csv"))
  }

  if(mode=="supervised"){
    groups<-unlist(strsplit(contrasts,"-"))
    if(test_method=="ttest"){
      statistics=apply(cf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_col]]);df<-df[df$group %in% groups,];test=t.test(value~group,data=df,paired=paired);return(test$p.value)})
    }
    if(test_method=='limma'){
      if(!log2_transform_cf){
        cf=log2(cf+background)
      }
      formula<-as.formula(paste("~0+",group_col,sep=""))
      design<-model.matrix(formula,data=sample_info)
      colnames(design)<-levels(factor(sample_info[[group_col]]))
      fit<-lmFit(cf,design)
      contrast.matrix<-makeContrasts(contrasts=contrasts,levels = design)
      fit <- contrasts.fit(fit,contrast.matrix)
      fit<-eBayes(fit)
      statistics<-topTable(fit,number=nrow(cf))
    }

    mean_<-list()
    for (group in groups){
      mean_[[paste(group,"_mean",sep="")]]<-rowMeans(cf[,sample_info[[group_col]]==group],na.rm=T)
    }
    mean_<-as.data.frame(mean_)
    if(test_method=="ttest"){
      results<-data.frame(mean_[names(statistics),],P.Value=statistics)
      results[["adj.P.Val"]]<-p.adjust(results[["P.Value"]],method='BH')
    }
    if(test_method=="limma"){
      results<-cbind(mean_[rownames(statistics),],statistics)
    }
    write.csv(results,file.path(project_dir,paste(group_col,test_method,"statistics.csv",sep="_")))
    stopifnot("rownames of cellular fraction should be same as rownames of statistics results"=all(rownames(results)==rownames(cf)))
  }
  filtered_cf<-cf[rowSums(cellular_fraction>background)>=(ncol(cf)/2),]
  if(mode=="supervised"){
    neg_log_pvalue=-log10(results[rownames(filtered_cf),][["P.Value"]])
    row_anno<-rowAnnotation(neg_log_pvalue=row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pvalue_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pvalue_cutoff)))+0.5)),width=unit(3,"cm"),annotation_name_side="top",annotation_name_rot=0)
    h<-Heatmap(t(scale(t(filtered_cf))),name="cellular_fraction",cluster_columns = cluster_within_group(t(scale(t(filtered_cf))),sample_info[[group_col]]),show_row_names = T,show_column_dend = F,show_column_names = T,right_annotation = row_anno,top_annotation = top_anno,...)
    tiff(filename = file.path(project_dir,"immun_deconvolution_heatmap.tiff"),width = 10,height=8,units = "in",res=300,compression = "lzw")
    draw(h)
    xpos<-0.95/(max(c(neg_log_pvalue,-log10(pvalue_cutoff)))+0.5)*(-log10(pvalue_cutoff))
    decorate_annotation("neg_log_pvalue",{
      grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
    })
    dev.off()
    return(list(cf=cf,statistics=results))
  } else {
    h<-Heatmap(t(scale(t(filtered_cf))),name="cellular_fraction",show_row_names = T,show_column_names = T,top_annotation = top_anno,...)
    tiff(filename = file.path(project_dir,"immun_deconvolution_heatmap.tiff"),width = 10,height=8,units = "in",res=300,compression = "lzw")
    draw(h)
    dev.off()
    return(cf)
  }

}
