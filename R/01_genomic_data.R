# =============================================================================
# GENOMIC DATA PROCESSING FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   1. isar_transform         - Copy number segment transformation (purity/ploidy)
#   2. df2granges             - Data frame to GRanges conversion
#   3. df2grangelist          - Data frame to GRangesList (by sample)
#   4. bingranges             - Bin genome into fixed windows
#   5. gr2linear              - GRanges to linear coordinates
#   6. column2namedVector     - Data frame column to named vector
#   7. align_df               - Align data frame to reference order
#   8. chisq_dist             - Chi-squared distance matrix
#   9. ggcopynumber           - ggplot2 copy number visualization
#  10. ggscores               - Genomic scores plot (Amp/Del)
#  11. ggarms                 - Chromosomal arm CNV plot
#  12. gggenome               - Genome-wide CNV visualization (multi-sample)
#
# SOURCE PROVENANCE:
#   functions.R: woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R (2026-03-02)
#   myscripts.R: woodman_lab.XLi23/myscripts.R (2023-06-29)
#   DNAfuncs.R:  IBC/DNAfuncs.R (2024-03-18)
# =============================================================================

# --- SECTION 1: COORDINATE UTILITIES ---

# SOURCE: myscripts.R
#' Adjust copy-number segment means for purity and ploidy.
#'
#' Adjust copy-number segment means for purity and ploidy.
#' @param segmean Function argument documented from the legacy interface.
#' @param purity Function argument documented from the legacy interface.
#' @param ploidy Function argument documented from the legacy interface.
#' @return A numeric vector of transformed segment means.
#' @details Source provenance: myscripts.R.
#'
#' @examples
#' segmean <- c(-0.3, 0, 0.2)
#' isar_transform(segmean, purity = 0.7, ploidy = 2)
#' @export
isar_transform<-function(segmean,purity,ploidy){
  #ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3966983/#R75
  # the formula R'(x) = q(x)/τ = R(x)/α − 2(1 − α)/(ατ) is wrong.
  # suggested parameters: noise threshold of 0.3, a broad length cutoff of 0.5 chromosome arms, a confidence level of 95%, and a copy-ratio cap of 1.5
  ratio<-2^segmean
  ratio_=(purity*ploidy*ratio+2*(1-purity)*ratio-2*(1-purity))/(purity*ploidy)
  segmean_<-log2(ratio_)
  segmean_[is.na(segmean_)]<--10 #very low coverage set as -10
  return(segmean_)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Convert a genomic interval data frame to a `GRanges` object.
#'
#' Convert a genomic interval data frame to a `GRanges` object.
#' @param df Input data frame or matrix.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param seqnames_col Column name used by the existing implementation.
#' @param start_col Column name used by the existing implementation.
#' @param end_col Column name used by the existing implementation.
#' @param strand_col Column name used by the existing implementation.
#' @param meta_cols Function argument documented from the legacy interface.
#' @return A `GenomicRanges::GRanges` object.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' df2granges(df = ..., genome = ...)
#' }
#' @export
df2granges<-function(df,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=TRUE,xy=FALSE,seqnames_col="chromosome",start_col="start",end_col="end",strand_col=NULL,meta_cols=NULL){
  genome=match.arg(genome)
  if(genome=='hg19'){
    txdb<-BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    txdb<-BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  seqinfo_<-GenomeInfoDb::seqinfo(txdb)
  seqs=GenomeInfoDb::seqlevels(txdb)
  if(simplified){
    seqs=GenomeInfoDb::seqlevels(txdb)[1:24]
    if(!xy){
      seqs=GenomeInfoDb::seqlevels(txdb)[1:22]
    }
    seqinfo_<-seqinfo_[seqs]
  }
  df[[seqnames_col]]<-paste("chr",df[[seqnames_col]],sep="")
  df<-df[df[[seqnames_col]] %in% seqs,]
  if(is.null(strand_col)){
    strand="*"
  } else {
    strand=df[[strand_col]]
  }
  if(is.null(meta_cols)) {
    mcols=NULL
  } else {
    mcols=df[,meta_cols,drop=F]
  }
  granges<-GenomicRanges::GRanges(seqnames=df[[seqnames_col]],ranges=IRanges::IRanges(start=df[[start_col]],end=df[[end_col]]),strand = strand,mcols=mcols,seqinfo=seqinfo_)
  colnames(S4Vectors::mcols(granges))=meta_cols
  return(granges)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Split a genomic interval data frame into a sample-wise `GRangesList`.
#'
#' Split a genomic interval data frame into a sample-wise `GRangesList`.
#' @param df Input data frame or matrix.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param sample_col Column name used by the existing implementation.
#' @param seqnames_col Column name used by the existing implementation.
#' @param start_col Column name used by the existing implementation.
#' @param end_col Column name used by the existing implementation.
#' @param strand_col Column name used by the existing implementation.
#' @param meta_cols Function argument documented from the legacy interface.
#' @return A named `GenomicRanges::GRangesList` object.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' df2grangelist(df = ..., genome = ...)
#' }
#' @export
df2grangelist<-function(df,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=TRUE,xy=FALSE,sample_col="sample",seqnames_col="chromosome",start_col="start",end_col="end",strand_col=NULL,meta_cols=NULL){
  genome=match.arg(genome)
  seqlevelsStyle=match.arg(seqlevelsStyle)
  samples<-unique(df[[sample_col]])
  grl<-list()
  grl<-lapply(samples,function(samp){print(samp);sub_df=df[df[[sample_col]]==samp,];
  gr=df2granges(sub_df,genome = genome,seqlevelsStyle = seqlevelsStyle,simplified = simplified,xy=xy,seqnames_col = seqnames_col,start_col = start_col,end_col = end_col,strand_col = strand_col,meta_cols = meta_cols);
  #  gr$REF<-DNAStringSet(mcols(gr)[[meta_cols[1]]],start=1,width=nchar(mcols(gr)[[meta_cols[1]]]));
  #  ALT_l<-lapply(mcols(gr)[[meta_cols[2]]],function(alt){DNAStringSet(x=alt,start=1,width=nchar(alt))});
  #  gr$ALT<-DNAStringSetList(ALT_l);
  return(gr)})
  grl<-GenomicRanges::GRangesList(grl)
  names(grl)<-samples
  return(grl)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Tile a genome into fixed windows and transfer overlapping metadata.
#'
#' Tile a genome into fixed windows and transfer overlapping metadata.
#' @param gr Genomic ranges input used by the function.
#' @param window Numeric tuning parameter used by the existing implementation.
#' @return A binned `GenomicRanges::GRanges` object with metadata columns copied from overlaps.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' bingranges(gr = ..., window = ...)
#' }
#' @export
bingranges<-function(gr,window=500000){
  seqinfo_<-GenomeInfoDb::seqinfo(gr)
  bingenome<-lapply(1:length(seqinfo_),function(i){chromosome=seqinfo_@seqnames[i];
    start=seq(1,seqinfo_@seqlengths[i],by=window);
    end=start+window-1;
    end[length(end)]=seqinfo_@seqlengths[i];
    return(GenomicRanges::GRanges(seqnames=chromosome,ranges = IRanges::IRanges(start=start,end=end),strand="*",seqinfo = seqinfo_))}
    )
  bingenome<-do.call(c,bingenome)
  default_dict<-list("character"="","logical"=NA,"integer"=0L,"numeric"=0)
  for (col in colnames(S4Vectors::mcols(gr))){
    S4Vectors::mcols(bingenome)[[col]]<-default_dict[[class(S4Vectors::mcols(gr)[[col]])]]
  }
  ind<-GenomicRanges::findOverlaps(bingenome,gr,select="first")
  S4Vectors::mcols(bingenome)[!is.na(ind),]<-S4Vectors::mcols(gr)[ind[!is.na(ind)],]
  return(bingenome)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Convert genomic coordinates to linear genome positions.
#'
#' Convert genomic coordinates to linear genome positions.
#' @param gr Genomic ranges input used by the function.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @return A data frame with chromosome, start, end, and metadata columns on a linear genome scale.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' gr2linear(gr = ..., genome = ...)
#' }
#' @export
gr2linear<-function(gr,genome=c("hg19",'hg38'),seqlevelsStyle=c("NCBI","UCSC"),simplified=T,xy=F){
  genome_=match.arg(genome)
  seqlevelsStyle_=match.arg(seqlevelsStyle)
  if((!is.null(gr)) & (!any(is.na(GenomeInfoDb::seqlengths(gr))))){
    seqinfo_=GenomeInfoDb::seqinfo(gr)
  } else {
    if(genome=='hg19'){
      txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else {
      txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    GenomeInfoDb::seqlevelsStyle(txdb)<-seqlevelsStyle
    seqinfo_<-GenomeInfoDb::seqinfo(txdb)
    seqs=GenomeInfoDb::seqlevels(txdb)
    if(simplified){
      seqs=GenomeInfoDb::seqlevels(txdb)[1:24]
      if(!xy){
        seqs=GenomeInfoDb::seqlevels(txdb)[1:22]
      }
      seqinfo_<-seqinfo_[seqs]
    }
  }
  linear_genome<-as.data.frame(seqinfo_)[,"seqlengths",drop=F]
  linear_genome$End<-cumsum(as.numeric(linear_genome$seqlengths))
  linear_genome$Start<-c(1,linear_genome[["End"]][-nrow(linear_genome)]+1)
  linear_genome<-linear_genome[,c("Start","End")]

  linear_df<-cbind(data.frame(Chromosome=GenomeInfoDb::seqnames(gr),Start=BiocGenerics::start(gr),End=BiocGenerics::end(gr)),S4Vectors::mcols(gr))
  linear_df$Start<-linear_df$Start+linear_genome$Start[match(linear_df$Chromosome,rownames(linear_genome))]-1
  linear_df$End<-linear_df$End+linear_genome$Start[match(linear_df$Chromosome,rownames(linear_genome))]-1
  return(linear_df)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Convert one data-frame column to a named vector.
#'
#' Convert one data-frame column to a named vector.
#' @param df Input data frame or matrix.
#' @param column Function argument documented from the legacy interface.
#' @return A named vector built from the requested column.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' df <- data.frame(score = c(1, 2), row.names = c("A", "B"))
#' column2namedVector(df, "score")
#' @export
column2namedVector<-function(df,column){
  vec<-df[[column]]
  names(vec)<-rownames(df)
  return(vec)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Align rows or columns to the order defined by a reference object.
#'
#' Align rows or columns to the order defined by a reference object.
#' @param aligned_df Function argument documented from the legacy interface.
#' @param aligned_name Function argument documented from the legacy interface.
#' @param keep_names Function argument documented from the legacy interface.
#' @param refered_df Function argument documented from the legacy interface.
#' @param refered_name Function argument documented from the legacy interface.
#' @return A reordered data frame.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' align_df(aligned_df = ..., aligned_name = ...)
#' }
#' @export
align_df<-function(aligned_df,aligned_name,keep_names=NULL,refered_df,refered_name){
  if(refered_name=="rownames"){
    refered_=rownames(refered_df)
  } else  if(refered_name=="colnames"){
    refered_=colnames(refered_df)
  } else{
    refered_=refered_df[[refered_name]]
  }
  stopifnot(length(refered_)==length(unique(refered_)))
  if(aligned_name=='rownames'){
    if(is.numeric(keep_names)){
      keep_df<-aligned_df[keep_names,,drop=F]
      aligned_df<-aligned_df[-keep_names,,drop=F]
    } else if(is.character(keep_names)){
      keep_df<-aligned_df[rownames(aligned_df) %in% keep_names,,drop=F]
      aligned_df<-aligned_df[!rownames(aligned_df) %in% keep_names,,drop=F]
    } else {
      keep_df<-NULL
    }
    aligned_<-rbind(aligned_df[refered_,,drop=F],keep_df)
  } else if(aligned_name=="colnames"){
    if(is.numeric(keep_names)){
      keep_df<-aligned_df[,keep_names,drop=F]
      aligned_df<-aligned_df[,-keep_names,drop=F]
    } else if(is.character(keep_names)){
      keep_df<-aligned_df[,colnames(aligned_df) %in% keep_names,drop=F]
      aligned_df<-aligned_df[,!colnames(aligned_df) %in% keep_names,drop=F]
    } else {
      keep_df<-NULL
    }
    aligned_<-cbind(aligned_df[,refered_],keep_df)
  } else {
    if(is.numeric(keep_names)){
      keep_df<-aligned_df[,keep_names,drop=F]
      aligned_df<-aligned_df[,-keep_names,drop=F]
    } else if(is.character(keep_names)){
      keep_df<-aligned_df[,aligned_df[[aligned_name]] %in% keep_names,drop=F]
      aligned_df<-aligned_df[,!aligned_df[[aligned_name]] %in% keep_names,drop=F]
    } else {
      keep_df<-NULL
    }
    aligned_<-rbind(aligned_df[match(refered_,aligned_df[[aligned_name]]),,drop=F],keep_df)
  }
  return(aligned_)
}

# SOURCE: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R)
#' Compute a pairwise chi-squared distance matrix.
#'
#' Compute a pairwise chi-squared distance matrix.
#' @param df Input data frame or matrix.
#' @return A `stats::dist` object.
#' @details Source provenance: functions.R (woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R).
#'
#' @examples
#' \dontrun{
#' chisq_dist(df = ...)
#' }
#' @export
chisq_dist<-function(df){
  statistics<-apply(df,2,function(col_1){apply(df,2,function(col_2){p<-stats::chisq.test(col_1,col_2,correct = T);return(p$statistic)})})
  return(stats::as.dist(statistics))
}

# --- SECTION 2: COPY NUMBER VISUALIZATION ---

# SOURCE: myscripts.R
#' Plot segmented copy-number profiles across the genome.
#'
#' Plot segmented copy-number profiles across the genome.
#' @param segments Function argument documented from the legacy interface.
#' @param sample_order Function argument documented from the legacy interface.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param chromosome_colors Color specification used when plotting.
#' @param anno_df Function argument documented from the legacy interface.
#' @param anno_sample_col Column name used by the existing implementation.
#' @param anno_seqnames_col Column name used by the existing implementation.
#' @param anno_start_col Column name used by the existing implementation.
#' @param anno_end_col Column name used by the existing implementation.
#' @param anno_meta_cols Function argument documented from the legacy interface.
#' @param anno_gene_text_size Function argument documented from the legacy interface.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: myscripts.R.
#'
#' @examples
#' \dontrun{
#' ggcopynumber(segments = ..., sample_order = ...)
#' }
#' @export
ggcopynumber<-function(segments=NULL,sample_order=NULL,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=T,xy=F,chromosome_colors=c("grey","black"),
                       anno_df=NULL,anno_sample_col="Tumor_Sample_Barcode",anno_seqnames_col = "Chromosome",anno_start_col = "Start_position",anno_end_col = "End_Position",anno_meta_cols = c("Hugo_Symbol","Variant_Classification","Variant_Type"),
                       anno_gene_text_size=2){
  if((!is.factor(segments[["sample"]])) | (!is.ordered(segments[["sample"]]))){
    segments[["sample"]]<-factor(segments[["sample"]],levels=sort(unique(segments[["sample"]]),decreasing = T),ordered=T)
  }
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  genome_=match.arg(genome)
  seqlevelsStyle_=match.arg(seqlevelsStyle)
  if(genome=='hg19'){
    txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  seqlevelsStyle(txdb)<-seqlevelsStyle
  seqinfo_<-seqinfo(txdb)
  seqs=seqlevels(txdb)
  if(simplified){
    seqs=seqlevels(txdb)[1:24]
    if(!xy){
      seqs=seqlevels(txdb)[1:22]
    }
    seqinfo_<-seqinfo_[seqs]
  }
  linear_genome<-as.data.frame(seqinfo_)[,"seqlengths",drop=F]
  linear_genome$End<-cumsum(as.numeric(linear_genome$seqlengths))
  linear_genome$Start<-c(1,linear_genome[["End"]][-nrow(linear_genome)]+1)
  linear_genome<-linear_genome[,c("Start","End")]
  linear_genome$Mid<-(linear_genome$Start+linear_genome$End)/2
  linear_genome$nudge<-c(-0.1,-0.3)
  linear_genome$fill<-chromosome_colors

  segments_l<-lapply(unique(segments$sample),function(samp){
    seg<-segments[segments[["sample"]]==samp,]
    seg_gr<-df2granges(seg,genome="hg19",seqlevelsStyle = "NCBI",seqnames_col = "chromosome",start_col = "start",end_col = "end",meta_cols = c("sample","probes","segmean"))
    seg_ldf<-gr2linear(seg_gr)
    return(seg_ldf)
  })
  segments_ldf<-do.call(rbind,segments_l)
  if(is.null(sample_order)){
    sample_order=sort(unique(segments_ldf[["sample"]]),decreasing = T)
  }
  segments_ldf$sample<-factor(segments_ldf[["sample"]],levels=sample_order,ordered = T)

  if(!is.null(anno_df)){
    anno_ldf<-lapply(unique(anno_df[[anno_sample_col]]),function(samp){data=anno_df[anno_df[[anno_sample_col]]==samp,];
    gr=df2granges(df=data,genome = "hg19",seqlevelsStyle = "NCBI",seqnames_col = anno_seqnames_col,start_col = anno_start_col,end_col = anno_end_col,meta_cols = anno_meta_cols);
    ldf<-gr2linear(gr,genome = "hg19",seqlevelsStyle = "NCBI")
    ldf[["sample"]]<-samp
    return(ldf)})
    anno_ldf<-do.call(rbind,anno_ldf)
  }
  anno_ldf$sample<-factor(anno_ldf[["sample"]],levels=sample_order,ordered = T)
  anno_ldf<-as.data.frame(anno_ldf %>% group_by(sample,Hugo_Symbol) %>% mutate(hits=n()))
  anno_ldf<-anno_ldf[!duplicated(anno_ldf[,c("sample","Hugo_Symbol")]),]
  anno_ldf[["Variant_Classification"]]<-ifelse(anno_ldf[["hits"]]==1,anno_ldf[["Variant_Classification"]],"Multi_Hits")

  gene_anno_ldf<-as.data.frame(anno_ldf %>% group_by(Hugo_Symbol) %>% mutate(mut_rate=round(n_distinct(sample)/length(unique(segments_ldf[["sample"]]))*100)))
  gene_anno_ldf<-gene_anno_ldf[!duplicated(gene_anno_ldf[["Hugo_Symbol"]]),]
  gene_anno_ldf[["Label"]]<-paste(gene_anno_ldf[["Hugo_Symbol"]],paste(as.character(gene_anno_ldf[["mut_rate"]]),"%)",sep=""),sep="(")

  p<-ggplot()+geom_rect(data=segments_ldf,aes(ymin=Start,ymax=End,xmin=as.numeric(sample)-1,xmax=as.numeric(sample),fill=segmean))+
    scale_fill_gradientn(colours = c("blue","#FBFCFE","white","#FEFBFB","red"),values = c(0,0.4,0.5,0.6,1),limits=c(-1.5,1.5),oob=squish)+
    coord_cartesian(ylim=c(max(linear_genome$End),1),xlim=c(-3,max(as.numeric(segments_ldf[["sample"]]))+3),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(-3,max(as.numeric(segments_ldf[["sample"]]))+3),expand=c(0,0))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=max(as.numeric(segments_ldf[["sample"]]))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(ymin=Start,ymax=End,xmin=-1.1,xmax=-0.1),fill=linear_genome[["fill"]],alpha=0.2)+
    geom_text(data=linear_genome,aes(y=Mid,x=nudge-1.5,label=rownames(linear_genome)))+
    geom_jitter(data=anno_ldf,aes(x=as.numeric(sample)-0.5,y=Start),position=position_jitter(width=0,height=0.1),color="black",shape=4,size=1)+
    geom_text_repel(data=gene_anno_ldf,aes(y=Start,label=Label),size=anno_gene_text_size,segment.inflect=T,x=ceiling(max(as.numeric(segments_ldf[["sample"]])))+0.1,xlim=c(ceiling(max(as.numeric(segments_ldf[["sample"]])))+1,NA),force=0.9,direction = "y",hjust="left",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom")

  return(p)
}

# SOURCE: myscripts.R
#' Plot genome-wide amplification and deletion scores.
#'
#' Plot genome-wide amplification and deletion scores.
#' @param scores_gr Genomic ranges input used by the function.
#' @param type_colors Color specification used when plotting.
#' @param threshold Numeric tuning parameter used by the existing implementation.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param chromosome_colors Color specification used when plotting.
#' @param scores_anno_gr Genomic ranges input used by the function.
#' @param text_size Function argument documented from the legacy interface.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: myscripts.R.
#'
#' @examples
#' \dontrun{
#' ggscores(scores_gr = ..., type_colors = ...)
#' }
#' @export
ggscores<-function(scores_gr,type_colors=c("Amp"="red","Del"="blue"),threshold=1,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=T,xy=F,chromosome_colors=c("grey","black"),scores_anno_gr,text_size=3){
  # prepare linear genome data
  stopifnot(all(seqlengths(scores_gr)==seqlengths(scores_anno_gr)))
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  genome_=match.arg(genome)
  seqlevelsStyle_=match.arg(seqlevelsStyle)
  if((!is.null(scores_gr)) & (!any(is.na(seqlengths(scores_gr))))){
    seqinfo_=seqinfo(scores_gr)
  } else {
    if(genome=='hg19'){
      txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
    } else {
      txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
    }
    seqlevelsStyle(txdb)<-seqlevelsStyle
    seqinfo_<-seqinfo(txdb)
    seqs=seqlevels(txdb)
    if(simplified){
      seqs=seqlevels(txdb)[1:24]
      if(!xy){
        seqs=seqlevels(txdb)[1:22]
      }
      seqinfo_<-seqinfo_[seqs]
    }
  }
  linear_genome<-as.data.frame(seqinfo_)[,"seqlengths",drop=F]
  linear_genome$End<-cumsum(as.numeric(linear_genome$seqlengths))
  linear_genome$Start<-c(1,linear_genome[["End"]][-nrow(linear_genome)]+1)
  linear_genome<-linear_genome[,c("Start","End")]
  linear_genome$Mid<-(linear_genome$Start+linear_genome$End)/2
  linear_genome$nudge<-c(-0.5,-0.8)
  linear_genome$fill<-chromosome_colors

  #prepare Amp data
  Amp_scores_gr<-scores_gr[scores_gr$Type=="Amp"]
  Amp_scores_ldf<-gr2linear(Amp_scores_gr)
  Amp_scores_ldf$color<-ifelse(Amp_scores_ldf$neg_log10_q_value>threshold,type_colors["Amp"],"grey")
  #prepare cosmic genes located in sequence with significant q value with Amp linear genome
  Amp_sig_scores_gr<-Amp_scores_gr[Amp_scores_gr$neg_log10_q_value>=threshold]
  Amp_sig_scores_anno_gr<-scores_anno_gr[unique(queryHits(findOverlaps(scores_anno_gr,Amp_sig_scores_gr)))]
  Amp_sig_anno_ldf<-gr2linear(Amp_sig_scores_anno_gr)
  # plot Amp scores with cosmic genes
  Amp_p<-ggplot()+
    geom_segment(data=Amp_scores_ldf,aes(y=Start,x=0,yend=End,xend=neg_log10_q_value,color=I(color)))+
    coord_cartesian(xlim=c(0,max(Amp_scores_ldf$neg_log10_q_value)+4),ylim=c(max(linear_genome$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(0,ceiling(max(Amp_scores_ldf$neg_log10_q_value))),expand=c(0,0),breaks=c(0,1,2,5,10))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=-Inf,xmax=ceiling(max(Amp_scores_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(xmin=0,xmax=ceiling(max(Amp_scores_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Amp"])+
    geom_text_repel(data=Amp_sig_anno_ldf,aes(y=(Start+End)/2,label=Symbol),size=text_size,segment.inflect=T,x=ceiling(max(Amp_scores_ldf$neg_log10_q_value))+0.1,xlim=c(ceiling(max(Amp_scores_ldf$neg_log10_q_value))+1,NA),force=0.9,direction = "y",hjust="left",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  #prepare Del data
  Del_scores_gr<-scores_gr[scores_gr$Type=="Del"]
  Del_scores_ldf<-gr2linear(Del_scores_gr)
  Del_scores_ldf$color<-ifelse(Del_scores_ldf$neg_log10_q_value>threshold,type_colors["Del"],"grey")
  #prepare cosmic genes located in sequence with significant q value with Del linear genome
  Del_sig_scores_gr<-Del_scores_gr[Del_scores_gr$neg_log10_q_value>=threshold]
  Del_sig_scores_anno_gr<-scores_anno_gr[unique(queryHits(findOverlaps(scores_anno_gr,Del_sig_scores_gr)))]
  Del_sig_anno_ldf<-gr2linear(Del_sig_scores_anno_gr)
  #plot Del scores with cosmic gene
  Del_p<-ggplot()+
    geom_segment(data=Del_scores_ldf,aes(y=Start,x=0,yend=End,xend=neg_log10_q_value,color=I(color)))+
    coord_cartesian(xlim=c(max(Del_scores_ldf$neg_log10_q_value)+3,-1),ylim=c(max(linear_genome$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_reverse()+
    scale_x_continuous(limits=c(-1,ceiling(max(Del_scores_ldf$neg_log10_q_value))),expand=c(0,0),breaks=c(0,1,2,5,10))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=ceiling(max(Del_scores_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(xmin=0,xmax=ceiling(max(Del_scores_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Del"])+
    geom_text(data=linear_genome,aes(y=Mid,x=nudge,label=rownames(linear_genome)))+
    geom_text_repel(data=Del_sig_anno_ldf,aes(y=(Start+End)/2,label=Symbol),size=text_size,segment.inflect=T,x=ceiling(max(Del_scores_ldf$neg_log10_q_value))+0.1,xlim=c(NA,ceiling(max(Del_scores_ldf$neg_log10_q_value))+1),force=0.9,direction = "y",hjust="right",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  return(ggarrange(Del_p,Amp_p,nrow = 1,align = "h",widths = c(1.2,1),label.x="neg_log10_p_value"))
}

# SOURCE: myscripts.R
#' Visualize chromosome-arm level copy-number events.
#'
#' Visualize chromosome-arm level copy-number events.
#' @param arms Input data frame or matrix.
#' @param arm_significance Function argument documented from the legacy interface.
#' @param sample_order Function argument documented from the legacy interface.
#' @param threshold Numeric tuning parameter used by the existing implementation.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param type_colors Color specification used when plotting.
#' @param arm_colors Color specification used when plotting.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: myscripts.R.
#'
#' @examples
#' \dontrun{
#' ggarms(arms = ..., arm_significance = ...)
#' }
#' @export
ggarms<-function(arms,arm_significance,sample_order=NULL,threshold=1,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=T,xy=F,type_colors=c("Amp"="red","Del"="blue"),arm_colors=c("grey","black")){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  library(GenVisR)
  genome_=match.arg(genome)
  seqlevelsStyle_=match.arg(seqlevelsStyle)
  arms_<-cytoGeno %>% dplyr::filter(genome==genome_)
  arms_[["arm"]]<-paste(gsub("chr","",arms_[["chrom"]]),substr(arms_[["name"]],1,1),sep="")
  arms_<-as.data.frame(arms_ %>% group_by(arm) %>% summarise(chromosome=unique(gsub("chr","",chrom)),start=min(chromStart)+1,end=max(chromEnd)))
  arms_gr=sort(df2granges(arms_,genome = genome_,seqlevelsStyle = seqlevelsStyle_,simplified = simplified,xy=xy,seqnames_col = "chromosome",start_col = "start",end_col = "end",meta_cols = "arm"))
  linear_arms<-gr2linear(arms_gr)
  linear_arms$Mid<-(linear_arms$Start+linear_arms$End)/2
  linear_arms$nudge<-c(-0.5,-2)
  linear_arms$fill<-arm_colors

  arms<-pivot_longer(arms,cols=!arm,names_to = "sample",values_to = "copynumber")
  arms_ldf<-merge(arms,linear_arms,by="arm")

  if(is.null(sample_order)){
    sample_order=sort(unique(arms_ldf[["sample"]]),decreasing = T)
  }
  arms_ldf$sample<-factor(arms_ldf[["sample"]],levels=sample_order,ordered = T)

  arm_p<-ggplot()+geom_rect(data=arms_ldf,aes(ymin=Start,ymax=End,xmin=as.numeric(sample)-1,xmax=as.numeric(sample),fill=copynumber))+
    scale_fill_gradientn(colours = c("blue","#FBFCFE","white","#FEFBFB","red"),values = c(0,0.4,0.5,0.6,1),limits=c(-1.5,1.5),oob=squish)+
    coord_cartesian(ylim=c(max(linear_arms$End),1),xlim=c(-4,max(as.numeric(arms_ldf[["sample"]]))),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(-4,max(as.numeric(arms_ldf[["sample"]]))),expand=c(0,0))+
    scale_y_continuous(limits=c(1,max(linear_arms$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=max(as.numeric(arms_ldf[["sample"]]))),fill=NA,color='black')+
    geom_rect(data=linear_arms,aes(ymin=Start,ymax=End,xmin=-1.1,xmax=-0.1),fill=linear_arms[["fill"]],alpha=0.2)+
    geom_text(data=linear_arms,aes(y=Mid,x=nudge-1.5,label=arm))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          legend.position = "left")


  Amp_arm_significance=data.frame(arm=arm_significance[["Arm"]],neg_log10_q_value=-log10(arm_significance[["Amp q-value"]]+1e-20))
  Amp_arm_significance_ldf<-merge(Amp_arm_significance,linear_arms,by="arm")
  Amp_arm_significance_ldf$sig_color<-ifelse(Amp_arm_significance_ldf$neg_log10_q_value>threshold,type_colors["Amp"],"black")
  Amp_arm_significance_p<-ggplot()+
    geom_rect(data=Amp_arm_significance_ldf,aes(ymin=Start,xmin=0,ymax=End,xmax=neg_log10_q_value,fill=I(sig_color)))+
    coord_cartesian(xlim=c(0,max(Amp_arm_significance_ldf$neg_log10_q_value)+4),ylim=c(max(linear_arms$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(0,ceiling(max(Amp_arm_significance_ldf$neg_log10_q_value))),expand=c(0,0),breaks=c(0,1,2,3))+
    scale_y_continuous(limits=c(1,max(linear_arms$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=-Inf,xmax=ceiling(max(Amp_arm_significance_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_arms,aes(xmin=0,xmax=ceiling(max(Amp_arm_significance_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Amp"])+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  Del_arm_significance=data.frame(arm=arm_significance[["Arm"]],neg_log10_q_value=-log10(arm_significance[["Del q-value"]]+1e-20))
  Del_arm_significance_ldf<-merge(Del_arm_significance,linear_arms,by="arm")
  Del_arm_significance_ldf$sig_color<-ifelse(Del_arm_significance_ldf$neg_log10_q_value>threshold,type_colors["Del"],"black")
  linear_arms$nudge<-c(-0.3,-0.7)
  Del_arm_significance_p<-ggplot()+
    geom_rect(data=Del_arm_significance_ldf,aes(ymin=Start,xmin=0,ymax=End,xmax=neg_log10_q_value,fill=I(sig_color)))+
    coord_cartesian(xlim=c(max(Del_arm_significance_ldf$neg_log10_q_value)+0.1,-1),ylim=c(max(linear_arms$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_reverse()+
    scale_x_continuous(limits=c(-1,ceiling(max(Del_arm_significance_ldf$neg_log10_q_value))+0.1),expand=c(0,0),breaks=c(0,1,2,3))+
    scale_y_continuous(limits=c(1,max(linear_arms$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=ceiling(max(Del_arm_significance_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_arms,aes(xmin=0,xmax=ceiling(max(Del_arm_significance_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Del"])+
    geom_text(data=linear_arms,aes(y=Mid,x=nudge,label=arm))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  return(ggarrange(arm_p,Del_arm_significance_p,Amp_arm_significance_p,nrow = 1,align = "h",widths = c(2.5,1,1)))
}

# SOURCE: myscripts.R
# NOTE: gggenome in DNAfuncs.R is identical in logic; the myscripts.R version is used here as it
# includes the integrated scores panel (CN_p + Del_p + Amp_p).
# SOURCE: myscripts.R
#' Render multi-sample genome-wide copy-number plots.
#'
#' Render multi-sample genome-wide copy-number plots.
#' @param segments Function argument documented from the legacy interface.
#' @param sample_col Column name used by the existing implementation.
#' @param seqnames_col Column name used by the existing implementation.
#' @param start_col Column name used by the existing implementation.
#' @param end_col Column name used by the existing implementation.
#' @param meta_cols Function argument documented from the legacy interface.
#' @param sample_order Function argument documented from the legacy interface.
#' @param segments_anno_df Function argument documented from the legacy interface.
#' @param anno_sample_col Column name used by the existing implementation.
#' @param anno_seqnames_col Column name used by the existing implementation.
#' @param anno_start_col Column name used by the existing implementation.
#' @param anno_end_col Column name used by the existing implementation.
#' @param anno_meta_cols Function argument documented from the legacy interface.
#' @param anno_gene_text_size Function argument documented from the legacy interface.
#' @param scores_df Function argument documented from the legacy interface.
#' @param score_seqnames_col Column name used by the existing implementation.
#' @param score_start_col Column name used by the existing implementation.
#' @param score_end_col Column name used by the existing implementation.
#' @param score_meta_cols Function argument documented from the legacy interface.
#' @param type_colors Color specification used when plotting.
#' @param threshold Numeric tuning parameter used by the existing implementation.
#' @param scores_anno_df Function argument documented from the legacy interface.
#' @param score_anno_seqnames_col Column name used by the existing implementation.
#' @param score_anno_start_col Column name used by the existing implementation.
#' @param score_anno_end_col Column name used by the existing implementation.
#' @param score_anno_meta_cols Function argument documented from the legacy interface.
#' @param anno_score_text_size Function argument documented from the legacy interface.
#' @param genome Option controlling how the function runs.
#' @param seqlevelsStyle Option controlling how the function runs.
#' @param simplified Logical flag controlling optional behavior.
#' @param xy Logical flag controlling optional behavior.
#' @param chromosome_colors Color specification used when plotting.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: myscripts.R.
#'
#' @examples
#' \dontrun{
#' gggenome(segments = ..., sample_col = ...)
#' }
#' @export
gggenome<-function(segments,sample_col="sample",seqnames_col = "chromosome",start_col = "start",end_col = "end",meta_cols = c("sample","probes","segmean"),sample_order=NULL,
                   segments_anno_df=NULL,anno_sample_col="Tumor_Sample_Barcode",anno_seqnames_col = "Chromosome",anno_start_col = "Start_position",anno_end_col = "End_Position",anno_meta_cols = c("Hugo_Symbol","Variant_Classification","Variant_Type"),anno_gene_text_size=2,
                   scores_df=NULL,score_seqnames_col = "Chromosome",score_start_col = "Start",score_end_col = "End",score_meta_cols = c("Type","neg_log10_q_value"),type_colors=c("Amp"="red","Del"="blue"),threshold=1,
                   scores_anno_df=NULL,score_anno_seqnames_col = "Chromosome",score_anno_start_col = "Start",score_anno_end_col = "End",score_anno_meta_cols = "Symbol",anno_score_text_size=2,
                   genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=T,xy=F,chromosome_colors=c("grey","black")){
  # prepare linear genome data
  #stopifnot(all(seqlengths(scores_gr)==seqlengths(scores_anno_gr)))
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  genome_=match.arg(genome)
  seqlevelsStyle_=match.arg(seqlevelsStyle)
  if(genome=='hg19'){
    txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  seqlevelsStyle(txdb)<-seqlevelsStyle
  seqinfo_<-seqinfo(txdb)
  seqs=seqlevels(txdb)
  if(simplified){
    seqs=seqlevels(txdb)[1:24]
    if(!xy){
      seqs=seqlevels(txdb)[1:22]
    }
    seqinfo_<-seqinfo_[seqs]
  }

  linear_genome<-as.data.frame(seqinfo_)[,"seqlengths",drop=F]
  linear_genome$End<-cumsum(as.numeric(linear_genome$seqlengths))
  linear_genome$Start<-c(1,linear_genome[["End"]][-nrow(linear_genome)]+1)
  linear_genome<-linear_genome[,c("Start","End")]
  linear_genome$Mid<-(linear_genome$Start+linear_genome$End)/2
  linear_genome$nudge<-c(-0.5,-0.8)
  linear_genome$fill<-chromosome_colors

  segments_l<-lapply(unique(segments[[sample_col]]),function(samp){
    seg<-segments[segments[[sample_col]]==samp,]
    seg_gr<-df2granges(seg,genome=genome,seqlevelsStyle = seqlevelsStyle,seqnames_col = seqnames_col,start_col = start_col,end_col = end_col,meta_cols = meta_cols)
    seg_ldf<-gr2linear(seg_gr)
    return(seg_ldf)
  })
  segments_ldf<-do.call(rbind,segments_l)
  colnames(segments_ldf)[colnames(segments_ldf)=="segmean"]="copynumber"
  if(is.null(sample_order)){
    sample_order=sort(unique(segments_ldf[[sample_col]]),decreasing = T)
  }
  if((!is.factor(segments[[sample_col]])) | (!is.ordered(segments[[sample_col]]))){
    segments[[sample_col]]<-factor(segments[[sample_col]],levels=sort(unique(segments[[sample_col]]),decreasing = T),ordered=T)
  }
  segments_ldf$sample<-factor(segments_ldf[[sample_col]],levels=sample_order,ordered = T)

  if(!is.null(segments_anno_df)){
    segments_anno_ldf<-lapply(unique(segments_anno_df[[anno_sample_col]]),function(samp){data=segments_anno_df[segments_anno_df[[anno_sample_col]]==samp,];
    gr=df2granges(df=data,genome = genome,seqlevelsStyle = seqlevelsStyle,seqnames_col = anno_seqnames_col,start_col = anno_start_col,end_col = anno_end_col,meta_cols = anno_meta_cols);
    ldf<-gr2linear(gr)
    ldf[[anno_sample_col]]<-samp
    return(ldf)})
    segments_anno_ldf<-do.call(rbind,segments_anno_ldf)


  segments_anno_ldf$sample<-factor(segments_anno_ldf[[anno_sample_col]],levels=sample_order,ordered = T)
  segments_anno_ldf<-as.data.frame(segments_anno_ldf %>% group_by(sample,Hugo_Symbol) %>% mutate(hits=n()))
  segments_anno_ldf<-segments_anno_ldf[!duplicated(segments_anno_ldf[,c("sample","Hugo_Symbol")]),]
  segments_anno_ldf[["Variant_Classification"]]<-ifelse(segments_anno_ldf[["hits"]]==1,segments_anno_ldf[["Variant_Classification"]],"Multi_Hits")

  gene_segments_anno_ldf<-as.data.frame(segments_anno_ldf %>% group_by(Hugo_Symbol) %>% mutate(mut_rate=round(n_distinct(sample)/length(unique(segments_ldf[["sample"]]))*100)))
  gene_segments_anno_ldf<-gene_segments_anno_ldf[!duplicated(gene_segments_anno_ldf[["Hugo_Symbol"]]),]
  gene_segments_anno_ldf[["Label"]]<-paste(gene_segments_anno_ldf[["Hugo_Symbol"]],paste(as.character(gene_segments_anno_ldf[["mut_rate"]]),"%)",sep=""),sep="(")
  gene_segments_anno_ldf$sample<-factor(gene_segments_anno_ldf$sample,levels=sample_order,ordered = T)

  CN_p<-ggplot()+geom_rect(data=segments_ldf,aes(ymin=Start,ymax=End,xmin=as.numeric(sample)-1,xmax=as.numeric(sample),fill=copynumber))+
    scale_fill_gradientn(colours = c("blue","#FBFCFE","white","#FEFBFB","red"),values = c(0,0.4,0.5,0.6,1),limits=c(-1.5,1.5),oob=squish)+
    coord_cartesian(ylim=c(max(linear_genome$End),1),xlim=c(-3,max(as.numeric(segments_ldf[["sample"]]))+4),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(-3,max(as.numeric(segments_ldf[["sample"]]))+4),expand=c(0,0))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=max(as.numeric(segments_ldf[["sample"]]))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(ymin=Start,ymax=End,xmin=-1.1,xmax=-0.1),fill=linear_genome[["fill"]],alpha=0.2)+
    geom_text(data=linear_genome,aes(y=Mid,x=nudge-1.5,label=rownames(linear_genome)))+
    geom_jitter(data=segments_anno_ldf,aes(x=as.numeric(sample)-0.5,y=Start),position=position_jitter(width=0,height=0.1),color="black",shape=4,size=1)+
    geom_text_repel(data=gene_segments_anno_ldf,aes(y=Start,label=Label),size=anno_gene_text_size,segment.inflect=T,x=ceiling(max(as.numeric(segments_ldf[["sample"]])))+0.1,xlim=c(ceiling(max(as.numeric(segments_ldf[["sample"]])))+1,NA),force=0.9,direction = "y",hjust="left",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank(),
          legend.position = "left")

  } else {
    CN_p<-ggplot()+geom_rect(data=segments_ldf,aes(ymin=Start,ymax=End,xmin=as.numeric(sample)-1,xmax=as.numeric(sample),fill=copynumber))+
      scale_fill_gradientn(colours = c("blue","#FBFCFE","white","#FEFBFB","red"),values = c(0,0.4,0.5,0.6,1),limits=c(-1.5,1.5),oob=squish)+
      coord_cartesian(ylim=c(max(linear_genome$End),1),xlim=c(-3,max(as.numeric(segments_ldf[["sample"]]))),clip="on")+
      scale_y_reverse()+
      scale_x_continuous(limits=c(-3,max(as.numeric(segments_ldf[["sample"]]))),expand=c(0,0))+
      scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
      geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=max(as.numeric(segments_ldf[["sample"]]))),fill=NA,color='black')+
      geom_rect(data=linear_genome,aes(ymin=Start,ymax=End,xmin=-1.1,xmax=-0.1),fill=linear_genome[["fill"]],alpha=0.2)+
      geom_text(data=linear_genome,aes(y=Mid,x=nudge-1.5,label=rownames(linear_genome)))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x.bottom = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y = element_blank(),
            legend.position = "left")
  }




  #prepare Amp data
  colnames(scores_df)[5]<-"neg_log10_q_value"
  scores_gr<-df2granges(scores_df,genome=genome,seqlevelsStyle = seqlevelsStyle,seqnames_col = score_seqnames_col,start_col = score_start_col,end_col = score_end_col,meta_cols = score_meta_cols)
  scores_anno_gr<-df2granges(scores_anno_df,genome=genome,seqlevelsStyle = seqlevelsStyle,seqnames_col = score_anno_seqnames_col,start_col = score_anno_start_col,end_col = score_anno_end_col,meta_cols = score_anno_meta_cols)

  Amp_scores_gr<-scores_gr[scores_gr$Type=="Amp"]
  Amp_scores_ldf<-gr2linear(Amp_scores_gr)
  Amp_scores_ldf$color<-ifelse(Amp_scores_ldf$neg_log10_q_value>threshold,type_colors["Amp"],"grey")
  #prepare cosmic genes located in sequence with significant q value with Amp linear genome
  Amp_sig_scores_gr<-Amp_scores_gr[Amp_scores_gr$neg_log10_q_value>=threshold]
  Amp_sig_scores_anno_gr<-scores_anno_gr[unique(queryHits(findOverlaps(scores_anno_gr,Amp_sig_scores_gr)))]
  Amp_sig_anno_ldf<-gr2linear(Amp_sig_scores_anno_gr)
  # plot Amp scores with cosmic genes
  Amp_p<-ggplot()+
    geom_segment(data=Amp_scores_ldf,aes(y=Start,x=0,yend=End,xend=neg_log10_q_value,color=I(color)))+
    coord_cartesian(xlim=c(0,max(Amp_scores_ldf$neg_log10_q_value)+4),ylim=c(max(linear_genome$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_continuous(limits=c(0,ceiling(max(Amp_scores_ldf$neg_log10_q_value))),expand=c(0,0),breaks=c(0,1,2,5,10))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=-Inf,xmax=ceiling(max(Amp_scores_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(xmin=0,xmax=ceiling(max(Amp_scores_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Amp"])+
    geom_text_repel(data=Amp_sig_anno_ldf,aes(y=(Start+End)/2,label=Symbol),size=anno_score_text_size,segment.inflect=T,x=ceiling(max(Amp_scores_ldf$neg_log10_q_value))+0.1,xlim=c(ceiling(max(Amp_scores_ldf$neg_log10_q_value))+1,NA),force=0.9,direction = "y",hjust="left",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  #prepare Del data
  Del_scores_gr<-scores_gr[scores_gr$Type=="Del"]
  Del_scores_ldf<-gr2linear(Del_scores_gr)
  Del_scores_ldf$color<-ifelse(Del_scores_ldf$neg_log10_q_value>threshold,type_colors["Del"],"grey")
  #prepare cosmic genes located in sequence with significant q value with Del linear genome
  Del_sig_scores_gr<-Del_scores_gr[Del_scores_gr$neg_log10_q_value>=threshold]
  Del_sig_scores_anno_gr<-scores_anno_gr[unique(queryHits(findOverlaps(scores_anno_gr,Del_sig_scores_gr)))]
  Del_sig_anno_ldf<-gr2linear(Del_sig_scores_anno_gr)
  #plot Del scores with cosmic gene
  Del_p<-ggplot()+
    geom_segment(data=Del_scores_ldf,aes(y=Start,x=0,yend=End,xend=neg_log10_q_value,color=I(color)))+
    coord_cartesian(xlim=c(max(Del_scores_ldf$neg_log10_q_value)+3,-1),ylim=c(max(linear_genome$End),1),clip="on")+
    scale_y_reverse()+
    scale_x_reverse()+
    scale_x_continuous(limits=c(-1,ceiling(max(Del_scores_ldf$neg_log10_q_value))),expand=c(0,0),breaks=c(0,1,2,5,10))+
    scale_y_continuous(limits=c(1,max(linear_genome$End)),expand=c(0,0))+
    geom_rect(data=NULL,aes(ymin=-Inf,ymax=Inf,xmin=0,xmax=ceiling(max(Del_scores_ldf$neg_log10_q_value))),fill=NA,color='black')+
    geom_rect(data=linear_genome,aes(xmin=0,xmax=ceiling(max(Del_scores_ldf$neg_log10_q_value)),ymin=Start,ymax=End,fill=I(fill)),alpha=0.2)+
    geom_vline(xintercept=1,linetype="dashed",color=type_colors["Del"])+
    geom_text(data=linear_genome,aes(y=Mid,x=nudge,label=rownames(linear_genome)))+
    geom_text_repel(data=Del_sig_anno_ldf,aes(y=(Start+End)/2,label=Symbol),size=anno_score_text_size,segment.inflect=T,x=ceiling(max(Del_scores_ldf$neg_log10_q_value))+0.1,xlim=c(NA,ceiling(max(Del_scores_ldf$neg_log10_q_value))+1),force=0.9,direction = "y",hjust="right",segment.size =0.2,segment.curvature = 0,segment.ncp=3,segment.angle=20,max.iter = 100000,max.overlaps = Inf)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y = element_blank())

  return(ggarrange(CN_p,Del_p,Amp_p,nrow = 1,align = "h",widths = c(1.5,1.2,1)))
}
