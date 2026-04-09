# =============================================================================
# MUTATION ANALYSIS FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   get_combined_tcga_mutation_table    - Combine tumor-specific TCGA mutations with background
#   filter_mutects                      - SNV filtering pipeline (13+ filters)
#   filter_pindels                      - Indel filtering pipeline (16+ filters)
#   convert_mutations                   - Convert mutation format (mutsig2cv/maftools/pynbs)
#   complex_oncoplot                    - Complex oncoprint visualization
#   annotate_genomic_mutations_with_panel_info - Annotate mutations with panel intervals
#   annotate_mdl_mutations_with_timepoints     - Add timepoint info to mutations
#   annotate_mdl_mutations_with_mutationdb     - Annotate with OncoKB/COSMIC
#   get_min_multihit_distance           - Minimum distance between multi-hit codons
#   harmonize_mutations                 - Harmonize mutations across assay types
#   oncogenicpathway_tableplot          - Pathway-level mutation table visualization
#   cancerhallmark_heatmap              - Cancer hallmark gene heatmap
#
# SOURCE PROVENANCE:
#   woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mutations.R (2026-03-02)
#   This is the CANONICAL version (most recent). IBC/DNAfuncs.R (2024-03-18) has
#   an older filter_mutects/filter_pindels with fewer filter options.
#
# DUPLICATE NOTE: filter_mutects and filter_pindels also exist in:
#   - Pembro/funcsInPembro.R (2025-07-02) — same API as here
#   - IBC/DNAfuncs.R (2024-03-18) — older, missing treat_noisy_samples param
#   - myscripts.R (2023-06-29) — older version
#   RECOMMEND: Use this file (2026-03-02) as it is most recent.
# =============================================================================

#' get_combined_tcga_mutation_table
#' combine tumor specific tcga mutation simplified table with pre-existing non-tumor specific tcga mutation simplified table
#'
#' @param tcga_tumor_abbreviation character. abbreviation for tcga tumor.
#' @param tcga_studies data.frame. TCGA study list.
#' @param tcga_snp_table data.frame. simplified TCGA snp frequency table.
#' @param tcga_indel_table data.frame. simplified TCGA ins/del frequency table.
#'
#' @return list. contain combined simplified tcga snp/indel frequency table
#' @export
#'
get_combined_tcga_mutation_table<-function(tcga_tumor_abbreviation="",tcga_studies,tcga_snp_table,tcga_indel_table){
  stopifnot("tumor type is not in tcga tumor type list!"=tcga_tumor_abbreviation %in% tcga_studies$Study_Abbreviation)
  tcga_project<-paste("TCGA",tcga_tumor_abbreviation,sep="-")
  query <- TCGAbiolinks::GDCquery(
    project = tcga_project,
    data.category = "Simple Nucleotide Variation",
    access = "open",
    legacy = FALSE,
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  TCGAbiolinks::GDCdownload(query)
  mutations <- TCGAbiolinks::GDCprepare(query)
  n_sample<-length(unique(mutations$Tumor_Sample_Barcode))
  subset_cosmic_tcga_snps<-mutations[mutations$Variant_Type=="SNP",]
  simplified_subset_cosmic_tcga_snps<-subset_cosmic_tcga_snps %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarize(n_mutect=dplyr::n())
  simplified_subset_cosmic_tcga_snps<-data.frame(gene=simplified_subset_cosmic_tcga_snps$Hugo_Symbol,n_mutect=simplified_subset_cosmic_tcga_snps$n_mutect,n_sample=n_sample)
  combined_tcga_snp_table<-rbind(tcga_snp_table[!(tcga_snp_table$gene %in% simplified_subset_cosmic_tcga_snps$gene),],simplified_subset_cosmic_tcga_snps)
  assertthat::are_equal(length(unique(combined_tcga_snp_table$gene)),nrow(combined_tcga_snp_table))


  subset_cosmic_tcga_indels<-mutations[mutations$Variant_Type %in% c("INS","DEL"),]
  simplified_subset_cosmic_tcga_indels<-subset_cosmic_tcga_indels %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarize(n_pindel=dplyr::n())
  simplified_subset_cosmic_tcga_indels<-data.frame(gene=simplified_subset_cosmic_tcga_indels$Hugo_Symbol,n_pindel=simplified_subset_cosmic_tcga_indels$n_pindel,n_sample=n_sample)
  combined_tcga_indel_table<-rbind(tcga_indel_table[!(tcga_indel_table$gene %in% simplified_subset_cosmic_tcga_indels$gene),],simplified_subset_cosmic_tcga_indels)
  assertthat::are_equal(length(unique(combined_tcga_indel_table$gene)),nrow(combined_tcga_indel_table))
  return(list(combined_tcga_snp_table=combined_tcga_snp_table,combined_tcga_indel_table=combined_tcga_indel_table))
}


#' filter_mutects
#'
#' filter mutects using various filter sets
#'
#' @param mutects  dataframe. containing snv annotation.

#' @param traced_columns  character. default c("Sample_ID","gene","type","chr","start","end"),used for removing dupicated snv annotations.

#' @param low_purity_filter  logic. default TRUE, remove low purity sample or keep.
#' @param low_purity_samples character. list of samples with low purity.

#' @param treat_noisy_samples logic. default TRUE, remove specially treat noisy samples (remove all mutects not present in clean/rest samples).
#' @param noisy_samples character. list of samples with too much noise.

#' @param exonutronly_filter logic. default TRUE, keep only snv located in exon, utr or not.
#' @param keep character. default c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"), list snv locus (column func.knowngene) to keep.

#' @param off_target_filter logic. default TRUE, remove off-target snv or not.
#' @param intervals GenomicRanges.providing the targeted regions information.

#' @param germline_filter logic. default TRUE, remove germline snv or not.
#' @param cutoff_n_vaf numeric. defalut 0.01, cutoff for normal variant allele fraction.

#' @param treat_known_cancer_gene_specially logic. defalut TRUE, split known cancer genes and treat differently.
#' @param known_cancer_genes character. list of cancer gene names.
#' @param special_mutect_genes character. list of genes whose snv not through filter sets.

#' @param mapping_quality_filter logic. default TRUE, remove low mapping quality snv or not.
#' @param cutoff_mq numeric. default 30, cuttoff for maximal mapping quality of reads.

#' @param statistical_filter logic.default TRUE, remove snv not passing proportional test or not.
#' @param cutoff_prop_test_p numeric. default 0.01, cutoff for proportial test pvalue.

#' @param tumor_f_filter logic. default TRUE, remove snv with low tumor fraction or not.
#' @param cutoff_tumor_f numeric. default 0.0, cufoff for tumor fraction.

#' @param tumor_coverage_filter logic. default TRUE, remove snv with low tumor coverage or not.
#' @param cutoff_t_coverage numeric. default 10, cutoff for tumor coverage (column t_alt+ column t_ref).

#' @param normal_coverage_filter logic. default TRUE, remove snv with low normal coverage or not.
#' @param cutoff_n_coverage numeric. default 10, cutoff for normal coverage (column n_alt+ column n_ref).

#' @param lodt_filter logic. default TRUE, remove snv with low log likelyhood of tumor event or not.
#' @param cutoff_lodt numeric. default 6.3, cutoff for number of log likelyhood of tumor event.

#' @param consistent_mutect_statistic_filter logic. default TRUE, remove snv failed to pass consistency statistical test (binom test).
#' @param mutect_for_consistent_statistic data.frame. containing gene,number of snv, number of samples (columns: gene,n_pindel,n_sample), used for binom test with targeted snv.
#' @param cutoff_binom_pval numeric. default 0.05, cutoff for binom test pvalue.

#' @param description_file  character. default description.txt, file name for deposit of information of each operation (filter sets).
#' @param save_traced_mutect logic. default TRUE, save traced_mutect_file.
#' @param traced_mutect_file character. default traced_pindel.txt, file name for deposit results after each operation (filter sets) if save_traced_mutect=TRUE.
#'
#' @return data.frame. filtered mutects(snp).
#'


filter_mutects<-function(mutects, traced_columns=c("Sample_ID","gene","type","chr","start","end"),
                         low_purity_filter=T,low_purity_samples,
                         treat_noisy_samples=F,noisy_samples=c(),
                         exonutronly_filter=T,keep=c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),
                         off_target_filter=T,intervals=intervals,
                         germline_filter=T,cutoff_n_vaf=0.01,
                         treat_known_cancer_gene_specially=T,known_cancer_genes=c("TP53","GNA11","GNAQ"),special_mutect_genes=c(),
                         mapping_quality_filter=T,cutoff_mq=30,
                         statistical_filter=T,cutoff_prop_test_p=0.01,
                         tumor_f_filter=T,cutoff_tumor_f=0,
                         tumor_coverage_filter=T,cutoff_t_coverage=10,
                         normal_coverage_filter=T,cutoff_n_coverage=10,
                         lodt_filter=T,cutoff_lodt=6.3,
                         consistent_mutect_statistic_filter=T,mutect_for_consistent_statistic=NA,cutoff_binom_pval=0.05,
                         description_file="description.txt",
                         save_traced_mutect=TRUE,traced_mutect_file="traced_mutect.txt"){
  options(warn=-1)
  mutects<-mutects[!duplicated(do.call(paste,mutects[,traced_columns])),]
  fileConn<-description_file
  write(paste(paste("Starting from ",nrow(mutects),sep="")," mutects.",sep=""),fileConn,append=F)
  mutects_traced<-list()
  traced_names<-c()
  mutects_traced[["original_mutects"]]<-mutects[,traced_columns]
  traced_names<-c(traced_names,"original_mutects")

  #low purity filter
  if(low_purity_filter){
    mutects<-mutects[!mutects[["Sample_ID"]] %in% low_purity_samples,]
    write("After filter out low purity (0.5);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_low_purity"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_low_purity")
  }
  #special treatment of very noisy samples (too many calls)
  if(treat_noisy_samples){
    clean_mutects<-mutects[!mutects[["Sample_ID"]] %in% noisy_samples,]
    noisy_mutects<-mutects[mutects[["Sample_ID"]] %in% noisy_samples,]
    denoisy_mutects<-noisy_mutects[noisy_mutects[["gene"]] %in% unique(clean_mutects[['gene']]),]
    denoisy_mutects <- denoisy_mutects %>% dplyr::group_by(Sample_ID,gene) %>% dplyr::mutate(hits=dplyr::n(),max_t_alt_count=max(t_alt_count))
    denoisy_mutects <- as.data.frame(denoisy_mutects)
    denoisy_mutects<-denoisy_mutects[denoisy_mutects$hits==1 | (denoisy_mutects$hits>1 & denoisy_mutects$t_alt_count==denoisy_mutects$max_t_alt_count),]
    denoisy_mutects<-denoisy_mutects[,colnames(clean_mutects)]
    mutects<-rbind(denoisy_mutects,clean_mutects)
    write("After specially treatment of noisy samples;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["treat noisy samples"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"treat_noisy_samples")
  }
  #exon UTRs only

  if(exonutronly_filter){
    mutects<-mutects[mutects$func.knowngene %in% c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),]
    write("After filter out non exonic or UTRs", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_nonexonutr"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #offtarget filter
  if(off_target_filter){

    mutects_gr<-df2granges(mutects,genome = "hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,seqnames_col = traced_columns[4],start_col = traced_columns[5],end_col =traced_columns[6],meta_cols = colnames(mutects))
    GenomeInfoDb::seqlevels(intervals)<-GenomeInfoDb::seqlevels(mutects_gr)
    GenomeInfoDb::seqinfo(intervals)<-GenomeInfoDb::seqinfo(mutects_gr)
    targets<-plyranges::filter_by_overlaps(mutects_gr,intervals)
    mutects<-as.data.frame(S4Vectors::mcols(targets))
    write("After filter out off-targets mutects;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_off_target"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_off_target")
  }

  #germline cutoff 0.01
  if(germline_filter){
    mutects<-mutects[(mutects$n_alt_count/(mutects$n_alt_count+mutects$n_ref_count))<=cutoff_n_vaf,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out germline mutects with cutoff_n_vaf ",cutoff_n_vaf,sep=" "),fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_germline"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_germline")
  }

  #treat known cancer gene specially
  if(treat_known_cancer_gene_specially){
    if(length(known_cancer_genes)>0){
      known_cancer_gene_mutects<-mutects[mutects$gene %in% setdiff(known_cancer_genes,special_mutect_genes),]
    }
    if(length(special_mutect_genes)>0){
      special_cancer_gene_mutects<-mutects[mutects$gene %in% special_mutect_genes,,drop=F]
    }

    mutects<-mutects[!(mutects$gene %in% known_cancer_genes),]
    write("After split known cancer gene (cosmic or cancer_genetic);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," non known cancer gene mutects.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(known_cancer_gene_mutects),sep="")," known cancer gene mutects.",sep=""), fileConn,append = T)
    mutects_traced[["non_known_cancer_gene_mutect"]]<-mutects[,traced_columns]
    mutects_traced[["known_cancer_gene_mutect"]]<-known_cancer_gene_mutects[,traced_columns]
    traced_names<-c(traced_names,"non_known_cancer_gene_mutect")
    traced_names<-c(traced_names,"known_cancer_gene_mutect")
    if(length(special_mutect_genes)>0){
      write(paste(paste("There are still ",nrow(special_cancer_gene_mutects),sep="")," special cancer gene mutects.",sep=""), fileConn,append = T)
      mutects_traced[["special_cancer_gene_mutect"]]<-special_cancer_gene_mutects[,traced_columns]
      traced_names<-c(traced_names,"special_cancer_gene_mutect")
    }
  }

  #mapping quality filter
  if(mapping_quality_filter){
    mutects<-mutects[(mutects$t_ref_max_mapq>=cutoff_mq & mutects$t_alt_max_mapq>=cutoff_mq),]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out mutects with low average mapping quality with cutoff_mq ",cutoff_mq,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_mapping_quality"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_mapping_quality")
  }

  #statistical filter
  #browser()
  if(statistical_filter){
    pval<-apply(mutects[,c("t_ref_count","t_alt_count",'n_ref_count',"n_alt_count")],1,function(row){stats::prop.test(x=c(as.integer(row["t_alt_count"]),as.integer(row["n_alt_count"])),n=c(as.integer(row["t_alt_count"])+as.integer(row["t_ref_count"])+1,as.integer(row["n_alt_count"])+as.integer(row["n_ref_count"])+1))$p.value})
    mutects<-mutects[pval<=cutoff_prop_test_p,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out prop test failed mutects with cutoff_prop_test_p ",cutoff_prop_test_p,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_prop_statistic"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_prop_statistic")
  }

  #tumor_f filter
  if(tumor_f_filter){
    mutects<-mutects[mutects$tumor_f>=cutoff_tumor_f,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low tumor fraction mutects with  cutoff_tumor_f ",cutoff_tumor_f,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_tumor_f"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_f")
  }

  #tumor depth coverage filter
  if(tumor_coverage_filter){
    mutects<-mutects[(mutects$t_alt_count+mutects$t_ref_count)>=cutoff_t_coverage,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low tumor coverage mutects with cutoff_t_coverage ",cutoff_t_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_tumor_coverage"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_coverage")
  }

  #normal depth coverage filter
  if(normal_coverage_filter){
    mutects<-mutects[(mutects$n_alt_count+mutects$n_ref_count)>=cutoff_n_coverage,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low normal coverage mutects with cutoff_n_coverage ",cutoff_n_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_normal_coverage"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_normal_coverage")
  }

  #lodt filter
  if(lodt_filter){
    mutects<-mutects[as.numeric(mutects$t_lod_fstar)>=cutoff_lodt,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out low t_lod_fstar mutects with cutoff_lodt ",cutoff_lodt,sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["filter_lodt"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_lodt")
  }

  #combine known_cancer gene mutects
  if(treat_known_cancer_gene_specially & length(known_cancer_genes)>0){
    mutects<-rbind(mutects,known_cancer_gene_mutects)
    write("After combine known cancer gene mutects;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["combined_known_mutects"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"combined_known_mutects")
  }

  #browser()
  if(consistent_mutect_statistic_filter){
    n_sample_subject<-length(unique(mutects[["Sample_ID"]]))
    mutects_subject<-mutects %>% dplyr::group_by(gene) %>% dplyr::summarise(n_mutect_subject=dplyr::n())
    mutects_subject[["n_sample_subject"]]<-n_sample_subject

    mutects_for_statistic<-merge(mutects_subject,mutect_for_consistent_statistic,by=c("gene"),all.x=T)
    mutects_for_statistic[["n_sample"]][is.na(mutects_for_statistic[["n_sample"]])]=max(mutect_for_consistent_statistic$n_sample)
    mutects_for_statistic[["n_mutect"]][is.na(mutects_for_statistic[["n_mutect"]])]<-0
    mutects_for_statistic<-set_column_as_rownames(mutects_for_statistic,"gene")
    mutects_for_statistic$n_mutect_subject<-apply(mutects_for_statistic[,c("n_mutect_subject","n_sample_subject")],1,min)
    mutects_for_statistic$n_mutect<-apply(mutects_for_statistic[,c("n_mutect","n_sample")],1,min)
    pvals<-apply(mutects_for_statistic,1,function(r){test<-stats::binom.test(r["n_mutect_subject"],r["n_sample_subject"],p=r["n_mutect"]/r["n_sample"],alternative="greater");return(test$p.value)})
    adjpvals<-stats::p.adjust(pvals,method="BH")
    keep_genes<-names(adjpvals)[adjpvals>cutoff_binom_pval]
    mutects<-mutects[mutects$gene %in% keep_genes,]
    mutects<-mutects[rownames(mutects)!="NA",]
    write(paste("After filter out non consitent mutect with tcga mutect with cutoff_binom_pval ",cutoff_binom_pval,sep=" "), fileConn,append = T)
    write(paste(paste("There are still totally ",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    if(treat_known_cancer_gene_specially){
      write(paste(paste("There are still ",sum(mutects$gene %in% known_cancer_gene_mutects$gene),sep="")," known cancer gene mutects.",sep=""), fileConn,append = T)
      write(paste(paste("There are ",nrow(known_cancer_gene_mutects)-sum(mutects$gene %in% known_cancer_genes),sep="")," known cancer gene mutects removed.",sep=""), fileConn,append = T)
    }
    mutects_traced[["filter_nonconsistent_mutect_statistic"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"filter_nonconsistent_mutect_statistic")
  }

  #combine special cancer gene mutects
  if(treat_known_cancer_gene_specially & length(special_mutect_genes)>0){
    mutects<-rbind(mutects,special_cancer_gene_mutects)
    write("After combine special cancer gene mutects;", fileConn,append = T)
    write(paste(paste("There are still totally",nrow(mutects),sep="")," mutects.",sep=""), fileConn,append = T)
    mutects_traced[["combined_final_mutects"]]<-mutects[,traced_columns]
    traced_names<-c(traced_names,"combined_final_mutects")
  }



  if(save_traced_mutect){
    mutect_order=NA
    for(traced_name in rev(traced_names)){
      if(traced_name %in% c("special_cancer_gene_mutect","known_cancer_gene_mutect")){
        next
      } else{
        mutects_temp=mutects_traced[[traced_name]]
        rownames(mutects_temp)<-do.call(paste,mutects_temp)
        mutects_temp<-mutects_temp[c(intersect(mutect_order,rownames(mutects_temp)),setdiff(rownames(mutects_temp),mutect_order)),]
        mutect_order=rownames(mutects_temp)
        rownames(mutects_temp)<-NULL
        mutects_traced[[traced_name]]=mutects_temp
      }
    }
    for(traced_name in traced_names){
      mutects_temp<-mutects_traced[[traced_name]]
      colnames(mutects_temp)<-paste(traced_name,colnames(mutects_temp),sep=":")
      utils::write.table(mutects_temp,file = traced_mutect_file,append = T,row.names = F,sep="\t")
    }
  }
  return(mutects)
}

#' filter_pindels
#'
#' filter pindels using various filter sets
#'
#' @param pindels dataframe. containing insertion and deletion annotation.

#' @param traced_columns character. default c("Sample_ID","gene","type","chr","start","end"),used for removing dupicated pindel annotations.

#' @param low_purity_filter logic. default TRUE, remove low purity sample or keep.
#' @param low_purity_samples character. list of samples with low purity.

#' @param treat_noisy_samples logic. default TRUE, remove specially treat noisy samples (remove all pindels not present in clean/rest samples).
#' @param noisy_samples character. list of samples with too much noise.

#' @param exonutronly_filter logic. default TRUE, keep only pindels located in exon, utr or not.
#' @param keep character. default c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"), list pindel locus (column func.knowngene) to keep.

#' @param remove_td_filter logic. default TRUE, remove tandem duplicated or not.

#' @param td_coverage_filter logic. default TRUE, remove low coverage tandom duplicate.
#' @param cutoff_td_coverage numeric. default 1.5, cutoff for ratio of tandem duplicate to average genome coverage.

#' @param off_target_filter logic. default TRUE, remove off-target pindels or not.
#' @param intervals GenomicRanges.providing the targeted regions information.

#' @param germline_filter logic. default TRUE, remove germline pindels or not.
#' @param cutoff_n_vaf numeric. defalut 0.01, cutoff for normal variant allele fraction.

#' @param treat_known_cancer_gene_specially logic. defalut TRUE, split known cancer genes and treat differently.
#' @param known_cancer_genes character. list of cancer gene names.
#' @param special_pindel_genes character. list of genes whose pindel no through filter sets.

#' @param mapping_quality_filter logic. default TRUE, remove low mapping quality pindels or not.
#' @param cutoff_mq numeric. default 25, cuttoff for averaged mapping quality of anchor reads.

#' @param statistical_filter logic.default TRUE, remove pindels not passing proportional test or not.
#' @param cutoff_prop_test_p numeric. default 0.01, cutoff for proportial test pvalue.

#' @param t_vaf_filter logic. default TRUE, remove pindels with low tumor variant allele fraction or not.
#' @param cutoff_t_vaf numeric. default 0.0, cutoff for tumor variant allele fraction.

#' @param tumor_f_filter logic. default TRUE, remove pindels with low tumor fraction or not.
#' @param cutoff_tumor_f numeric. default 0.0, cufoff for tumor fraction.

#' @param deletion_length_filter logic. defalut TRUE, remove too long deletions or not.
#' @param cutoff_d_length numeric. default 1000, cutoff for deletion length.

#' @param insertion_length_filter logic. defalut TRUE, remove too short insertions or not.
#' @param cutoff_i_length numeric. default 6, cutoff for insertiontion length.

#' @param tumor_coverage_filter logic. default TRUE, remove pindels with low tumor coverage or not.
#' @param cutoff_t_coverage numeric. default 10, cutoff for tumor coverage (column t_alt+ column t_ref).

#' @param normal_coverage_filter logic. default TRUE, remove pindels with low normal coverage or not.
#' @param cutoff_n_coverage numeric. default 10, cutoff for normal coverage (column n_alt+ column n_ref).

#' @param support_filter logic default TRUE, remove pindels with low support reads or not.
#' @param cutoff_support numeric. default 3, cutoff for number of support reads.

#' @param homopolymer_filter logic. default TRUE, remove pindels with too many homopolymer or not.
#' @param maxn_homopolymer numeric. default 6, allowed maximal numer of homopolymer.

#' @param repeat_pindel_filter logic. default TRUE, remove repeated pindels or not.
#' @param cutoff_repeat numeric. default 2, cutoff for repeating frequency.

#' @param consistent_pindel_statistic_filter logic. default TRUE, remove pindels failed to pass consistency statistical test (binom test).
#' @param pindel_for_consistent_statistic data.frame. containing gene,number of pindels, number of samples (columns: gene,n_pindel,n_sample), used for binom test with targeted pindels.
#' @param cutoff_binom_pval numeric. default 0.05, cutoff for binom test pvalue.

#' @param description_file character. default description.txt, file name for deposit of information of each operation (filter sets).
#' @param save_traced_pindel logic. default TRUE, save traced_pindel_file.
#' @param traced_pindel_file character. default traced_pindel.txt, file name for deposit results after each operation (filter sets) if save_traced_pindel=TRUE.

#'
#' @return data.frame. filtered pindels.
#'


filter_pindels<-function(pindels,traced_columns=c("Sample_ID","gene","type","chr","start","end"),
                         low_purity_filter=T,low_purity_samples,
                         treat_noisy_samples=F,noisy_samples=c(),
                         exonutronly_filter=T,keep=c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),
                         remove_td_filter=T,
                         td_coverage_filter=T,cutoff_td_coverage=1.5,
                         off_target_filter=T,intervals=intervals,
                         germline_filter=T,cutoff_n_vaf=0.01,
                         treat_known_cancer_gene_specially=T,known_cancer_genes=c("TP53","GNA11","GNAQ"),special_pindel_genes=c("BAP1","PRKDC"),
                         mapping_quality_filter=T,cutoff_mq=25,
                         statistical_filter=T,cutoff_prop_test_p=0.01,
                         t_vaf_filter=T,cutoff_t_vaf=0.0,
                         tumor_f_filter=T,cutoff_tumor_f=0,
                         deletion_length_filter=T,cutoff_d_length=1000,
                         insertion_length_filter=T,cutoff_i_length=6,
                         tumor_coverage_filter=T,cutoff_t_coverage=10,
                         normal_coverage_filter=T,cutoff_n_coverage=10,
                         support_filter=T,cutoff_support=3,
                         homopolymer_filter=T,maxn_homopolymer=6,
                         repeat_pindel_filter=T,cutoff_repeat=2,
                         consistent_pindel_statistic_filter=T,pindel_for_consistent_statistic=NA,cutoff_binom_pval=0.05,
                         description_file="description.txt",
                         save_traced_pindel=T,traced_pindel_file="traced_pindel.txt"){
  #filter derived from SVFilter,gatk,and github.com/genome/pindel
  options(warn=-1)
  pindels<-pindels[!duplicated(do.call(paste,pindels[,traced_columns])),]
  fileConn<-description_file
  write(paste(paste("Starting from ",nrow(pindels),sep="")," indels.",sep=""),fileConn,append=F)
  genome_coverage_tf<- pindels %>% dplyr::group_by(Sample_ID) %>% dplyr::summarise(genome_coverage=sum(t_ref_count+t_alt_count)/dplyr::n())
  genome_coverage<-genome_coverage_tf[["genome_coverage"]]
  names(genome_coverage)<-genome_coverage_tf[["Sample_ID"]]
  pindels_traced<-list()
  traced_names<-c()
  pindels_traced[["original_pindels"]]<-pindels[,traced_columns]
  traced_names<-c(traced_names,"original_pindels")

  #low purity filter
  if(low_purity_filter){
    pindels<-pindels[!pindels[["Sample_ID"]] %in% low_purity_samples,]
    write("After filter out low purity (0.5);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_low_purity"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_low_purity")
  }
  # special treatment of very noisy samples (too many calls)
  if(treat_noisy_samples){
    clean_pindels<-pindels[!pindels[["Sample_ID"]] %in% noisy_samples,]
    noisy_pindels<-pindels[pindels[["Sample_ID"]] %in% noisy_samples,]
    denoisy_pindels<-noisy_pindels[noisy_pindels[["gene"]] %in% unique(clean_pindels[['gene']]),]
    denoisy_pindels <- denoisy_pindels %>% dplyr::group_by(Sample_ID,gene) %>% dplyr::mutate(hits=dplyr::n(),max_t_alt_count=max(t_alt_count))
    denoisy_pindels <- as.data.frame(denoisy_pindels)
    denoisy_pindels<-denoisy_pindels[denoisy_pindels$hits==1 | (denoisy_pindels$hits>1 & denoisy_pindels$t_alt_count==denoisy_pindels$max_t_alt_count),]
    denoisy_pindels<-denoisy_pindels[,colnames(clean_pindels)]
    pindels<-rbind(denoisy_pindels,clean_pindels)
    write("After specially treatment of noisy samples;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," pindels.",sep=""), fileConn,append = T)
    pindels_traced[["treat noisy samples"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"treat_noisy_samples")
  }
  #exon UTRs only
  if(exonutronly_filter){
    pindels<-pindels[pindels$func.knowngene %in% c("exonic","exonic;splicing","splicing","UTR3","UTR5","UTR5;UTR3"),]
    write("After filter out non exonic or UTRs", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_nonexonutr"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #remove TD
  if(remove_td_filter){
    pindels<-pindels[pindels$type!="TD",]
    write("After filter out tandem duplicate indels;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_TD"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  } else if(td_coverage_filter){
    #sequencing relative depth (ratio between the average of sequencing depth in the duplicated region and that over entire genome) filter >=1.5
    genome_coverage<-genome_coverage[pindels[['Sample_ID']]]
    pindels<-pindels[!(pindels$type=="TD" & ((pindels$t_ref_count+pindels$t_alt_count)>=(cutoff_td_coverage*genome_coverage))),]
    write(paste("After filter out low coverage tandem duplicate indels with cutoff_td_coverage ",cutoff_td_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_TD_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonexonutr")
  }

  #offtarget filter
  if(off_target_filter){
    pindels_gr<-df2granges(pindels,genome = "hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,seqnames_col = traced_columns[4],start_col = traced_columns[5],end_col = traced_columns[6],meta_cols = colnames(pindels))
    GenomeInfoDb::seqlevels(intervals)<-GenomeInfoDb::seqlevels(pindels_gr)
    GenomeInfoDb::seqinfo(intervals)<-GenomeInfoDb::seqinfo(pindels_gr)
    targets<-plyranges::filter_by_overlaps(pindels_gr,intervals)
    pindels<-as.data.frame(S4Vectors::mcols(targets))
    write("After filter out off-targets indel;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_off_target"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_off_target")
  }

  #germline cutoff 0.01
  if(germline_filter){
    pindels<-pindels[pindels$n_vaf<=cutoff_n_vaf,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out germline indel with cutoff_n_vaf ",cutoff_n_vaf,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_germline"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_germline")
  }

  #treat known cancer gene specially
  if(treat_known_cancer_gene_specially){
    if(length(known_cancer_genes)>0){
      known_cancer_gene_pindels<-pindels[pindels$gene %in% setdiff(known_cancer_genes,special_pindel_genes),]
    }
    if(length(special_pindel_genes)>0){
      special_cancer_gene_pindels<-pindels[pindels$gene %in% special_pindel_genes,,drop=F]
    }
    pindels<-pindels[!(pindels$gene %in% known_cancer_genes),]
    write("After split known cancer gene (cosmic or cancer_genetic);", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," non known cancer gene pindels.",sep=""), fileConn,append = T)
    write(paste(paste("There are still ",nrow(known_cancer_gene_pindels),sep="")," known cancer gene pindels.",sep=""), fileConn,append = T)
    pindels_traced[["non_known_cancer_gene_pindel"]]<-pindels[,traced_columns]
    pindels_traced[["known_cancer_gene_pindel"]]<-known_cancer_gene_pindels[,traced_columns]
    traced_names<-c(traced_names,"non_known_cancer_gene_pindel")
    traced_names<-c(traced_names,"known_cancer_gene_pindel")
    if(length(special_pindel_genes)>0){
      write(paste(paste("There are still ",nrow(special_cancer_gene_pindels),sep="")," special cancer gene pindels.",sep=""), fileConn,append = T)
      pindels_traced[["special_cancer_gene_pindel"]]<-special_cancer_gene_pindels[,traced_columns]
      traced_names<-c(traced_names,"special_cancer_gene_pindel")
    }
  }

  #mapping quality filter
  if(mapping_quality_filter){
    pindels<-pindels[(pindels$mq/pindels$support)>=cutoff_mq,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low mapping quality indels with cutoff_mq ",cutoff_mq,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_mapping_quality"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_mapping_quality")
  }

  #statistical filter
  if(statistical_filter){
    pval<-apply(pindels[,c("t_ref_count","t_alt_count",'n_ref_count',"n_alt_count")],1,function(row){stats::prop.test(x=c(as.integer(row[2]),as.integer(row[4])),n=c(as.integer(row[1])+as.integer(row[2])+1,as.integer(row[3])+as.integer(row[4])+1))$p.value})
    pindels<-pindels[pval<=cutoff_prop_test_p,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out prop test failed indels with cufoff_prop_test_p ",cutoff_prop_test_p,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_prop_statistic"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_prop_statistic")
  }

  # tumor_f filter
  if(tumor_f_filter){
    pindels<-pindels[pindels$tumor_f>=cutoff_tumor_f,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out lwo tumor fraction indels with cutoff_tumor_f ",cutoff_tumor_f,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_tumor_f"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_f")
  }

  # t_vaf filter
  if(t_vaf_filter){
    pindels<-pindels[pindels$t_vaf>=cutoff_t_vaf,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low t_vaf indels with cutoff_t_vaf ",cutoff_t_vaf,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_t_vaf"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_t_vaf")
  }

  #deletion length filter
  if(deletion_length_filter){
    pindels<-pindels[pindels$length<=cutoff_d_length,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out  too long deletions with cutoff_d_length",cutoff_d_length,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_deletion_length"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_deletion_length")
  }

  #insertion length filter
  if(insertion_length_filter){
    pindels<-pindels[pindels$insertlen<=cutoff_i_length,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out too long insertion with cutoff_i_length ",cutoff_i_length,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_insertion_length"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_insertion_length")
  }

  #tumor depth coverage filter
  if(tumor_coverage_filter){
    pindels<-pindels[(pindels$t_alt_count+pindels$t_ref_count)>=cutoff_t_coverage,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low tumor coverage indels with cutoff_t_coverage ",cutoff_t_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_tumor_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_tumor_coverage")
  }

  #normal depth coverage filter
  if(normal_coverage_filter){
    pindels<-pindels[(pindels$n_alt_count+pindels$n_ref_count)>=cutoff_n_coverage,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low normal coverage indels with cutoff_n_coverage ",cutoff_n_coverage,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_normal_coverage"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_normal_coverage")
  }

  #support filter
  if(support_filter){
    pindels<-pindels[pindels$support>=cutoff_support,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out low support indels with cutoff_support",cutoff_support,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_low_support"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_low_support")
  }

  #homopolymer filter
  if(homopolymer_filter){
    #homopolymer(many copies of a single repeating unit) filer maxlength=6
    homopolymer_pat<-paste(c(paste("A{",maxn_homopolymer,",}",sep=""),paste("T{",maxn_homopolymer,",}",sep=""),paste("C{",maxn_homopolymer,",}",sep=""),paste("G{",maxn_homopolymer,",}",sep="")),collapse="|")
    pindels<-pindels[!grepl(homopolymer_pat,pindels$refs,perl = T) | grepl(homopolymer_pat,pindels$sams,perl = T),]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out too much homopolymer indels with cutoff_maxn_homopolymer ",maxn_homopolymer,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["filter_homopolymer"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_homopolymer")
  }

  #combine known_cancer gene pindels
  if(treat_known_cancer_gene_specially & length(known_cancer_genes)>0){
    pindels<-rbind(pindels,known_cancer_gene_pindels)
    write("After combine known cancer gene pindels;", fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["combined_known_pindels"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"combined_known_pindels")
  }

  #reapeat pindel filter
  if(repeat_pindel_filter){
    pindels_<-pindels
    pindels_[["Pindel_ID"]]<-do.call(paste,pindels_[,c("gene","chr","start")])
    repeat_pindel_IDs<-names(table(pindels_[["Pindel_ID"]]))[table(pindels_[["Pindel_ID"]])>=cutoff_repeat]
    pindels_<-pindels_[!(pindels_[["Pindel_ID"]] %in% repeat_pindel_IDs),]
    pindels<-pindels_[,-ncol(pindels_)]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out high repeated indels with cutoff_repeat ",cutoff_repeat,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    if(treat_known_cancer_gene_specially){
      write(paste(paste("There are still ",sum(pindels$gene %in% known_cancer_gene_pindels$gene),sep="")," known cancer gene indels.",sep=""), fileConn,append = T)
      write(paste(paste("There are ",nrow(known_cancer_gene_pindels)-sum(pindels$gene %in% known_cancer_genes),sep="")," known cancer gene indels removed.",sep=""), fileConn,append = T)
    }
    pindels_traced[["filter_repeat_pindel"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_repeat_pindel")
  }

  if(consistent_pindel_statistic_filter){
    n_sample_subject<-length(unique(pindels[["Sample_ID"]]))
    pindels_subject<-pindels %>% dplyr::group_by(gene) %>% dplyr::summarise(n_pindel_subject=dplyr::n())
    pindels_subject[["n_sample_subject"]]<-n_sample_subject
    pindels_for_statistic<-merge(pindels_subject,pindel_for_consistent_statistic,by="gene",all.x=T)
    pindels_for_statistic[["n_sample"]][is.na(pindels_for_statistic[["n_sample"]])]=max(pindel_for_consistent_statistic$n_sample)
    pindels_for_statistic[["n_pindel"]][is.na(pindels_for_statistic[["n_pindel"]])]<-0
    pindels_for_statistic<-set_column_as_rownames(pindels_for_statistic,"gene")
    pindels_for_statistic<-pindels_for_statistic[pindels_for_statistic[["n_pindel_subject"]]<=pindels_for_statistic[["n_sample_subject"]],]
    pindels_for_statistic$n_pindel_subject<-apply(pindels_for_statistic[,c("n_pindel_subject","n_sample_subject")],1,min)
    pindels_for_statistic$n_pindel<-apply(pindels_for_statistic[,c("n_pindel","n_sample")],1,min)
    pvals<-apply(pindels_for_statistic,1,function(r){test<-stats::binom.test(r["n_pindel_subject"],r["n_sample_subject"],p=r["n_pindel"]/r["n_sample"],alternative="greater");return(test$p.value)})
    adjpvals<-stats::p.adjust(pvals)
    keep_genes<-names(adjpvals)[adjpvals>cutoff_binom_pval]
    pindels<-pindels[pindels[["gene"]] %in% keep_genes,]
    pindels<-pindels[rownames(pindels)!="NA",]
    write(paste("After filter out non consitent pindel with tcga pindel with cutoff_binom_pval ",cutoff_binom_pval,sep=" "), fileConn,append = T)
    write(paste(paste("There are still ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    if(treat_known_cancer_gene_specially){
      write(paste(paste("There are still ",sum(pindels$gene %in% known_cancer_gene_pindels$gene),sep="")," known cancer gene indels.",sep=""), fileConn,append = T)
      write(paste(paste("There are ",nrow(known_cancer_gene_pindels)-sum(pindels$gene %in% known_cancer_genes),sep="")," known cancer gene indels removed.",sep=""), fileConn,append = T)

    }
    pindels_traced[["filter_nonconsistent_pindel_statistic"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"filter_nonconsistent_pindel_statistic")
  }




  #combine special cancer gene pindels
  if(treat_known_cancer_gene_specially & length(special_pindel_genes)>0){
    pindels<-rbind(pindels,special_cancer_gene_pindels)
    write("After combine special cancer gene pindels;", fileConn,append = T)
    write(paste(paste("There are still totally ",nrow(pindels),sep="")," indels.",sep=""), fileConn,append = T)
    pindels_traced[["combined_final_pindels"]]<-pindels[,traced_columns]
    traced_names<-c(traced_names,"combined_final_pindels")
  }


  if(save_traced_pindel){
    pindel_order=NA
    for(traced_name in rev(traced_names)){
      if(traced_name %in% c("special_cancer_gene_pindel","known_cancer_gene_pindel")){
        next
      } else{
        pindels_temp=pindels_traced[[traced_name]]
        rownames(pindels_temp)<-do.call(paste,pindels_temp)
        pindels_temp<-pindels_temp[c(intersect(pindel_order,rownames(pindels_temp)),setdiff(rownames(pindels_temp),pindel_order)),]
        pindel_order=rownames(pindels_temp)
        rownames(pindels_temp)<-NULL
        pindels_traced[[traced_name]]=pindels_temp
      }
    }
    for(traced_name in traced_names){
      pindels_temp<-pindels_traced[[traced_name]]
      colnames(pindels_temp)<-paste(traced_name,colnames(pindels_temp),sep=":")
      utils::write.table(pindels_temp,file = traced_pindel_file,append = T,row.names = F,sep="\t")
    }
  }
  return(pindels)
}


#' convert_mutations
#'
#' convert mutations used for mutsig2cv, maftools and pynbs
#'
#' @param mutations data.frame. mutations data.frame including SNPs and Indels.
#' @param target character. targeted algorithm.
#' @param msc_snp_conversion_dict character. a named character mapping vector between cgl snp and mutsig2cv snp.
#' @param msc_indel_conversion_dict character. a named character mapping vector between cgl indel and mutsig2cv indel.
#' @param maftools_mutation_conversion_dict character. a named character mapping vector between cgl snp/indel and maftools snp/indel.
#' @param pynbs_nonSilent_variant_classication character. nonSilent variant classification for pyNBS.
#'
#' @return data.frame. converted data.frame fit specifid algorithm
#' @export
#'

convert_mutations<-function(mutations,target=c("mutsig2cv","maftools","pynbs"),msc_snp_conversion_dict,msc_indel_conversion_dict,maftools_mutation_conversion_dict,pynbs_nonSilent_variant_classication){
  target=match.arg(target)
  snps<-mutations[mutations$type=="SNP",]
  indels<-mutations[mutations$type!="SNP",]

  #SNP
  Hugo_Symbol<-snps$gene.knowngene
  Entrez_Gene_Id<-snps$entrez_gene_id
  Center="IPCT"
  NCBI_Build<-"hg19"
  Chromosome<-snps$chrom
  Start_position<-snps$start
  End_position<-snps$end
  Strand="*"
  Reference_Allele<-snps$ref_allele
  Tumor_Seq_Allele2<-snps$alt_allele
  Variant_Classification<-ifelse(snps$exonicfunc.knowngene==".",snps$func.knowngene,snps$exonicfunc.knowngene)     #nonsynonymous SNV -- missense, stopgain -- nonsense, synonymous SNV -- quiet
  Variant_Classification<-unname(msc_snp_conversion_dict[Variant_Classification])
  Variant_Type<-"SNP"
  Tumor_Sample_Barcode<-snps[["Sample_ID"]]
  i_TumorVAF_WU<-as.numeric(snps$tumor_f)*100
  cdna=snps$cdna
  Protein_Change<-snps$aaannotation
  snps<-data.frame(Hugo_Symbol=Hugo_Symbol,
                   Entrez_Gene_Id=Entrez_Gene_Id,
                   Center=Center,
                   NCBI_Build=NCBI_Build,
                   Chromosome=Chromosome,
                   Start_position=Start_position,
                   End_position=End_position,
                   Strand=Strand,
                   Variant_Classification=Variant_Classification,
                   Variant_Type=Variant_Type,
                   Reference_Allele=Reference_Allele,
                   Tumor_Seq_Allele1=Reference_Allele,
                   Tumor_Seq_Allele2=Tumor_Seq_Allele2,
                   Tumor_Sample_Barcode=Tumor_Sample_Barcode,
                   Protein_Change=Protein_Change,
                   i_TumorVAF_WU=i_TumorVAF_WU,
                   cdna=cdna,
                   stringsAsFactors = F)
  snps<-snps[snps$Variant_Classification %in% msc_snp_conversion_dict,]
  #INDEL
  Hugo_Symbol<-indels$gene.knowngene
  Entrez_Gene_Id<-indels$entrez_gene_id
  Center="IPCT"
  NCBI_Build<-"hg19"
  Chromosome<-indels$chrom
  Start_position<-indels$start
  End_position<-indels$end
  Strand="*"
  Reference_Allele<-indels$ref_allele
  Tumor_Seq_Allele2<-indels$alt_allele
  Variant_Classification<-ifelse(indels$exonicfunc.knowngene==".",
                                 ifelse(indels$type=="INS",paste(indels$func.knowngene,"_Ins",sep=""),ifelse(indels$type=="DEL",paste(indels$func.knowngene,"_Del",sep=""),paste(indels$func.knowngene,"_TD",sep=""))),
                                 ifelse(indels$exonicfunc.knowngene=="stopgain",ifelse(indels$type=="I",paste("Frame_shift_",indels$exonicfunc.knowngene,"_Ins",sep=""),paste("Frame_shift_",indels$exonicfunc.knowngene,"_Del",sep="")),indels$exonicfunc.knowngene)
  ) # INS --insertion, DEL -- Deletion, TD -- tandem duplicate
  Variant_Classification<-unname(msc_indel_conversion_dict[Variant_Classification])
  Variant_Type<-indels$type
  Tumor_Sample_Barcode<-indels[["Sample_ID"]]
  i_TumorVAF_WU<-as.numeric(indels$tumor_f)*100
  cdna=indels$cdna
  Protein_Change<-indels$aaannotation
  indels<-data.frame(Hugo_Symbol=Hugo_Symbol,
                     Entrez_Gene_Id=Entrez_Gene_Id,
                     Center=Center,
                     NCBI_Build=NCBI_Build,
                     Chromosome=Chromosome,
                     Start_position=Start_position,
                     End_position=End_position,
                     Strand=Strand,
                     Variant_Classification=Variant_Classification,
                     Variant_Type=Variant_Type,
                     Reference_Allele=Reference_Allele,
                     Tumor_Seq_Allele1=Reference_Allele,
                     Tumor_Seq_Allele2=Tumor_Seq_Allele2,
                     Tumor_Sample_Barcode=Tumor_Sample_Barcode,
                     Protein_Change=Protein_Change,
                     i_TumorVAF_WU=i_TumorVAF_WU,
                     cdna=cdna,
                     stringsAsFactors = F)
  indels<-indels[indels$Variant_Classification %in% msc_indel_conversion_dict,]
  mutations<-rbind(snps,indels)
  mutations<-mutations[!is.na(mutations$Variant_Classification),]

  if(target=="mutsig2cv"){
    return(mutations)
  }
  if(target=="maftools"){
    mutations$Variant_Classification<-unname(maftools_mutation_conversion_dict[mutations$Variant_Classification])
    mutations<-mutations[mutations$Variant_Classification %in% maftools_mutation_conversion_dict,]
    mutations<-mutations[!is.na(mutations$Variant_Classification),]
    return(mutations)
  }
  if(target=="pynbs"){
    mutations<-mutations[mutations$Variant_Classification %in% pynbs_nonSilent_variant_classication,c("Tumor_Sample_Barcode","Hugo_Symbol")]
    mutations<-mutations[!duplicated(mutations),]
    return(mutations)
  }
}


#' complex_oncoplot
#' oncoplot of mutation and cnv
#'
#' @param snp_indels data.frame. longer_style snp and indels (only nonSilent), and also cnv data, should be provided without default, column of sample id should be factor.
#' @param for_row character. col used for row names of oncoplot data,default value: Hugo_Symbol.
#' @param for_column character. col used for column names of oncoplot data,default value: Tumor_Sample_Barcode.
#' @param for_value character. col used for value of oncoplot data,default value: Variant_Classification.
#' @param cnv logic. include cnv into the oncoplot,default TRUE.
#' @param cnv_types character vector. if cnv is TRUE, what kinds of alteration would be regarded as cnv
#' @param annotate_source logic. hether to mark mutation source.
#' @param source_column character. column for source.
#' @param source_label character. label for source.
#' @param sample_order character. sample order.
#' @param selected_genes character. selected genes that will be shown.
#' @param gene_order character. gene order.
#' @param mutsigCV2_sig_genes data.frame. significant genes from mutsig2CV.
#' @param mutsigCV2_sig_p_threshold numeric. threshold used to filter significant driver genes.
#' @param multi_hit logic. combine various mutation and cnv as multi_hit, default FALSE.
#' @param multi_hit_col character. color for multi_hit, default black.
#' @param top_filter numeric. number of top frequent genes would be plotted,default 20.
#' @param col character.color pallete for oncoplot, default NA that oncoplot will use its pallete.
#' @param bg_color character.background color, default grey (#' dcdedc).
#' @param remove_empty_sample logic. remove empty samples, default FALSE.
#' @param remove_macromolecular_gene logic. remove macromolecular genes that detemined by macromolecular_threshold, default FALSE.
#' @param mark_macromolecular_gene logic. mark macromolecular genes with macromolecular_mark.
#' @param macromolecular_mark character. used to mark macromolecular gene if set remove_macromolecular_gene as FALSE, default *.
#' @param macromolecular_threshold numeric. macromolecular threshold that determine whether macromolecular, default 10000bp.
#' @param top_annotation HeatmapAnotation. top annotation for oncoplot, default NULL, if NULL, total alteration burden gona be used as top annotaion.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param output_dir character. output plot directory.
#' @param ... list. other params used for oncoprint.
#'
#' @return data.frame. mutations.
#' @export
#'




complex_oncoplot<-function(snp_indels,
                           for_row="Hugo_Symbol",
                           for_column="Tumor_Sample_Barcode",
                           for_value="Variant_Classification",
                           cnv=TRUE,
                           cnv_types=c("Amplification","Deletion"),
                           annotate_source=F,
                           source_column="Source",
                           source_label=c("MDL"="*"),
                           sample_order=NULL,
                           selected_genes,
                           gene_order,
                           mutsigCV2_sig_genes=NA,
                           mutsigCV2_sig_p_threshold=0.05,
                           multi_hit=FALSE,
                           multi_hit_col="black",
                           top_filter=20L,
                           col,
                           bg_color="#dcdedc",
                           remove_empty_sample=F,
                           remove_macromolecular_gene=F,
                           mark_macromolecular_gene=F,
                           macromolecular_mark="*",
                           macromolecular_threshold=10000,
                           top_annotation=NULL,
                           top_anno_sample_order=NULL,
                           output_dir="./",
                           ...){

  f<-function(l,multi_hit) {
    if(length(l)==1){
      return(as.character(l))
    } else{
      if(multi_hit){
        cn_term<-l[l %in% c("Amp","Gain","Loss","Del")]
        mutation_terms<-l[!l %in% c("Amp","Gain","Loss","Del")]
        if(length(mutation_terms)>1){
          mutation_term="Multi_Hit"
        } else{
          mutation_term=mutation_terms
        }
        return (paste(c(mutation_term,cn_term),collapse = ";"))
      } else{
        return(paste(l,collapse=";"))
      }
    }}
  #extract and prepare snp indel data
  snp_indels<-snp_indels[,c(for_row,for_column,for_value,source_column)]
  colnames(snp_indels)[1:3]<-c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification")
  snp_indels[["Tumor_Sample_Barcode"]]<-factor(snp_indels[["Tumor_Sample_Barcode"]])
  if(!cnv){
    snp_indels<-snp_indels[! snp_indels$Variant_Classification %in% cnv_types,]
    col=col[!names(col) %in% cnv_types]
  }

  DNA_samples<-top_anno_sample_order
  #only mutsigcv significant gene mutation allowed for following analysis
  snp_indels<-snp_indels[snp_indels$Variant_Classification!="Silent" & (!is.na(snp_indels$Variant_Classification)),]
  snp_indels_o<-snp_indels
  #prepare tumor mutation burden data
  # make sure rownames of tmb should be consistent with colnames top_snp_indels
  ctmb=snp_indels %>% dplyr::group_by(Tumor_Sample_Barcode,Variant_Classification) %>% dplyr::summarise(freq=dplyr::n()) %>% tidyr::pivot_wider(id_cols = "Tumor_Sample_Barcode", names_from="Variant_Classification",values_from="freq",values_fill=0)
  ctmb<-as.data.frame(ctmb)
  rownames(ctmb)<-ctmb[["Tumor_Sample_Barcode"]]
  ctmb<-ctmb[,-1]
  ctmb<-ctmb[,colnames(ctmb)!="NA"]
  ctmb<-ctmb[DNA_samples,]
  tmb<-ctmb[,!colnames(ctmb) %in% cnv_types]
  #tmb<-tmb[colnames(top_snp_indels),]


  # prepare top n snp indel data and then data for oncoprint
  snp_indels_p<-as.data.frame(tidyr::pivot_wider(snp_indels,id_cols = "Hugo_Symbol",names_from="Tumor_Sample_Barcode",values_from="Variant_Classification",values_fn=pryr::partial(f,multi_hit = multi_hit)))
  rownames(snp_indels_p)<-snp_indels_p[[1]]
  snp_indels_p<-as.matrix(snp_indels_p[,-1,drop=F])
  snp_indels_p[is.na(snp_indels_p)]<-""
  if(!is.na(match("None",rownames(snp_indels_p)))){
    snp_indels_p<-snp_indels_p[-(match("None",rownames(snp_indels_p))),]
  }
  snp_indels_p<-snp_indels_p[,DNA_samples,drop=F]

  if(missing(selected_genes)){
    selected_genes<-mutsigCV2_sig_genes[["gene"]][mutsigCV2_sig_genes[["p"]]<=mutsigCV2_sig_p_threshold]
    macromolecular_genes<-mutsigCV2_sig_genes[["gene"]][mutsigCV2_sig_genes[["codelen"]]>=macromolecular_threshold]
  } else{
    remove_macromolecular_gene=F
    mark_macromolecular_gene=F
  }
  snp_indels_p<-snp_indels_p[rownames(snp_indels_p) %in% selected_genes,,drop=F]

  top_filter<-min(c(nrow(snp_indels_p),top_filter))
  if(is.integer(top_filter)){
    top_mutated_genes<-rownames(snp_indels_p)[order(rowSums(snp_indels_p!=""),decreasing = T)][1:top_filter]
  } else{
    top_mutated_genes<-rownames(snp_indels_p)[(rowSums(snp_indels_p!="")/ncol(snp_indels_p))>=top_filter]
  }
  snp_indels_p<-snp_indels_p[top_mutated_genes,]

  if(remove_macromolecular_gene){snp_indels_p<-snp_indels_p[setdiff(rownames(snp_indels_p),macromolecular_genes),]}
  if(mark_macromolecular_gene) {rownames(snp_indels_p)[rownames(snp_indels_p) %in% macromolecular_genes] <-paste(macromolecular_mark,rownames(snp_indels_p)[rownames(snp_indels_p) %in% macromolecular_genes],sep="")}



  # if remove empty sample
  if(remove_empty_sample){
    empty_samples<-colnames(snp_indels_p)[colSums(snp_indels_p!="")==0]
    DNA_samples_<-setdiff(DNA_samples,empty_samples)
    anno_DNA<-anno_DNA[match(DNA_samples_,DNA_samples)]
    tmb<-tmb[DNA_samples_,,drop=F]
    snp_indels_p<-snp_indels_p[,DNA_samples_,drop=F]
    sample_order<-sample_order[-match(empty_samples,sample_order)]
  }


  if(annotate_source){
    source_p<-as.data.frame(tidyr::pivot_wider(snp_indels_o,id_cols = "Hugo_Symbol",names_from="Tumor_Sample_Barcode",values_from=source_column,values_fn=pryr::partial(f,multi_hit = FALSE)))
    rownames(source_p)<-source_p[[1]]
    source_p<-as.matrix(source_p[,-1,drop=F])
    source_p[is.na(source_p)]<-""
    if(!is.na(match("None",rownames(source_p)))){
      source_p<-source_p[-(match("None",rownames(source_p))),]
    }

    if(remove_macromolecular_gene){
      source_p<-source_p[setdiff(rownames(source_p),macromolecular_genes),]
    }
    if(mark_macromolecular_gene){
      rownames(source_p)[rownames(source_p) %in% macromolecular_genes] <-paste(macromolecular_mark,rownames(source_p)[rownames(source_p) %in% macromolecular_genes],sep="")
    }
    source_p<-source_p[rownames(snp_indels_p),colnames(snp_indels_p)]
  }




  # if remove macromolecular genes
  # extract macrocolecular genes
  #hs<-org.Hs.eg.db
  #hg19<-TxDb.Hsapiens.UCSC.hg19.knownGene
  #cds<-transcriptLengths(hg19,with.cds_len = T)
  #cds<-cds[(cds$cds_len!=0) & (!is.na(cds$gene_id)),]
  #max_cds<-cds %>% dplyr::group_by(gene_id) %>% dplyr::summarise(max_cds_len=max(cds_len)) %>% arrange(desc(max_cds_len))
  #macromolecular_genes<-max_cds$gene_id[max_cds$max_cds_len>=macromolecular_theshold]
  #macromolecular_genes<-AnnotationDbi::select(hs,keys=macromolecular_genes,columns=c("ENTREZID", "SYMBOL"),keytype = "ENTREZID")$SYMBOL

  #if col is na
  if(missing(col)){
    palette19<-c("#800000","#9A6324","#808000","#469990","#000075","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","#fabed4","#ffd8b1","#fffac8","#aaffc3","#dcbeff")
    col<-sample(palette19,ncol(ctmb),replace = F)
    names(col)<-colnames(ctmb)
  } else{
    col<-col[colnames(ctmb)]
  }
  tmb_col=col[!names(col) %in% cnv_types]

  testit::assert("row names of tmb should be same as col names of snp_indels_p ", all(rownames(ctmb)==colnames(snp_indels_p)))
  testit::assert("names of color should be same as col names as tmb",all(names(col)==colnames(ctmb)))

  # plot heatmap annotation of tmb

  tmb_anno<-ComplexHeatmap::HeatmapAnnotation(tmb=ComplexHeatmap::anno_barplot(tmb,gp=grid::gpar(fill=tmb_col,col=tmb_col),border=F,bar_width = 0.85,height=grid::unit(4,"cm")),annotation_name_side="left")
  #plot oncoplot
  if(multi_hit){
    col=c(col,"Multi_Hit"=multi_hit_col)
  }
  if(cnv){multi_hit=F}
  onco_p<-ComplexHeatmap::oncoPrint(snp_indels_p,
                                    alter_fun=function(x,y,w,h,v){
                                      if(multi_hit){
                                        grid::grid.rect(x,y,w*0.9,h*0.9,gp=grid::gpar(fill=bg_color,col=NA))
                                        for(variant_c in names(v)){
                                          if(v[variant_c]){
                                            grid::grid.rect(x, y, w*0.9, h*0.9,gp = grid::gpar(fill = col[variant_c],col=NA))
                                          }
                                        }
                                      } else{
                                        if(cnv){
                                          mutation_v<-v[! names(v) %in% cnv_types ]
                                          cnv_v<-v[names(v) %in% cnv_types]
                                        } else{
                                          mutation_v<-v
                                        }
                                        cell_height=grid::unit(1/nrow(snp_indels_p),"npc")
                                        cell_bar_margin=0.05*cell_height
                                        stacks_in_cell=ifelse(sum(mutation_v)==0,1,sum(mutation_v))
                                        y_=y-cell_height*0.5+cell_bar_margin+grid::unit(1/nrow(snp_indels_p)/stacks_in_cell*0.9/2,"npc")
                                        grid::grid.rect(x,y,w*0.9,h*0.9,gp=grid::gpar(fill=bg_color,col=NA))
                                        for(variant_c in names(mutation_v)){
                                          if(mutation_v[variant_c]){
                                            grid::grid.rect(x, y_, w*0.9, h*0.9*(1/stacks_in_cell),gp = grid::gpar(fill = col[variant_c], col = NA))
                                            y_=y_+h*0.9*(1/stacks_in_cell)
                                          }
                                        }
                                        if(cnv){
                                          for(cnv_c in names(cnv_v)){
                                            if(cnv_v[cnv_c]){
                                              grid::grid.rect(x, y, w*0.75, h*0.9,gp = grid::gpar(fill = NA, col = col[cnv_c],lwd=3))
                                            }
                                          }
                                        }
                                      }
                                    },
                                    top_annotation = c(tmb_anno,top_annotation),
                                    pct_side = "right",
                                    row_names_side = "left",
                                    col=col,
                                    show_column_names = T,
                                    column_order=sample_order,...)
  method_lgd = ComplexHeatmap::Legend(labels = "MDL", title = "Support", type = "points", pch = 8, legend_gp = grid::gpar(col = "black"), background = bg_color)
  grDevices::svg(file.path(output_dir,"complex_oncoplot.svg"),width=8,height=12)
  ComplexHeatmap::draw(onco_p,annotation_legend_list=list(method_lgd))
  if(annotate_source){
    snp_indels_p_<-snp_indels_p[ComplexHeatmap::row_order(onco_p),ComplexHeatmap::column_order(onco_p)]
    source_p<-source_p[ComplexHeatmap::row_order(onco_p),ComplexHeatmap::column_order(onco_p)]
    heatmap_oncoprint_body<-gsub("heatmap_","",ComplexHeatmap::list_components()[grepl("heatmap_oncoPrint",ComplexHeatmap::list_components())])
    x_step=1/ncol(source_p)
    x_start=1/ncol(source_p)/2
    y_step=1/nrow(source_p)
    y_start=1/nrow(source_p)/2

    for(i in 1:nrow(source_p)){
      for(j in 1:ncol(source_p)){
        #if(top_source[i,j]==names(source_label)){
        #  ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.rect(x_start+(j-1)*x_step,1-y_start-(i-1)*y_step,x_step*0.9,y_step*0.9,gp=grid::gpar(fill=bg_color,col=col[top_snp_indels_[i,j]],lwd=2)))
        #}
        if(all(grepl("MDL",unlist(strsplit(source_p[i,j],";"))))){
          ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.rect(x_start+(j-1)*x_step,1-y_start-(i-1)*y_step,x_step*0.9,y_step*0.9,gp=grid::gpar(fill=bg_color,col=col[snp_indels_p_[i,j]],lwd=2)))
        }
        if(grepl(names(source_label),source_p[i,j])){
          ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.text(source_label,x=x_start+(j-1)*x_step,y=1-y_start*1.3-(i-1)*y_step))
        }
      }
    }
  }
  grDevices::dev.off()
  return(snp_indels_p)
}

#' complex_oncoplot_o
#' oncoplot of mutation and cnv
#'
#' @param snp_indels data.frame. longer_style snp and indels (only nonSilent), and also cnv data, should be provided without default, column of sample id should be factor.
#' @param for_row character. col used for row names of oncoplot data,default value: Hugo_Symbol.
#' @param for_column character. col used for column names of oncoplot data,default value: Tumor_Sample_Barcode.
#' @param for_value character. col used for value of oncoplot data,default value: Variant_Classification.
#' @param cnv logic. include cnv into the oncoplot,default TRUE.
#' @param cnv_types character vector. if cnv is TRUE, what kinds of alteration would be regarded as cnv
#' @param annotate_source logic. hether to mark mutation source.
#' @param source_column character. column for source.
#' @param source_label character. label for source.
#' @param sample_order character. sample order.
#' @param selected_genes character. selected genes that will be shown.
#' @param gene_order character. gene order.
#' @param mutsigCV2_sig_genes data.frame. significant genes from mutsig2CV.
#' @param mutsigCV2_sig_p_threshold numeric. threshold used to filter significant driver genes.
#' @param multi_hit logic. combine various mutation and cnv as multi_hit, default FALSE.
#' @param multi_hit_col character. color for multi_hit, default black.
#' @param top_filter numeric. number of top frequent genes would be plotted,default 20.
#' @param col character.color pallete for oncoplot, default NA that oncoplot will use its pallete.
#' @param bg_color character.background color, default grey (#' dcdedc).
#' @param remove_empty_sample logic. remove empty samples, default FALSE.
#' @param remove_macromolecular_gene logic. remove macromolecular genes that detemined by macromolecular_threshold, default FALSE.
#' @param mark_macromolecular_gene logic. mark macromolecular genes with macromolecular_mark.
#' @param macromolecular_mark character. used to mark macromolecular gene if set remove_macromolecular_gene as FALSE, default *.
#' @param macromolecular_threshold numeric. macromolecular threshold that determine whether macromolecular, default 10000bp.
#' @param top_annotation HeatmapAnotation. top annotation for oncoplot, default NULL, if NULL, total alteration burden gona be used as top annotaion.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param output_dir character. output plot directory.
#' @param ... list. other params used for oncoprint.
#'
#' @return data.frame. mutations.
#' @export
#'




complex_oncoplot_o<-function(snp_indels,
                           for_row="Hugo_Symbol",
                           for_column="Tumor_Sample_Barcode",
                           for_value="Variant_Classification",
                           cnv=TRUE,
                           cnv_types=c("Amplification","Deletion"),
                           annotate_source=F,
                           source_column="Source",
                           source_label=c("MDL"="*"),
                           sample_order=NULL,
                           selected_genes,
                           gene_order,
                           mutsigCV2_sig_genes=NA,
                           mutsigCV2_sig_p_threshold=0.05,
                           multi_hit=FALSE,
                           multi_hit_col="black",
                           top_filter=20L,
                           col,
                           bg_color="#dcdedc",
                           remove_empty_sample=F,
                           remove_macromolecular_gene=F,
                           mark_macromolecular_gene=F,
                           macromolecular_mark="*",
                           macromolecular_threshold=10000,
                           top_annotation=NULL,
                           top_anno_sample_order=NULL,
                           output_dir="./",
                           ...){

  f<-function(l,multi_hit) {
    if(length(l)==1){
      return(as.character(l))
    } else{
      if(multi_hit){
        return ("Multi_Hit")
      } else{
        return(paste(l,collapse=";"))
      }
    }}
  DNA_samples<-top_anno_sample_order
  #only mutsigcv significant gene mutation allowed for following analysis
  snp_indels<-snp_indels[snp_indels$Variant_Classification!="Silent" & (!is.na(snp_indels$Variant_Classification)),]
  snp_indels_o<-snp_indels
  #prepare tumor mutation burden data
  # make sure rownames of tmb should be consistent with colnames top_snp_indels

  tmb=snp_indels %>% dplyr::group_by(Tumor_Sample_Barcode,Variant_Classification) %>% dplyr::summarise(freq=dplyr::n()) %>% tidyr::pivot_wider(id_cols = "Tumor_Sample_Barcode", names_from="Variant_Classification",values_from="freq",values_fill=0)
  tmb<-as.data.frame(tmb)
  rownames(tmb)<-tmb[["Tumor_Sample_Barcode"]]
  tmb<-tmb[,-1]
  tmb<-tmb[,colnames(tmb)!="NA"]
  tmb<-tmb[DNA_samples,]
  #tmb<-tmb[colnames(top_snp_indels),]

  #extract and prepare snp indel data
  snp_indels<-snp_indels[,c(for_row,for_column,for_value)]
  colnames(snp_indels)<-c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification")
  snp_indels[["Tumor_Sample_Barcode"]]<-factor(snp_indels[["Tumor_Sample_Barcode"]])
  if(!cnv){
    snp_indels<-snp_indels[! snp_indels$Variant_Classification %in% cnv_types,]
    col=col[!names(col) %in% cnv_types]
  }

  # prepare top n snp indel data and then data for oncoprint
  snp_indels_p<-as.data.frame(tidyr::pivot_wider(snp_indels,id_cols = "Hugo_Symbol",names_from="Tumor_Sample_Barcode",values_from="Variant_Classification",values_fn=pryr::partial(f,multi_hit = multi_hit)))
  rownames(snp_indels_p)<-snp_indels_p[[1]]
  snp_indels_p<-as.matrix(snp_indels_p[,-1,drop=F])
  snp_indels_p[is.na(snp_indels_p)]<-""
  if(!is.na(match("None",rownames(snp_indels_p)))){
    snp_indels_p<-snp_indels_p[-(match("None",rownames(snp_indels_p))),]
  }
  snp_indels_p<-snp_indels_p[,DNA_samples,drop=F]

  if(missing(selected_genes)){
    selected_genes<-mutsigCV2_sig_genes[["gene"]][mutsigCV2_sig_genes[["p"]]<=mutsigCV2_sig_p_threshold]
    macromolecular_genes<-mutsigCV2_sig_genes[["gene"]][mutsigCV2_sig_genes[["codelen"]]>=macromolecular_threshold]
  } else{
    remove_macromolecular_gene=F
    mark_macromolecular_gene=F
  }
  snp_indels_p<-snp_indels_p[rownames(snp_indels_p) %in% selected_genes,,drop=F]

  top_filter<-min(c(nrow(snp_indels_p),top_filter))
  if(is.integer(top_filter)){
    top_mutated_genes<-rownames(snp_indels_p)[order(rowSums(snp_indels_p!=""),decreasing = T)][1:top_filter]
  } else{
    top_mutated_genes<-rownames(snp_indels_p)[(rowSums(snp_indels_p!="")/ncol(snp_indels_p))>=top_filter]
  }
  snp_indels_p<-snp_indels_p[top_mutated_genes,]

  if(remove_macromolecular_gene){snp_indels_p<-snp_indels_p[setdiff(rownames(snp_indels_p),macromolecular_genes),]}
  if(mark_macromolecular_gene) {rownames(snp_indels_p)[rownames(snp_indels_p) %in% macromolecular_genes] <-paste(macromolecular_mark,rownames(snp_indels_p)[rownames(snp_indels_p) %in% macromolecular_genes],sep="")}



  # if remove empty sample
  if(remove_empty_sample){
    empty_samples<-colnames(snp_indels_p)[colSums(snp_indels_p!="")==0]
    DNA_samples_<-setdiff(DNA_samples,empty_samples)
    anno_DNA<-anno_DNA[match(DNA_samples_,DNA_samples)]
    tmb<-tmb[DNA_samples_,,drop=F]
    snp_indels_p<-snp_indels_p[,DNA_samples_,drop=F]
    sample_order<-sample_order[-match(empty_samples,sample_order)]
  }


  if(annotate_source){
    source_p<-as.data.frame(tidyr::pivot_wider(snp_indels_o,id_cols = "Hugo_Symbol",names_from="Tumor_Sample_Barcode",values_from=source_column,values_fn=pryr::partial(f,multi_hit = FALSE)))
    rownames(source_p)<-source_p[[1]]
    source_p<-as.matrix(source_p[,-1,drop=F])
    source_p[is.na(source_p)]<-""
    if(!is.na(match("None",rownames(source_p)))){
      source_p<-source_p[-(match("None",rownames(source_p))),]
    }

    if(remove_macromolecular_gene){
      source_p<-source_p[setdiff(rownames(source_p),macromolecular_genes),]
    }
    if(mark_macromolecular_gene){
      rownames(source_p)[rownames(source_p) %in% macromolecular_genes] <-paste(macromolecular_mark,rownames(source_p)[rownames(source_p) %in% macromolecular_genes],sep="")
    }
    source_p<-source_p[rownames(snp_indels_p),colnames(snp_indels_p)]
  }




  # if remove macromolecular genes
  # extract macrocolecular genes
  #hs<-org.Hs.eg.db
  #hg19<-TxDb.Hsapiens.UCSC.hg19.knownGene
  #cds<-transcriptLengths(hg19,with.cds_len = T)
  #cds<-cds[(cds$cds_len!=0) & (!is.na(cds$gene_id)),]
  #max_cds<-cds %>% dplyr::group_by(gene_id) %>% dplyr::summarise(max_cds_len=max(cds_len)) %>% arrange(desc(max_cds_len))
  #macromolecular_genes<-max_cds$gene_id[max_cds$max_cds_len>=macromolecular_theshold]
  #macromolecular_genes<-AnnotationDbi::select(hs,keys=macromolecular_genes,columns=c("ENTREZID", "SYMBOL"),keytype = "ENTREZID")$SYMBOL

  #if col is na
  if(missing(col)){
    palette19<-c("#800000","#9A6324","#808000","#469990","#000075","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#42d4f4","#4363d8","#911eb4","#f032e6","#fabed4","#ffd8b1","#fffac8","#aaffc3","#dcbeff")
    col<-sample(palette19,ncol(tmb),replace = F)
    names(col)<-colnames(tmb)
  } else{
    col<-col[colnames(tmb)]
  }

  testit::assert("row names of tmb should be same as col names of snp_indels_p ", all(rownames(tmb)==colnames(snp_indels_p)))
  testit::assert("names of color should be same as col names as tmb",all(names(col)==colnames(tmb)))

  # plot heatmap annotation of tmb

  tmb_anno<-ComplexHeatmap::HeatmapAnnotation(tmb=ComplexHeatmap::anno_barplot(tmb,gp=grid::gpar(fill=col,col=col),border=F,bar_width = 0.85,height=grid::unit(4,"cm")),annotation_name_side="left")
  #plot oncoplot
  if(multi_hit){
    col=c(col,"Multi_Hit"=multi_hit_col)
  }

  onco_p<-ComplexHeatmap::oncoPrint(snp_indels_p,
                                    alter_fun=function(x,y,w,h,v){
                                      if(multi_hit){
                                        grid::grid.rect(x,y,w*0.9,h*0.9,gp=grid::gpar(fill=bg_color,col=NA))
                                        for(variant_c in names(v)){
                                          if(v[variant_c]){
                                            grid::grid.rect(x, y, w*0.9, h*0.9,gp = grid::gpar(fill = col[variant_c],col=NA))
                                          }
                                        }
                                      } else{
                                        if(cnv){
                                          mutation_v<-v[! names(v) %in% cnv_types ]
                                          cnv_v<-v[names(v) %in% cnv_types]
                                        } else{
                                          mutation_v<-v
                                        }
                                        cell_height=grid::unit(1/nrow(snp_indels_p),"npc")
                                        cell_bar_margin=0.05*cell_height
                                        stacks_in_cell=ifelse(sum(mutation_v)==0,1,sum(mutation_v))
                                        y_=y-cell_height*0.5+cell_bar_margin+grid::unit(1/nrow(snp_indels_p)/stacks_in_cell*0.9/2,"npc")
                                        grid::grid.rect(x,y,w*0.9,h*0.9,gp=grid::gpar(fill=bg_color,col=NA))
                                        for(variant_c in names(mutation_v)){
                                          if(mutation_v[variant_c]){
                                            grid::grid.rect(x, y_, w*0.9, h*0.9*(1/stacks_in_cell),gp = grid::gpar(fill = col[variant_c], col = NA))
                                            y_=y_+h*0.9*(1/stacks_in_cell)
                                          }
                                        }
                                        if(cnv){
                                          for(cnv_c in names(cnv_v)){
                                            if(cnv_v[cnv_c]){
                                              grid::grid.rect(x, y, w*0.5, h*0.9,gp = grid::gpar(fill = col[cnv_c], col = NA))
                                            }
                                          }
                                        }
                                      }
                                    },
                                    top_annotation = c(tmb_anno,top_annotation),
                                    pct_side = "right",
                                    row_names_side = "left",
                                    col=col,
                                    show_column_names = T,
                                    column_order=sample_order,...)
  method_lgd = ComplexHeatmap::Legend(labels = "MDL", title = "Support", type = "points", pch = 8, legend_gp = grid::gpar(col = "black"), background = bg_color)
  grDevices::svg(file.path(output_dir,"complex_oncoplot.svg"),width=8,height=12)
  ComplexHeatmap::draw(onco_p,annotation_legend_list=list(method_lgd))
  if(annotate_source){
    snp_indels_p_<-snp_indels_p[ComplexHeatmap::row_order(onco_p),ComplexHeatmap::column_order(onco_p)]
    source_p<-source_p[ComplexHeatmap::row_order(onco_p),ComplexHeatmap::column_order(onco_p)]
    heatmap_oncoprint_body<-gsub("heatmap_","",ComplexHeatmap::list_components()[grepl("heatmap_oncoPrint",ComplexHeatmap::list_components())])
    x_step=1/ncol(source_p)
    x_start=1/ncol(source_p)/2
    y_step=1/nrow(source_p)
    y_start=1/nrow(source_p)/2

    for(i in 1:nrow(source_p)){
      for(j in 1:ncol(source_p)){
        #if(top_source[i,j]==names(source_label)){
        #  ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.rect(x_start+(j-1)*x_step,1-y_start-(i-1)*y_step,x_step*0.9,y_step*0.9,gp=grid::gpar(fill=bg_color,col=col[top_snp_indels_[i,j]],lwd=2)))
        #}
        if(all(grepl("MDL",unlist(strsplit(source_p[i,j],";"))))){
          ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.rect(x_start+(j-1)*x_step,1-y_start-(i-1)*y_step,x_step*0.9,y_step*0.9,gp=grid::gpar(fill=bg_color,col=col[snp_indels_p_[i,j]],lwd=2)))
        }
        if(grepl(names(source_label),source_p[i,j])){
          ComplexHeatmap::decorate_heatmap_body(heatmap_oncoprint_body,grid::grid.text(source_label,x=x_start+(j-1)*x_step,y=1-y_start*1.3-(i-1)*y_step))
        }
      }
    }
  }
  grDevices::dev.off()
  return(snp_indels_p)
}


#' annotate_genomic_mutations_with_panel_info
#'
#' annotate genomic mutations (WES/WGS) with provided panel information.
#'
#' @param genomic_mutations data.frame. wes mutations data.frame.
#' @param mutation_chromosome_col chararcter. genomic mutation column for chromosome.
#' @param mutation_start_col chararcter. genomic mutation column for start.
#' @param mutation_end_col chararcter. genomic mutation column for end.
#' @param mutation_gene_col chararcter. genomic mutation column for genename.
#' @param panel_info data.frame. information  of panel used for  MDL mutation.
#' @param panel_id_col character. column for panel id
#' @param panel_chromosome_col chararcter. panel information column for chromosome.
#' @param panel_start_col chararcter. panel information column for start.
#' @param panel_end_col chararcter. panel information column for end.
#' @param panel_gene_col chararcter. panel information column for genename.
#' @param used_panels character. panel ids used for current cohort
#'
#' @return data.frame.
#' @export
#'
annotate_genomic_mutations_with_panel_info<-function(genomic_mutations,mutation_chromosome_col="Chromosome",mutation_start_col="Start_position",mutation_end_col="End_position",mutation_gene_col="gene",
                                                     panel_info,panel_id_col="document_type",panel_chromosome_col="chromosome",panel_start_col="start_position",panel_end_col="end_position",panel_gene_col="gene",
                                                     used_panels=c("HP MD Solid Tumor Genomic Assay V1 Report","HP MD Solid Tumor Genomic Assay-DNA 2018 Report","HP MD Liquid Biopsy Panel V1 Report")){
  library(plyranges)
  genomic_mutations_gr=df2granges(df=genomic_mutations,genome="hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,
                                  seqnames_col = mutation_chromosome_col,start_col = mutation_start_col,end_col = mutation_end_col,
                                  meta_cols = setdiff(colnames(genomic_mutations),c(mutation_chromosome_col,mutation_start_col,mutation_end_col)))

  for(panel in used_panels){
    if(panel %in% unique(panel_info[[panel_id_col]])){
      sub_panel_info<-panel_info[panel_info[[panel_id_col]]==panel,]
      in_panel_genelist_ind<-which(genomic_mutations[[mutation_gene_col]] %in% unique(sub_panel_info[[panel_gene_col]]))
      genomic_mutations[[panel]]<-"out of panel genelist"
      genomic_mutations[[panel]][in_panel_genelist_ind]<-"in panel genelist"
      sub_panel_info=sub_panel_info[!(is.na(sub_panel_info[[panel_chromosome_col]]) | is.na(sub_panel_info[[panel_start_col]]) | is.na(sub_panel_info[[panel_end_col]])),]
      if(nrow(sub_panel_info)==0){
        genomic_mutations[[panel]]<-paste(genomic_mutations[[panel]],"unlucky no postion information available",sep=";")
      } else {
        sub_panel_gr<-df2granges(df=sub_panel_info,genome="hg19",seqlevelsStyle = "NCBI",simplified = T,xy=T,
                                 seqnames_col = panel_chromosome_col,start_col = panel_start_col,end_col = panel_end_col,
                                 meta_cols = setdiff(colnames(sub_panel_info),c(panel_chromosome_col,panel_start_col,panel_end_col)))

        covered_ind<-queryHits(findOverlaps(genomic_mutations_gr, sub_panel_gr))
        noncovered_ind<-setdiff(in_panel_genelist_ind,covered_ind)
        genomic_mutations[[panel]][noncovered_ind]<-paste(genomic_mutations[[panel]][noncovered_ind],"; but not covered by panel sequencing.",sep="")
        genomic_mutations[[panel]][covered_ind]<-paste(genomic_mutations[[panel]][covered_ind],"; and covered by panel sequencing.",sep="")
      }
    } else {
      genomic_mutations[[panel]]<-"no panel information provided now!"
    }
  }
  return(genomic_mutations)
}


#' annotate_mdl_mutations_with_timepoints
#'
#' annotate mdl mutations with timepoints data. mark sample as well as mutations with time points(like "Baseline","TP2",...) according to the minimal gap between collect time.
#'
#' @param mdl_mutations data.frame. mdl mutation data.frame.
#' @param mdl_mutation_patient_id_col character. mdl mutation column for patient id
#' @param mdl_mutation_collect_time_col character. mdl mutation column for collect time
#' @param mdl_mutation_collect_date_format character. mdl mutation collect time column date format.
#' @param timepoint_info data.frame. mdl sample timepoint information
#' @param timepoint_info_patient_id_col character. mdl sample timepoint information column for patient id
#' @param timepoint_info_timepoint_col character. mdl sample timepoint information column for timepoint
#' @param timepoint_info_collect_time_col character. mdl sample timepoint information column for collect time.
#' @param timepoint_info_collect_date_format character. time format for timepoint info collect time.
#'
#' @return data.frame. timepoint annotated mdl mutations
#' @export
#'

annotate_mdl_mutations_with_timepoints<-function(mdl_mutations,mdl_mutation_patient_id_col="Patient_ID",mdl_mutation_collect_time_col,mdl_mutation_collect_date_format="%Y-%m-%d",timepoint_info,timepoint_info_patient_id_col="Patient_ID",timepoint_info_timepoint_col,timepoint_info_collect_time_col,timepoint_info_collect_date_format="%m/%d/%Y"){
  stopifnot("Patient ID not found in mutations"=(mdl_mutation_patient_id_col %in% colnames(mdl_mutations)))
  stopifnot("sample info contains no collect date"=(!any(is.na(timepoint_info[[timepoint_info_collect_time_col]]))))
  if(!c("Timepoint") %in% colnames(mdl_mutations)){mdl_mutations[["Timepoint"]]=""}
  if(!c("Comment_Timepoint") %in% colnames(mdl_mutations)){mdl_mutations[["Comment_Timepoint"]]=paste(mdl_mutation_collect_time_col,mdl_mutations[[mdl_mutation_collect_time_col]],sep=":")}
  results<-lapply(unique(mdl_mutations[[mdl_mutation_patient_id_col]]),function(patient){
    print(patient)
    patient_mutations<-mdl_mutations[mdl_mutations[[mdl_mutation_patient_id_col]]==patient,,drop=F]
    if(patient %in% timepoint_info[[timepoint_info_patient_id_col]]) {
      collect_time=timepoint_info[timepoint_info[[timepoint_info_patient_id_col]]==patient,c(timepoint_info_timepoint_col,timepoint_info_collect_time_col),drop=F]
      collect_time[[timepoint_info_timepoint_col]][is.na(collect_time[[timepoint_info_timepoint_col]])]<-"Unknown"
      collect_time<-collect_time[!duplicated(collect_time),,drop=F]
      collect_time_<-collect_time[[timepoint_info_collect_time_col]]
      names(collect_time_)<-collect_time[[timepoint_info_timepoint_col]]
      collect_time_<-collect_time_[order(names(collect_time_))]
      collect_time_<-gsub(" [AP]M.*","",collect_time_,perl=T)
      patient_mutations$Timepoint<-names(collect_time_)[apply(abs(outer(as.Date(collect_time_,format=timepoint_info_collect_date_format),as.Date(gsub(" [AP]M.*","",patient_mutations[[mdl_mutation_collect_time_col]],perl=T),format=mdl_mutation_collect_date_format),FUN="-")),2,which.min)]
      patient_mutations$Comment_Timepoint<-paste(patient_mutations$Comment_Timepoint,paste(paste(names(collect_time_),unname(collect_time_),sep=":"),collapse = ";"),sep=" | ")
    } else {
      patient_mutations$Comment_Timepoint<-paste(patient_mutations$Comment_Timepoint,paste(patient,"not listed in DNA sample info maybe not yet sequenced!",sep=":"),sep=" | ")
    }
    return(patient_mutations)
  })

  mdl_mutations<-as.data.frame(do.call(rbind,results))
  mdl_sample_info<-mdl_mutations[,c(mdl_mutation_patient_id_col,"Timepoint",mdl_mutation_collect_time_col)]
  #mdl_sample_info$Sample_ID<-paste(mdl_sample_info[[mdl_mutation_patient_id_col]],mdl_sample_info[["Timepoint"]],sep="_")
  mdl_sample_info<-mdl_sample_info[!duplicated(mdl_sample_info),]
  Comment_Timepoint<-apply(mdl_mutations,1,function(mdl_mutations_r){
    sub_mdl_sample_info<-mdl_sample_info[mdl_sample_info[[mdl_mutation_patient_id_col]]==mdl_mutations_r[mdl_mutation_patient_id_col] & mdl_sample_info[["Timepoint"]]==mdl_mutations_r["Timepoint"],,drop=F]
    if(nrow(sub_mdl_sample_info)!=1){
      return(paste(mdl_mutations_r["Comment_Timepoint"],"\n* Caution: multiple collect_time mapped to same timepoint \n",paste(sub_mdl_sample_info[[mdl_mutation_collect_time_col]],collapse=";"),sep=""))
    } else {
      return(mdl_mutations_r["Comment_Timepoint"])
    }
  })
  mdl_mutations$Comment_Timepoint<-Comment_Timepoint
  return(mdl_mutations)
}


#' annotate_mdl_mutations_with_mutationdb
#'
#' annotat mdl mutation with known mutation database.
#'
#' @param mdl_mutations data.frame. mdl mutations data.frame.
#' @param mdl_mutations_gene_col chararcter. mdl mutation column for genename.
#' @param mdl_mutations_aa_col chararcter. mdl mutation column for amino acid.
#' @param with_oncokb logic. annotate with oncoKB or not.
#' @param oncokb_db data.frame. oncokb database
#' @param oncokb_gene_col character. oncoKB column for genename
#' @param with_cosmic logic. annotate with cosmic or not.
#' @param cosmic_db data.frame. cosmic database.
#' @param cosmic_gene_col character. cosmic column for genename
#' @param cosmic_aa_col character. cosmic column for amino acid
#'
#' @return data.frame. mdl mutation with known mutation database.
#' @export
#'

annotate_mdl_mutations_with_mutationdb<-function(mdl_mutations,mdl_mutations_gene_col="gene_name",mdl_mutations_aa_col="protein_description",with_oncokb=T,oncokb_db,oncokb_gene_col="Hugo Symbol",with_cosmic=T,cosmic_db,cosmic_gene_col="GENE_NAME",cosmic_aa_col="Mutation AA"){
  if(with_oncokb){
    mdl_mutations[["oncokb"]]<-ifelse(mdl_mutations[[mdl_mutations_gene_col]] %in% oncokb_db[[oncokb_gene_col]],"Yes","No")
  }
  if(with_cosmic){
    mdl_mutations$mdl_mutations_protein_codon<-paste(mdl_mutations[[mdl_mutations_gene_col]],gsub("\\*.*","",mdl_mutations[[mdl_mutations_aa_col]]),sep=":")
    cosmic_protein_codon<-unique(paste(cosmic_db[[cosmic_gene_col]],gsub("\\*.*","",cosmic_db[[cosmic_aa_col]]),sep=":"))
    mdl_mutations[["cosmic"]]<-ifelse(mdl_mutations$mdl_mutations_protein_codon %in% cosmic_protein_codon,"Yes","No")
  }
  return(mdl_mutations)
}


#' get_min_multihit_distance
#'
#' calculte the minimal distance for multiple mutations
#'
#' @param codons numeric. a numeric vector describ the protein codons
#'
#' @return numeric. minimal distance among codons.
#' @export
#'
get_min_multihit_distance<-function(codons){
  if(length(codons)==1){
    return(NA)
  } else {
    return(sapply(codons,function(codon){min(abs(codon-codons)[-which.min(abs(codon-codons))],na.rm=T)}))
  }
}
#' harmonize_mutations
#'
#' annotate genomic mutations(WES/WGS) and mdl mutations
#' merge them if they are consistent
#'
#' @param wes_sample_info data.frame. WES sample information.
#' @param wes_sample_id_col character. WES sample information column for sample id.
#' @param wes_patient_id_col character. WES sample information column for patient id.
#' @param wes_timepoint_col character. WES sample information column for timepoint.
#' @param wes_sample_collect_col character. WES sample information column for collect time.
#' @param wes_mutations data.frame. WES mutations.
#' @param wes_mutation_sample_id_col character. WES mutation column for sample id.
#' @param wes_gene_col  character. WES mutation column for genename.
#' @param wes_chromosome_col  character. WES mutation column for chromosome.
#' @param wes_start_col  character. WES mutation column for start.
#' @param wes_end_col  character. WES mutation column for end.
#' @param wes_keep_columns  character. WES mutation column kept for the output.
#' @param wes_seq_intervals grangeslist. intervals for WES sequencing.
#' @param wgs_sample_info data.frame. WGS sample information.
#' @param wgs_sample_id_col character. WGS sample information column for sample id.
#' @param wgs_patient_id_col character. WGS sample information column for patient id.
#' @param wgs_timepoint_col character. WGS sample information column for timepoint.
#' @param wgs_sample_collect_col character. WGS sample information column for collect time.
#' @param wgs_mutations data.frame. WGS mutations.
#' @param wgs_mutation_sample_id_col character. WGS mutation column for sample id.
#' @param wgs_gene_col  character. WGS mutation column for genename.
#' @param wgs_chromosome_col character. WGS mutation column for chromosome.
#' @param wgs_start_col character. WGS mutation column for start.
#' @param wgs_end_col character. WGS mutation column for end.
#' @param wgs_keep_columns character. WGS mutation column kept for the output.
#' @param mdl_sample_info data.frame. mdl sample information.
#' @param mdl_sample_id_col character. mdl sample information column for sample id.
#' @param mdl_patient_id_col character. mdl sample information column for patient id.
#' @param mdl_timepoint_col character. mdl sample information column for timepoint.
#' @param mdl_sample_collect_time_col character. mdl sample information column for collect time.
#' @param mdl_sample_collect_date_format character. timepoint information collect time column date format.
#' @param mdl_mutations data.frame. mdl mutations.
#' @param mdl_mutation_sample_id_col character. mdl mutation column for sample id.
#' @param mdl_mutation_patient_id_col character. mdl mutation column for patient id.
#' @param mdl_mutation_collect_time_col character. mdl sample information column for collect time.
#' @param mdl_mutation_collect_date_format character. mdl mutation collect time column date format.
#' @param mdl_mutation_gene_col character. mdl sample information column for genename.
#' @param mdl_mutation_indication_col character. mdl sample information column for indication(like "wildtype","mutation").
#' @param mdl_mutations_aa_col character. mdl sample information column for amino acid
#' @param mdl_mutation_panel_col character. mdl sample information column for panel id.
#' @param mdl_mutation_keep_columns character. mdl mutation column kept for the output.
#' @param panel_info data.frame. panel information.
#' @param panel_id_col character. panel information column for panel id.
#' @param panel_chromosome_col character. panel information column for chromosome.
#' @param panel_start_col character. panel information column for start.
#' @param panel_end_col character. panel information column for end.
#' @param panel_gene_col character. panel information column for genename.
#' @param with_oncokb logic. annotated with oncoKB or not.
#' @param oncokb_db data.frame. oncoKB database.
#' @param oncokb_gene_col character. oncoKB column for genename
#' @param with_cosmic logic. annotated with cosmic or not.
#' @param cosmic_db data.frame. cosmic database.
#' @param cosmic_gene_col character. cosmic column for genename
#' @param cosmic_aa_col character. cosmic column for amino acid
#' @param mdl_columns2wes_columns character. a named character used for mapping column names of mdl mutations and wes(wgs) mutations
#' @param combined logic. if TRUE, consistent mutations between mdl mutations and wes/wgs mutations will be combined into one mutation.
#'
#' @return data.frame. annotated and merged mutations(WES/WGS,MDL)
#' @export
#'
harmonize_mutations<-function(wes_sample_info,wes_sample_id_col,wes_patient_id_col,wes_timepoint_col,wes_sample_collect_col,
                              wes_mutations=NULL,wes_mutation_sample_id_col,wes_gene_col='gene',wes_chromosome_col="Chromosome",wes_start_col="Start_Position",wes_end_col="End_Position",wes_keep_columns=c("Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Protein_Change","i_TumorVAF_WU","cdna"),
                              wes_seq_intervals,
                              wgs_sample_info,wgs_sample_id_col,wgs_patient_id_col,wgs_timepoint_col,wgs_sample_collect_col,
                              wgs_mutations=NULL,wgs_mutation_sample_id_col,wgs_gene_col='gene',wgs_chromosome_col="Chromosome",wgs_start_col="Start_Position",wgs_end_col="End_Position",wgs_keep_columns=c("transcript","aaannotation","func.knowngene","exonicfunc.knowngene","cdna","chrom","start","end","ref_allele","alt_allele","t_ref_count","t_alt_count","n_ref_count","n_alt_count"),
                              mdl_sample_info,mdl_sample_id_col="Sample_ID",mdl_patient_id_col="Patient_ID",mdl_timepoint_col,mdl_sample_collect_time_col,mdl_sample_collect_date_format="%Y-%m-%d",
                              mdl_mutations=NULL,mdl_mutation_sample_id_col=NULL,mdl_mutation_patient_id_col="Patient_ID",mdl_mutation_collect_time_col,mdl_mutation_collect_date_format="%Y-%m-%d",mdl_mutation_gene_col="gene",mdl_mutation_indication_col="indication_name",mdl_mutations_aa_col="protein_description",mdl_mutation_panel_col="document_subtype_description",mdl_mutation_keep_columns=c("aa_change","indication_name","nucleotide_change" ,"reference_nucleotide","alternate_nucleotide","nucleotide_change_type" ,"nucleotide_change_subtype"),
                              panel_info,panel_id_col="document_type",panel_chromosome_col="Chromosome",panel_start_col="Start_Position",panel_end_col="End_Position",panel_gene_col="gene",
                              with_oncokb=T,oncokb_db,oncokb_gene_col="Hugo Symbol",with_cosmic=T,cosmic_db,cosmic_gene_col="GENE_NAME",cosmic_aa_col="Mutation AA",
                              mdl_columns2wes_columns=c("aa_change"="Protein_Change","nucleotide_change"= "cdna","reference_nucleotide"="Reference_Allele","alternate_nucleotide"="Tumor_Seq_Allele2"),
                              combined=F
) {
  panels<-unique(mdl_mutations[[mdl_mutation_panel_col]])

  #annotate mdl_mutations with timepoint
  if(!is.null(mdl_mutations)){
    if(is.null(mdl_mutation_sample_id_col)){
      mdl_mutations[["Sample_ID"]]<-paste(mdl_mutations[[mdl_mutation_patient_id_col]],mdl_mutations[[mdl_mutation_collect_time_col]],sep="_")
      mdl_mutation_sample_id_col<-"Sample_ID"
    }
    #multi_hit annotation
    #mdl_mutations_wt<-mdl_mutations[mdl_mutations[[mdl_mutation_indication_col]]=="wildtype",]
    #mdl_mutations_wt[["mdl_multihits"]]<-NA
    #mdl_mutations<-mdl_mutations[mdl_mutations[[mdl_mutation_indication_col]]!="wildtype",]
    mdl_mutation_ids<-paste(mdl_mutations[[mdl_mutation_sample_id_col]],mdl_mutations[[mdl_mutation_gene_col]],sep=":")
    multihit_mdl_mutation_ids<-names(table(mdl_mutation_ids)[table(mdl_mutation_ids)>1])
    mdl_mutations[["mdl_multihits"]]<-ifelse(mdl_mutation_ids %in% multihit_mdl_mutation_ids,"Yes","No")
    #mdl_mutations<-rbind(mdl_mutations_wt,mdl_mutations_mut)
    mdl_mutation_keep_columns<-c(mdl_mutation_keep_columns,"mdl_multihits")

    #annotation with mutation database
    mdl_mutations<-annotate_mdl_mutations_with_mutationdb(mdl_mutations,mdl_mutations_gene_col="gene_name",mdl_mutations_aa_col="protein_description",with_oncokb=T,oncokb_db=oncokb_cancer_genes_info,oncokb_gene_col="Hugo Symbol",with_cosmic=T,cosmic_db=cosmic_mutations,cosmic_gene_col="GENE_NAME",cosmic_aa_col="Mutation AA")
    mdl_mutation_keep_columns<-c(mdl_mutation_keep_columns,c("oncokb","cosmic"))
    mdl_mutations<-annotate_mdl_mutations_with_timepoints(mdl_mutations = mdl_mutations,mdl_mutation_patient_id_col = mdl_mutation_patient_id_col,mdl_mutation_collect_time_col = mdl_mutation_collect_time_col,mdl_mutation_collect_date_format=mdl_mutation_collect_date_format,
                                                          timepoint_info = mdl_sample_info,timepoint_info_patient_id_col = mdl_patient_id_col,timepoint_info_timepoint_col = mdl_timepoint_col,timepoint_info_collect_time_col = mdl_sample_collect_time_col,timepoint_info_collect_date_format=mdl_sample_collect_date_format)


    mdl_mutations<-data.frame(Source=mdl_mutations[[mdl_mutation_panel_col]],
                              Sample_ID=mdl_mutations[[mdl_mutation_sample_id_col]],
                              Patient_ID=mdl_mutations[[mdl_mutation_patient_id_col]],
                              Timepoint=mdl_mutations[["Timepoint"]],
                              Gene=mdl_mutations[[mdl_mutation_gene_col]],
                              Variant_Classification=mdl_mutations[["Variant_Classification"]],
                              Variant_Type=mdl_mutations[["Variant_Type"]],
                              Collect_Time=as.character(mdl_mutations[[mdl_mutation_collect_time_col]]),
                              mdl_mutations[,mdl_mutation_keep_columns],
                              Comment_Timepoint=mdl_mutations[["Comment_Timepoint"]]
    )
    colnames(mdl_mutations)[match(names(mdl_columns2wes_columns),colnames(mdl_mutations))]<-unname(mdl_columns2wes_columns)
  } else{
    mdl_mutations<-NULL
  }
  #annotate wes_mutations with panel_info
  if(!is.null(wes_mutations)){
    wes_mutations<-annotate_genomic_mutations_with_panel_info(genomic_mutations = wes_mutations,mutation_chromosome_col = wes_chromosome_col,mutation_start_col = wes_start_col,mutation_end_col = wes_end_col,mutation_gene_col=wes_gene_col,
                                                              panel_info = panel_info,panel_id_col = panel_id_col,panel_chromosome_col = panel_chromosome_col,panel_start_col = panel_start_col,panel_end_col = panel_end_col,panel_gene_col=panel_gene_col,
                                                              used_panels = panels)
    wes_mutations<-merge(wes_mutations,wes_sample_info[,c(wes_sample_id_col,wes_patient_id_col,wes_sample_collect_col,wes_timepoint_col)],by.x=wes_mutation_sample_id_col,by.y=wes_sample_id_col,all.x=T)
    wes_mutations[["Source"]]<-"WES"
    wes_mutations<-data.frame(Source=wes_mutations[["Source"]],
                              Sample_ID=wes_mutations[[wes_mutation_sample_id_col]],
                              Patient_ID=wes_mutations[[wes_patient_id_col]],
                              Timepoint=wes_mutations[[wes_timepoint_col]],
                              Gene=wes_mutations[[wes_gene_col]],
                              Variant_Classification=wes_mutations[["Variant_Classification"]],
                              Variant_Type=wes_mutations[["Variant_Type"]],
                              Collect_Time=wes_mutations[[wes_sample_collect_col]],
                              wes_mutations[,c(wes_keep_columns,panels)]
    )
    #wes_mutations<-wes_mutations[,c("Source",wes_mutation_sample_id_col,wes_gene_col,wes_patient_id_col,wes_timepoint_col,panels,wes_keep_columns)]
  } else {
    wes_mutations<-NULL
  }

  #annotate wgs_mutations with panel_info
  if(!is.null(wgs_mutations)){
    wgs_mutations<-annotate_genomic_mutations_with_panel_info(genomic_mutations = wgs_mutations,mutation_chromosome_col = wgs_chromosome_col,mutation_start_col = wgs_start_col,mutation_end_col = wgs_end_col,
                                                              panel_info = panel_info,panel_id_col = panel_id_col,panel_chromosome_col = panel_chromosome_col,panel_start_col = panel_start_col,panel_end_col = panel_end_col,
                                                              used_panels = panels)
    wgs_mutations<-merge(wgs_mutations,wgs_sample_info[,c(wgs_sample_id_col,wgs_patient_id_col,wgs_timepoint_col)],by.x=wes_mutation_sample_id_col,by.y=wgs_sample_id_col,all.x=T)
    wgs_mutaions[["Source"]]<-"WGS"
    wgs_mutations<-data.frame(Source=wgs_mutations[["Source"]],
                              Sample_ID=wgs_mutations[[wgs_mutation_sample_id_col]],
                              Patient_ID=wgs_mutations[[wgs_patient_id_col]],
                              Timepoint=wgs_mutations[[wgs_timepoint_col]],
                              Gene=wgs_mutations[[wgs_mutation_gene_col]],
                              Variant_Classification=wgs_mutations[["Variant_Classification"]],
                              Variant_Type=wgs_mutations[["Variant_Type"]],
                              Collect_Time=wgs_mutations[[wgs_sample_collect_col]],
                              wgs_mutations[,c(wgs_mutation_keep_columns,panels)],
    )
    #wgs_mutations<-wgs_mutations[,c("Source",wgs_mutation_sample_id_col,wgs_gene_col,wgs_patient_id_col,wgs_timepoint_col,panels,wgs_keep_columns)]
  } else {
    wgs_mutations<-NULL
  }

  cat("before run the harmonize_muations, please check colnames of wes_mutations, wgs_mutations, mdl_mutations and various informations are consistent!")
  #patient_info<-data.frame(Patient_ID=patient_info[[patient_id_col]],
  #                         patient_info[,setdiff(colnames(patient_info),patient_id_col)])
  combined_mutations<-dplyr::bind_rows(mdl_mutations,wes_mutations,wgs_mutations)
  combined_mutations$codon<-as.numeric(sapply(combined_mutations[["Protein_Change"]],function(a){regmatches(a,regexpr("[0-9]+",a,perl=T))}))
  combined_mutations<-combined_mutations %>% dplyr::group_by(Patient_ID,Timepoint,Gene) %>% dplyr::mutate(GeneHits=dplyr::n(),min_multihit_dist=get_min_multihit_distance(codon))
  combined_mutations<-combined_mutations %>% dplyr::group_by(Patient_ID,Timepoint,Gene,codon) %>% dplyr::mutate(CodonHits=dplyr::n()) %>% dplyr::arrange(Patient_ID,Timepoint,dplyr::desc(GeneHits),dplyr::desc(CodonHits),codon)

  if(combined){
    #browser()
    combined_mutations<- combined_mutations %>% dplyr::group_by(Patient_ID,Timepoint,Gene,codon) %>% dplyr::mutate(N_codoncallplatform=length(unique(Source)))
    combined_mutations<- combined_mutations %>% dplyr::group_by(Patient_ID,Timepoint,Gene,cdna) %>% dplyr::mutate(N_cdnacallplatform=length(unique(Source)))
    singlecall_combined_mutations<-combined_mutations %>% dplyr::filter(N_codoncallplatform==1) %>% dplyr::mutate_if(is.numeric,as.character)
    singlecall_combined_mutations$Combined=FALSE
    multicall_combined_mutations<-combined_mutations %>% dplyr::filter(N_codoncallplatform>1) %>% dplyr::mutate_if(is.numeric,as.character) %>% dplyr::group_by(Patient_ID,Timepoint,Gene,codon) %>% dplyr::summarise(dplyr::across(dplyr::everything(),~paste(.x,collapse=";")))
    multicall_combined_mutations$Combined=TRUE
    if(nrow(multicall_combined_mutations)>0){
      combined_mutations<-dplyr::bind_rows(multicall_combined_mutations,singlecall_combined_mutations)
    }

  }

  combined_mutations$Variant_Type<-sapply(combined_mutations$Variant_Type,function(item){unlist(strsplit(item,";"))[1]})
  combined_mutations$Variant_Classification<-sapply(combined_mutations$Variant_Classification,function(item){unlist(strsplit(item,";"))[1]})
  combined_mutations$Support<- sapply(combined_mutations$Source,function(item){items<-unlist(strsplit(item,";"));ifelse(all(grepl("WES",items)),"WES",ifelse(all(!grepl("WES",items)),"MDL","MDL;WES"))})
  return(combined_mutations) # need manual curation
}


#' oncogenicpathway_tableplot
#'
#' oncogenic pathway table plot
#'
#' @param pathways list. no default, genesets for gsea, refer fgsea::gmtPathways
#' @param mutations data.frame. mutations data
#' @param sample_info data.frame. sample information.
#' @param gene_column character. mutation column for genename
#' @param sample_id_column character. sample information scolumn for sample id
#' @param pathway_order character. ordered pathway name.
#' @param sample_order character. ordered sample names.
#' @param top_anno HeatmapAnnotation. top annoation for oncogenic pathway table plot
#' @param top_anno_sample_order character. sample names for top annotation, consistent order with top annotaion
#' @param show_mutations logic. show mutation or not.
#' @param ... list. extra parameter passed to ComplexHeatmap::Heatmap
#'
#' @return Heatmap. oncogenic pathway table plot
#' @export
#'

oncogenicpathway_tableplot<-function(pathways,mutations,sample_info,gene_column="",sample_id_column="Sample_ID",pathway_order,sample_order=NULL,top_anno,top_anno_sample_order,show_mutations=F,...){
  stopifnot("gene_column is not found in columns of mutations"=gene_column %in% colnames(mutations))
  mutation_sample_order<-colnames(mutations)[-1]
  if(!missing(sample_info)){
    stopifnot("sample_id_column is not found in columns of sample_info"=sample_id_column %in% colnames(sample_info))
    stopifnot("mutations colnames don't match sample_info sample_id"=all(mutation_sample_order==sample_info[[sample_id_column]]))
  }
  pathway_sample_count_list<-lapply(pathways,function(pathway){
    sapply(mutations[,mutation_sample_order],function(sample_mutations){
      ifelse(length(intersect(mutations[[gene_column]][(!is.na(sample_mutations)) & (sample_mutations!="Silent")],pathway))==0,NA,length(intersect(mutations[[gene_column]][(!is.na(sample_mutations)) & (sample_mutations!="Silent")],pathway)))
    })
  })

  pathway_sample_mutation_list<-lapply(pathways,function(pathway){
    sapply(mutations[,mutation_sample_order],function(sample_mutations){
      paste(intersect(mutations[[gene_column]][(!is.na(sample_mutations)) & (sample_mutations!="Silent")],pathway),collapse=",\n")})
  })

  names(pathway_sample_count_list)<-sapply(names(pathway_sample_count_list),function(pathway){return(paste(pathway,"(",as.character(length(pathways[[pathway]])),")",sep=""))})
  names(pathway_sample_mutation_list)<-names(pathway_sample_count_list)
  pathway_sample_count_table<-as.data.frame(t(as.data.frame(pathway_sample_count_list,check.names=F)))[,mutation_sample_order]
  pathway_sample_mutation_table<-as.data.frame(t(as.data.frame(pathway_sample_mutation_list,check.names=F)))[,mutation_sample_order]
  top_anno<-top_anno[match(colnames(pathway_sample_count_table),top_anno_sample_order),]

  if(missing(pathway_order)){
    pathway_order<-rownames(pathway_sample_count_table)[order(rowSums(!is.na(as.matrix(pathway_sample_count_table))),decreasing = T)]
  }

  if(!is.null(sample_order)){
    if((!missing(sample_info)) & any(sample_order %in% colnames(sample_info))){
      sample_order_<-colnames(pathway_sample_count_table)[order(sample_info[[sample_order]],colSums(!is.na(as.matrix(pathway_sample_count_table))),decreasing=T)]
    } else {
      sample_order_<-sample_order
    }
  } else {
    sample_order_<-colnames(pathway_sample_count_table)[order(colSums(!is.na(as.matrix(pathway_sample_count_table))),decreasing=T)]
  }
  col=circlize::colorRamp2(c(0,5),c("#ffe6e6","#cc0000"))
  if(show_mutations){
    heatmap<-ComplexHeatmap::Heatmap(pathway_sample_count_table,name="Mutation",col=col,cluster_columns = F,cluster_rows = F,row_order = pathway_order,column_order = sample_order_,top_annotation = top_anno,
                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                       grid::grid.text(ifelse(is.na(pathway_sample_count_table[i, j]),"",paste(pathway_sample_count_table[i, j],pathway_sample_mutation_table[i, j],sep=":\n")), x, y,gp=grid::gpar(fontsize=5))
                     },
                     ...)
  } else {
    heatmap<-ComplexHeatmap::Heatmap(pathway_sample_count_table,name="Mutation",col=col,cluster_columns = F,cluster_rows = F,row_order = pathway_order,column_order = sample_order_,top_annotation = top_anno,
                     cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                       grid::grid.text(ifelse(is.na(pathway_sample_count_table[i, j]),"",pathway_sample_count_table[i, j]), x, y)
                     },
                     ...)
  }

  return(heatmap)
}


#' cancerhallmark_heatmap
#' heatmap of cancer hallmarks
#' @param assay data.frame. assay data, such as mutations.
#' @param assay_format character. format of assay, like "l" for long format
#' @param gene_id_column character. column for gene id.
#' @param sample_id_col character. column for sample id.
#' @param value_col character. column for value.
#' @param cancerhallmarks data.frame. cancer hallmarks.
#' @param top_anno HeatmapAnnotation. top annotation.
#' @param top_anno_sample_order character. sample order of top annotation.
#' @param hallmark_colors character. colors for hallmarks.
#' @param impact_colors character. colors for impact.
#' @param show_assay_color logic. whether show assay color individually.
#' @param assay_colors character. colors for assay value.
#' @param ...  list. parameter passed to Heatmap.
#'
#' @return list. contains hallmark_heatmap and underlying matrix.
#' @export
#'
cancerhallmark_heatmap<-function(assay,assay_format="l",gene_id_column="Hugo_Symbol",sample_id_col="Sample_ID",value_col,cancerhallmarks=CosmicCancerGeneCensusHallmarks,top_anno,top_anno_sample_order,hallmark_colors=NULL,impact_colors=NULL,show_assay_color=F,assay_colors=NULL,...){
  if(assay_format=="l"){
    assay<-tidyr::pivot_wider(assay,id_cols = gene_id_column,names_from = sample_id_col,values_from = value_col,values_fn = function(items){ifelse(length(items)>1,"Multihits",items)})
  } else {
    assay[[gene_id_column]]<-rownames(assay)
  }
  assay<-assay[assay[[gene_id_column]] %in% unique(cancerhallmarks$GENE_SYMBOL),]
  assay<-assay[,c(gene_id_column,top_anno_sample_order)]
  assay<-merge(cancerhallmarks,assay,by.x="GENE_SYMBOL",by.y=gene_id_column)

  n_samples<-rowSums(!is.na(assay[,(ncol(cancerhallmarks)+1):(ncol(assay))]))
  assay<-assay[order(assay$HALLMARK,n_samples,decreasing=T),]
  if(is.null(hallmark_colors)){
    hallmark_colors<-sample(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075'),length(unique(cancerhallmarks$HALLMARK)))
    names(hallmark_colors)<-unique(cancerhallmarks$HALLMARK)
  }
  if(is.null(impact_colors)){
    impact_colors<-c("suppresses"="red","promotes"="blue","promotes;suppresses"="yellow","unknown"="grey")
  }
  hallmark_row_anno<-ComplexHeatmap::rowAnnotation(Hallmark=assay$HALLMARK,Impact=assay$IMPACT,col=list("Hallmark"=hallmark_colors,"Impact"=impact_colors))
  heatmap_matrix<-as.matrix(assay[,(ncol(cancerhallmarks)+1):(ncol(assay))])
  rownames(heatmap_matrix)<-assay$GENE_SYMBOL
  if(show_assay_color){
    if(is.null(assay_colors)){
      assay_colors<-sample(c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075'),length(stats::na.omit(unique(as.vector(heatmap_matrix)))))
      names(assay_colors)<-stats::na.omit(unique(as.vector(heatmap_matrix)))
      if("Silent" %in% names(assay_colors)) {assay_colors["Silent"]="white"}
    }
    heatmap_plot<-ComplexHeatmap::Heatmap(heatmap_matrix,name="Variant_Impact",cluster_rows = F,cluster_columns = F,col=assay_colors,na_col = "grey",top_annotation = top_anno,left_annotation = hallmark_row_anno,row_split = assay$HALLMARK,row_title_gp = grid::gpar(fontsize=0),row_names_side = "left",row_names_gp = grid::gpar(fontsize=3),...)

  } else{
    assay_colors<-rep("red",length(stats::na.omit(unique(as.vector(heatmap_matrix)))))
    names(assay_colors)<-stats::na.omit(unique(as.vector(heatmap_matrix)))
    if("Silent" %in% names(assay_colors)) {assay_colors["Silent"]="white"}
    heatmap_plot<-ComplexHeatmap::Heatmap(heatmap_matrix,name="Variant_Impact",cluster_rows = F,cluster_columns = F,col=assay_colors,na_col = "grey",top_annotation = top_anno,left_annotation = hallmark_row_anno,row_split = assay$HALLMARK,row_title_gp = grid::gpar(fontsize=0),row_names_side = "left",row_names_gp = grid::gpar(fontsize=3),show_heatmap_legend = F,...)
  }
  return(list(heatmap_matrix=assay,heatmap_plot=heatmap_plot))
}
