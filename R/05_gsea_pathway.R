# =============================================================================
# GSEA AND PATHWAY ANALYSIS FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
# CONTENTS:
#   1. gsea                - fGSEA with enrichment maps and table output
#   2. pathway_gene_table  - Pathway-to-gene summary table
#   3. fgseatableplot      - GSEA result table visualization (ComplexHeatmap)
#   4. sankey_heatmap      - Sankey+heatmap pathway-gene visualization
#   5. read_pathways       - Read pathway gene sets from file
#   6. goEnrich            - Gene ontology enrichment with multi-plot output
#   7. geneset_activity    - ssGSEA/GSVA/zscore geneset scoring
# SOURCE: Pembro/funcsInPembro.R (2025-07-02) | myscripts.R (2023-06-29)
# =============================================================================


# SOURCE: Pembro/funcsInPembro.R (2025-07-02)
#' Gsea.
#'
#' Gsea.
#' @param statistics Function argument documented from the legacy interface.
#' @param output_dir Output plotting or file parameter used by the existing implementation.
#' @param logFC_col Column name used by the existing implementation.
#' @param pval_col Column name used by the existing implementation.
#' @param pval_cutoff Numeric tuning parameter used by the existing implementation.
#' @param pathways Function argument documented from the legacy interface.
#' @param enrichment_padj_cutoff Numeric tuning parameter used by the existing implementation.
#' @param pathway_gene_table_file Output plotting or file parameter used by the existing implementation.
#' @param fgseatableplot_file Output plotting or file parameter used by the existing implementation.
#' @param enrich_pathways Function argument documented from the legacy interface.
#' @param enrichment_map Function argument documented from the legacy interface.
#' @param FC_cutoff Numeric tuning parameter used by the existing implementation.
#' @param kappa_cutoff Numeric tuning parameter used by the existing implementation.
#' @param ... Additional arguments passed through to downstream functions.
#' @return The result object produced by the analysis.
#' @details Source provenance: Pembro/funcsInPembro.R (2025-07-02).
#'
#' @examples
#' \dontrun{
#' gsea(...)
#' }
#' @export
gsea<-function(statistics,output_dir=".",logFC_col="logFC",pval_col="P.Value",pval_cutoff=0.05,pathways,enrichment_padj_cutoff=0.05,pathway_gene_table_file,fgseatableplot_file="fgseatableplot.tiff",enrich_pathways="up2down2",enrichment_map=T,FC_cutoff=1.5,kappa_cutoff=0.3,...){
  if(!dir.exists(output_dir)){dir.create(output_dir)}
  if(class(statistics)=="data.frame"){
    stats<-statistics[[logFC_col]]
    names(stats)<-rownames(statistics)
  }
  fgseaRes <- fgsea(pathways = pathways, stats = stats,minSize=15,maxSize=500,nperm=10000)
  fgseaRes<-fgseaRes[order(fgseaRes$NES,decreasing=T)]
  sig_pathways<-pathways[fgseaRes$pathway[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj))]]
  if(length(sig_pathways)==0) {
    cat("no significant signaling pathway was found.")
    return(NULL)
  }
  if(sum(fgseaRes$padj<=enrichment_padj_cutoff  & (!is.na(fgseaRes$padj)))>30){
    sig_pathways<-pathways[fgseaRes$pathway[order(fgseaRes$padj)][1:30]]
  }
  pathway_gene_table=pathway_gene_table(sig_pathways,stats)
  if(missing(pathway_gene_table_file)){
    pathway_gene_table_file="pathway_gene_table.csv"
  }
  write.csv(pathway_gene_table,file.path(output_dir,pathway_gene_table_file))

  tiff(filename = file.path(output_dir,fgseatableplot_file),units="in", width=15, height=5, res=300, compression = 'lzw')
  fgseatableplot(pathways = sig_pathways,stats = stats,fgseaRes = fgseaRes)
  dev.off()

  if(!is.na(enrich_pathways)){
    up_enrich_pathways<-NULL
    down_enrich_pathways<-NULL
    enrich_pathways_<-NULL
    if(grepl("up",enrich_pathways)){
      n_up_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=up)[0-9]+",enrich_pathways,perl=T))))
      up_enrich_pathways<-head(fgseaRes$pathway,n_up_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(grepl("down",enrich_pathways)){
      n_down_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=down)[0-9]+",enrich_pathways,perl=T))))
      down_enrich_pathways<-tail(fgseaRes$pathway,n_down_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(is.null(enrich_pathways_)){enrich_pathways_<-enrich_pathways}
  }
  enrichpathwaysplot<-list()
  for(enrich_pathway in enrich_pathways_){
    enrich_filename= paste(enrich_pathway,".tiff",sep="")
    tiff(filename = file.path(output_dir,enrich_filename),units="in",width=6,height=4,res=300,compression='lzw')
    plot_enrichment<-plotEnrichment(pathways[[enrich_pathway]],stats)+labs(title=enrich_pathway)
    print(plot_enrichment)
    dev.off()
    enrichpathwaysplot[[enrich_pathway]]<-paste(enrich_pathway,".tiff",sep="")
  }
  enrichment_map_data<-NULL
  if(enrichment_map){
    sig_genes<-rownames(statistics)[statistics[[pval_col]]<=pval_cutoff]
    enrichment_map_data<-fgseaRes[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj)),]
    enrichment_map_data<-data.frame(ID=paste("pathway",1:nrow(enrichment_map_data),sep="_"),
                                    Term_Description=enrichment_map_data[["pathway"]],
                                    lowest_p=enrichment_map_data[["padj"]])
    changed_genes<-sapply(enrichment_map_data[["Term_Description"]],function(path){up_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]>=FC_cutoff]);
    up_regulated<-up_regulated[up_regulated %in% sig_genes];
    down_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]<=(-FC_cutoff)]);
    down_regulated<-down_regulated[down_regulated %in% sig_genes];
    return(c("Up_regulated"=paste(up_regulated,collapse = ", "),"Down_regulated"=paste(down_regulated,collapse = ", ")))
    })
    changed_genes<-as.data.frame(t(changed_genes))
    assertthat::are_equal(enrichment_map_data$Term_Description,rownames(changed_genes))
    enrichment_map_data<-cbind(enrichment_map_data,changed_genes)
    enrichment_map_data<-enrichment_map_data[(enrichment_map_data[["Up_regulated"]]!="") | (enrichment_map_data[["Down_regulated"]]!=""),]
    if(nrow(enrichment_map_data)>=2){
      tiff(filename = file.path(output_dir,"enrichment_map.tiff"),units="in", width=8, height=8, res=300, compression = 'lzw')
      enrichment_map_data<-cluster_enriched_terms(enrichment_map_data,use_description = T,plot_clusters_graph = T,plot_dend=T,kappa_threshold=kappa_cutoff)
      dev.off()
    } else{
      cat("rows of enrichment_map_data is less than 2. therefore, no enriched terms plot will be produced.")
    }
  }

  return(list(fgseatable=fgseaRes,sig_pathways=sig_pathways,pathway_gene_table=pathway_gene_table,enrichment_map_data=enrichment_map_data,pathway_gene_table_file=pathway_gene_table_file,fgseatableplot_file=fgseatableplot_file,enrichpathwaysplot=enrichpathwaysplot))
}


# SOURCE: Pembro/funcsInPembro.R (2025-07-02)
#' Pathway gene table.
#'
#' Pathway gene table.
#' @param pathways Function argument documented from the legacy interface.
#' @param stats Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: Pembro/funcsInPembro.R (2025-07-02).
#'
#' @examples
#' \dontrun{
#' pathway_gene_table(pathways = ..., stats = ...)
#' }
#' @export
pathway_gene_table<-function(pathways,stats){
  stats<-data.frame(Gene_ID=names(stats),logFC=stats)
  genes_ordered_by_freq<-names(genes_freq<-sort(table(Reduce(c,pathways))[intersect(stats[["Gene_ID"]],unique(Reduce(c,pathways)))],decreasing = T))
  results<-data.frame(Gene=genes_ordered_by_freq)
  for(pathway in names(pathways)){
    genes_of_pathway<-intersect(genes_ordered_by_freq,pathways[[pathway]])
    results[[pathway]]=NA
    results[[pathway]][match(genes_of_pathway,genes_ordered_by_freq)]<-stats[["logFC"]][match(genes_of_pathway,stats[["Gene_ID"]])]
  }
  return(results)
}


# SOURCE: Pembro/funcsInPembro.R (2025-07-02)
fgseatableplot<-function (pathways, stats, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 3, 1.2, 1.2), render = TRUE)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  pathways <- pathways[sapply(pathways, length) > 0]
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(pn, just = "right", x = unit(0.95, "npc")),
         ggplot() +
           geom_segment(aes(x = p, xend = p, y = 0, yend = statsAdj[p]), size = 0.2) +
           scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         ggplot()+geom_rect(aes(xmin=0,xmax=fgseaRes$NES[fgseaRes$pathway==pn],ymin=-0.6,ymax=0.6),fill=ifelse(sign(fgseaRes$NES[fgseaRes$pathway==pn])==1,"green","red")) +
           scale_x_continuous(limits = c(min(fgseaRes$NES,na.rm=T)*1.1, max(max(fgseaRes$NES,na.rm=T)*1.1,0)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         #         textGrob(sprintf("%.2f", annotation$NES)),
         textGrob(sprintf("%.1e", annotation$pval)),
         textGrob(sprintf("%.1e", annotation$padj)))
  })
  rankPlot_segment <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  rankPlot_rect <- ggplot() + geom_blank() + scale_x_continuous(limits = c(min(fgseaRes$NES)*1.1, max(fgseaRes$NES)*1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(), rankPlot_segment, rankPlot_rect, nullGrob(), nullGrob()))

  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))

  p <- arrangeGrob(grobs = grobs[grobsToDraw], ncol = sum(colwidths != 0), widths = colwidths[colwidths != 0])
  if (render) {
    grid.draw(p)
  }
  else {
    p
  }
}


# SOURCE: Pembro/funcsInPembro.R (2025-07-02)
#' Sankey heatmap.
#'
#' Sankey heatmap.
#' @param sig_expressions Function argument documented from the legacy interface.
#' @param scale Option controlling how the function runs.
#' @param pathways Function argument documented from the legacy interface.
#' @param keep_other Function argument documented from the legacy interface.
#' @param other_color Color specification used when plotting.
#' @param pathway_colors Color specification used when plotting.
#' @param sankey_width Function argument documented from the legacy interface.
#' @param line_size Function argument documented from the legacy interface.
#' @param text_size Function argument documented from the legacy interface.
#' @param ... Additional arguments passed through to downstream functions.
#' @return A heatmap-like plotting object or rendered plot, depending on the code path.
#' @details Source provenance: Pembro/funcsInPembro.R (2025-07-02).
#'
#' @examples
#' \dontrun{
#' sankey_heatmap(...)
#' }
#' @export
sankey_heatmap<-function(sig_expressions,scale=TRUE,pathways,keep_other=TRUE,other_color="gray",pathway_colors,sankey_width=unit(10, "cm"),line_size=1,text_size=6,...){
  if(scale){
    sig_expressions<-t(scale(t(sig_expressions)))
  }
  if(missing(pathway_colors)){
    pathway_colors<-sample(c("#006400","#00008b","#b03060","#ff4500","#ffd700","#7fff00","#00ffff","#ff00ff","#6495ed","#ffdab9"),length(pathways))
    names(pathway_colors)<-names(pathways)
  }
  pathways<-data.frame(gene=unname(unlist(pathways)),pathway=rep(names(pathways),times=sapply(pathways,length)))
  sankey_df<-pathways[pathways[["gene"]] %in% rownames(sig_expressions),]
  sankey_df<-sankey_df[!duplicated(sankey_df),]
  if(keep_other){
    sankey_df_other<-data.frame(gene=setdiff(rownames(sig_expressions),sankey_df[["gene"]]),pathway="Other")
    sankey_df<-rbind(sankey_df,sankey_df_other)
    pathway_colors<-c(pathway_colors,"Other"=other_color)
  }
  sankey_anno = HeatmapAnnotation(sankey= anno_empty(border = F, width = sankey_width),which = "row")
  ht = Heatmap(sig_expressions, right_annotation = sankey_anno,...)
  ht = draw(ht)
  h_genes = data.frame(gene=rownames(sig_expressions)[row_order(ht)])
  h_genes[["h_order"]]=nrow(h_genes):1
  sankey_df[["h_order"]]<-h_genes[["h_order"]][match(sankey_df[["gene"]],h_genes[["gene"]])]
  sankey_df<-sankey_df[order(sankey_df[["pathway"]],sankey_df[["h_order"]],decreasing=T),]
  sankey_df[["point_x"]]=0.02
  sankey_df[["point_y"]]=sankey_df[["h_order"]]
  sankey_df[["point_color"]]<-pathway_colors[sankey_df[["pathway"]]]
  sankey_df[['point_color']][which(sankey_df[["gene"]] %in% sankey_df[["gene"]][duplicated(sankey_df[["gene"]])])]<-"black"
    sankey_df[["rect_x"]]=0.95
    sankey_df[["rect_y"]]=nrow(sankey_df):1/nrow(sankey_df)*nrow(sig_expressions)
    sankey_df[["rect_width"]]=0.05
    sankey_df[["rect_height"]]=1*nrow(sig_expressions)/nrow(sankey_df)
    sankey_df[["curve_x_start"]]<-0.05
    sankey_df[["curve_x_end"]]<-0.95
    sankey_df[["curve_y_start"]]<-sankey_df[["h_order"]]
    sankey_df[["curve_y_end"]]<-sankey_df[["rect_y"]]
    sankey_df[["text_x"]]=0.94
    sankey_df<-as.data.frame(sankey_df %>% group_by(pathway) %>% mutate(text_y=mean(rect_y)))
    sankey_df[["color"]]<-pathway_colors[sankey_df[["pathway"]]]
    decorate_annotation("sankey", {
      pushViewport(viewport(xscale = c(0, 1),yscale = c(0.5, nrow(sig_expressions)+0.5)))
      sigmoid=1 / (1 + exp(-seq(-5,5,length=100)))
      sigmoid<-(sigmoid-sigmoid[1])/(sigmoid[100]-sigmoid[1])
      for(i in 1:nrow(sankey_df)){
        grid.circle(x=sankey_df[i,"point_x"], y=sankey_df[i,"point_y"], r=unit(1,"mm"),gp = gpar(col="black",fill = sankey_df[i,"point_color"],alpha=0.5), default.units = "native")
        grid.lines(x=seq(sankey_df[i,"curve_x_start"],sankey_df[i,"curve_x_end"],length=100),y=(sankey_df[i,"curve_y_start"]+sigmoid)+sigmoid*(sankey_df[i,"curve_y_end"]-sankey_df[i,"curve_y_start"]-1),gp=gpar(col= sankey_df[i,"color"],alpha=0.2,lwd=line_size),default.units = "native")
        grid.rect(x=sankey_df[i,"rect_x"],y=sankey_df[i,"rect_y"],width=sankey_df[i,"rect_width"],height=sankey_df[i,"rect_height"],gp=gpar(fill=sankey_df[i,"color"],col=NA,alpha=1),just="left",default.units = "native")
        grid.text(label=sankey_df[i,"pathway"],x=sankey_df[i,"text_x"],y=sankey_df[i,"text_y"],just = "right",gp=gpar(fontsize=text_size),default.units = "native")
      }
      popViewport()
    })
}


# SOURCE: myscripts.R (2023-06-29)
#' Read pathways.
#'
#' Read pathways.
#' @param pathway_file Output plotting or file parameter used by the existing implementation.
#' @param splitter Function argument documented from the legacy interface.
#' @param pathwayname_ind Function argument documented from the legacy interface.
#' @param gene_start_ind Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: myscripts.R (2023-06-29).
#'
#' @examples
#' \dontrun{
#' read_pathways(pathway_file = ..., splitter = ...)
#' }
#' @export
read_pathways<-function(pathway_file,splitter="\t",pathwayname_ind=1,gene_start_ind=3){
  pathways=list()
  con=file(pathway_file,"r")
  while(TRUE){
    line=readLines(con,n=1)
    if(length(line)==0){break}
    words<-unlist(strsplit(line,splitter))
    pathways[[words[pathwayname_ind]]]<-words[gene_start_ind:length(words)]
  }
  close(con)
  return(pathways)
}


# SOURCE: myscripts.R (2023-06-29)
#' Go Enrich.
#'
#' Go Enrich.
#' @param genelist Function argument documented from the legacy interface.
#' @param ont Option controlling how the function runs.
#' @param enrichGO_params Function argument documented from the legacy interface.
#' @param rendergodagplot Function argument documented from the legacy interface.
#' @param dagplot_params Function argument documented from the legacy interface.
#' @param renderdotplot Function argument documented from the legacy interface.
#' @param dotplot_params Function argument documented from the legacy interface.
#' @param rendersemanticemapplot Function argument documented from the legacy interface.
#' @param emapplot_params Function argument documented from the legacy interface.
#' @param rendercnetplot Function argument documented from the legacy interface.
#' @param cnetplot_params Function argument documented from the legacy interface.
#' @param qvalue_cutoff Numeric tuning parameter used by the existing implementation.
#' @param similarity_obj Function argument documented from the legacy interface.
#' @param term_similarity_params Function argument documented from the legacy interface.
#' @param simplifygo Function argument documented from the legacy interface.
#' @param simplifygo_params Function argument documented from the legacy interface.
#' @param enrichmap Function argument documented from the legacy interface.
#' @param termsimilarity_cutoff Numeric tuning parameter used by the existing implementation.
#' @param enrichmap_params Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: myscripts.R (2023-06-29).
#'
#' @examples
#' \dontrun{
#' goEnrich(genelist = ..., ont = ...)
#' }
#' @export
goEnrich<-function(genelist,ont="BP",enrichGO_params=NULL,rendergodagplot=T,dagplot_params=NULL,renderdotplot=T,dotplot_params=NULL,rendersemanticemapplot=T,emapplot_params=NULL,rendercnetplot=T,cnetplot_params=NULL,qvalue_cutoff=0.05,similarity_obj=c("gene","term"),term_similarity_params=NULL,simplifygo=T,simplifygo_params=NULL,enrichmap=T,termsimilarity_cutoff=0.2,enrichmap_params=NULL){
  similarity_obj<-match.arg(similarity_obj)
  go<-do.call(clusterProfiler::enrichGO, c(list(gene=genelist,OrgDb = org.Hs.eg.db, keyType = "SYMBOL",ont= ont),enrichGO_params))
  write.csv(go@result,file = "goEnrich_statistics.csv",row.names = F)
  if(rendergodagplot){
    tiff(filename = "godagplot.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    goPlot<-do.call(enrichplot::goplot,c(list(x=go),dagplot_params))
    print(goPlot)
    dev.off()
  }
  if(renderdotplot){
    tiff(filename = "godotplot.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    dotPlot<-do.call(enrichplot::dotplot,c(list(object=go),dotplot_params))
    print(dotPlot)
    dev.off()
  }
  if(rendersemanticemapplot){
    semanticsimilarities=pairwise_termsim(go)
    tiff(filename = "gosemanticenrichmap.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    emapPlot<-do.call(enrichplot::emapplot,c(list(x=semanticsimilarities),emapplot_params))
    print(emapPlot)
    dev.off()
  }
  if(rendercnetplot){
    tiff(filename = "gocnetplot.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    cnetPlot<-do.call(enrichplot::cnetplot,c(list(x=go),cnetplot_params))
    print(cnetPlot)
    dev.off()
  }
  sig_goid<-go@result$ID[go@result$qvalue<=qvalue_cutoff]
  if(similarity_obj=="gene"){
    gotable<-GOSemSim::godata(OrgDb=org.Hs.eg.db,keytype='SYMBOL',ont="BP")@geneAnno
    gotable<-gotable[(gotable$SYMBOL %in% genelist) & (gotable$GO %in% sig_goid),]
    gol<-lapply(unique(gotable$GO),function(id){return(gotable$SYMBOL[gotable$GO==id])})
    names(gol)<-unique(gotable$GO)
    termsimilarities<-do.call(term_similarity,c(list(gl=gol),term_similarity_params))
  } else {
    termsimilarities<-do.call(GO_similarity,c(list(go_id=sig_go,ont=ont),term_similarity_params))
  }
  if(simplifygo){
    tiff(filename = "goclusterplot.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    df = do.call(simplifyGO,c(list(mat=termsimilarities),simplifygo_params))
    dev.off()
  }
  if(enrichmap){
    sig_go_result<-go@result[go@result$qvalue<=qvalue_cutoff,]
    termsimilarities_<-termsimilarities
    diag(termsimilarities_)<-0
    graph<-igraph::graph.adjacency(termsimilarities_,weighted=T,mode="lower")
    graph <- igraph::delete.edges(graph, E(graph)[ weight < termsimilarity_cutoff ])
    tiff(filename = "enrichmapplot.tiff",width = 12,height=8,units = "in",res=300,compression = "lzw")
    enrichmapPlot=do.call(plot,c(list(x=graph,vertex.color=rainbow(n=max(df$cluster))[df$cluster],
                                      vertex.size=-log10(sig_go_result$qvalue),
                                      vertex.label=sig_go_result$Description),enrichmap_params))
    print(enrichmapPlot)
    dev.off()
  }
  return(list(go=go,termsimilarities=termsimilarities,gocluster=df,enrichmap_graph=graph))
}

# Shared helper `geneset_activity` is defined canonically in 03_expression_analysis.R.
