# =============================================================================
# SPATIAL IMMUNE PROFILING FUNCTIONS (mIF/Vectra)
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
# CONTENTS (Stats):
#   defineCell, reassignNA, getCtsPcDens, compute_all_nearest_distance,
#   find_nearest_distance_dist, unique_phenotypes, distance_matrix,
#   select_rows, getMeanNearestNeighbourDistance, count_within_many,
#   count_within, getCountWithin, performTtestsAllRows,
#   performTtestsAllClassesOneVsRest, robustscale, uniCoxPh, multiCoxPh,
#   getSigPair, bi_ripleys_k_mod, getMeanCountWithin
# CONTENTS (Plots):
#   PieDonut, plotPieDonutBySlides, plot_immunoflo, trimMif4plots,
#   plotImmuneHighlights
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#          woodman_lab.XLi23/mIFvectra/R/vectraPlotsFun.R (2024-04-09)
# NOTE: Use mIFvectra version (3 panels) for general work;
#       Updated_Panel_10/vectraStatsFun.R (2024-08-05) has expanded Panel 10 cell defs
# =============================================================================


# =============================================================================
# STATS FUNCTIONS
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
# =============================================================================

#' define Cell
#'
#' @param x Function argument documented from the legacy interface.
#' @param Cpheno Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' defineCell(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
defineCell=function(x,Cpheno){
  apply(x[,Cpheno],1,function(x) paste(na.omit(Cpheno[x]),collapse = ""))
}


Panels=list(
  `panel 10`=c("CD8+","CK+CD3+CD68+","Foxp3+","Ki67+","PD1+","PDL1+"),
  `panel 1`=c("CD8+","CK+CD3+CD68+","PD1+","PDL1+"),
  `panel 2`=c("CD8+","CD45RO+","CK+CD3+","GB+")
)

Panels.EssentialCellDefs=list(
  `panel 10`=list(
    c("CK"),
    c("CK","PDL1"),
    c("CK","Ki67"),
    c("CK","PDL1","Ki67"),
    c("CD3"),
    c("CD3","PDL1"),
    c("CD3","PD1"),
    c("CD3","PD1","PDL1"),
    c("CD3","KI67"),
    c("CD3","PDL1","Ki67"),
    c("CD3","PD1","Ki67"),
    c("CD3","PD1","PDL1","Ki67"),
    c("CD3","Foxp3"),
    c("CD3","CD8"),
    c("CD3","CD8","PDL1"),
    c("CD3","CD8","PD1"),
    c("CD3","CD8","PD1","PDL1"),
    c("CD3","CD8","Ki67"),
    c("CD3","CD8","PDL1","Ki67"),
    c("CD3","CD8","PD1","Ki67"),
    c("CD3","CD8","PD1","PDL1","Ki67"),
    c("CD68"),
    c("CD68","PDL1")#
  )
)



#' reassign NA in a vector
#'
#' @param x Function argument documented from the legacy interface.
#' @param y Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' reassignNA(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
reassignNA=function(x,y){
  x[is.na(x)]=y
  return(x)}




#### counts, percentage and density ####

#' get counts, percentage and density
#'
#' @param dfs a data frame object that contains slide ID, ROI ID and cell ID column
#' @param ROIIDs a vector of ROI IDs
#' @param cellDefs a vector of cell definition -- positive markers sorted and combined.
#' @param TissueCategoryAreas Function argument documented from the legacy interface.
#' @param mutually_exclusive Function argument documented from the legacy interface.
#' @param byROI Function argument documented from the legacy interface.
#' @param bySlide Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' getCtsPcDens(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
getCtsPcDens=function(dfs,ROIIDs,cellDefs,TissueCategoryAreas,mutually_exclusive=T,byROI=T,bySlide=T,selectedROIs){
  tissCat=rownames(TissueCategoryAreas)
  results=list()
  if(!mutually_exclusive) consideredCombns=lapply(strsplit(cellDefs,"\\+"), function(x) paste(x,"+",sep = ""))
  FeatureDf=data.frame(
    `Marker Combn`=rep(cellDefs,3),
    `Tissue Category`=rep(tissCat,each=length(cellDefs)),
    check.names = F
  )
  if(byROI){results[["ctsByROI"]]=results[["pcByROI"]]=results[["densByROI"]]=FeatureDf}
  if(bySlide){results[["ctsBySlide"]]=results[["pcBySlide"]]=results[["densBySlide"]]=FeatureDf}
  for(ROI in ROIIDs){
    TissueCategoryArea=t(TissueCategoryAreas[,ROI])

# counts_percentage_density -----------------------------------------------


    df=dfs[[ROI]]
    if(mutually_exclusive){
      tmp=df[,c("Cell ID","Tissue Category","cellDef")] %>% mutate(bin=T) %>%
        tidyr::spread(key = cellDef,value = bin) %>%
        tibble::add_column(!!!setNames(rep(FALSE,length(setdiff(cellDefs,colnames(.)))),setdiff(cellDefs,colnames(.)))) %>%
        select(match(c("Tissue Category",cellDefs),colnames(.))) %>% reassignNA(F)

    }else{
      tmp=data.frame(df[,c("Tissue Category"),drop=F],check.names = F)

      for(i in 1:length(cellDefs)){
        consideredCombn = consideredCombns[[i]]
        tmp[[cellDefs[i]]]=apply(df[,consideredCombn,drop=F],1,all)
      }
    }

    df2bind=data.frame(`Mutually Exclusive`=mutually_exclusive,check.names = F)
    # binary combination variables
    results[["binByROI"]][[ROI]]=cbind(df2bind,tmp)
    results[["TotalCtsByROI"]][[ROI]]=totalCts=nrow(df)
    tmp=sapply(tissCat[tissCat!="all"],function(x) tmp%>%filter(`Tissue Category`==x)%>%select(grep("\\+",colnames(tmp)))%>%colSums)
    tmp=cbind(tmp,all=rowSums(tmp))
    if(byROI){
      # counts
      results[["ctsByROI"]][[ROI]]=reshape2::melt(tmp)[,"value"]
      # percentage
      temp=tmp/totalCts
      results[["pcByROI"]][[ROI]]=reshape2::melt(temp)[,"value"]
      # density
      temp=tmp/TissueCategoryArea[rep(1,nrow(tmp))]
      results[["densByROI"]][[ROI]]=reshape2::melt(temp)[,"value"]
    }
  }

  if(byROI){
    # Add df2bind and feature type
    by="ROI"
    tmp=setNames(paste(c("cts","pc","dens"),by,sep = "By"),c("counts","percentage","density"))
    for(i in 1:length(tmp)){
      results[[tmp[i]]]=cbind(df2bind,`feature type`=names(tmp)[i],results[[tmp[i]]])
    }
  }

  if(bySlide){
  selectedROIs=lapply(selectedROIs,function(x) x[x%in%ROIIDs])
  selectedROIs=selectedROIs[sapply(selectedROIs, function(x) length(x)>0)]
    for(slide in names(selectedROIs)){
      df2bind=data.frame(`Mutually Exclusive`=mutually_exclusive,check.names = F)
      ROIs=selectedROIs[[slide]]
      TissueCategoryArea=t(rowSums(TissueCategoryAreas[,ROIs,drop=F]))
      totalCts=sum(unlist(results[["TotalCtsByROI"]][ROIs]))


      results[["binBySlide"]][[slide]] = tmp= Reduce(rbind,results[["binByROI"]][ROIs])
      tmp=tmp[,(ncol(df2bind)+1):ncol(tmp)]
      tmp=sapply(tissCat[tissCat!="all"],function(x) tmp%>%filter(`Tissue Category`==x)%>%select(grep("\\+",colnames(tmp)))%>%colSums)
      tmp=cbind(tmp,all=rowSums(tmp))

      # counts
      results[["ctsBySlide"]][[slide]]=reshape2::melt(tmp)[,"value"]
      # percentage
      temp=tmp/totalCts
      results[["pcBySlide"]][[slide]]=reshape2::melt(temp)[,"value"]
      # density
      temp=tmp/TissueCategoryArea[rep(1,nrow(tmp))]
      results[["densBySlide"]][[slide]]=reshape2::melt(temp)[,"value"]
    }

    # Add df2bind and feature type
    by="Slide"
    tmp=setNames(paste(c("cts","pc","dens"),by,sep = "By"),c("counts","percentage","density"))
    for(i in 1:length(tmp)){
      results[[tmp[i]]]=cbind(df2bind,`feature type`=names(tmp)[i],results[[tmp[i]]])
    }
  }
  results[["TotalCtsByROI"]]=NULL
  # if(!byROI){results[["binByROI"]]=NULL}
  return(results)
}


#### Nearest neighbors distance ####
#' Nearest neighbors from a file.
#'
#' Compute nearest distance to each phenotype for each cell in a
#' (possibly merged) inForm cell seg table. Add `Distance to <phenotype>`
#' columns.
#' Write the result to a new file.
#'
#' NOTE: The input file is read using [read_cell_seg_data] so the conversions
#' and cleanup it does will be applied to the output data.
#'
#' @param cell_table_path Path to an inForm cell seg data file, or NULL
#' to prompt for the path.
#' @param out_path Path to the output file, or NULL to create a path from the
#' input file path.
#' @importFrom magrittr "%>%"
#' @export
#' @seealso [find_nearest_distance] which performs the distance calculation.
#' @family distance functions
#' @md
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
compute_all_nearest_distance <- function(cell_table_path=NULL, out_path=NULL) {
  # Get the path to the cell seg table and check it
  if (is.null(cell_table_path))
    cell_table_path = file.choose()

  # Read the table
  cat('Reading', cell_table_path, '\n')
  csd = read_cell_seg_data(cell_table_path)

  # Compute the distances
  cat('Computing distances\n')
  result = NULL
  phenos = unique_phenotypes(csd)
  result = csd %>%
    dplyr::group_by(!!rlang::sym(field_column(csd))) %>%
    dplyr::do(dplyr::bind_cols(., find_nearest_distance(., phenos)))

  if (is.null(out_path))
    out_path = sub('\\.txt$', '_dist.txt', cell_table_path)
  cat('Writing', out_path, '\n')
  readr::write_tsv(result, out_path, na='#N/A')
}

#' Distance-matrix implementation of `find_nearest_distance`.
#' @param csd A data frame with `Cell X Position`,
#'        `Cell Y Position` and `Phenotype` columns,
#'        such as the result of calling
#'        [read_cell_seg_data].
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' `unique_phenotypes(csd)` will be used.
#' @param dst Optional distance matrix. If provided, this should be
#' `distance_matrix(csd)`.
#' @return A `tibble` containing a `Distance to <phenotype>` column
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' and `Cell ID <phenotype>` column for each phenotype.
#' Columns will contain `NA` values where there is no other cell
#' of the phenotype.
#' @seealso find_nearest_distance
#' @md
#' @keywords internal
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
find_nearest_distance_dist = function(csd, phenotypes=NULL, dst=NULL) {

  phenotypes = validate_phenotypes(phenotypes, csd)

  if (is.null(dst))
    dst = distance_matrix(csd)

  # Removing dimnames from dst gives unnamed columns in results
  # The names are clutter, take up space and confuse the tests
  dimnames(dst) = NULL

  result = purrr::map2_dfc(
    names(phenotypes), phenotypes,
   function(name, phenotype) {
     # Which cells are in the target phenotype?
     phenotype_cells = select_rows(csd, phenotype)
     if (sum(phenotype_cells)>0) {
       # Subset columns of the distance matrix to just phenotype cells
       phenotype_dist = dst[, phenotype_cells, drop=FALSE]

       # Find the minimum distance > 0; i.e. from cells to not-self cells
       dist_col = apply(phenotype_dist, 1, row_min)

       # Find the index of the minimum distance > 0
       # and use this to index the Cell IDs of the target phenotypes
       which_dist_col = apply(phenotype_dist, 1, which_row_min)
       cell_id_col = csd$`Cell ID`[phenotype_cells][which_dist_col]
       pheno_cols = tibble::tibble(dist_col, cell_id_col)
     }
     else {
       # No cells of the selected phenotype
       na_col = rep(NA_integer_, nrow(csd))
       pheno_cols = tibble::tibble(dist_col=na_col, cell_id_col=na_col)
     }
     pheno_cols  %>%
       rlang::set_names(paste(c('Distance to', 'Cell ID'), name))
   })
}

unique_phenotypes = function (csd) {
  if ("Phenotype" %in% names(csd))
    return(purrr::discard(sort(unique(csd$Phenotype)), ~.x ==
                            ""))
  phenos = names(csd) %>% stringr::str_subset("Phenotype ") %>%
    stringr::str_remove("Phenotype ") %>% stringr::str_c("+")
  if (length(phenos) == 0)
    stop("Cell seg table does not have a phenotype column.")
  rlang::set_names(phenos)
}

distance_matrix = function (csd) {
  stopifnot("Cell X Position" %in% names(csd), "Cell Y Position" %in%
              names(csd))
  as.matrix(stats::dist(csd[, c("Cell X Position", "Cell Y Position")]))
}

#' Flexibly select rows of a data frame.
#'
#' Select rows of a data frame based on phenotypes or other
#' expressions.
#'
#' `select_rows` implements a flexible mechanism for selecting cells (rows)
#' from a cell segmentation table. Cells may be selected by single or
#' multiple phenotype, by expression level, or combinations of both.
#'
#' See the tutorial
#' [Selecting cells within a cell segmentation table](https://akoyabio.github.io/phenoptr/articles/selecting_cells.html)
#'for extensive documentation and examples.
#'
#' @param csd A data frame
#' @param sel May be a character vector, a one-sided formula, a list
#'   containing such or `NA`. A character vector is interpreted as
#'   the name(s) of one or
#'   more phenotypes and selects any matching phenotype. A formula is
#'   interpreted as an expression on the columns of `csd`.
#'   Multiple list items are joined with AND. `NA` is interpreted
#'   as "select all". It is convenient for lists of selection criteria.
#' @return A logical vector of length `nrow(csd)` which selects rows
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'   according to `sel`.
#' @export
#' @examples
#' csd <- sample_cell_seg_data
#'
#' # Select tumor cells with PDL1 expression > 3
#' selector <- list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3)
#' pdl1_pos_tumor <- csd[select_rows(csd, selector),]
#' range(pdl1_pos_tumor$`Entire Cell PDL1 (Opal 520) Mean`)
#'
#' # Select all T-cells. Note: Use c() to combine phenotypes, not list()
#' selector <- c('CD8+', 'FoxP3+')
#' tcells <- csd[select_rows(csd, selector),]
#' table(tcells$Phenotype)
#' @md
#' @seealso [parse_phenotypes] for a convenient way to create selectors
#' for most common phenotypes.
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
select_rows <- function(csd, sel) {
  stopifnot(is.data.frame(csd))

  # Evaluate a single phenotype in a per-marker file
  evaluate_per_marker = function(s) {
    if (!stringr::str_detect(s, '[+-]$'))
      stop(paste0(s, ' is not a valid per-marker phenotype name.'))
    column_name = paste('Phenotype', stringr::str_remove(s, '[+-]$'))
    if (!column_name %in% names(csd))
      stop(paste0("No '", column_name, "' column in data."))
    csd[[column_name]] == s
  }

  # Evaluate a single selector
  select_one = function(s) {
    if (length(s)==1 && is.na(s)) {
      # NA means select all
      rep(TRUE, nrow(csd))
    } else if (is.character(s)) {
      # Selector is one or more phenotype names,
      # look for match with phenotype column
      # Any match qualifies
      if ('Phenotype' %in% names(csd)) {
        csd[['Phenotype']] %in% s
      }
      else {
        # Phenotype per-marker has multiple columns
        col_selections = purrr::map(s, evaluate_per_marker)
        purrr::reduce(col_selections, `|`)
      }
    } else {
      # Selector is a function, evaluate it on csd
      col_selections = lazyeval::f_eval(s, csd)

      # Check for valid result
      if (class(col_selections) != 'logical' ||
          length(col_selections) != nrow(csd))
        stop('Invalid expression in select_rows: ~', lazyeval::f_text(s))
      col_selections
    }
  }

  # Everything is selected by default
  result = rep(TRUE, nrow(csd))
  if (!is.list(sel)) sel = list(sel)
  for (s in sel)
    result = result & select_one(s)

  # Don't return NA values, treat them as false
  result %>% tidyr::replace_na(FALSE)
}



#' get Mean Nearest Neighbour Distance
#'
#' @param dfs Function argument documented from the legacy interface.
#' @param ROIIDs Function argument documented from the legacy interface.
#' @param cellDefs Function argument documented from the legacy interface.
#' @param Cpheno Function argument documented from the legacy interface.
#' @param byROI Function argument documented from the legacy interface.
#' @param bySlide Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param mutually_exclusive Function argument documented from the legacy interface.
#' @param cellDefColName Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' getMeanNearestNeighbourDistance(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
getMeanNearestNeighbourDistance=function(dfs,ROIIDs,cellDefs,Cpheno,byROI=F,bySlide=T,selectedROIs,mutually_exclusive,cellDefColName="cellDef"){
  results=list()
  results[["NNmxByROI"]]=results[["NNDPhenoByROI"]]=results[["NNDPhenoBySlide"]]=list()

  CTs= sort(cellDefs)
  if (mutually_exclusive){
    phenoColName=cellDefColName
  }else{
    phenoColName=c(cellDefColName,Cpheno)
    phenoList=lapply(strsplit(cellDefs,"\\+"), function(x) as.list(paste(x,"+",sep="")))
    names(phenoList)=cellDefs
  }
  df2bind=data.frame(Phenotype=rep(CTs,each=length(CTs)),`Distance to`=rep(CTs,length(CTs)),check.names = F)
  df2bind=data.frame(`Tissue Category`=rep(c("stroma","tumor","all"),each=nrow(df2bind)),
                     df2bind[rep(1:nrow(df2bind),3),],check.names = F)
  results[["meanNNDByROI"]]=results[["meanNNDBySlide"]]=cbind(df2bind,`feature type`="nearest neighbour distance")
  headers=c("Cell ID","Tissue Category",phenoColName,"Cell X Position","Cell Y Position")

  for (ROI in ROIIDs){
    csd=dfs[[ROI]][,headers]
    #
    if(!mutually_exclusive){
      csd=csd[csd[,cellDefColName]!="",]
      csd[,Cpheno]=sapply(Cpheno, function(x) ifelse(csd[,x],x,sub("\\+","\\-",x)))
      csd=expr::changeColNames(csd,ind =match(Cpheno,colnames(csd)),newNames = paste("Phenotype",gsub("\\+","",Cpheno)))
      distances <- find_nearest_distance(csd,phenotypes =phenoList)
      colnames(csd)[colnames(csd)==cellDefColName]="Phenotype"
    }else{
      colnames(csd)[colnames(csd)==cellDefColName]="Phenotype"
      csd$Phenotype[!csd$Phenotype%in%CTs]="other"
      csd=csd %>% filter(Phenotype!="other")
      distances <- find_nearest_distance(csd)
    }
    results[["NNmxByROI"]][[ROI]]=tmp= csd[,c("Cell ID","Tissue Category","Phenotype")] %>%
      bind_cols(data.frame(ID=ROI),.,distances)
    tmp= tmp[,-grep("Cell ID .*\\+",colnames(tmp))]
    tmp=tmp %>% reshape2::melt(
      .,value.name = "Nearest Neighbor Distance",
      variable.name = "Distance to",id.vars=grep("Distance to",colnames(tmp),invert = T,value = T))
    tmp$`Distance to`=gsub("Distance to ","",tmp$`Distance to`)
    results[["NNDPhenoByROI"]][[ROI]]=tmp
    if(byROI){
      tmp1=tmp%>%
        group_by(`Tissue Category`, Phenotype,`Distance to`)%>%
        summarise(meanNNdist=mean(`Nearest Neighbor Distance`,na.rm=T))
      tmp2=tmp%>%
        group_by(Phenotype,`Distance to`)%>%
        summarise(meanNNdist=mean(`Nearest Neighbor Distance`))%>%
        bind_cols(`Tissue Category`="all",.)
      tmp=left_join(df2bind,bind_rows(tmp1,tmp2),by=c("Tissue Category","Phenotype","Distance to"),check.names=F)
      results[["meanNNDByROI"]][[ROI]]=tmp[["meanNNdist"]]
    }
  }
  if(bySlide){
    selectedROIs=lapply(selectedROIs,function(x) x[x%in%ROIIDs])
    selectedROIs=selectedROIs[sapply(selectedROIs, function(x) length(x)>0)]
    for(slide in names(selectedROIs)){
      ROIs=selectedROIs[[slide]]
      results[["NNDPhenoBySlide"]][[slide]]=tmp=Reduce(rbind,results[["NNDPhenoByROI"]][ROIs])

      tmp1=tmp%>%
        group_by(`Tissue Category`, Phenotype,`Distance to`)%>%
        summarise(meanNNdist=mean(`Nearest Neighbor Distance`,na.rm=T))
      tmp2=tmp%>%
        group_by(Phenotype,`Distance to`)%>%
        summarise(meanNNdist=mean(`Nearest Neighbor Distance`))%>%
        bind_cols(`Tissue Category`="all",.)
      tmp=left_join(df2bind,bind_rows(tmp1,tmp2),by=c("Tissue Category","Phenotype","Distance to"),check.names=F)
      results[["meanNNDBySlide"]][[slide]]=tmp[["meanNNdist"]]
    }
  }
  results[["NNDPhenoBySlide"]]=results[["NNDPhenoByROI"]]=NULL
  return(results)
}

#### count within ####

#' Count cells within a radius for multiple tissue categories and phenotypes
#' in a single field.
#'
#' This is a wrapper around [count_within()] which supports counting
#' multiple phenotype pairs and tissue categories within a single field.
#' For each given tissue category, pair of
#' 'from' phenotype and 'to' phenotype, and radius, it counts the number of
#'  'from' cells
#' having a 'to' cell within `radius` microns.
#'
#' The `category` parameter may be a single category or a list of categories.
#'
#' See the tutorial
#' [Selecting cells within a cell segmentation table](https://akoyabio.github.io/phenoptr/articles/selecting_cells.html)
#'  for more on
#' the use of `pairs` and `phenotype_rules`.
#'
#' @param csd A cell seg data table.
#' @param pairs A list of pairs of phenotypes. Each entry is a two-element
#'   vector. The result will contain values for each pair.
#' @param radius The radius or radii to search within.
#' @param phenotype_rules (Optional) A named list.
#'   Item names are phenotype names and must match entries in `pairs`.
#'   Item values are selectors for [select_rows()].
#' @param category Optional tissue categories to restrict both `from` and
#' `to` phenotypes.
#' @param verbose If TRUE, display progress.
#' @return A `tibble` containing these columns:
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'   \describe{
#'    \item{\code{slide_id}}{Slide ID from the data, if available.}
#'    \item{\code{source}}{Source field name.}
#'    \item{\code{field}}{Name of the individual field, if available.}
#'    \item{\code{category}}{Tissue category, if provided as a parameter,
#'    or "all".}
#'    \item{\code{from}}{From phenotype.}
#'    \item{\code{to}}{To phenotype.}
#'    \item{\code{radius}, \code{from_count}, \code{to_count},
#'    \code{from_with}, \code{within_mean}}{Results from [count_within]
#'    for this data file and tissue category.}
#'  }
#' @examples
#' csd <- sample_cell_seg_data
#'
#' # Count tumor cells near macrophages, and tumor cells near CD8 separately,
#' # in tumor and stroma tissue categories separately.
#' pairs <- list(c('CK+', 'CD68+'),
#'              c('CK+', 'CD8+'))
#' radius <- c(10, 25)
#' category <- list('Tumor', 'Stroma')
#' count_within_many(csd, pairs, radius, category)
#'
#' # Count tumor cells near any T cell in all tissue categories.
#' # Use `phenotype_rules` to define the T cell phenotype
#' pairs <- c('CK+', 'T cell')
#' rules <- list(
#' 'T cell'=c('CD8+', 'FoxP3+'))
#' count_within_many(csd, pairs, radius, phenotype_rules=rules)
#' @md
#' @export
#' @family distance functions
#' @importFrom magrittr "%>%"
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
count_within_many <- function(csd, pairs, radius, category=NA,
                              phenotype_rules=NULL, verbose=TRUE) {
  # count_within_many_impl_rtree gives incorrect results if category
  # is a list than includes both NA and named categories.
  # This is not something we need to support, just disallow it
  if (any(is.na(category)) && !all(is.na(category)))
    stop('Category argument cannot include both NA and named categories.')

  pairs = clean_pairs(pairs)

  all_phenotypes = unique(do.call(c, pairs))
  phenotype_rules = make_phenotype_rules(all_phenotypes, phenotype_rules)

  combos = purrr::cross(list(pair=pairs, category=category))

  # Try to get a name for this field
  field_col = dplyr::if_else('Annotation ID' %in% names(csd),
                             'Annotation ID', 'Sample Name')
  name = ifelse(field_col %in% names(csd),
                csd[[1, field_col]], NA_character_)

  if (verbose) cat('Processing', name, '\n')

  count_within_many_impl(csd, name, combos, radius, phenotype_rules)
}

# Helper functions for count_within_batch and count_within_many
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Clean pairs.
#'
#' Clean pairs.
#' @param pairs Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' clean_pairs(pairs = ...)
#' }
#' @export
clean_pairs = function(pairs) {
  # Allow a single pair to be specified as a plain vector
  if (is.character(pairs) && length(pairs)==2)
    pairs = list(pairs)

  stopifnot(is.list(pairs), length(pairs) > 0)
  pairs
}

#' Helper function for count_within_batch and count_within_many.
#' This does the actual work of calling count_within multiple times and
#' accumulating the result.
#' @param csd Cell seg data for a single field.
#' @param name Name associated with `csd`, for example the basename of the
#' image file.
#' @param combos List of pairs of (from phenotype name, to phenotype name)
#' and tissue category.
#' @param radii Vector of radii.
#' @param phenotype_rules Named list of phenotype rules.
#' @keywords internal
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
count_within_many_impl <- function(csd, name, combos, radii, phenotype_rules) {
  # if (getOption('use.rtree.if.available') &&
  #     requireNamespace('rtree', quietly=TRUE))
  #   counts = count_within_many_impl_rtree(
  #     csd, name, combos, radii, phenotype_rules)
  # else
    counts = count_within_many_impl_dist(
      csd, name, combos, radii, phenotype_rules)

  # Add columns for slide and source
  counts = counts %>%
    tibble::add_column(source=name, .before=1)

  if ('Slide ID' %in% names(csd)) {
    slide = as.character(csd[1, 'Slide ID'])
    counts = counts %>% tibble::add_column(slide_id=slide, .before=1)
  }

  counts
}

#' Distance matrix implementation of count_within_many_impl
#' @param csd Cell seg data for a single field.
#' @param name Name associated with `csd`, for example the basename of the
#' image file.
#' @param combos List of pairs of (from phenotype name, to phenotype name)
#' and tissue category.
#' @param radii Vector of radii.
#' @param phenotype_rules Named list of phenotype rules.
#' @seealso count_within_many_impl
#' @md
#' @keywords internal
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
count_within_many_impl_dist <- function(
    csd, name, combos, radii, phenotype_rules) {
  category = combos %>% purrr::map_chr('category') %>% unique()

  # Subset to what we care about, for faster distance calculation
  if (!anyNA(category))
    csd = csd %>% dplyr::filter(`Tissue Category` %in% category)

  # Compute the distance matrix for these cells
  dst = distance_matrix(csd)

  # Compute counts for each from, to, and category in combos
  row_count = purrr::map_df(combos, function(row) {
    # Call count_within for each item in combos
    # count_within handles multiple radii
    from = row$pair[1]
    from_sel = phenotype_rules[[from]]
    to = row$pair[2]
    to_sel = phenotype_rules[[to]]
    count_within(csd=csd, from=from_sel, to=to_sel,
                 category=row$category,
                 radius=radii, dst=dst) %>%
      # Add columns for from, to, category
      tibble::add_column(
        category = ifelse(is.na(row$category), 'all', row$category),
        from=from,
        to=to,
        .before=1)
  })
}

#' Count cells within a radius for a single field.
#'
#' Count the number of \code{from} cells having a \code{to} cell within
#' \code{radius} microns in tissue category \code{category}.
#' Compute the average number of \code{to} cells
#' within \code{radius} of \code{from} cells.
#'
#' For each \code{from} cell, count the number of \code{to} cells within
#' \code{radius} microns. Report the number of \code{from} cells containing
#' at least \emph{one} \code{to} cell within \code{radius} as \code{from_with}.
#' Report the \emph{average} number of \code{to} cells per
#' \code{from} cell as \code{within_mean}.
#'
#' \code{count_within} counts cells within a single field. It will give an
#' error if run on a merged cell seg data file. To count cells in a merged file,
#' use \code{\link[dplyr]{group_by}} and \code{\link[dplyr]{do}} to call
#' \code{count_within} for each sample in the merged file. See the Examples.
#'
#' There are some subtleties to the count calculation.
#' \itemize{
#'   \item It is not symmetric in \code{from} and \code{to}.
#'   For example the number of tumor cells with a
#'   macrophage within 25 microns is not the same as the number of macrophages
#'   with a tumor cell within 25 microns.
#'   \item \code{from_count*within_mean} is \emph{not} the number of
#'   \code{to} cells within \code{radius} of a \code{from} cell, it may
#'   count \code{to} cells multiple times.
#'   \item Surprisingly, \code{from_count*within_mean} is symmetric in
#'   \code{from} and \code{to}. The double-counting works out.
#' }
#'
#' To aggregate \code{within_mean} across multiple samples (e.g. by Slide ID)
#' see the examples below.
#'
#' If \code{category} is specified, all reported values are for cells within
#' the given tissue category. If \code{category} is NA, values are reported
#' for the entire data set.
#'
#' \code{radius} may be a vector with multiple values.
#'
#' @param csd A data frame with \code{Cell X Position},
#'        \code{Cell Y Position} and \code{Phenotype} columns,
#'        such as the result of calling \code{\link{read_cell_seg_data}}.
#' @param from,to Selection criteria for the
#' rows and columns. Accepts all formats accepted by \code{\link{select_rows}}.
#' @param radius The radius or radii to search within.
#' @param category Optional tissue category to restrict both \code{from} and
#' \code{to}.
#' @param dst Optional distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#'
#' @return A \code{\link{tibble}} with five columns and one row for each
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'   value in \code{radius}:
#'   \describe{
#'    \item{\code{radius}}{The value of \code{radius} for this row.}
#'    \item{\code{from_count}}{The number of \code{from} cells found in
#'     \code{csd}.}
#'    \item{\code{to_count}}{The number of \code{to} cells found in \code{csd}.}
#'    \item{\code{from_with}}{The number of \code{from} cells with a
#'    \code{to} cell within \code{radius}.}
#'    \item{\code{within_mean}}{The average number of \code{to} cells found
#'    within \code{radius} microns of each \code{from} cell.}
#'  }
#' @export
#' @family distance functions
#' @examples
#' library(dplyr)
#' csd <- sample_cell_seg_data
#'
#' # Find the number of macrophages with a tumor cell within 10 or 25 microns
#' count_within(csd, from='CD68+', to='CK+', radius=c(10, 25))
#'
#' # Find the number of tumor cells with a macrophage within 10 or 25 microns
#' count_within(csd, from='CK+', to='CD68+', radius=c(10, 25))
#'
#' \dontrun{
#' # If 'merged' is a merged cell seg file, this will run count_within for
#' # each field:
#' distances = merged %>% group_by(`Slide ID`, `Sample Name`) %>%
#'   do(count_within(., from='CK+', to='CD68+', radius=c(10, 25)))
#'
#' # This will aggregate the fields by Slide ID:
#' distances %>% group_by(`Slide ID`, radius) %>%
#'   summarize(within=sum(from_count*within_mean, na.rm=TRUE),
#'             from_count=sum(from_count),
#'             to_count=sum(to_count),
#'             from_with=sum(from_with),
#'             within_mean=within/from_count) %>%
#'   select(-within)
#' }

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
count_within <- function(csd, from, to, radius, category=NA, dst=NULL) {
  stop_if_multiple_fields(csd)
  stopifnot(length(radius) > 0, all(radius>0))

  # If a category is provided, subset now
  if (!is.na(category)) {
    category_cells = csd$`Tissue Category`==category
    csd = csd[category_cells, ]
    if (!is.null(dst))
      dst = dst[category_cells, category_cells, drop=FALSE]
  }

  if (is.null(dst))
    dst = distance_matrix(csd)

  dst = subset_distance_matrix(csd, dst, from, to)
  if (prod(dim(dst))>0) {
    purrr::map_df(radius, function(rad) {
      within = apply(dst, 1, function(r) sum(r>0 & r<=rad))
      tibble::tibble(
        radius = rad,
        from_count = dim(dst)[1], # Number of from cells
        to_count = dim(dst)[2],   # Number of to cells
        from_with = sum(within>0), # Number of from cells having a
        # to cell within radius
        within_mean = mean(within) # Mean number of to cells within
        # radius of a from cell
      )}
    )
  } else {
    tibble::tibble(
      radius = radius,
      from_count = dim(dst)[1],
      to_count = dim(dst)[2],
      from_with = 0L,
      within_mean = NA
    )
  }
}

make_phenotype_rules = function (phenotypes, existing_rules = NULL)
{
  if (is.null(existing_rules))
    existing_rules = list()
  else if (!is.list(existing_rules) || (length(existing_rules) >
                                        0 && is.null(names(existing_rules))))
    stop("existing_rules must be a named list.")
  existing_names = names(existing_rules)
  extra_names = setdiff(existing_names, phenotypes)
  if (length(extra_names) > 0)
    stop("A rule was given for an unused phenotype: ", paste(extra_names,
                                                             sep = ", "))
  missing_names = setdiff(phenotypes, existing_names)
  new_rules = purrr::set_names(as.list(missing_names))
  c(existing_rules, new_rules)
}

stop_if_multiple_fields=function (csd) {
  col = field_column_(csd)
  if (col %in% names(csd) && length(unique(csd[[col]])) > 1)
    stop("Data contains multiple samples, ", "please select one or set whole_slide=TRUE.")
}

field_column_=function (csd) {
  dplyr::if_else("Annotation ID" %in% names(csd), "Annotation ID",
                 "Sample Name")
}

#' Subset the rows and columns of a distance matrix.
#' @param csd A data frame containing cell segmentation data,
#'        such as the result of
#'        \code{\link{read_cell_seg_data}}.
#' @param dst The distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#' @param row_selection,col_selection Selection criteria for the
#' rows and columns. Accepts all formats accepted by
#' \code{\link{select_rows}}.
#' @return The input matrix \code{dst} subsetted to include only the
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' rows corresponding to \code{row_selection} and columns
#' corresponding to \code{col_selection}.
#' @family distance functions
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
subset_distance_matrix <- function(csd, dst, row_selection, col_selection) {
  # Check for pre-0.1.0.9002 parameter order
  if (is.matrix(csd) && is.data.frame(dst))
    stop(
      'csd and dst parameters to subset_distance_matrix are in the wrong order')

  rows = select_rows(csd, row_selection)
  cols = select_rows(csd, col_selection)
  dst[rows, cols, drop=FALSE]
}

#' get Counts Within
#'
#' @param dfs Function argument documented from the legacy interface.
#' @param ROIIDs Function argument documented from the legacy interface.
#' @param cellDefs Function argument documented from the legacy interface.
#' @param byROI Function argument documented from the legacy interface.
#' @param bySlide Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param mutually_exclusive Function argument documented from the legacy interface.
#' @param cellDefColName Function argument documented from the legacy interface.
#' @param ExclusiveMarkers Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' getCountWithin(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
getCountWithin=function(dfs,ROIIDs,cellDefs,byROI=F,bySlide=T,selectedROIs,mutually_exclusive,cellDefColName="cellDef",ExclusiveMarkers){
  results=list()

  if (mutually_exclusive){
    CTs= sort(cellDefs) # all cell types
    phenoColName=cellDefColName
  }else{
    CTs = sort(ExclusiveMarkers) # all cell types
    phenoColName=ExclusiveMarkers
  }

  headers=c("Sample Name","Cell ID","Tissue Category",phenoColName,"Cell X Position","Cell Y Position")
  pairs <- as.list(as.data.frame(combn(CTs,2)))
  df2bind=Reduce(rbind,pairs);colnames(df2bind)=c("from","to")
  df2bind=data.frame(
    df2bind[rep(1:nrow(df2bind),each=length(radius)),],
    radius= rep(radius,nrow(df2bind)),
    row.names=NULL,check.names = F)
  df2bind=data.frame(category=rep(c("stroma","tumor","all"),each=nrow(df2bind)),
                     df2bind[rep(1:nrow(df2bind),3),],check.names = F)
  results[["ctsWtByROI"]]=results[["ctsWtBySlide"]]=cbind(df2bind,`feature type`="counts within")

  for (ROI in ROIIDs){
    csd=dfs[[ROI]][,headers]
    if(!mutually_exclusive){
      csd[["CT"]]=defineCell(csd,Cpheno = ExclusiveMarkers)
      csd[["CT"]][csd[["CT"]]==""]="other"
      phenoColName="CT"
      csd=csd[,-which(colnames(csd)%in%ExclusiveMarkers)]
    }
    colnames(csd)[colnames(csd)==phenoColName]="Phenotype"
    csd$Phenotype[!csd$Phenotype%in%CTs]="other"
    csd=csd %>% filter(Phenotype!="other")
    tictoc::tic(ROI)
    tmp1=count_within_many(csd, pairs, radius, category)
    tmp2=count_within_many(csd, pairs, radius) %>% mutate(category="all")
    tictoc::toc()
    results[["ctsWtRaw"]][[ROI]]=distances =rbind(tmp1,tmp2)
    if(byROI){
      results[["ctsWtByROI"]][[ROI]]=left_join(df2bind,distances,by=c("category","from","to","radius"))[["within_mean"]]
    }
  }
  if(bySlide){
    selectedROIs=lapply(selectedROIs,function(x) x[x%in%ROIIDs])
    selectedROIs=selectedROIs[sapply(selectedROIs, function(x) length(x)>0)]
    for (slide in names(selectedROIs)){
      ROIs=selectedROIs[[slide]]
      distances = Reduce(rbind,results[["ctsWtRaw"]][ROIs]) %>% select(-`source`)
      distances = distances %>% group_by(category, from,to,radius) %>%
        summarize(within=sum(from_count*within_mean, na.rm=TRUE),
                  from_count=sum(from_count),
                  to_count=sum(to_count),
                  from_with=sum(from_with),
                  within_mean=within/from_count) %>%
        select(-within)
      results[["ctsWtBySlide"]][[slide]]=left_join(df2bind,distances,by=c("category","from","to","radius"))[["within_mean"]]
    }
  }
  return(results)
}

#' get t-test result
#'
#' perform profile-wise t-test by row
#' @param dataGroup1 The group 1 data matrix. Feature by row, sample by column.
#' @param dataGroup2 The group 2 data matrix. Feature by row, sample by column.
#' @return A data frame of mean differences and p-values
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @examples
#' t.result = performTtestsAllClassesOneVsRest(dataMatrix=data_mx,classVector=c("treatment1","treatment2","treatment1","control"));
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
performTtestsAllRows = function(dataGroup1,dataGroup2,...){
  nGroup1 = ncol(dataGroup1)
  nGroup2 = ncol(dataGroup2)
  dataAll = cbind(dataGroup1,dataGroup2)
  tTestWithErrorHandling = function(x){
    if(nGroup1==1){
      testResult = try(t.test(mu=unlist(x[1:nGroup1]),x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE,...);
      if(is.character(testResult)){
        warning(testResult)
        c(NA,NA,NA)
      }else{
        c(testResult$p.value,unlist(x[1:nGroup1]),testResult$estimate)
      }

    }else if(nGroup2==1){
      testResult = try(t.test(x[1:nGroup1],mu=unlist(x[(nGroup1+1):(nGroup1+nGroup2)])),silent=TRUE,...);
      if(is.character(testResult)){
        warning(testResult)
        c(NA,NA,NA)
      }else{
        c(testResult$p.value,testResult$estimate,unlist(x[(nGroup1+1):(nGroup1+nGroup2)]))
      }
    }else{
      testResult = try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)]),silent=TRUE,...);
      if(is.character(testResult)){
        warning(testResult)
        c(NA,NA,NA)
      }else{
        c(testResult$p.value,testResult$estimate)
      }
    }
  }

  results = matrix(unlist(apply(dataAll,1,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
  colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
  rownames(results) = rownames(dataGroup1)
  results
}

#' get paired t-test result by row
#'
#' perform profile-wise t-test by row
#' @param dataGroup1 The group 1 data matrix. Feature by row, sample by column.
#' @param dataGroup2 The group 2 data matrix. Feature by row, sample by column.
#' @return A data frame of mean differences and p-values
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @examples
#' t.result = performTtestsAllClassesOneVsRest(dataMatrix=data_mx,classVector=c("treatment1","treatment2","treatment1","control"));
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
performPairedTtestsAllRows = function(dataGroup1,dataGroup2,...){
  nGroup1 = ncol(dataGroup1)
  nGroup2 = ncol(dataGroup2)
  stopifnot(nGroup1==nGroup2 & nGroup1>1)

  dataAll = cbind(dataGroup1,dataGroup2)
  tTestWithErrorHandling = function(x){
    testResult = try(t.test(x[1:nGroup1],x[(nGroup1+1):(nGroup1+nGroup2)],paired=T),silent=TRUE,...);
    if(is.character(testResult)){
      warning(testResult)
      c(NA,NA,NA)
    }else{
      c(testResult$p.value,testResult$estimate,0)
    }
  }

  results = matrix(unlist(apply(dataAll,1,tTestWithErrorHandling)),ncol=3,byrow=TRUE)
  colnames(results) = c("P.value","Mean.group.1","Mean.group.2")
  rownames(results) = rownames(dataGroup1)
  results
}

#' get t-test result
#'
#' get t-test result by group, one versus the rest
#' @param dataMatrix The data matrix
#' @param classVector vector matching column names to categories
#' @examples
#' t.result = performTtestsAllClassesOneVsRest(dataMatrix,classVector)
#' @return A data frame of mean differences and p-values
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
performTtestsAllClassesOneVsRest = function(dataMatrix,classVector,paired=FALSE,...){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)

  allPvalues = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  allDiffMeans = matrix(NA,nrow=nrow(dataMatrix),ncol=nClasses)
  colnames(allPvalues) = possibleClasses
  rownames(allPvalues) = rownames(dataMatrix)
  colnames(allDiffMeans) = possibleClasses
  rownames(allDiffMeans) = rownames(dataMatrix)

  for(i in 1:nClasses){
    class = possibleClasses[i]
    if(paired){
      resultTest = performPairedTtestsAllRows(dataMatrix[,classVector==class],dataMatrix[,classVector!=class],...)
    }else{
      resultTest = performTtestsAllRows(dataMatrix[,classVector==class],dataMatrix[,classVector!=class],...)
    }
    allPvalues[,i] = resultTest[,1]
    allDiffMeans[,i] = resultTest[,2]-resultTest[,3]
  }
  result = list(allPvalues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}


#' get t-test result by group, pairwise
#'
#' get t-test result by group, pairwise
#' @param dataMatrix The data matrix
#' @param classVector vector matching column names to categories
#' @examples
#' t.result = performTtestsAllClassesEachPair(dataMatrix,classVector)
#' @return A data frame of mean differences and p-values
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
performTtestsAllClassesEachPair = function(dataMatrix,classVector,paired=F,...){
  if(ncol(dataMatrix)!=length(classVector)){
    stop("Number of columns of data matrix must be equal to the length of the class vector")
  }
  possibleClasses = unique(classVector)
  nClasses = length(possibleClasses)

  allPValues = NULL
  allDiffMeans = NULL
  names = NULL
  for(i in 1:(nClasses-1)){
    for(j in (i+1):nClasses){
      class1 = possibleClasses[i]
      class2 = possibleClasses[j]
      names = c(names,paste(class1,class2,sep="."))
      if(paired){
        result = performPairedTtestsAllRows(dataMatrix[,classVector==class1],dataMatrix[,classVector==class2],...)
      }else{
        result = performTtestsAllRows(dataMatrix[,classVector==class1],dataMatrix[,classVector==class2],...)
      }
      allPValues = cbind(allPValues,result[,1])
      allDiffMeans = cbind(allDiffMeans,result[,2] - result[,3])
    }
  }
  colnames(allPValues) = names
  colnames(allDiffMeans) = names
  result = list(allPValues,allDiffMeans)
  names(result) = c("P.Values","Difference.Between.Means")
  return(result)
}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Get Association.
#'
#' Get Association.
#' @param D Function argument documented from the legacy interface.
#' @param ress Function argument documented from the legacy interface.
#' @param preds Function argument documented from the legacy interface.
#' @param ctrlVs Function argument documented from the legacy interface.
#' @param padjMethod Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' getAssociation(D = ..., ress = ...)
#' }
#' @export
getAssociation=function(D,ress,preds,ctrlVs=NULL,padjMethod="bonferroni"){
  assc=data.frame(outcome=character(),exposure=character(),coef=numeric(),pVal=numeric(),padj=numeric())
  r=0
  # tmp=apply(D[,preds,drop=F],2,function(x) length(unique(na.exclude(x[!is.infinite(x)]))))
  tmp=apply(D[,preds,drop=F],2,function(x) sd(x,na.rm = T))
  preds=preds[!tmp%in%c(0,1) & !is.na(tmp)]
  for (res in ress){
    for (pred in preds){
      r=r+1
      formu = paste0(res, " ~ ", paste(c(sprintf("`%s`",pred),ctrlVs),collapse = " + "))
      mod=lm(data = D, as.formula(formu))
      assc[r,c("outcome","exposure","coef","pVal","n")]=
        c(res,pred,
          tryCatch(summary(mod)$coefficients[2,c(1,4)],error = function(e) c(NA,NA)),
          tryCatch(summary(mod)$df[2],error = function(e) c(NA,NA)))
    }
    assc$coef=as.numeric(assc$coef)
    assc$pVal=as.numeric(assc$pVal)
    assc[,"padj"]=p.adjust(assc[,"pVal"],method = padjMethod)#
  }
  return(assc)
}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Get Cor.
#'
#' Get Cor.
#' @param mx Function argument documented from the legacy interface.
#' @param resdf Function argument documented from the legacy interface.
#' @param binSize Function argument documented from the legacy interface.
#' @param method Option controlling how the function runs.
#' @param alternative Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' getCor(mx = ..., resdf = ...)
#' }
#' @export
getCor=function(mx,resdf,binSize,method="pearson",alternative = "two.sided"){
  is=ncol(mx) %/% binSize +1
  corResult=list()
  for(i in 1:is){
    tictoc::tic(i)
    r1=(i-1)*binSize+1
    r2=min(i*binSize,ncol(mx))
    corResult[[i]]=rstatix::cor_test(
      cbind(mx[,r1:r2],resdf),
      vars = colnames(resdf),
      vars2 = colnames(mx)[r1:r2],
      alternative = alternative,
      method = method,
      conf.level = 0.95,
      use = "pairwise.complete.obs"
    )
    tictoc::toc()
  }
  corResult=Reduce(rbind,corResult) %>% slice(order(var1,p))
  return(corResult)
}

#' robust scaling
#' uses median an mad instead of mean and row
#' applies the scaling to the columns (samples) by default
#' @export
#' @param data matrix or data.frame
#' @param dim should rows (1) or columns (2:default) be scaled
#' @param center subract median (default:TRUE)
#' @param scale scale by mad  (default:TRUE)
#' @param preserveScale default TRUE , equalize scales but do not change them
#' @examples
#' library(quantable)
#' tmp = matrix(rep((1:100),times = 4) + rnorm(100*4,0,3),ncol=4)
#' mean = c(20,30,10,40)
#' sd = c(4,3,4,5)
#' tmp = sweep(tmp,2,sd,"*")
#' tmp = sweep(tmp,2,mean,"+")
#' boxplot(tmp)
#' tmp = robustscale(tmp)
#' boxplot(tmp$data)
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraPlotsFun.R (2024-04-09).
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraPlotsFun.R (2024-04-09)
robustscale <- function(data, dim=1, center=TRUE, scale=TRUE, scale_max=Inf,
                        preserveScale = TRUE){
  medians = NULL
  if(center){
    medians <- apply(data,dim,median,na.rm=TRUE)
    data = sweep(data,dim,medians,"-")
  }
  mads=NULL
  if(scale){
    mads <- apply(data,dim, mad,na.rm =TRUE)
    if(preserveScale){
      mads <- mads/mean(mads)
    }
    data = (sweep(data,dim,mads,"/"))
    if(scale_max != Inf){
      data[data > scale_max] <- scale_max
      data[data < scale_max] <- -scale_max
    }
  }
  return(list(data=data,medians=medians,mads=mads))
}

# Shared helpers `uniCoxPh` and `multiCoxPh`
# are defined canonically in 06_survival_analysis.R.

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Get Sig Pair.
#'
#' Get Sig Pair.
#' @param corResult.sigs Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' getSigPair(corResult.sigs = ...)
#' }
#' @export
getSigPair=function(corResult.sigs){
  sigPairs=data.frame()
  sigSingles=data.frame()
  for(resV in names(corResult.sigs)){
    tmpdf=corResult.sigs[[resV]]
    tmpdf%>%filter(`feature type` %in% c("nearest neighbour distance","degree of clustering CW","percentage ratio"))%>%filter()
    tmp=strsplit(gsub(" @.*$","",tmpdf$`Marker Combn`),split = " TO | To | R ")
    ind=sapply(tmp, function(x) length(unique(x)))==2
    tmpdf.pairs=tmpdf[ind,]
    tmpdf.singles=tmpdf[!ind,]
    tmp=sapply(tmp,function(x) paste(sort(unique(x)),collapse = ";"))
    tmp.pairs=tmp[ind]
    tmp.singles=tmp[!ind]
    sigPairs=rbind(sigPairs,data.frame(trait=resV,p=-log10(tmpdf.pairs$p),cor=tmpdf.pairs$cor,comb=tmpdf.pairs$`Marker Combn`,tcat=tmpdf.pairs$`Tissue Category`,pair=tmp.pairs)%>%arrange(desc(p)))
    sigSingles=rbind(sigSingles,data.frame(trait=resV,p=-log10(tmpdf.singles$p),cor=tmpdf.singles$cor,comb=tmpdf.singles$`Marker Combn`,tcat=tmpdf.singles$`Tissue Category`,single=tmp.singles)%>%arrange(desc(p)))
  }
  return(list(sigPairs,sigSingles))
}

bi_ripleys_k_mod=function (mif, mnames, r_range = 0:100, edge_correction = "translation",
                           num_permutations = 50, permute = FALSE, keep_permutation_distribution = FALSE,
                           overwrite = TRUE, workers = 6, big = 1000, nlarge = 1000,
                           xloc = NULL, yloc = NULL)
{
  require(spatialTIME)
  Label = Anchor = Counted = `Exact CSR` = NULL
  if (!inherits(mif, "mif")) {
    stop("Please use a mIF object for mif")
  }
  if (!inherits(mnames, "character") & !inherits(mnames, "data.frame")) {
    stop("Please use either a character vector or data frame of marker combinations for mnames")
  }
  if (!(0 %in% r_range)) {
    r_range = c(0, r_range)
  }
  # out = parallel::mclapply(names(mif$spatial), function(spatial_name){
  out=list()
  n=1
  for (spatial_name in names(mif$spatial)[!names(mif$spatial)%in%out$`ROI ID`]){#
    spat = mif$spatial[[spatial_name]]
    if (is.null(xloc) & is.null(yloc)) {
      spat = spat %>% dplyr::mutate(xloc = (XMin + XMax)/2,
                                    yloc = (YMin + YMax)/2)
    }else {
      spat = spat %>% dplyr::rename(`:=`("xloc", xloc),
                                    `:=`("yloc", yloc))
    }
    win = spatstat.geom::convexhull.xy(
      c(spat$xloc,min(spat$xloc)-10,min(spat$xloc)-10,max(spat$xloc)+10,max(spat$xloc)+10),
      c(spat$yloc,min(spat$yloc)-10,max(spat$yloc)+10,min(spat$yloc)-10,max(spat$yloc)+10))
    # win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)

    area = spatstat.geom::area(win)
    spat = as.matrix(spat[, c("xloc", "yloc", mnames)])
    if (inherits(mnames, "data.frame")) {
      m_combos = mnames
    }
    if (inherits(mnames, "character")) {
      m_combos = expand.grid(anchor = mnames, counted = mnames) %>%
        dplyr::filter(anchor != counted)
    }
    if (!permute) {
      exact_K = spatialTIME:::calculateK(i_dat = spat, j_dat = spat,
                                         anchor = "anchor", counted = "counted", area = area,
                                         win = win, big = big, r_range = r_range, edge_correction = edge_correction,
                                         cores = workers)
    }
    res = parallel::mclapply(1:nrow(m_combos), function(combo) {
      anchor = m_combos[combo, ] %>% dplyr::pull(anchor) %>%
        as.character()
      counted = m_combos[combo, ] %>% dplyr::pull(counted) %>%
        as.character()
      cat(spatial_name, "\t", combo, "\t", anchor, "\t",
          counted, "\n")
      spat_tmp = spat[!(spat[, anchor] == 1 & spat[, counted] ==
                          1), c("xloc", "yloc", anchor, counted)]
      i_dat = spat_tmp[spat_tmp[, 3] == 1, ]
      j_dat = spat_tmp[spat_tmp[, 4] == 1, ]
      if (sum(spat_tmp[, 3]) < 2 | sum(spat_tmp[, 4]) <
          2) {
        final = data.frame(Label = spatial_name, r = r_range,
                           Anchor = anchor, Counted = counted, `Theoretical CSR` = pi *
                             r_range^2, `Observed K` = NA, `Permuted CSR` = NA,
                           `Exact CSR` = NA, check.names = FALSE)
        if (permute) {
          final = final %>% dplyr::full_join(expand.grid(r = r_range,
                                                         iter = seq(num_permutations)), by = "r")
        }
        else {
          final$iter = 1
        }
        return(final)
      }
      K_obs = data.frame(r = r_range, `Theoretical CSR` = pi *
                           r_range^2, check.names = FALSE)
      K_obs$`Observed K` = spatialTIME:::calculateK(i_dat = i_dat, j_dat = j_dat,
                                                    anchor = anchor, counted = counted, area = area,
                                                    win = win, big = big, r_range = r_range, edge_correction = edge_correction,
                                                    cores = workers)
      K_obs$Anchor = anchor
      K_obs$Counted = counted
      if (permute) {
        perm_rows = lapply(seq(num_permutations), function(x) {
          sample(1:nrow(spat), sum(nrow(i_dat), nrow(j_dat)),
                 replace = FALSE)
        })
        kpermed = parallel::mclapply(seq(perm_rows),
                                     function(perm_n) {
                                       cat(perm_n)
                                       perm = perm_rows[[perm_n]]
                                       dat = spat[perm, 1:2]
                                       label = c(rep(anchor, nrow(i_dat)), rep(counted,
                                                                               nrow(j_dat)))
                                       i_dat = dat[label == anchor, ]
                                       j_dat = dat[label == counted, ]
                                       permed = data.frame(r = r_range, `Theoretical CSR` = pi *
                                                             r_range^2, iter = perm_n, check.names = FALSE)
                                       permed$`Permuted CSR` = spatialTIME:::calculateK(i_dat = i_dat,
                                                                                        j_dat = j_dat, anchor = anchor, counted = counted,
                                                                                        area = area, win = win, big = big, r_range = r_range,
                                                                                        edge_correction = edge_correction, cores = workers)
                                       return(permed)
                                     }, mc.preschedule = FALSE, mc.allow.recursive = TRUE) %>%
          do.call(dplyr::bind_rows, .)
        kpermed$`Exact CSR` = NA
      }
      else {
        kpermed = data.frame(r = r_range, `Theoretical CSR` = pi *
                               r_range^2, iter = 1, check.names = FALSE)
        kpermed$`Permuted CSR` = NA
        kpermed$`Exact CSR` = exact_K
      }
      final = dplyr::full_join(K_obs, kpermed, by = c("r",
                                                      "Theoretical CSR")) %>% dplyr::mutate(Label = spatial_name,
                                                                                            .before = 1)
      return(final)
    }) %>% do.call(dplyr::bind_rows, .)
    res = res[, c(1, 2, 7, 5, 6, 3, 4, 8, 9)]
    # return(res)
    out[[spatial_name]]=res
    print(spatial_name)
    print(n)
    n=n+1
  }
  # , mc.cores = workers, mc.preschedule = FALSE, mc.allow.recursive = TRUE) %>%
  out=out %>%do.call(dplyr::bind_rows, .) %>% dplyr::rename(`:=`(!!mif$sample_id, Label))
  save(out,file = sprintf("%s/Mut%s.bi.out.RData",params$OPdir,mutually_exclusive))
  if (!keep_permutation_distribution & permute) {
    out = out %>% dplyr::select(-iter) %>% dplyr::group_by(dplyr::across(mif$sample_id),
                                                           r, Anchor, Counted) %>% dplyr::summarise_all(~mean(.,
                                                                                                              na.rm = TRUE)) %>% dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` -
                                                                                                                                                 `Permuted CSR`, `Degree of Clustering Theoretical` = `Observed K` -
                                                                                                                                                 `Theoretical CSR`, `Exact CSR` = NA) %>% dplyr::mutate(iter = num_permutations,
                                                                                                                                                                                                        .before = Anchor)
  }
  if (overwrite) {
    mif$derived$bivariate_Count = out %>% dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` -
                                                          `Theoretical CSR`, `Degree of Clustering Permutation` = `Observed K` -
                                                          `Permuted CSR`, `Degree of Clustering Exact` = `Observed K` -
                                                          `Exact CSR`) %>% dplyr::mutate(Run = 1)
  }
  if (!overwrite) {
    mif$derived$bivariate_Count = mif$derived$bivariate_Count %>%
      dplyr::bind_rows(out %>% dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` -
                                               `Theoretical CSR`, `Degree of Clustering Permutation` = `Observed K` -
                                               `Permuted CSR`, `Degree of Clustering Exact` = `Observed K` -
                                               `Exact CSR`) %>% dplyr::mutate(Run = ifelse(exists("bivariate_Count",
                                                                                                  mif$derived), max(mif$derived$bivariate_Count$Run) +
                                                                                             1, 1)))
  }
  return(mif)
}


# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Get Mean Count Within.
#'
#' Get Mean Count Within.
#' @param radii Function argument documented from the legacy interface.
#' @param dfs Function argument documented from the legacy interface.
#' @param ROIIDs Function argument documented from the legacy interface.
#' @param cellDefs Function argument documented from the legacy interface.
#' @param Cpheno Function argument documented from the legacy interface.
#' @param byROI Function argument documented from the legacy interface.
#' @param bySlide Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param mutually_exclusive Function argument documented from the legacy interface.
#' @param cellDefColName Function argument documented from the legacy interface.
#' @param OP Function argument documented from the legacy interface.
#' @param aggMethod Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' getMeanCountWithin(radii = ..., dfs = ...)
#' }
#' @export
getMeanCountWithin=function(radii,dfs,ROIIDs,cellDefs,Cpheno,byROI=F,bySlide=T,selectedROIs,mutually_exclusive,cellDefColName="cellDef",OP="OPmIF",aggMethod){
  require(phenoptr)
  require(phenoptrReports)
  results=list()

  CTs= sort(cellDefs)
  if (mutually_exclusive){
    phenoColName=cellDefColName
  }else{
    phenoColName=c(cellDefColName,Cpheno)
    phenoList=lapply(strsplit(cellDefs,"\\+"), function(x) as.list(paste(x,"+",sep="")))
    names(phenoList)=cellDefs
  }
  # results[["meanCWByROI"]]=results[["meanCWBySlide"]]=cbind(df2bind,`feature type`="count within")
  headers=c("Cell ID","Tissue Category",phenoColName,"Cell X Position","Cell Y Position")

  if(byROI){
    csd=Reduce(rbind,dfs)
    df2bind=data.frame(From=rep(CTs,each=length(CTs)),To=rep(CTs,length(CTs)),check.names = F)
    df2bind=data.frame(`Tissue Category`=rep(c("tumor","stroma","all"),each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),3),],check.names = F)
    df2bind=data.frame(`Radius`=rep(radii,each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),length(radii)),],check.names = F)
    df2bind=data.frame(`Sample Name`=rep(names(dfs),each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),length(names(dfs))),],check.names = F)
    if(!mutually_exclusive){
      csd=csd[csd[,cellDefColName]!=""&(!is.na(csd[,cellDefColName])),]
      csd[,Cpheno]=sapply(Cpheno, function(x) ifelse(csd[,x],x,sub("\\+","\\-",x)))
      csd=expr::changeColNames(csd,ind =match(Cpheno,colnames(csd)),newNames = paste("Phenotype",gsub("\\+","",Cpheno)))
      CWs.byROI <- phenoptrReports::count_within_summary(csd,radii = radii,phenotypes =phenoList,.by = "Sample Name",categories=c("tumor","stroma")) %>%
        mutate(`Tissue Category`=tolower(`Tissue Category`))
    }else{
      colnames(csd)[colnames(csd)==cellDefColName]="Phenotype"
      csd$Phenotype[!csd$Phenotype%in%CTs]="other"
      csd=csd %>% filter(Phenotype!="other")
      CWs.byROI <-phenoptrReports::count_within_summary(csd,radii = radii,.by = "Sample Name",categories=c("tumor","stroma")) %>%
        mutate(`Tissue Category`=tolower(`Tissue Category`))
    }
    results[["meanCWByROI.raw"]]=CWs.byROI
    results[["meanCWByROI"]]=left_join(df2bind,CWs.byROI)%>%reshape2::dcast(Radius+`Tissue Category`+From+To~`Sample Name`,value.var = "Within mean")
  }

  if(bySlide){
    csd=Reduce(rbind,dfs[unlist(selectedROIs)])

    df2bind=data.frame(From=rep(CTs,each=length(CTs)),To=rep(CTs,length(CTs)),check.names = F)
    df2bind=data.frame(`Tissue Category`=rep(c("tumor","stroma","all"),each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),3),],check.names = F)
    df2bind=data.frame(`Radius`=rep(radii,each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),length(radii)),],check.names = F)
    df2bind=data.frame(`Slide ID`=rep(unique(csd$`Slide ID`),each=nrow(df2bind)),
                       df2bind[rep(1:nrow(df2bind),length(unique(csd$`Slide ID`))),],check.names = F)


    if(!mutually_exclusive){
      csd=csd[csd[,cellDefColName]!=""&(!is.na(csd[,cellDefColName])),]
      csd[,Cpheno]=sapply(Cpheno, function(x) ifelse(csd[,x],x,sub("\\+","\\-",x)))
      csd=expr::changeColNames(csd,ind =match(Cpheno,colnames(csd)),newNames = paste("Phenotype",gsub("\\+","",Cpheno)))
      CWs.bySlide <- phenoptrReports::count_within_summary(csd,radii = radii,phenotypes =phenoList,.by = "Slide ID",categories=c("tumor","stroma")) %>%
        mutate(`Tissue Category`=tolower(`Tissue Category`))
    }else{
      colnames(csd)[colnames(csd)==cellDefColName]="Phenotype"
      csd$Phenotype[!csd$Phenotype%in%CTs]="other"
      csd=csd %>% filter(Phenotype!="other")
      CWs.bySlide <-phenoptrReports::count_within_summary(csd,radii = radii,.by = "Slide ID",categories=c("tumor","stroma")) %>%
        mutate(`Tissue Category`=tolower(`Tissue Category`))
    }

    results[["meanCWBySlide.raw"]]=CWs.bySlide
    save(CWs.bySlide,file=sprintf("%s/%s_mutuallyExclusive_%s_aggMethod_%s_byROI.RData",OP,"CountWithin.raw",mutually_exclusive,aggMethod))
    results[["meanCWBySlide"]]=left_join(df2bind,CWs.bySlide)%>%reshape2::dcast(Radius+`Tissue Category`+From+To~`Slide ID`,value.var = "Within mean")
  }

  return(results)
}


# =============================================================================
# PLOTS FUNCTIONS
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraPlotsFun.R (2024-04-09)
# =============================================================================


# ggforce
# grid
# gtable::gtable_add_grob
# gridExtra::grid.arrange
# tidyr

#' PieDonut
#'
#' @param data Function argument documented from the legacy interface.
#' @param mapping Function argument documented from the legacy interface.
#' @param mainCol Function argument documented from the legacy interface.
#' @param start Function argument documented from the legacy interface.
#' @param addPieLabel Function argument documented from the legacy interface.
#' @param addDonutLabel Function argument documented from the legacy interface.
#' @param showRatioDonut Function argument documented from the legacy interface.
#' @param showRatioPie Function argument documented from the legacy interface.
#' @param ratioByGroup Function argument documented from the legacy interface.
#' @param showRatioThreshold Function argument documented from the legacy interface.
#' @param labelposition Function argument documented from the legacy interface.
#' @param labelpositionThreshold Function argument documented from the legacy interface.
#' @param r0 Function argument documented from the legacy interface.
#' @param r1 Function argument documented from the legacy interface.
#' @param r2 Function argument documented from the legacy interface.
#' @param explode Function argument documented from the legacy interface.
#' @param selected Function argument documented from the legacy interface.
#' @param explodePos Function argument documented from the legacy interface.
#' @param color Function argument documented from the legacy interface.
#' @param pieAlpha Function argument documented from the legacy interface.
#' @param donutAlpha Function argument documented from the legacy interface.
#' @param maxx Function argument documented from the legacy interface.
#' @param showPieName Function argument documented from the legacy interface.
#' @param showDonutName Function argument documented from the legacy interface.
#' @param title Function argument documented from the legacy interface.
#' @param pieLabelSize Function argument documented from the legacy interface.
#' @param donutLabelSize Function argument documented from the legacy interface.
#' @param titlesize Function argument documented from the legacy interface.
#' @param explodePie Function argument documented from the legacy interface.
#' @param explodeDonut Function argument documented from the legacy interface.
#' @param use.label Function argument documented from the legacy interface.
#' @param use.labels Function argument documented from the legacy interface.
#' @param family Function argument documented from the legacy interface.
#' @param draw Function argument documented from the legacy interface.
#'
#' @import ggplot2
#' @import grid
#' @import tidyr
#' @importFrom gridExtra grid.arrange
#' @importFrom gtable gtable_add_grob
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' PieDonut(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
PieDonut=function(
    data, mapping,
    mainCol=NULL,
    start = getOption("PieDonut.start", 0),
    addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE,showRatioPie = TRUE, ratioByGroup = TRUE,
    showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.02),
    labelposition = getOption("PieDonut.labelposition", 2),
    labelpositionThreshold = 0.1,
    r0 = getOption("PieDonut.r0",0.3),
    r1 = getOption("PieDonut.r1", 1),
    r2 = getOption("PieDonut.r2", 1.2),
    explode = NULL, selected = NULL, explodePos = 0.1,
    color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL,
    showPieName = TRUE, showDonutName = FALSE, title = NULL,
    pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE,
    explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE,
    family = getOption("PieDonut.family", ""),
    draw=F)
{
  require(ggplot2)
  require(ggforce)
  require(grid)
  (cols = colnames(data))
  if (use.labels)
    data = addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping))
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies))
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies))
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts))
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts))
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie)
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label,
                                                             "\n(", scales::percent(df$ratio), ")"), as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]]))
    df[[pies]] <- factor(df[[pies]])
  df
  if(missing(mainCol)|is.null(mainCol)){
    mainCol = gg_color_hue(nrow(df))
  }
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus !=
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 *
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "tidyr::complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]]))
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]]))
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]),
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) -
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2)
        df3$label = paste0(df3$label, "(", df3$ratio,
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio,
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 *
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut)
        df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) +
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) +
        df3$y
      if (labelposition == 2)
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut)
        df3$radius[df3$focus != 0] = df3$radius[df3$focus !=
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0)
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0),
                                    r = as.character(r1), start = "start1", end = "end1",
                                    fill = pies), alpha = pieAlpha, color = color, data = df) +
    transparent() + scale_fill_manual(values = mainCol) +
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy",
                                       xend = "segxend", yend = "segyend"), data = df) +
      geom_text(aes_string(x = "segxend", y = "segyend",
                           label = "label", hjust = "hjust", vjust = "vjust"),
                size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy",
                                       xend = "segxend", yend = "segyend"), data = df[df$ratio <
                                                                                        labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend",
                                                                                                                                          y = "segyend", label = "label", hjust = "hjust",
                                                                                                                                          vjust = "vjust"), size = pieLabelSize, data = df[df$ratio <
                                                                                                                                                                                             labelpositionThreshold, ], family = family) + geom_text(aes_string(x = "labelx",
                                                                                                                                                                                                                                                                y = "labely", label = "label"), size = pieLabelSize,
                                                                                                                                                                                                                                                     data = df[df$ratio >= labelpositionThreshold, ],
                                                                                                                                                                                                                                                     family = family)
  }
  else {
    p1 <- p1 +
      geom_text(aes_string(x = "labelx", y = "labely",
                           label = "label"), size = pieLabelSize, data = df,
                family = family)
    # geom_text_repel(aes_string(x = "labelx", y = "labely",
    #                            label = "label"), size = pieLabelSize, data = df,
    #                 family = family)
  }
  if (showPieName)
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies,
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family),plot.margin = unit(c(0, 0, 0, 0), "pt"))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y",
                                        r0 = as.character(r1), r = as.character(r2),
                                        start = "start1", end = "end1", fill = "no",
                                        explode = "focus"), alpha = donutAlpha, color = color,
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y",
                                        r0 = as.character(r1), r = as.character(r2),
                                        start = "start1", end = "end1", fill = "no"),
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) +
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy",
                                         xend = "segxend", yend = "segyend"), data = df3) +
        geom_text(aes_string(x = "segxend", y = "segyend",
                             label = "label", hjust = "hjust", vjust = "vjust"),
                  size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely",
                                      label = "label"), size = donutLabelSize, data = df3,
                           family = family)
    }
    else {
      p3 <- p3 +
        geom_segment(aes_string(x = "segx", y = "segy",  xend = "segxend", yend = "segyend"),
                     data = df3[df3$ratio1 < labelpositionThreshold, ]) +
        geom_text(aes_string(x = "segxend", y = "segyend", label = "label", hjust = "hjust", vjust = "vjust"),
                  size = donutLabelSize, data = df3[df3$ratio1 < labelpositionThreshold, ], family = family) +
        # geom_text_repel(aes_string(x = "segxend", y = "segyend", label = "label", hjust = "hjust", vjust = "vjust"),
        #           size = donutLabelSize, data = df3[df3$ratio1 < labelpositionThreshold, ], family = family) +
        geom_text(aes_string(x = "labelx", y = "labely", label = "label"),
                  size = donutLabelSize, data = df3[df3$ratio1 >= labelpositionThreshold, ], family = family)
      # geom_text_repel(aes_string(x = "labelx", y = "labely", label = "label"),
      #           size = donutLabelSize, data = df3[df3$ratio1 >= labelpositionThreshold, ], family = family)
    }
    if (!is.null(title))
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title,
                          size = titlesize, family = family)
    else if (showDonutName)
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3,
                          label = donuts, hjust = 0, size = titlesize,
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family),plot.margin = unit(c(0, 0, 0, 0), "pt"))

    # overlay two plots
    g1 <- ggplot_gtable(ggplot_build(p1))
    g2 <- ggplot_gtable(ggplot_build(p3))
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g <- gtable::gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t,pp$l, pp$b, pp$l) # add g2 panel onto g1
    g = gridExtra::arrangeGrob(g)
    if(draw){
      grid.draw(g)
    }
    return(g)
  }else {
    if(draw){p1}
    return(p1)
  }
}

addLabelDf = function (data, mapping = NULL) {
  if (!is.null(mapping)) {
    (mapnames = names(mapping))
    cols = c()
    for (i in 1:length(mapnames)) {
      temp = getMapping(mapping, mapnames[i])
      cols = c(cols, temp)
    }
    cols = unique(cols)
    data[cols] = lapply(data[cols], function(x) to_label(x,
                                                         add.non.labelled = TRUE))
  }
  else {
    data[] = lapply(data, function(x) to_label(x, add.non.labelled = TRUE))
  }
  data
}

getMapping=function (mapping, varname) {
  require(sjmisc)
  if (is.null(mapping))
    return(NULL)
  result = paste(mapping[varname])
  if (result == "NULL")
    result <- NULL
  if (!is.null(result)) {
    if (packageVersion("ggplot2") > "2.2.1") {
      result = stringr::str_replace_all(result, "~", "")
    }
    result = stringr::str_replace_all(result, stringr::fixed("c("),
                                      "")
    result = stringr::str_replace_all(result, stringr::fixed(")"),
                                      "")
    result = stringr::str_replace_all(result, " ", "")
    if (stringr::str_detect(result, ",")) {
      result = unlist(stringr::str_split(result, ","))
    }
  }
  result
}

makeSubColor=function (main, no = 3) {
  result = c()
  for (i in 1:length(main)) {
    temp = ztable::gradientColor(main[i], n = no + 2)[2:(no + 1)]
    result = c(result, temp)
  }
  result
}

transparent=function (size = 0) {
  temp = theme(rect = element_rect(fill = "transparent", size = size),
               panel.background = element_rect(fill = "transparent"),
               panel.border = element_rect(size = size), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())
  temp
}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Get Pie Donut Legend.
#'
#' Get Pie Donut Legend.
#' @param mainCol Function argument documented from the legacy interface.
#' @return An object returned by the function based on the requested query or input.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' getPieDonutLegend(mainCol = ...)
#' }
#' @export
getPieDonutLegend=function(mainCol){
  require(ggplot2)
  tmp=ggplot(data=data.frame(Marker=names(mainCol),Value=1), aes(x=Marker, fill=Marker)) +
    geom_bar() +
    scale_fill_manual(values=mainCol)+
    theme(legend.position = "bottom")
  leg=cowplot::get_legend(tmp)
  return(leg)
}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Plot Pie Donut By Slides.
#'
#' Plot Pie Donut By Slides.
#' @param slides2plot Function argument documented from the legacy interface.
#' @param slideMapSample.dict Function argument documented from the legacy interface.
#' @return A plot object, grob, or side-effect plot generated by the function.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' plotPieDonutBySlides(...)
#' }
#' @export
plotPieDonutBySlides=function(slides2plot,slideMapSample.dict){
  ROIs2plot=unlist(slideMapSample.dict[slides2plot])
  nmax=max(sapply(slideMapSample.dict[slides2plot], length))

  rmax=0
  lays=matrix(NA,ncol = nmax,nrow = 0)
  for(x in slideMapSample.dict[slides2plot]){
    y= rep(NA,nmax)
    y[1:length(x)]=1:length(x)+rmax
    rmax=max(y,na.rm = T)
    lays=rbind(lays,y)
  }
  rownames(lays)=slides2plot
  lays=cbind(1:nrow(lays)+rmax,lays)
  tmp=gridExtra::arrangeGrob(grobs = c(PieDonutPlots[ROIs2plot],lapply(slides2plot, grid::textGrob)), layout_matrix = lays)
  leg=gridExtra::arrangeGrob(getPieDonutLegend(mainCol))
  p=gridExtra::arrangeGrob(tmp,leg,heights=c(10,1))
  return(p)
}

#' plot_immunoflo
#'
#' Plot multiplex immunofluorescence marker summaries.
#' @param mif Function argument documented from the legacy interface.
#' @param plot_title Function argument documented from the legacy interface.
#' @param mnames Function argument documented from the legacy interface.
#' @param mcolors Function argument documented from the legacy interface.
#' @param cell_type Function argument documented from the legacy interface.
#' @param filename Function argument documented from the legacy interface.
#' @param path Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' plot_immunoflo(...)
#' }
plot_immunoflo=function (mif, plot_title, mnames, mcolors = NULL, cell_type = NULL, base_size=7, dot_size=2,
                         filename = NULL, path = NULL)
{
  if (missing(mif))
    stop("MIF is missing; please provide the appropriate data")
  if (!is(mif, "mif"))
    stop("Please use a mif object")
  pb <- dplyr::progress_estimated(length(mif[["spatial"]]))
  plot <- lapply(mif[["spatial"]], function(x) {
    pb$tick()$print()
    plot_data <- x %>% dplyr::select(plot_title, .data$`Cell X Position`,
                                     .data$`Cell Y Position`, !!mnames, cell_type) %>%
      tidyr::pivot_longer(cols = !!mnames, names_to = "marker",
                          values_to = "indicator") %>% dplyr::mutate(xloc = .data$`Cell X Position`, yloc =  .data$`Cell Y Position`,
                                                                     marker = factor(.data$marker, levels = mnames))
    plot_title <- if (length(plot_title) == 1) {
      paste0("ID: ", unique(x[[plot_title]]))
    }
    else {
      paste0("ID: ", paste(unique(x[, plot_title]), collapse = ", "))
    }
    if (is.null(mcolors)) {
      mcolors = RColorBrewer::brewer.pal(length(mnames),
                                         "Paired")
    }
    if (is.null(cell_type)) {
      basic_plot <- plot_data %>% dplyr::filter(.data$indicator ==
                                                  1) %>% ggplot2::ggplot(ggplot2::aes(x = .data$xloc,
                                                                                      y = .data$yloc, color = .data$marker)) + ggplot2::geom_point(data = plot_data[plot_data$indicator ==
                                                                                                                                                                      0, ], color = "gray70") + ggplot2::geom_point(size = dot_size) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(5)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
        ggplot2::ggtitle(plot_title) + ggplot2::scale_color_manual(NULL, values = mcolors, drop = FALSE) + ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    }
    else {
      basic_plot <- plot_data %>% dplyr::filter(.data$indicator == 1) %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$xloc, y = .data$yloc, color = .data$marker, shape = .data[[cell_type]])) +
        ggplot2::geom_point(data = plot_data[plot_data$indicator == 0, ], color = "gray70") + ggplot2::geom_point(size = dot_size) +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(5)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
        ggplot2::ggtitle(plot_title) + ggplot2::scale_color_manual(NULL,
                                                                   values = mcolors, drop = FALSE) + ggplot2::scale_shape_manual(NULL, values = c(3, 16), drop = FALSE) + ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    }
    basic_plot = basic_plot + ggplot2::scale_y_reverse()
    return(basic_plot)
  })
  if (!is.null(filename)) {
    grDevices::pdf(sprintf("%s.pdf", filename), height = 10,
                   width = 10)
    on.exit(dev.off())
    invisible(lapply(seq_along(plot), function(x) {
      print(plot[[x]])
    }))
    grDevices::dev.off()
  }
  mif$derived$spatial_plots = plot
  return(mif)
}

#' subset mif object to the ones having response/outcome data only
#'
#' @param mif Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param resVar Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' hasResVarOnly(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
hasResVarOnly=function(mif,selectedROIs,resVar) {
  tmp=mif$clinical%>%
    filter(mrn %in% mif$clinical$mrn[apply(mif$clinical[resVar],1,function(x) any(!is.na(x)))]) %>%
    select(mrn) %>% tibble::deframe()
  tmp=mif$sample %>% filter(mrn %in% tmp) %>% select(`Slide ID`) %>% tibble::deframe() %>% unique
  selectedROIs=selectedROIs[tmp[tmp %in% names(selectedROIs)]]
  return(selectedROIs)
}

#' Trim mif objects
#'
#' @param mif Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param resVar Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' trimMif4plots(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
trimMif4plots=function(mif,selectedROIs,resVar=NULL){
  if(!is.null(resVar)){
    selectedROIs=hasResVarOnly(mif,selectedROIs,resVar)
  }
  selectedROIs=lapply(selectedROIs, function(x) x[1])
  mif$sample=mif$sample[mif$sample[[mif$sample_id]]%in% unlist(selectedROIs),]
  mif$clinical=mif$clinical[mif$clinical[[mif$patient_id]] %in% mif$sample[[mif$patient_id]],]
  mif$spatial=mif$spatial[unlist(selectedROIs)]
  return(mif)
}


#' plot paired slides
#'
#' @param pairedSlides Function argument documented from the legacy interface.
#' @param mif Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' pairedSlidePlots(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
pairedSlidePlots=function(pairedSlides,mif,selectedROIs){
  selectedROIs=lapply(selectedROIs, function(i) i[i%in%names(mif$spatial)])
  slideCode=pairedSlides_%>%unlist%>%setNames(rep(1,length(.)),.)
  slideCode[is.na(names(slideCode))]=NA
  tmp=c()
  for(i in 1:length(slideCode)){
    tmp[i]=sum(slideCode[1:i],na.rm = T)
  }
  tmp[is.na(names(slideCode))]=NA
  slideCode=setNames(tmp,names(slideCode))
  layout=matrix(slideCode,nrow=nrow(pairedSlides_),ncol = ncol(pairedSlides_),byrow=F)
  slideCode=na.exclude(slideCode)
  tmp=mif$derived$spatial_plots[unlist(selectedROIs[names(slideCode)])]
  legend=cowplot::get_legend(tmp[[1]])
  tmp=lapply(tmp, function(x) x+theme(legend.position = "none",title = element_text(size = 10))+ coord_fixed())
  a=eval(parse(text=sprintf("arrangeGrob(arrangeGrob(%s,layout_matrix =%s),legend, ncol=2, widths=c(9, 0.8))",paste(paste("tmp[[",1:length(tmp),"]]",sep =""),collapse = ","),"layout")))
  return(a)
}

#' plot Immune Highlights
#'
#' @param mif Function argument documented from the legacy interface.
#' @param selectedROIs Function argument documented from the legacy interface.
#' @param featureName Function argument documented from the legacy interface.
#' @param featureType Function argument documented from the legacy interface.
#' @param trait Function argument documented from the legacy interface.
#' @param compartment Function argument documented from the legacy interface.
#' @param FeatureTable Function argument documented from the legacy interface.
#' @param marker Function argument documented from the legacy interface.
#'
#' @return The value returned by the current implementation.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#' @export
#'
#' @examples
#' \dontrun{
#' plotImmuneHighlights(...)
#' }
# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
plotImmuneHighlights=function(mif,selectedROIs,featureName,featureType,trait,compartment,FeatureTable,marker=NULL){
  mif=trimMif4plots(mif,selectedROIs = selectedROIs)
  tmp=sapply(c(featureName,featureType,compartment), is.null)
  fontSize.0=10
  if(any(tmp)){
    message(sprintf("Feature information not complete. %s not provided.",
                    paste(c("featureName","featureType","compartment")[tmp],collapse = ", ")))
    if(is.null(marker)) {stop()}
    else{message("Highlighting provided markers..")}
    orderedSlides=names(selectedROIs)
    title=sprintf(orderedSlides[i])
    fontSize=fontSize.0
  }else{
    marker=unlist(strsplit(trimws(gsub("@.*$","",featureName)),split = " TO | To | R "))

    orderedSlides=dplyr::inner_join(mif$clinical,mif$sample) %>%
      filter(`Slide ID` %in% names(selectedROIs)) %>%
      slice(order(.[[trait]]))%>%
      select(!!c(mif$patient_id),`Slide ID`)%>%tibble::deframe()
    # get trait measurement by slides order
    orderedTrait=(mif$clinical%>%slice(match(names(orderedSlides),mif$clinical$mrn)))[[trait]]
    # get feature value by slides orer
    values=FeatureTable%>%filter(`Marker Combn`==featureName,`feature type`==featureType,`Tissue Category`==compartment)%>%
      select(!!orderedSlides)%>%as.list%>%unlist()
    title=sapply(1:length(orderedSlides),function(i)sprintf(
      "%s\n%s: %s; value:%s",
      orderedSlides[i],trait,format(orderedTrait[i],digit=3),format(values[i],digit=3)))
    fontSize=round(0.75*(fontSize.0))
  }
  mif=plot_immunoflo(mif, plot_title = "ROI ID",  mnames = marker,
                     mcolors = setNames(cols4all::c4a("dark24",length(marker)),marker), base_size = fontSize.0, dot_size = 1.2,
                     cell_type = "Tissue Category")

  slidePlots=mif$derived$spatial_plots

  tmp=slidePlots[unlist(selectedROIs[orderedSlides])]
  legend=cowplot::get_legend(tmp[[1]])
  for(i in 1:length(tmp)){
    a=tmp[[i]]+
      labs( title = title)+
      theme(legend.position = "none",title = element_text(size = fontSize))
    tmp[[i]]=a
  }
  plots=eval(parse(text=sprintf("gridExtra::arrangeGrob(gridExtra::arrangeGrob(%s,nrow=%s),legend, ncol=2, widths=c(12, 1))",paste(paste("tmp[[",1:length(tmp),"]]",sep =""),collapse = ","),round(sqrt(length(tmp))))))
  return(plots)
}
