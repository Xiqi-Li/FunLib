# =============================================================================
# STATISTICAL AND DATA UTILITY FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS — Data utilities:
#   is.empty.data.frame       - Check if data frame has zero rows
#   getmode                   - Return mode (most frequent) value
#   chisq_dist                - Chi-squared distance matrix
#   removeoutlier             - Remove outliers by IQR
#   set_column_as_rownames    - Set one or more columns as row names
#   discrete_quantile         - Quantile normalization within class/batch groups
#   align_df                  - Align data frame rows/columns to reference order
#   column2namedVector        - Convert a column to a named vector
#   filter_assay_info         - Filter assay info by sample list
#   format_infotable          - Standardize info table formatting (columns, values)
#   extract_paired_sample_info - Extract longitudinally paired samples
#   read_info                 - Read + format info table (caches formatted version)
#   read_foundry_datasets     - Read datasets from MD Anderson Foundry API
#   write_foundry_datasets    - Write a dataset to MD Anderson Foundry API
#
# CONTENTS — Statistical functions:
#   uni_cox                   - Univariate Cox regression (p-value or HR)
#   cal_riskscore             - Calculate risk score from glmnet betas
#   survivalSignatures        - Survival-associated signature extraction (glmnet)
#   svmsurvivalSignatures     - Survival-associated signature extraction (SVM)
#   predictGroup              - Nearest-neighbour group prediction
#   compare_two_groups        - Two-group comparison (t-test or Wilcoxon)
#   safe_compare_two_groups   - Two-group comparison with error handling
#   add_p_adjust              - Add adjusted p-values (BH or other) to a table
#   title_case_smart          - Smart title-case formatting (preserves small words)
#   compareFeatures           - Association/correlation testing across features
#
# SOURCE PROVENANCE:
#   functions.R:  woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R (2026-03-02)
#   funcsInPembro_codebase.R: Pembro/Pembro_codebase/funcsInPembro_codebase.R (2026-01-23)
# =============================================================================

# --- SECTION 1: DATA UTILITIES -----------------------------------------------

# SOURCE: functions.R (2026-03-02)
is.empty.data.frame<-function(df){
  return (nrow(df)==0)
}

# SOURCE: functions.R (2026-03-02)
#' Return mode (most frequent) value
#'
#' Return mode (most frequent) value
#' @param v Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' getmode(c("A", "B", "A"))
#' @export
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# SOURCE: functions.R (2026-03-02)
#' Remove outliers by IQR
#'
#' Remove outliers by IQR
#' @param data Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' removeoutlier(data = ...)
#' }
#' @export
removeoutlier<-function(data){
  quartiles<-stats::quantile(data,probs=c(0.25,0.75),na.rm=T)
  IQR<-stats::IQR(data)
  Lower <- quartiles[1] - 1.5*IQR
  Upper <- quartiles[2] + 1.5*IQR
  data_wo_outlier <- subset(data, data > Lower & data < Upper)
  return(data_wo_outlier)
}

# SOURCE: functions.R (2026-03-02)
#' Set one or more columns as row names
#'
#' Set one or more columns as row names
#' @param df Input data object used by the function.
#' @param columns Column name or column selection used by the existing implementation.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' set_column_as_rownames(df = ..., columns = ...)
#' }
#' @export
set_column_as_rownames<-function(df,columns){
  if(length(columns)>1){
    rn=apply(df[,columns],1,function(row){paste(row,collapse = "_")})
  } else {
    rn=df[,columns]
  }
  target_df<-df[,-match(columns,colnames(df)),drop=F]
  rownames(target_df)<-rn
  return(target_df)
}

# SOURCE: functions.R (2026-03-02)
#' Quantile normalization within class/batch groups
#'
#' Quantile normalization within class/batch groups
#' @param counts Input data object used by the function.
#' @param sample_info Input data object used by the function.
#' @param class_factors Function argument documented from the legacy interface.
#' @param batch_factors Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' discrete_quantile(counts = ..., sample_info = ...)
#' }
#' @export
discrete_quantile<-function(counts,sample_info,class_factors,batch_factors){
  assertthat::are_equal(colnames(counts),sample_info[["Sample_ID"]])
  if(length(batch_factors)>1){
    batchfactors<-apply(sample_info[,batch_factors],1,function(row){paste(row,collapse = "_")})
  } else {
    batchfactors=sample_info[,batch_factors]
  }
  if(length(class_factors)>1){
    classfactors<-apply(sample_info[,class_factors],1,function(row){paste(row,collapse = "_")})
  } else {
    classfactors<-sample_info[,class_factors]
  }
  groups=paste(classfactors,batch_factors,sep="_")
  results<-lapply(unique(groups),function(group){sub_counts<-counts[,groups==group];return(limma::normalizeQuantiles(sub_counts))})
  results<-do.call(cbind,results)
  results<-results[,colnames(counts)]
  return(results)
}

# Shared helpers `chisq_dist`, `align_df`, and `column2namedVector`
# are defined canonically in 01_genomic_data.R.

# SOURCE: functions.R (2026-03-02)
#' Filter assay info by sample list
#'
#' Filter assay info by sample list
#' @param assay_info Input data object used by the function.
#' @param sample_info Input data object used by the function.
#' @return A transformed data object returned by the function.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' filter_assay_info(assay_info = ..., sample_info = ...)
#' }
#' @export
filter_assay_info<-function(assay_info,sample_info){
  return(assay_info[assay_info$Sample_ID %in% sample_info$Sample_ID,])
}

# SOURCE: functions.R (2026-03-02)
#' Standardize info table formatting (columns, values)
#'
#' Standardize info table formatting (columns, values)
#' @param infotable Input data object used by the function.
#' @param nchar_ABBV Function argument documented from the legacy interface.
#' @param feature_dictionary Function argument documented from the legacy interface.
#' @param convert_terms Function argument documented from the legacy interface.
#' @param column_content_operation Function argument documented from the legacy interface.
#' @param check_rownames Function argument documented from the legacy interface.
#' @return A transformed data object returned by the function.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' format_infotable(infotable = ..., nchar_ABBV = ...)
#' }
#' @export
format_infotable<-function(infotable,nchar_ABBV=5,feature_dictionary=NULL,convert_terms=F,column_content_operation=list(),check_rownames=F){
  column_names<-colnames(infotable)
  if(!is.null(feature_dictionary)){
    column_names[column_names %in% feature_dictionary[["Other_Term"]]]<-feature_dictionary[["Official_Term"]][match(column_names[column_names %in% feature_dictionary[["Other_Term"]]],feature_dictionary[["Other_Term"]])]
  }
  column_name_dict<-data.frame(original_column_name=column_names,squeezed_column_name=toupper(gsub("[^A-Za-z0-9]","",column_names,perl = T)))
  if(any(duplicated(column_name_dict[["squeezed_column_name"]]))){
    message("duplicated columns detected!")
    duplicated_squeeze_column_names=names(table(column_name_dict[["squeezed_column_name"]])[table(column_name_dict[["squeezed_column_name"]])>1])
    for(duplicated_squeeze_column_name in duplicated_squeeze_column_names){
      cat("Columns:",paste(paste(column_name_dict[["original_column_name"]][which(column_name_dict[["squeezed_column_name"]]==duplicated_squeeze_column_name)],collapse=", ")," ---may be same columns!\n\n",sep=""))
    }
    cat("Please verify the column title!\n")
    stop()
  }
  column_names<-gsub("[^A-Za-z0-9]"," ",column_names,perl = T)
  column_names<-stringr::str_trim(column_names)
  column_names<-gsub(" +","_",column_names)
  column_names<-sapply(column_names,function(column_name){
    if(grepl("_",column_name)){
      words<-unlist(stringr::str_split(column_name,"_"))
      words<-ifelse(nchar(words)<=nchar_ABBV & useful::upper.case(words),words,snakecase::to_upper_camel_case(words))
      paste(words,collapse = "_")
    } else{
      if(nchar(column_name)<=nchar_ABBV){
        ifelse(useful::upper.case(column_name),column_name,snakecase::to_upper_camel_case(column_name))
      } else{
        snakecase::to_upper_camel_case(column_name)
      }
    }
  })
  colnames(infotable)<-column_names
  character_factor_column_names<-column_names[sapply(infotable,class) %in% c("character","factor")]
  for(col in character_factor_column_names){
    column_contents<-infotable[[col]]
    column_content_dict<-data.frame(original_column_contents=column_contents,squeezed_column_contents=toupper(gsub("[^A-Za-z0-9]","",column_contents,perl = T)))
    column_content_table<-table(column_content_dict[["original_column_contents"]],column_content_dict[["squeezed_column_contents"]])
    if(any(colSums(column_content_table!=0)>1)){
      cat("Column ",col," may have wrong labels.\n")
      for(j in (1:ncol(column_content_table))[colSums(column_content_table!=0)>1]){
        wrong_labels<-names(column_content_table[,j][column_content_table[,j]!=0])
        cat(paste(wrong_labels,collapse=", ")," may be same but wrongly labeled.\n")
        base_label<-names(which.max(column_content_table[,j][column_content_table[,j]!=0]))
        if(convert_terms){
          cat("convert terms",paste(setdiff(wrong_labels,base_label),collapse=', '), "into ", base_label,".\n")
          for(other_label in setdiff(wrong_labels,base_label)){column_contents[column_contents==other_label]<-base_label}
        }
      }
    }
    infotable[[col]]<-column_contents
    if(col %in% names(column_content_operation)){
      connector=ifelse(column_content_operation[[col]][2]=="None","_",column_content_operation[[col]][2])
      word_format=column_content_operation[[col]][1]
      transform_word<-c("original"=as.character,"camel"=snakecase::to_upper_camel_case,"upper"=toupper,"lower"=tolower)
      column_contents<-gsub("[^A-Za-z0-9]"," ",column_contents,perl = T)
      column_contents<-stringr::str_trim(column_contents)
      column_contents<-gsub(" +",connector,column_contents)
      column_contents[!is.na(column_contents)]<-sapply(column_contents[!is.na(column_contents)],function(column_content){
        if(grepl(connector,column_content)){
          words<-unlist(stringr::str_split(column_content,connector))
          words<-transform_word[[word_format]](words)
          paste(words,collapse = connector)
        } else{
          transform_word[[word_format]](column_content)
        }
      })
      infotable[[col]]<-column_contents
    }
  }
  if(check_rownames){
    row_names<-rownames(infotable)
    row_name_dict<-data.frame(original_row_name=row_names,squeezed_row_name=toupper(gsub("[^A-Za-z0-9]","",row_names,perl = T)))
    if(any(duplicated(row_name_dict[["squeezed_row_name"]]))){
      message("duplicated columns detected!")
      duplicated_squeeze_row_names=names(table(row_name_dict[["squeezed_row_name"]])[table(row_name_dict[["squeezed_row_name"]])>1])
      for(duplicated_squeeze_row_name in duplicated_squeeze_row_names){
        cat("Columns:",paste(paste(row_name_dict[["original_row_name"]][which(row_name_dict[["squeezed_row_name"]]==duplicated_squeeze_row_name)],collapse=", ")," ---may be same columns!\n",sep=""))
      }
    }
    row_names<-gsub("[^A-Za-z0-9]"," ",row_names,perl = T)
    row_names<-stringr::str_trim(row_names)
    row_names<-gsub(" +","_",row_names)
    row_names<-sapply(row_names,function(row_name){
      if(grepl("_",row_name)){
        words<-unlist(stringr::str_split(row_name,"_"))
        words<-ifelse(nchar(words)<=nchar_ABBV & useful::upper.case(words),words,snakecase::to_upper_camel_case(words))
        paste(words,collapse = "_")
      } else{
        if(nchar(row_name)<=nchar_ABBV){
          ifelse(useful::upper.case(row_name),row_name,snakecase::to_upper_camel_case(row_name))
        } else{
          snakecase::to_upper_camel_case(row_name)
        }
      }
    })
    rownames(infotable)<-row_names
  }
  return(infotable)
}

# SOURCE: functions.R (2026-03-02)
#' Extract longitudinally paired samples
#'
#' Extract longitudinally paired samples
#' @param sample_info Input data object used by the function.
#' @param assay_col Column name or column selection used by the existing implementation.
#' @param patient_id_col Column name or column selection used by the existing implementation.
#' @param sample_id_col Column name or column selection used by the existing implementation.
#' @param timepoint_col Column name or column selection used by the existing implementation.
#' @param timepoint_terms Function argument documented from the legacy interface.
#' @param compact Function argument documented from the legacy interface.
#' @return A transformed data object returned by the function.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' extract_paired_sample_info(sample_info = ..., assay_col = ...)
#' }
#' @export
extract_paired_sample_info<-function(sample_info,assay_col="Assay",patient_id_col="Patient_ID",sample_id_col="Sample_ID",timepoint_col="Timepoint",timepoint_terms=c("Baseline","TP2"),compact=F){
  sample_info<-sample_info[sample_info[[timepoint_col]] %in% timepoint_terms,]
  paired_samples<-NULL
  paired_patients<-NULL
  for(assay in unique(sample_info[[assay_col]])){
    assay_sample_info<-sample_info[sample_info[[assay_col]]==assay,]
    assay_sample_info_<-data.frame(Patient_ID=assay_sample_info[[patient_id_col]],Sample_ID=assay_sample_info[[sample_id_col]],Timepoint=assay_sample_info[[timepoint_col]])
    assay_sample_info_w<-tidyr::pivot_wider(assay_sample_info_,id_cols="Patient_ID",names_from = "Timepoint",values_from = "Sample_ID",values_fill = NA,values_fn = function(l){paste(l, collapse=";")})
    if(any(rowSums(!is.na(assay_sample_info_w[,-1]))==2)){
      assay_sample_info_w<-assay_sample_info_w[rowSums(!is.na(assay_sample_info_w[,-1]))==2,]
      assay_paired_samples<-unname(unlist(assay_sample_info_w[,timepoint_terms]))
      assay_paired_patients<-unique(assay_sample_info_w[[patient_id_col]])
      paired_samples<-c(paired_samples,assay_paired_samples)
      paired_patients<-c(paired_patients,assay_paired_patients)
    }
  }
  paired_sample_info<-sample_info[sample_info[[sample_id_col]] %in% paired_samples,]
  if(compact){
    paired_patients<-names(table(paired_patients)==length(unique(sample_info[[assay_col]])))
  } else {
    paired_patients<-unique(paired_patients)
  }
  return(paired_sample_info[paired_sample_info[[patient_id_col]] %in% paired_patients,])
}

# SOURCE: functions.R (2026-03-02)
#' Read + format info table (caches formatted version)
#'
#' Read + format info table (caches formatted version)
#' @param origin_dir Function argument documented from the legacy interface.
#' @param origin_info_file Function argument documented from the legacy interface.
#' @param origin_sep Function argument documented from the legacy interface.
#' @param destination_dir Function argument documented from the legacy interface.
#' @param destination_info_file Function argument documented from the legacy interface.
#' @param destination_sep Function argument documented from the legacy interface.
#' @param ... Additional arguments passed through to downstream functions.
#' @return The object returned after reading the requested input.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' read_info(origin_dir = ..., origin_info_file = ...)
#' }
#' @export
read_info<-function(origin_dir,origin_info_file,origin_sep,destination_dir,destination_info_file,destination_sep,...){
  if(!dir.exists(destination_dir)){
    message("No destination directory exists!\n")
    message(paste(destination_dir," will be created!",sep=""))
    dir.create(destination_dir,recursive = T)
  }
  if(!file.exists(file.path(destination_dir,destination_info_file))){
    if(missing(origin_sep)){
      if(endsWith(origin_info_file,".csv")){
        origin_info<-read.csv(file.path(origin_dir,origin_info_file),header=T,stringsAsFactors = F,check.names = F)
      } else{
        origin_info<-read.table(file.path(origin_dir,origin_info_file),sep=origin_sep,header = T,stringsAsFactors = F,check.names = F)
      }
    } else{
      if(endsWith(origin_info_file,".csv")){
        origin_info<-read.csv(file.path(origin_dir,origin_info_file),header=T,sep=origin_sep,stringsAsFactors = F,check.names = F)
      } else{
        origin_info<-read.table(file.path(origin_dir,origin_info_file),sep=origin_sep,header = T,stringsAsFactors = F,check.names = F)
      }
    }
    origin_info<-format_infotable(infotable = origin_info,...)
    if(missing(destination_sep)){
      write.csv(origin_info,file = file.path(destination_dir,destination_info_file),row.names = F,quote = F)
    } else {
      write.table(origin_info,file = file.path(destination_dir,destination_info_file),sep=destination_sep,row.names = F,quote = F)
    }
    destination_info<-origin_info
  } else {
    if(missing(destination_sep)){
      if(endsWith(destination_info_file,".csv")){
        destination_info<-read.csv(file.path(destination_dir,destination_info_file),header=T,stringsAsFactors = F,check.names = F)
      } else {
        destination_info<-read.table(file.path(destination_dir,destination_info_file),header=T,stringsAsFactors = F,sep=destination_sep,check.names = F)
      }
    } else{
      if(endsWith(destination_info_file,".csv")){
        destination_info<-read.csv(file.path(destination_dir,destination_info_file),header=T,sep=destination_sep,stringsAsFactors = F,check.names = F)
      } else {
        destination_info<-read.table(file.path(destination_dir,destination_info_file),header=T,stringsAsFactors = F,sep=destination_sep,check.names = F)
      }
    }
  }
  return(destination_info)
}

# SOURCE: functions.R (2026-03-02)
#' Read datasets from MD Anderson Foundry API
#'
#' Read datasets from MD Anderson Foundry API
#' @param foundry.hostname Function argument documented from the legacy interface.
#' @param foundry.token Function argument documented from the legacy interface.
#' @param rid Function argument documented from the legacy interface.
#' @param aliases.yml Function argument documented from the legacy interface.
#' @return The object returned after reading the requested input.
#' @details Source provenance: foundry utility addition (2026-04-09).
#'
#' @examples
#' \dontrun{
#' read_foundry_datasets(foundry.hostname = ..., foundry.token = ...)
#' }
#' @export
read_foundry_datasets<-function(foundry.hostname="foundry.mdanderson.edu",foundry.token,rid,aliases.yml="~/.foundry/aliases.yml"){
  stopifnot("token missing!"=!missing(foundry.token))
  options(foundry.hostname = foundry.hostname,foundry.token = foundry.token)
  if(!dir.exists("~/.foundry")){
    dir.create("~/.foundry",recursive = T)
    file.create("~/.foundry/aliases.yml")
  }
  if(!missing(rid)){
    content=paste("foundry_data:\n    rid: ",rid,sep="")
    cat(content,file="~/.foundry/aliases.yml")
  } else {
    if(!file.exists("~/.foundry/aliases.yml")){
      file.copy(aliases.yml,"~/.foundry/aliases.yml")
    }
  }
  aliases<-yaml::read_yaml("~/.foundry/aliases.yml")
  datasets<-lapply(names(aliases), function(aliase){
    df<-as.data.frame(foundry::datasets.read_table(aliase))
    df
  })
  names(datasets)<-names(aliases)
  return(datasets)
}

# SOURCE: foundry utility addition (2026-04-09)
#' Write a dataset to MD Anderson Foundry API
#'
#' Write a dataset to MD Anderson Foundry API
#' @param df Input data object used by the function.
#' @param foundry.hostname Function argument documented from the legacy interface.
#' @param foundry.token Function argument documented from the legacy interface.
#' @param alias Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' write_foundry_datasets(df = ..., foundry.hostname = ...)
#' }
#' @export
write_foundry_datasets<-function(df,foundry.hostname="foundry.mdanderson.edu",foundry.token,alias){
  stopifnot("token missing!"=!missing(foundry.token))
  options(foundry.hostname = foundry.hostname,foundry.token = foundry.token)
  foundry::datasets.write_table(df,alias)
}

# --- SECTION 2: SURVIVAL STATISTICS ------------------------------------------

# SOURCE: functions.R (2026-03-02)
#' Univariate Cox regression (p-value or HR)
#'
#' Univariate Cox regression (p-value or HR)
#' @param time Function argument documented from the legacy interface.
#' @param status Function argument documented from the legacy interface.
#' @param covariate Function argument documented from the legacy interface.
#' @param returnvalue Option controlling how the function runs.
#' @return The result object produced by the survival analysis.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' uni_cox(time = ..., status = ...)
#' }
#' @export
uni_cox<-function(time,status,covariate,returnvalue=c("pvalue","coef","both")){
  if(is.null(returnvalue)){
    returnvalue="pvalue"} else{
      returnvalue<-match.arg(returnvalue)
    }
  fit<-survival::coxph(survival::Surv(time,status)~covariate)
  if(returnvalue=="pvalue"){
    values=summary(fit)$logtest["pvalue"]
  } else  if(returnvalue=="coef"){
    values=summary(fit)$coefficients[1,"coef"]
  } else{
    values=c(summary(fit)$coefficients[1,"coef"],summary(fit)$logtest["pvalue"])
    names(values)=c("coef",'pvalue')
  }
  return(values)
}

# SOURCE: functions.R (2026-03-02)
#' Calculate risk score from glmnet betas
#'
#' Calculate risk score from glmnet betas
#' @param selected_betas Function argument documented from the legacy interface.
#' @param covariates Input data object used by the function.
#' @return The value returned by the current implementation.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' cal_riskscore(selected_betas = ..., covariates = ...)
#' }
#' @export
cal_riskscore<-function(selected_betas,covariates){
  selected_covariates<-covariates[names(selected_betas),]
  risk_scores<-colSums((selected_covariates-apply(selected_covariates,1,mean))*selected_betas)
  return(risk_scores)
}

# SOURCE: functions.R (2026-03-02)
#' Survival-associated signature extraction (glmnet)
#'
#' Survival-associated signature extraction (glmnet)
#' @param time Function argument documented from the legacy interface.
#' @param status Function argument documented from the legacy interface.
#' @param covariates Input data object used by the function.
#' @param uni_ret Function argument documented from the legacy interface.
#' @param beta_filter Function argument documented from the legacy interface.
#' @param uni_cox_pvalue_cutoff Function argument documented from the legacy interface.
#' @param uni_cox_coef_cutoff Function argument documented from the legacy interface.
#' @param glmnet_alpha Function argument documented from the legacy interface.
#' @param scale Function argument documented from the legacy interface.
#' @return The result object produced by the survival analysis.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' survivalSignatures(time = ..., status = ...)
#' }
#' @export
survivalSignatures<-function(time,status,covariates,uni_ret="pvalue",beta_filter=0,uni_cox_pvalue_cutoff=0.1,uni_cox_coef_cutoff=0,glmnet_alpha=0,scale=T){
  uni_cox_results_<-apply(covariates,1,function(row){uni_cox(time=time,status=status,covariate=row,returnvalue = uni_ret)})
  uni_cox_results<-as.data.frame(matrix(uni_cox_results_,nrow=nrow(covariates),byrow=T))
  rownames(uni_cox_results)<-rownames(covariates)
  if(uni_ret=="both"){
    coln_uni_cox_results<-c("coef","pvalue")
  } else{
    coln_uni_cox_results<-uni_ret
  }
  colnames(uni_cox_results)<-coln_uni_cox_results
  if("pvalue" %in% colnames(uni_cox_results)){
    uni_cox_results<-uni_cox_results[uni_cox_results[,"pvalue"]<=uni_cox_pvalue_cutoff,,drop=F]
  }
  if("coef" %in% colnames(uni_cox_results)){
    uni_cox_results<-uni_cox_results[abs(uni_cox_results[,"coef"])>uni_cox_coef_cutoff,,drop=F]
  }
  sub_covariates<-t(covariates[rownames(uni_cox_results),])
  survival_mat<-as.matrix(data.frame(time=time,status=status))
  fit<-glmnet::cv.glmnet(x=sub_covariates,y=survival_mat,family="cox",alpha=glmnet_alpha)
  betas<-stats::coef(fit,s="lambda.min")[,1]
  selected_betas<-betas[abs(betas)>beta_filter]
  names(selected_betas)<-names(betas)[abs(betas)>beta_filter]
  selected_covariates<-covariates[names(selected_betas),]
  risk_scores<-colSums((selected_covariates-apply(selected_covariates,1,mean))*selected_betas)
  require(survival)
  roc<-timeROC::timeROC(T=time,delta=status,marker=risk_scores,cause=1,weighting="marginal",times=quantile(time,seq(0.1,0.9,0.1)),iid=TRUE)
  sen_and_spe<-roc$TP+(1-roc$FP)
  optimal_cut<-c(sort(unique(risk_scores)),Inf)[which.max(rowSums(sen_and_spe))]
  risk_order<-order(risk_scores)
  beta_order<-order(selected_betas)
  top_anno<-ComplexHeatmap::HeatmapAnnotation(status=status[risk_order],time=ComplexHeatmap::anno_barplot(time[risk_order]),riskscore=ComplexHeatmap::anno_barplot(risk_scores[risk_order]))
  right_anno<-ComplexHeatmap::rowAnnotation(beta=ComplexHeatmap::anno_barplot(selected_betas[beta_order]))
  if(scale){
    heatmap_data<-t(scale(t(selected_covariates[beta_order,risk_order])))
  } else{
    heatmap_data<-selected_covariates[beta_order,risk_order]
  }
  h_map<-ComplexHeatmap::Heatmap(heatmap_data,name="expression",cluster_rows=F,cluster_columns=F,top_annotation=top_anno,right_annotation=right_anno)
  return(list(uni_cox_results=uni_cox_results_,selected_betas=selected_betas,selected_covariates=selected_covariates,risk_scores=risk_scores,roc=roc,optimal_cut=optimal_cut,heatmap=h_map))
}

# SOURCE: functions.R (2026-03-02)
#' Nearest-neighbour group prediction
#'
#' Nearest-neighbour group prediction
#' @param sample_feature Input data object used by the function.
#' @param group1_features Input data object used by the function.
#' @param group2_features Input data object used by the function.
#' @param sample_groups Function argument documented from the legacy interface.
#' @param dist_method Function argument documented from the legacy interface.
#' @param aggregate_method Function argument documented from the legacy interface.
#' @param permutation Function argument documented from the legacy interface.
#' @param resample_size Function argument documented from the legacy interface.
#' @param repeats Function argument documented from the legacy interface.
#' @return A transformed data object returned by the function.
#' @details Source provenance: functions.R (2026-03-02).
#'
#' @examples
#' \dontrun{
#' predictGroup(sample_feature = ..., group1_features = ...)
#' }
#' @export
predictGroup<-function(sample_feature,group1_features,group2_features,sample_groups=c("Group1","Group2"),dist_method="euclidean",aggregate_method="average",permutation=F,resample_size=0.8,repeats=100){
  cat("make sure predicted sample is not in the group1 nor the group2!\n\n")
  cat("make sure all data normalized at same scale!\n\n")
  groups=c()
  if(!permutation) {
    dist1<- dist(rbind(sample_feature,t(group1_features)),method = dist_method)
    dist1_<-switch (aggregate_method,"average"=mean(dist1[1:ncol(group1_features)]),"max"=max(dist1[1:ncol(group1_features)]),"min"=min(dist1[1:ncol(group1_features)]))
    dist2<- dist(rbind(sample_feature,t(group2_features)),method = dist_method)
    dist2_<-switch (aggregate_method,"average"=mean(dist2[1:ncol(group2_features)]),"max"=max(dist2[1:ncol(group2_features)]),"min"=min(dist2[1:ncol(group2_features)]))
    if(dist1_>=dist2_){
      groups<-c(groups,"Group2")
    }  else {
      groups<-c(groups,"Group1")
    }
    return(groups)
  } else {
    for(i in 1:repeats){
      n_sample<-min(round(ncol(group1_features)*resample_size,0),round(ncol(group2_features)*resample_size,0))
      if(ncol(group1_features)<=3){
        group1_features_<-group1_features
      } else {
        group1_features_<-group1_features[,sample(1:ncol(group1_features),n_sample,replace = T)]
      }
      if(ncol(group2_features)<=3){
        group2_features_<-group2_features
      } else {
        group2_features_<-group2_features[,sample(1:ncol(group2_features),n_sample,replace = T)]
      }
      dist1<- dist(rbind(sample_feature,t(group1_features_)),method = dist_method)
      dist1_<-switch (aggregate_method,"average"=mean(dist1[1:ncol(group1_features_)]),"max"=max(dist1[1:ncol(group1_features_)]),"min"=min(dist1[1:ncol(group1_features_)]))
      dist2<- dist(rbind(sample_feature,t(group2_features_)),method = dist_method)
      dist2_<-switch (aggregate_method,"average"=mean(dist2[1:ncol(group2_features_)]),"max"=max(dist2[1:ncol(group2_features_)]),"min"=min(dist2[1:ncol(group2_features_)]))
      if(dist1_>=dist2_){
        groups<-c(groups,sample_groups[2])
      }  else {
        groups<-c(groups,sample_groups[1])
      }
    }
    group=names(table(groups))[which.max(table(groups))]
    p=binom.test(x=table(groups)[group],n=repeats,p=0.5)$p.value
    return(list(Group=group,pval=unname(p)))
  }
}

# --- SECTION 3: PEMBRO-DERIVED STATISTICAL UTILITIES -------------------------
# SOURCE: Pembro/Pembro_codebase/funcsInPembro_codebase.R (2026-01-23)

# SOURCE: funcsInPembro_codebase.R (2026-01-23)
#' Smart title-case formatting (preserves small words)
#'
#' Smart title-case formatting (preserves small words)
#' @param x Function argument documented from the legacy interface.
#' @param small_words Function argument documented from the legacy interface.
#' @return The value returned by the current implementation.
#' @details Source provenance: Pembro/Pembro_codebase/funcsInPembro_codebase.R (2026-01-23).
#'
#' @examples
#' \dontrun{
#' title_case_smart(x = ..., small_words = ...)
#' }
#' @export
title_case_smart <- function(
    x,
    small_words = c("a","an","and","as","at","but","by","for","in","nor",
                    "of","on","or","per","so","the","to","via","with")
) {
  words <- strsplit(x, " ")[[1]]
  out <- mapply(function(w, i) {
    lw <- tolower(w)
    if (i == 1 || !lw %in% small_words) {
      paste0(toupper(substring(lw, 1, 1)), substring(lw, 2))
    } else {
      lw
    }
  }, words, seq_along(words))
  paste(out, collapse = " ")
}

# SOURCE: funcsInPembro_codebase.R (2026-01-23)
#' Add adjusted p-values (BH or other) to a table
#'
#' Add adjusted p-values (BH or other) to a table
#' @param df Input data object used by the function.
#' @param p_col Column name or column selection used by the existing implementation.
#' @param method Option controlling how the function runs.
#' @param group_cols Column name or column selection used by the existing implementation.
#' @return A transformed data object returned by the function.
#' @details Source provenance: funcsInPembro_codebase.R (2026-01-23).
#'
#' @examples
#' df <- data.frame(group = c("A", "A", "B"), p_value = c(0.01, 0.2, 0.03))
#' add_p_adjust(df, group_cols = "group")
#' @export
add_p_adjust <- function(df, p_col = "p_value", method = "BH", group_cols = NULL) {
  if (is.null(group_cols)) {
    df[[paste0(p_col, "_adj")]] <- p.adjust(df[[p_col]], method = method)
  } else {
    df <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::mutate("{p_col}_adj" := p.adjust(.data[[p_col]], method = method)) %>%
      dplyr::ungroup()
  }
  df
}

# SOURCE: funcsInPembro_codebase.R (2026-01-23)
#' Two-group comparison (t-test or Wilcoxon)
#'
#' Two-group comparison (t-test or Wilcoxon)
#' @param dat Input data object used by the function.
#' @param value_col Column name or column selection used by the existing implementation.
#' @param group_col Column name or column selection used by the existing implementation.
#' @param method Option controlling how the function runs.
#' @return A transformed data object returned by the function.
#' @details Source provenance: funcsInPembro_codebase.R (2026-01-23).
#'
#' @examples
#' \dontrun{
#' compare_two_groups(dat = ..., value_col = ...)
#' }
#' @export
compare_two_groups <- function(dat, value_col = "RankChange", group_col = "clinical_benefit",
                               method = c("t.test", "wilcox")) {
  method <- match.arg(method)
  groups <- unique(dat[[group_col]])
  if (length(groups) != 2) stop("Exactly 2 groups required.")
  g1 <- dat[[value_col]][dat[[group_col]] == groups[1]]
  g2 <- dat[[value_col]][dat[[group_col]] == groups[2]]
  if (method == "t.test") {
    res <- t.test(g1, g2)
  } else {
    res <- wilcox.test(g1, g2, exact = FALSE)
  }
  data.frame(
    group1 = groups[1], group2 = groups[2],
    mean1 = mean(g1, na.rm = TRUE), mean2 = mean(g2, na.rm = TRUE),
    n1 = sum(!is.na(g1)), n2 = sum(!is.na(g2)),
    statistic = unname(res$statistic),
    p_value = res$p.value,
    method = method,
    stringsAsFactors = FALSE
  )
}

# SOURCE: funcsInPembro_codebase.R (2026-01-23)
#' Two-group comparison with error handling
#'
#' Two-group comparison with error handling
#' @param dat Input data object used by the function.
#' @param value_col Column name or column selection used by the existing implementation.
#' @param group_col Column name or column selection used by the existing implementation.
#' @param method Option controlling how the function runs.
#' @return A transformed data object returned by the function.
#' @details Source provenance: funcsInPembro_codebase.R (2026-01-23).
#'
#' @examples
#' \dontrun{
#' safe_compare_two_groups(dat = ..., value_col = ...)
#' }
#' @export
safe_compare_two_groups <- function(dat, value_col, group_col, method = "wilcox") {
  tryCatch(
    compare_two_groups(dat, value_col = value_col, group_col = group_col, method = method),
    error = function(e) {
      data.frame(group1 = NA, group2 = NA, mean1 = NA, mean2 = NA,
                 n1 = NA, n2 = NA, statistic = NA, p_value = NA,
                 method = method, stringsAsFactors = FALSE)
    }
  )
}

# SOURCE: funcsInPembro_codebase.R (2026-01-23)
# Association testing (lm/logistic) or correlation (rstatix) across features x outcomes.
# Supports groupCol for stratified analysis and padjByOutcome for per-outcome FDR.
# SOURCE: funcsInPembro_codebase.R (2026-01-23)
#' Association/correlation testing across features
#'
#' Association/correlation testing across features
#' @param D Input data object used by the function.
#' @param features Function argument documented from the legacy interface.
#' @param outcomes Function argument documented from the legacy interface.
#' @param method Option controlling how the function runs.
#' @param ctrlVs Function argument documented from the legacy interface.
#' @param padjMethod Option controlling how the function runs.
#' @param padjByOutcome Function argument documented from the legacy interface.
#' @param corMethod Option controlling how the function runs.
#' @param alternative Function argument documented from the legacy interface.
#' @param binSize Function argument documented from the legacy interface.
#' @param groupCol Function argument documented from the legacy interface.
#' @return A transformed data object returned by the function.
#' @details Source provenance: funcsInPembro_codebase.R (2026-01-23).
#'
#' @examples
#' \dontrun{
#' compareFeatures(D = ..., features = ...)
#' }
#' @export
compareFeatures <- function(D,
                            features,
                            outcomes,
                            method = c("association","correlation"),
                            ctrlVs = NULL,
                            padjMethod = "BH",
                            padjByOutcome = FALSE,
                            corMethod = "pearson",
                            alternative = "two.sided",
                            binSize = 500,
                            groupCol = NULL) {

  method <- match.arg(method)

  # --- sanitize column names ---
  origNames <- colnames(D)
  safeNames <- make.names(origNames, unique = TRUE)
  nameMap <- setNames(origNames, safeNames)
  colnames(D) <- safeNames

  features <- make.names(features)
  outcomes <- make.names(outcomes)
  ctrlVs <- if (!is.null(ctrlVs)) make.names(ctrlVs) else NULL
  if (!is.null(groupCol)) groupCol <- make.names(groupCol)

  # helper to map back original col names in result columns
  nameMapFn <- function(x) ifelse(x %in% names(nameMap), unname(nameMap[x]), x)

  # -------- association helper --------
  run_assoc <- function(subD) {
    tmp_res <- data.frame()
    for (out in outcomes) {
      for (pred in features) {
        cols_to_use <- c(out, pred, ctrlVs)
        if (!all(cols_to_use %in% colnames(subD))) next

        df <- subD[, cols_to_use, drop = FALSE]

        df[[out]]  <- as.numeric(df[[out]])
        df[[pred]] <- as.numeric(df[[pred]])
        if (!is.null(ctrlVs)) {
          for (ctrl in ctrlVs) df[[ctrl]] <- as.numeric(df[[ctrl]])
        }

        df <- df[complete.cases(df[, cols_to_use, drop = FALSE]), , drop = FALSE]
        if (nrow(df) < 3) next

        if (all(df[[out]] %in% c(0, 1))) {
          mod <- glm(
            as.formula(paste0(out, " ~ ", paste(c(pred, ctrlVs), collapse = " + "))),
            data = df, family = binomial()
          )
          s <- summary(mod)
          coef_val <- s$coefficients[2, 1]
          p_val    <- s$coefficients[2, 4]
          ci <- tryCatch(stats::confint(mod, parm = 2), error = function(e) c(NA_real_, NA_real_))
          method_str <- "logistic"
          R2_val <- NA_real_
          SD_val <- NA_real_
        } else {
          mod <- lm(
            as.formula(paste0(out, " ~ ", paste(c(pred, ctrlVs), collapse = " + "))),
            data = df
          )
          s <- summary(mod)
          coef_val <- s$coefficients[2, 1]
          p_val    <- s$coefficients[2, 4]
          ci <- tryCatch(stats::confint(mod, level = 0.95)[2, ], error = function(e) c(NA_real_, NA_real_))
          method_str <- "lm"
          R2_val <- s$r.squared
          SD_val <- s$sigma
        }

        tmp_res <- rbind(tmp_res, data.frame(
          outcome   = out,
          exposure  = pred,
          n_total   = nrow(df),
          n         = nrow(df),
          coef      = coef_val,
          conf_low  = ci[1],
          conf_high = ci[2],
          pVal      = p_val,
          R2        = R2_val,
          SD_resid  = SD_val,
          method    = method_str,
          stringsAsFactors = FALSE
        ))
      }
    }
    tmp_res
  }

  # -------- correlation helper --------
  run_corr <- function(subD) {
    vars_out <- outcomes[outcomes %in% colnames(subD)]
    if (length(vars_out) == 0L) return(data.frame())

    nCols <- length(features)
    nBins <- ceiling(nCols / binSize)
    corRes <- vector("list", nBins)
    idx <- 0L

    feats_present <- features[features %in% colnames(subD)]
    if (length(feats_present) == 0L) return(data.frame())

    for (i in seq_len(nBins)) {
      r1 <- (i - 1L) * binSize + 1L
      r2 <- min(i * binSize, length(feats_present))

      vars2 <- feats_present[r1:r2]
      vars2 <- setdiff(vars2, vars_out)
      if (length(vars2) == 0L) next

      df <- subD[, c(vars_out, vars2), drop = FALSE]

      for (nm in colnames(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))

      cr <- tryCatch(
        rstatix::cor_test(
          df,
          vars        = vars_out,
          vars2       = vars2,
          method      = corMethod,
          alternative = alternative,
          conf.level  = 0.95,
          use         = "pairwise.complete.obs"
        ),
        error = function(e) data.frame()
      )

      if (nrow(cr) == 0) next

      if ("n" %in% names(cr)) {
        cr <- dplyr::rename(cr, n_total = n)
      } else if ("n.obs" %in% names(cr)) {
        cr <- dplyr::rename(cr, n_total = n.obs)
      } else if (!("n_total" %in% names(cr))) {
        cr$n_total <- NA_integer_
      }

      idx <- idx + 1L
      corRes[[idx]] <- cr
    }

    if (idx == 0L) return(data.frame())
    Reduce(rbind, corRes[seq_len(idx)])
  }

  # -------- main execution --------
  res <- data.frame()

  if (method == "association") {
    if (!is.null(groupCol) && groupCol %in% colnames(D)) {
      groups <- unique(D[[groupCol]])
      for (g in groups) {
        subD <- D[D[[groupCol]] == g, , drop = FALSE]
        grp_res <- run_assoc(subD)
        if (nrow(grp_res) > 0) grp_res[[groupCol]] <- g
        res <- rbind(res, grp_res)
      }
    } else {
      res <- run_assoc(D)
    }

    if (nrow(res) > 0) {
      if (padjByOutcome) {
        group_vars <- c(if (!is.null(groupCol) && groupCol %in% colnames(res)) groupCol else NULL, "outcome")
        res <- res %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
          dplyr::mutate(pAdj = p.adjust(pVal, method = padjMethod)) %>%
          dplyr::ungroup()
      } else {
        if (!is.null(groupCol) && groupCol %in% colnames(res)) {
          res <- res %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(groupCol))) %>%
            dplyr::mutate(pAdj = p.adjust(pVal, method = padjMethod)) %>%
            dplyr::ungroup()
        } else {
          res$pAdj <- p.adjust(res$pVal, method = padjMethod)
        }
      }
    } else {
      res$pAdj <- numeric()
    }

  } else if (method == "correlation") {
    if (!is.null(groupCol) && groupCol %in% colnames(D)) {
      groups <- unique(D[[groupCol]])
      for (g in groups) {
        subD <- D[D[[groupCol]] == g, , drop = FALSE]
        grp_res <- run_corr(subD)
        if (nrow(grp_res) > 0) grp_res[[groupCol]] <- g
        res <- rbind(res, grp_res)
      }
    } else {
      res <- run_corr(D)
    }

    if (nrow(res) > 0) {
      if (!("p" %in% names(res))) res$p <- NA_real_

      if (padjByOutcome) {
        group_vars <- c(if (!is.null(groupCol) && groupCol %in% colnames(res)) groupCol else NULL, "var1")
        res <- res %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
          dplyr::mutate(pAdj = p.adjust(p, method = padjMethod)) %>%
          dplyr::ungroup()
      } else {
        if (!is.null(groupCol) && groupCol %in% colnames(res)) {
          res <- res %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(groupCol))) %>%
            dplyr::mutate(pAdj = p.adjust(p, method = padjMethod)) %>%
            dplyr::ungroup()
        } else {
          res$pAdj <- p.adjust(res$p, method = padjMethod)
        }
      }
    } else {
      res$pAdj <- numeric()
    }
  }

  # -------- map back original column names --------
  if ("outcome" %in% colnames(res)) res$outcome <- nameMapFn(res$outcome)
  if ("exposure" %in% colnames(res)) res$exposure <- nameMapFn(res$exposure)
  if ("var1" %in% colnames(res))    res$var1    <- nameMapFn(res$var1)
  if ("var2" %in% colnames(res))    res$var2    <- nameMapFn(res$var2)

  return(res)
}
