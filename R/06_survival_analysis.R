# =============================================================================
# SURVIVAL ANALYSIS FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   plotKM            - Kaplan-Meier plot with log-rank + Cox statistics
#   getKM_medium      - Extract median survival time per group
#   uniCoxPh          - Univariate Cox proportional hazards
#   multiCoxPh        - Multivariate Cox proportional hazards
#
# NOTE: Additional survival utilities (uni_cox, cal_riskscore, survivalSignatures,
#       survivalSignatures (SVM), predictGroup) are in 10_statistics_utilities.R
#
# SOURCE PROVENANCE:
#   plotKM/getKM_medium: GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R (2025-08-11)
#   uniCoxPh/multiCoxPh: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#
# DUPLICATE NOTE: plotKM also exists in:
#   - GBM/funcsInGBM.R (2024-09-06) — same signature
#   - GBM/XL_figures/results/TCGA_Glass_EIS_IDH1wt_KMcurves_BAM.R — older simpler version
#   RECOMMEND: GBMnat/R/funcsGBMnat.R version (2025-08-11)
# =============================================================================

#' get KM medium survival time
#'
#' @param survdf data frame with survival data.
#' @param duration name of the column with the duration of the event.
#' @param variable name of the column with the variable to stratify the analysis.
#' @param censoring name of the column with the censoring information.
#' @param conversion "d2m" for days to months conversion.
#' @param timeUnit time unit for the duration column.
#'
#' @return a data frame with the medium survival time for each group.
#' @details Source provenance: GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R (2025-08-11).
#' @export getKM_medium
#'
#' @import dplyr survival

# SOURCE: GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R (2025-08-11)
getKM_medium=function(survdf,duration,variable=NULL,censoring, conversion="d2m",timeUnit="month"){
  require(dplyr)
  require(survival)
  survdf=as.data.frame(survdf)

  if(conversion=="d2m"){
    survdf[[duration]]=survdf[[duration]]*12/365.25
    timeUnit="month"
  }

  if(!is.null(variable)){
    survdf[[variable]]=as.factor(survdf[[variable]])
    colnames(survdf)[[which(colnames(survdf)==variable)]]="theVariable"
  }
  colnames(survdf)[[which(colnames(survdf)==duration)]]="theDuration"
  colnames(survdf)[[which(colnames(survdf)==censoring)]]="theStatus"

  if(!is.null(variable)){
    dff=summary(survfit(Surv(theDuration, theStatus) ~ theVariable, data = survdf))$table
    names(dff)=gsub("theVariable=","",names(dff))
    # dff=data.frame(variable=names(dff),medium=dff)
    # colnames(dff)=c(variable,"medium")
  }else{
    dff=summary(survfit(Surv(theDuration, theStatus) ~ 1, data = survdf))$table
  }

  return(dff)
}

#' generate Kaplan-Meier survival curve
#'
#' @param survdf \code{data.frame()}.
#' Should at least include columns of duration (eg. survival), censoring (eg. vital status), grouping varible.
#' @param duration  \code{character()}. The name of the column of survival.
#' @param censoring \code{character()}. The name of the column of vital status.
#' @param variable \code{character()}. The name of the column of grouping variable.
#' @param conversion if "d2m". Convert duration unit from days to months.
#' @param palette \code{character()}. color palette. character vector.
#' @param OrderedVarFactor \code{character()}. order of the factors in grouping variable.
#' @param plotTitle \code{character()}. plot title.
#' @param table_base_size annotation font size.
#'
#' @return a Kaplan-Meier plot.
#' @details Source provenance: GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R (2025-08-11).
#' @importFrom survminer surv_fit
#'
#' @export
#'
# SOURCE: GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R (2025-08-11)
plotKM=function(survdf,duration,censoring,variable=NULL,conversion="d2m",statOnly=F,palette=NULL,OrderedVarFactor=NULL,plotTitle=NULL,table_base_size=6,timeUnit="month",timeYLab="OS"){
  require(survminer)
  require("survival")
  require(gridExtra)
  survdf=as.data.frame(survdf)
  if(!is.null(variable)){
    survdf[[variable]]=as.factor(survdf[[variable]])
    if(is.null(palette)){
      palette=cols4all::c4a("dark24",n = nlevels(survdf[[variable]]))
    }

    if(!is.null(OrderedVarFactor)){
      survdf[[variable]]=factor(survdf[[variable]],levels = OrderedVarFactor)
      # names(palette)=paste(variable,"=",OrderedVarFactor,sep="")
      names(palette)=levels(survdf[[variable]])
    }
  }

  if(conversion=="d2m"){
    survdf[[duration]]=survdf[[duration]]*12/365.25
    timeUnit="month"
  }


  if(is.null(variable)){
    call=as.formula(sprintf("Surv(%s, %s) ~ 1",duration,censoring))
  }else{
    call=as.formula(sprintf("Surv(%s, %s) ~ %s",duration,censoring,variable))
  }
  fit=survminer::surv_fit(call,data = survdf)

  if(!is.null(variable)){
    logRankTest.oa=survdiff(call,data = survdf)$pvalue

    if(nlevels(survdf[[variable]])==2){
      tmp=summary(coxph(call, data = survdf))
      stats=setNames(
        c(logRankTest.oa,tmp$conf.int[,"exp(coef)"],tmp$conf.int[,"lower .95"],tmp$conf.int[,"upper .95"]),
        c("p.value","HR","lowerCI95","upperCI95"))
      stats.str=stats.str_=sprintf("p = %.3f \nHR = %.3f (%.3f - %.3f)",
                                   logRankTest.oa,tmp$conf.int[,"exp(coef)"],tmp$conf.int[,"lower .95"],tmp$conf.int[,"upper .95"])
    }
    if(nlevels(survdf[[variable]])>2){
      logRankTest.pw=pairwise_survdiff(call,data = survdf,p.adjust.method="none")$p.value
      HR=list()
      refU=unique(combn(levels(survdf[[variable]]),2)[1,])
      for (ref in refU){
        dataSurv=survdf[,c(duration,censoring,variable)]
        dataSurv[[variable]]=factor(dataSurv[[variable]],levels = c(ref,levels(dataSurv[[variable]])[levels(dataSurv[[variable]])!=ref]))
        tmp=summary(coxph(call, data =dataSurv))
        tmp=cbind(
          HR=tmp$coefficients[,"exp(coef)"] %>% format(digits=2),
          CI=paste(tmp$conf.int[,"lower .95"]%>% format(digits=2),tmp$conf.int[,"upper .95"]%>% format(digits=2),sep = "-"),
          `p.coxph`=tmp$coefficients[,"Pr(>|z|)"] %>% format(digits=2))
        rownames(tmp)=apply(combn(levels(dataSurv[[variable]]),2)[,1:(nlevels(dataSurv[[variable]])-1),drop=F], 2, function(x) paste(x[2],"vs", x[1],sep = ""))
        HR[[ref]]=tmp
      }
      if(nlevels(dataSurv[[variable]])>2){
        HR=Reduce(rbind,HR)
        HR=HR[!duplicated(lapply(strsplit(rownames(HR),"vs"), sort)),]
      }else{HR=HR[[1]]}
      stats=cbind(HR,`p.logrank`=na.omit(as.vector(logRankTest.pw)) %>% format(digits=2))
      stats=rbind(stats,overall=c(NA,NA,NA,logRankTest.oa%>% format(digits=2)))
      # stats.str=readr::format_delim(stats %>% as.data.frame%>% tibble::rownames_to_column("pair"),delim = "\t")
      stats.str=stats.str_=capture.output(stats %>% as.data.frame%>%knitr::kable("simple",align = 'c'))
      stats.str=paste(stats.str[stats.str!=""&!grepl("-----",stats.str)],collapse = "\n")
    }
    if(statOnly){
      # if(nlevels(survdf[[variable]])==2){
      #   return(logRankTest.oa)
      # }else if(nlevels(dataSurv[[variable]])>2){
      #   return(stats)
      return(stats)

    }else{
      ggsurv=ggsurvplot(fit, data = survdf,
                        palette = palette,
                        # pval = TRUE,
                        xlab = paste("Time",sprintf("(%s)",timeUnit)),
                        ylab = sprintf("%s (probability)",timeYLab),
                        conf.int = F,
                        risk.table = TRUE,
                        risk.table.col = "strata",
                        legend.title = variable,
                        risk.table.y.text = F,
                        legend.labs=levels(survdf[[variable]]),
                        risk.table.height = 0.25,
                        break.x.by = 10,
                        ggtheme = theme_bw())
      if(nlevels(survdf[[variable]]) %in% c(1,2)){
        ggsurv$plot=ggsurv$plot +
          annotate("text",
                   size=(table_base_size+0.4)/.pt,
                   x=range(survdf[[duration]],na.rm = T)[2]*1/4,y=0.2,
                   label=stats.str,alpha = .9,)

      }else if(nlevels(survdf[[variable]])>2){
        ggsurv$plot=ggsurv$plot +
          annotation_custom(gridExtra::tableGrob(stats,theme=ttheme_minimal(base_size = table_base_size)),
                            xmin=range(na.omit(survdf[[duration]]))[2]*2.5/5,
                            xmax=range(na.omit(survdf[[duration]]))[2]*2.5/5,
                            ymin=0.75,ymax=0.90)
      }
      if(!is.null(plotTitle)){
        ggsurv$plot=ggsurv$plot +
          labs(title = sprintf("Kaplan-Meier Curve on %s",duration))
      }

      ggsurv$plot=ggsurv$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position=c(.8,.8))
      ggsurv$table=ggsurv$table + theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())
    }
  }else{
    ggsurv=ggsurvplot(fit, data = survdf,
                      # palette = palette,
                      # pval = TRUE,
                      xlab = paste("Time",sprintf("(%s)",timeUnit)),
                      ylab = sprintf("%s (probability)",timeYLab),
                      conf.int = F,
                      risk.table = TRUE,
                      legend=c(.8,.8),
                      risk.table.col = "strata",
                      # legend.title = variable,
                      risk.table.y.text = F,
                      # legend.labs=levels(survdf[[variable]]),
                      risk.table.height = 0.25,
                      break.x.by = 10,
                      ggtheme = theme_bw())
  }

  return(ggsurv)

}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Univariate Cox proportional hazards
#'
#' Univariate Cox proportional hazards
#' @param time Function argument documented from the legacy interface.
#' @param status Function argument documented from the legacy interface.
#' @param covariateMx Function argument documented from the legacy interface.
#' @return The result object produced by the analysis.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' uniCoxPh(...)
#' }
#' @export
uniCoxPh=function(time,status,covariateMx){
  covariate_names=setNames(rownames(covariateMx),paste("feature",1:nrow(covariateMx),sep = ""))
  rownames(covariateMx)=names(covariate_names)
  temp2=cbind(time=resdf[,resV],status=status,t(covariateMx)) %>% as.data.frame

  univ_formulas <- sapply(names(covariate_names),
                          function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = temp2)})
  univ_results <- lapply(univ_models,
                         function(x){
                           x <- summary(x)
                           n=x$n
                           n.event=x$nevent
                           # p.value<-signif(x$wald["pvalue"], digits=2)
                           p.value.logLik=x$logtest["pvalue"]
                           p.value.sc=x$sctest["pvalue"]
                           p.value.wald<-x$waldtest["pvalue"]
                           # wald.test<-signif(x$wald["test"], digits=2)
                           # wald.test<-x$wald["test"]
                           # beta<-signif(x$coef[1], digits=2);# coeficient beta
                           beta<-x$coef[1]
                           HR <-signif(x$coef[2], digits=2); # exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (",
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, p.value.logLik,p.value.sc,p.value.wald,n,n.event)
                           names(res)<-c("beta", "HR (95% CI for HR)","p.value.logLik","p.value.logRank","p.value.Wald","n","n.event")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE)) %>%as.data.frame()
  res[,c("beta","p.value.logLik","p.value.logRank","p.value.Wald","n","n.event")]=apply(res[,c("beta","p.value.logLik","p.value.logRank","p.value.Wald","n","n.event")], 2, as.numeric)
  rownames(res)=covariate_names
  return(res)
}

# SOURCE: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02)
#' Multivariate Cox proportional hazards
#'
#' Multivariate Cox proportional hazards
#' @param time Function argument documented from the legacy interface.
#' @param status Function argument documented from the legacy interface.
#' @param covariateMx Function argument documented from the legacy interface.
#' @return The result object produced by the analysis.
#' @details Source provenance: woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R (2024-07-02).
#'
#' @examples
#' \dontrun{
#' multiCoxPh(...)
#' }
#' @export
multiCoxPh=function(time,status,covariateMx){
  covariate_names=setNames(rownames(covariateMx),paste("feature",1:nrow(covariateMx),sep = ""))
  rownames(covariateMx)=names(covariate_names)
  temp2=cbind(time=resdf[,resV],status=status,t(covariateMx)) %>% as.data.frame
  temp2%>%analyse_multivariate(
    vars(time,status),
    covariates = eval(parse(text = sprintf("vars(%s)",paste(rownames(covariateMx),collapse=",")))))
}
