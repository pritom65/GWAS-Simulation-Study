## Model Genotypic

ModelGenotypicCC <- function(){
  dataObj <- list2env(list(causalOddsRatio = NA, 
                           mar = NA,
                           n_case = NA, 
                           n_control = NA,
                           prevalance = NA,
                           phenotype = NA
  )
  )
  setParameters <- function(causalOddsRatio = NULL, mar = NULL, prevalance = NULL, n_case = NULL, n_control = NULL){
    dataObj$causalOddsRatio <- ifelse(is.null(causalOddsRatio),dataObj$causalOddsRatio,causalOddsRatio)
    dataObj$mar <- ifelse(is.null(mar),dataObj$mar,mar)
    dataObj$prevalance <- ifelse(is.null(prevalance),dataObj$prevalance,prevalance)
    dataObj$n_case <- ifelse(is.null(n_case),dataObj$n_case,n_case)
    dataObj$n_control <- ifelse(is.null(n_control),dataObj$n_control,n_control)
    dataObj$phenotype <- if(is.null(n_case) && is.null(n_control)){
      dataObj$phenotype
    } else {
      c(rep(x = 1, n_case), rep(x = 0, n_control))
    }
  }
  
  structure(.Data = list(data = dataObj,
                         setParms = setParameters), 
            class = "genotypic")
}

generateSNPsCausal.genotypic <- function(model){ }


encodeSNPs.genotypic <- function(model,SNPs){ SNPs }

decodeSNPs.genotypic <- function(model,SNPs){ SNPs }


print.genotypic <- function(model) {
  with(data = as.list(model$data),
       expr = glue::glue("The parameters sets for model {class(model)}\
                          Causal Odds: {causalOddsRatio}
                          Prevalance : {prevalance}
                          Sample Size: {n_case+n_control}
                              Case   : {n_case}
                              Control: {n_control}
                          MAF        : {mar}")) |> 
    print()
  invisible(model)
}

deepCopy.genotypic <- function(model){
  modelCopy <- ModelGenotypicCC()
  formalNames <- formalArgs(modelCopy$setParms)
  container <- vector(mode = "list")
  argument_list <- as.list(model$data)
  
  for(name in formalNames){
    if(!is.na(argument_list[[name]])){
      container[[name]] <- argument_list[[name]]
    }
  }
  
  do.call(what = modelCopy$setParms, args = container)
  return(modelCopy)
}