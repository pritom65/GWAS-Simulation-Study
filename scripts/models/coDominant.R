## Co-Dominant

ModelCoDominantCC <- function(){
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
            class = "co_dominant")
}

generateSNPsCausal.co_dominant <- function(model, seed = 1){
  set.seed(seed)
  p <- model$data$mar      # P(B)
  p <- 2*p*(1-p)           # P(AB)
  
  rho <- model$data$causalOddsRatio # Odds Ratio
  pi <- model$data$prevalance       # P(Case)
  
  if(rho == 1){
    return(rbinom(n = with(data = model$data,expr = n_case + n_control), size = 1, prob = pi))
  }
  
  a <- (1 - pi) * (rho - 1)
  b <- pi * rho + 1 - pi - p * (rho - 1)
  c <- -p
  s <- c(-b + sqrt(b^2-4*a*c), -b - sqrt(b^2-4*a*c))/(2*a)
  probOf_x_eq_1GivenControl <- s <- ifelse(0< s[1] & s[1] < 1,s[1],s[2])
  probOf_x_eq_1GivenCase <- r <- (rho*s)/(1+s*(rho-1))
  
  c(rbinom(n = model$data$n_case,    size = 1, prob = probOf_x_eq_1GivenCase), 
    rbinom(n = model$data$n_control, size = 1, prob = probOf_x_eq_1GivenControl))
}


encodeSNPs.co_dominant <- function(model,SNPs){
  p <- model$data$mar
  q <- 1-p
  index <- SNPs == 0
  n_OfZero <- sum(index)
  
  SNPs[index] <- sample(x = c(0,2),size = n_OfZero,replace = TRUE,prob = c(q^2,p^2)/(p^2+q^2) )
  SNPs
}

decodeSNPs.co_dominant <- function(model,SNPs){
  SNPs[SNPs == 2] <- 0
  SNPs
}

print.co_dominant <- function(model) {
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



deepCopy.co_dominant <- function(model){
  modelCopy <- ModelCoDominantCC()
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