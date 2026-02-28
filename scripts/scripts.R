library(tidyverse)
library(future)
library(furrr)
library(progressr)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = .5))

generateSNPsCausal <- function(model,...) UseMethod("generateSNPsCausal")
getCausalPvalue <- function(model,...) UseMethod("getCausalPvalue")
encodeSNPs <- function(model,...) UseMethod("encodeSNPs")
decodeSNPs <- function(model,...) UseMethod("decodeSNPs")
deepCopy <- function(model,...) UseMethod("deepCopy")

source(here::here("scripts","models","additive.R"))
source(here::here("scripts","models","coDominant.R"))
source(here::here("scripts","models","dominant.R"))
source(here::here("scripts","models","recissive.R"))
source(here::here("scripts","models","genotypicCC.R"))

# The main function factory that create an object framework to execute "ONE" experiment for generating 
# causal SNPs that passed into the argument modelRef and convert the data to the models that are in the 
# list object passing with the argument ModelComp. The output is just a vector that provide whether for 
# this experiment the models reject or accept those null hypothesis.

fn_NoOfDinamicIter_M <- function(theta, moe = 0.05,alpha = .1, max_iter = 500){
  fn_optimize <- function(n){
    a <- qbeta(p = alpha/2,    shape1 = theta*(n-1), shape2 = (1-theta)*(n-1))
    b <- qbeta(p = 1- alpha/2, shape1 = theta*(n-1), shape2 = (1-theta)*(n-1))
    b-a
  }
  
  suppressWarnings(
    M <- uniroot(f = \(n) fn_optimize(n) - moe, lower = 1,upper = 10e3)$root |> 
      ceiling()
    )
  min(M, max_iter)
}


objGetOR <- function(orMin = 1, orMax = 51, orStep = 12, maxRefOR = NULL, derivativeFilter = .01, noOfOR_evaluated = 45){
  # browser()
  orFilter <- function(x, orSep) {
    if(orSep == 1)
      return(x)
    
    repeat{
      if (orSep > 1)
        break
      orSep <- orSep * 2
      x <- x[x %% orSep != 0]
    }
    x
  }
  
  dfOR <- 
    tibble(orSep = 1/2^(0:orStep)) |>
    transmute(or = map(.x = orSep,.f = \(sep) seq(0,1,sep)),
              or = map2(.x = or, .y = orSep, .f = \(x,orSep) orFilter(x,orSep)),
              power = NA_real_
    ) |> 
    unnest(cols = or) |> 
    mutate(or = or*(orMax-orMin) + orMin) 
  
  data <- list2env(list(orMin = orMin, 
                        orMax = orMax, 
                        orStep = orStep, 
                        dfOR = dfOR, 
                        noOfOR_evaluated = noOfOR_evaluated, 
                        zero = 1e-3,
                        one = 1 - 1e-2,
                        derivativeFilter = derivativeFilter,
                        plainForwardOneStep = TRUE)
  )
  
  fn_getOR <- function(){
    with(data = data$dfOR, or[which.max(is.na(power))])
  }
  
  
  fn_SaveOR <- function(OR,power){
    # browser()
    # data$dfOR |> na.omit() |> as.matrix() |> print()
    
    if(data$zero < power & power < data$one){
      data$dfOR[data$dfOR$or == OR,"power"] <- power
    }else if(data$one < power){
      data$orMax <- min(OR,data$orMax)
      data$dfOR <- data$dfOR[data$dfOR$or < data$orMax,]
    }else if(power < data$zero){
      data$orMin <- max(OR,data$orMin)
      data$dfOR <- data$dfOR[data$dfOR$or > data$orMin,]
    }
    
    
    if(data$plainForwardOneStep & 10 <= nrow(na.omit(data$dfOR))){
      data$plainForwardOneStep <- FALSE
      
      data$dfOR <- na.omit(data$dfOR)
      
      # fn_exp <- function(x,a,b) pnorm(q = x,mean = a,sd = b)
      fn_exp <- function(x,a,b) data$one * pbeta(q = (x-data$orMin)/(data$orMax - data$orMin),shape1 = a,shape2 = b)
      parms <- with(data = data$dfOR, 
                    expr = optim(par = c(mean(or),2),
                                 fn = \(parms) abs(fn_exp(or,parms[1],parms[2]) - power) |> sum() )$par 
      )
      
      midOR <- uniroot(f = \(x) fn_exp(x,parms[1],parms[2]) - 0.5,interval = c(data$orMin,data$orMax))$root
      midOR <- (midOR-data$orMin)/(data$orMax-data$orMin)
      
      if(!is.null(maxRefOR)){
        if(10*maxRefOR < data$orMax){
          fn_exp <- function(x) data$one * pbeta(q = (x-data$orMin)/(data$orMax - data$orMin),shape1 = parms[1],shape2 = parms[2])
          fn_der <- function(f,h = .0001){
            function(x){
              (f(x+h) - f(x-h))/(2*h)
            }
          }
          
          fn_expPrime1 <- fn_der(fn_exp)
          fn_expPrime2 <- fn_der(fn_expPrime1)
          
          # browser()
          orForDer <- seq(data$orMin,data$orMax,length.out = 1000)
          
          prime1 <- sapply(X = orForDer, FUN = \(x) fn_expPrime1(x) )
          prime2 <- sapply(X = orForDer, FUN = \(x) fn_expPrime2(x) )
          
          prime1 <- prime1[prime2 < 0]
          orForDer <- orForDer[prime2 < 0]
          orMax <- orForDer[order(abs(prime1 - data$derivativeFilter)) |> head()] |> max()
          data$orMax <- orMax
        }
      }
      
      dfOR <- 
        tibble(orSep = 1/2^(0:data$orStep)) |> 
        transmute(or = map(.x = orSep,.f = \(sep) seq(0,1,sep)),
                  or = map2(.x = or, .y = orSep, .f = \(x,orSep) orFilter(x,orSep)),
                  or = map(.x = or,.f = \(x) x[order(abs(x-midOR))] ),
                  power = NA_real_
        ) |> 
        unnest(cols = or) |> 
        mutate(or = or*(data$orMax-data$orMin) + data$orMin) |>
        filter(data$orMin < or, or < data$orMax) |> 
        select(or,power)
      
      dfOR <- dfOR[!dfOR$or %in% data$dfOR$or,]
      # data$dfOR <- bind_rows(data$dfOR, dfOR)
      data$dfOR <- bind_rows(
        data$dfOR |> 
          filter(data$orMin < or, or < data$orMax),
        dfOR
      )
    }
  }
  
  fn_checkOR_LE_limit <- function(){
    # browser()
    noOfOR_evaluated <- 
      data$dfOR |> 
      na.omit() |> 
      nrow()
    
    if(nrow(data$dfOR) == 0 | noOfOR_evaluated  >= data$noOfOR_evaluated){
      data$dfOR <- 
        na.omit(data$dfOR) |> 
        bind_rows(
          tibble(or = with(data = data, expr = c(orMin, orMax)),
                 power = c(data$zero,data$one))
        )
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  
  return(list(data = data, fn_getOR = fn_getOR, fn_SaveOR = fn_SaveOR, fn_checkOR_LE_limit = fn_checkOR_LE_limit))
}


objGetORNonConverge <- function(orMin = 1, orMax = 51, orStep = 12, maxPower, maxRefOR = NULL, derivativeFilter = .005, noOfOR_evaluated = 45){
  # browser()
  orFilter <- function(x, orSep) {
    if(orSep == 1)
      return(x)
    
    repeat{
      if (orSep > 1)
        break
      orSep <- orSep * 2
      x <- x[x %% orSep != 0]
    }
    x
  }
  
  dfOR <- 
    tibble(orSep = 1/2^(0:orStep)) |>
    transmute(or = map(.x = orSep,.f = \(sep) seq(0,1,sep)),
              or = map2(.x = or, .y = orSep, .f = \(x,orSep) orFilter(x,orSep)),
              power = NA_real_
    ) |> 
    unnest(cols = or) |> 
    mutate(or = or*(orMax-orMin) + orMin) |> 
    filter(or < orMax) |> 
    bind_rows(tibble(or = orMax, power = maxPower))
  
  data <- list2env(list(orMin = orMin, 
                        orMax = orMax, 
                        orStep = orStep, 
                        dfOR = dfOR, 
                        noOfOR_evaluated = noOfOR_evaluated, 
                        zero = 1e-5,
                        maxPower = maxPower,
                        derivativeFilter = derivativeFilter,
                        plainForwardOneStep = TRUE)
  )
  
  fn_getOR <- function(){
    with(data = data$dfOR, or[which.max(is.na(power))])
  }
  
  
  fn_SaveOR <- function(OR,power){
    # browser()
    # data$dfOR |> na.omit() |> as.matrix() |> print()
    
    if(data$zero < power){
      data$dfOR[data$dfOR$or == OR,"power"] <- power
    }else if(power < data$zero){
      data$ORatPowerZero <- max(OR,data$ORatPowerZero)
      data$orMin <- max(OR,data$orMin)
      data$dfOR <- data$dfOR[data$dfOR$or > data$orMin,]
    }
    
    
    if(data$plainForwardOneStep & 10 <= nrow(na.omit(data$dfOR))){
      data$plainForwardOneStep <- FALSE
      # browser()
      
      data$dfOR <- na.omit(data$dfOR)
      
      model <- data$dfOR |> 
        bind_rows(tibble(or = data$orMin, power = data$zero)) |> 
        mutate(power = power/max(power)) |> 
        lm(formula = power ~ or + I(or^2))
      fn_exp <- function(x) {
        predict(object = model,newdata = tibble(or = x))
      }
      
      
      if(!is.null(maxRefOR)){
        if(10*maxRefOR < data$orMax){

          fn_der <- function(f,h = .0001){
            function(x){
              (f(x+h) - f(x-h))/(2*h)
            }
          }
          
          fn_expPrime1 <- fn_der(fn_exp)
          
          orForDer <- seq(data$orMin,data$orMax,length.out = 1000)
          
          prime1 <- sapply(X = orForDer, FUN = \(x) fn_expPrime1(x) )
          
          
          orMax <- orForDer[order(abs(prime1 - data$derivativeFilter)) |> head()] |> max()
          data$orMax <- max(orMax,min(data$dfOR$or))
          data$maxPower <- 
            data$dfOR |> 
            na.omit() |> 
            mutate(minAbsDev = abs(or - orMax)) |> 
            slice_min(minAbsDev) |> 
            pull(power)
        
          data$dfOR <- 
          data$dfOR |> 
          filter(or <= data$orMax)
        }
      }
      
      # midOR <- mean(c(data$orMin, data$orMax))
      # midOR <- (midOR - data$orMin) / (data$orMax - data$orMin)
      midOR <- .5
      # browser()
      # model <- lm(formula = power ~ or + I(or^2),data = data$dfOR) 
      # midOR <- predict(object = model,newdata = tibble(or = (data$orMin + data$orMax)/2 ))
      # midOR <- (midOR-data$orMin)/(data$orMax-data$orMin)
      
      dfOR <- 
        tibble(orSep = 1/2^(0:data$orStep)) |> 
        transmute(or = map(.x = orSep,.f = \(sep) seq(0,1,sep)),
                  or = map2(.x = or, .y = orSep, .f = \(x,orSep) orFilter(x,orSep)),
                  or = map(.x = or,.f = \(x) x[order(abs(x-midOR))] ),
                  power = NA_real_
        ) |> 
        unnest(cols = or) |> 
        mutate(or = or*(data$orMax-data$orMin) + data$orMin) |>
        filter(data$orMin < or, or < data$orMax) |> 
        select(or,power)
      
      dfOR <- dfOR[!dfOR$or %in% data$dfOR$or,]
      data$dfOR <- bind_rows(
        data$dfOR |> 
          filter(data$orMin < or),
        dfOR
        )
    }
  }
  
  fn_checkOR_LE_limit <- function(){
    # browser()
    noOfOR_evaluated <- 
      data$dfOR |> 
      na.omit() |> 
      nrow()
    
    # print(noOfOR_evaluated)
    
    if(data$maxPower < .1 | nrow(data$dfOR) == 0 | noOfOR_evaluated  >= data$noOfOR_evaluated){
      data$dfOR <- 
        na.omit(data$dfOR) |> 
        bind_rows(
          tibble(or = with(data = data, expr = orMin),
                 power = data$zero)
        )
      return(FALSE)
    }else {
      return(TRUE)
    }
  }
  
  return(list(data = data, fn_getOR = fn_getOR, fn_SaveOR = fn_SaveOR, fn_checkOR_LE_limit = fn_checkOR_LE_limit))
}

runSingleExperiment <- function(modelRef,modelComp,seed = 1234,sigLavel = 5e-8){
  objSeedSetter <- setSeedValue(seed)
  phenotype <- modelRef$data$phenotype
  
  obtainOneTest <- function(or){
    causalSNPs <- generateSNPsCausal(modelRef, seed = objSeedSetter(or))
    if(class(modelRef) != class(modelComp)){
      encodedSNPsAdditive <- encodeSNPs(model = modelRef, SNPs = causalSNPs)
      decodedSNPs <- decodeSNPs(model = modelComp, SNPs = encodedSNPsAdditive)
    } else {
      decodedSNPs <- causalSNPs
    }
    
    if(class(modelComp) == "additive"){
      # browser()
      additiveModelLogisticRegSummary <- glm(formula = y ~ x,
                                             family = "binomial",
                                             data = data.frame(x = decodedSNPs, y = phenotype)
      ) |> 
        summary()
      return(additiveModelLogisticRegSummary$coefficients[2,4] < sigLavel)
    }
    
    p_valueOfGeneticModel <- chisq.test(table(phenotype,decodedSNPs))$p.value
    p_valueOfGeneticModel < sigLavel
  }
  
  function(or){
    suppressWarnings({
      
      modelRef$setParms(causalOddsRatio = or)
      initialTestOfHypothesis <- replicate(n = 500,expr = obtainOneTest(or))
      theta <- mean(initialTestOfHypothesis)
      M <- fn_NoOfDinamicIter_M(theta,moe = .05,alpha = .1,max_iter = 500)
      finalTestOfHypothesis <- replicate(n = M, expr = obtainOneTest(or))
      mean(c(initialTestOfHypothesis, finalTestOfHypothesis))
      
    })
  }
}


runFullSingleExperiment <- function(mar, prev, modelRef, modelComp = NULL, maxOR = 50,maxRefOR = NULL, noOfOR_evaluated = 45, seed = 1234, sigLevel = 5e-8){
  # browser()
  modelRef <- deepCopy(modelRef)
  modelRef$setParms(mar = mar, prevalance = prev)
  
  if(!is.null(modelComp)){
    modelComp <- deepCopy(modelComp)
    expObject <- runSingleExperiment(modelRef = modelRef,modelComp = modelComp, sigLavel = sigLevel, seed = seed)
  } else {
    expObject <- runSingleExperimentRefClass(modelRef = modelRef, sigLavel = sigLevel, seed = seed)
  }
  
  power <- expObject(maxOR)
  # browser()
  if(is.null(modelComp) & power < .5){
    return( list(powerDf = tibble(or = c(1,maxOR),power = 0:1), orMax = NA_real_) ) 
  } 
  
  if(.99 <= power){
    if(!is.null(modelComp)){
      objOR <- objGetOR(orMax = maxOR, maxRefOR = maxRefOR, noOfOR_evaluated = noOfOR_evaluated)
    } else {
      objOR <- objGetOR(orMax = maxOR, noOfOR_evaluated = noOfOR_evaluated)
    }
    
  } else {
    if(!is.null(modelComp)){
      objOR <- objGetORNonConverge(orMax = maxOR, maxPower = power, maxRefOR = maxRefOR, noOfOR_evaluated = noOfOR_evaluated)
    } else {
      objOR <- objGetORNonConverge(orMax = maxOR, maxPower = power, noOfOR_evaluated = noOfOR_evaluated)
    }
  } 
  
  
  while(objOR$fn_checkOR_LE_limit()){
    or <- objOR$fn_getOR()
    if(is.null(modelComp)){
      expObject <- runSingleExperimentRefClass(modelRef = modelRef, sigLavel = sigLevel, seed = seed)
    } else {
      expObject <- runSingleExperiment(modelRef = modelRef, modelComp = modelComp, sigLavel = sigLevel, seed = seed)
    }
    
    power <- expObject(or)
    objOR$fn_SaveOR(or, power)
  }
  
  if(is.null(modelComp)){
    cat(sprintf("MAR: %.2f Prev: %.2f Model: %s \n", mar, prev, str_to_title(class(modelRef))) )
  } else {
    cat(sprintf("MAR: %.2f Prev: %.2f Model: %s \n", mar, prev, str_to_title(class(modelComp))) )
  }
  
  if(is.null(modelComp)){
    return( with(data = objOR$data, expr = list(powerDf = dfOR, orMax = getCeilingForExperiment(dfOR), maxRefOR = orMax)) )
  }
  
  objOR$data$dfOR
}


runSingleExperimentRefClass <- function(modelRef,seed=1234,sigLavel = 5e-8){
  modelComp <- deepCopy(modelRef)
  runSingleExperiment(modelRef = modelRef,modelComp = modelComp,seed = seed,sigLavel = sigLavel)
}

getCeilingForExperiment <- function(df){
  # browser()
  u <- NULL
  orMin <- min(df$or)
  
  fn_exp <- function(x,a,b) max(df$power) * pbeta(q = (x-min(df$or))/(3*max(df$or) + 20 - min(df$or)),shape1 = a,shape2 = b)
  suppressWarnings(
    par <- with(data = df, 
                expr = optim(par = c(1,1),
                             fn = \(parms) abs(fn_exp(or,parms[1],parms[2]) - power) |> sum()
                )$par
    )
  )
  
  fn_calArea <- function(parms){
    a <- parms[1]
    b <- parms[2]
    # browser()
    
    refCurve <- \(x) fn_exp(x = x,a = par[1],b = par[2])
    compCurve <- \(x) fn_exp(x = x,a = a,b = b)
    
    l <- optim(par = 0,fn = \(x) abs(refCurve(x) - 1e-4) )$par
    u <- optim(par = 10,fn = \(x) abs(compCurve(x) - .975))$par
    
    x <- seq(l,u,.01)
    # print(c(a,b,l,u))
    # (tibble(x,y1 = refCurve(x),y2 = compCurve(x)) |>
    #   ggplot(aes(x)) +
    #   geom_line(aes(y=y1),col = "red") +
    #   geom_line(aes(y=y2),col = "darkgreen") +
    #   geom_vline(xintercept = orMin)) |>
    #   print()
    u <<- u
    temp <- mean(abs(refCurve(x) - compCurve(x)))
    temp
  }
  
  optim(par = par + .1, fn = \(parms) (fn_calArea(parms) - .6)^2) |> 
    suppressWarnings()
  
  u |> ceiling()
} 

fn_NoOfDinamicIter_M <- function(theta, moe = 0.05,alpha = .1, max_iter = 500){
  fn_optimize <- function(n){
    a <- qbeta(p = alpha/2,    shape1 = theta*(n-1), shape2 = (1-theta)*(n-1))
    b <- qbeta(p = 1- alpha/2, shape1 = theta*(n-1), shape2 = (1-theta)*(n-1))
    b-a
  }
  
  suppressWarnings(
    M <- uniroot(f = \(n) fn_optimize(n) - moe, lower = 1,upper = 10e3)$root |> 
      ceiling()
  )
  min(M, max_iter)
}

setSeedValue <- function(seed){
  counter <- 0
  
  function(or){
    # browser()
    counter <<- counter + 1
    x <- glue::glue("{sprintf('%.3f', or) |> str_remove(pattern = '\\\\.')}{'001'}") |> 
      as.numeric() + 
      counter +
      seed
    x
  }
}




runFullExperiment <- function(modelRef,modelComp,parmsGrid,maxOR = 100,noOfOR_evaluated = 45, seed = 1234, sigLevel = 5e-8,leaveEmptyCore = 2){
  objs <- ls(envir = .GlobalEnv)
  funs <- objs[sapply(objs, function(x) is.function(get(x, envir = .GlobalEnv)))]
  timeStarted <- Sys.time()
  
  
  if(availableCores() > leaveEmptyCore){
    suppressMessages(plan(multisession, workers = availableCores() - leaveEmptyCore))
  } else {
    plan(sequential)
  }
  
  # handlers("cli")
  # options(cli.progress_bar_width = 40)
  handlers("txtprogressbar")
  options(progressr.txtprogressbar.width = 30)
  glue::glue("The reference Model : {stringr::str_to_title(class(modelRef))}")
  with_progress({
    
    pb <- progressor(steps = nrow(parmsGrid))
    # browser()
    refModelOrDf <- 
      parmsGrid |> 
      mutate(modelRef = list(modelRef),
             result = future_pmap(.l = list(mar,prevalance,modelRef),
                                  .f = \(m,p,modelRef) {
                                    # pb(sprintf(" MAR: %.2f Prev: %.2f Model: %s\n", mar, prevalance, stringr::str_to_title(class(modelRef))) )
                                    pb()
                                    suppressWarnings(suppressPackageStartupMessages({
                                      library("tidyverse")
                                      runFullSingleExperiment(mar = m,
                                                              prev = p,
                                                              modelRef = modelRef,
                                                              maxOR = maxOR,
                                                              noOfOR_evaluated = noOfOR_evaluated,
                                                              seed = seed,
                                                              sigLevel = sigLevel)
                                    })
                                    )
                                  },
                                  .options = furrr_options(globals = funs)
             ),
             maxOR = map_dbl(result,.f = \(l) l$orMax),
             maxRefOR = map_dbl(.x = result,
                                .f = \(l) ifelse("maxRefOR" %in% names(l),l$maxRefOR,NA)
                                ),
             powerDf = map(result,.f = \(l) l$powerDf),
      ) |> 
      select(mar,prevalance,powerDf,maxOR,maxRefOR) |> 
      filter(!is.na(maxOR))
  })
  
  glue::glue("The comparison Model: {paste(stringr::str_to_title(sapply(X = modelComp,FUN = class)), collapse = ', ')}")
  # handlers("cli")
  # options(cli.progress_bar_width = 40)
  with_progress({
    
    pb <- progressor(steps = nrow(parmsGrid)*length(modelComp))
    compModelDf <- 
      parmsGrid |> 
      right_join(
        refModelOrDf |> 
          select(mar,prevalance,maxOR,maxRefOR),
        by = join_by(mar, prevalance)
      ) |> 
      mutate(modelRef = list(modelRef),
             modelComp = list(modelComp)) |> 
      unnest(modelComp) |> 
      mutate(powerDf = future_pmap(.l = list(mar,prevalance,modelRef,modelComp,maxOR,maxRefOR),
                                   .f = \(m,p,modelRef,modelComp,maxOR,maxRefOR){
                                     # pb(sprintf(" MAR: %.2f Prev: %.2f Model: %s", mar,
                                     #            prevalance, stringr::str_to_title(class(modelComp))) )
                                     pb()
                                     suppressWarnings(suppressPackageStartupMessages({
                                       library("tidyverse")
                                       runFullSingleExperiment(mar = m,
                                                               prev = p,
                                                               modelRef = modelRef,
                                                               modelComp = modelComp,
                                                               maxOR = maxOR,
                                                               maxRefOR = maxRefOR,
                                                               noOfOR_evaluated = noOfOR_evaluated,
                                                               seed = seed,
                                                               sigLevel = sigLevel)
                                     })
                                     )
                                   },
                                   .options = furrr_options(globals = funs)
      )
      )
    
  })
  
  cat(sprintf("Time Required: %.2f Mins",Sys.time() - timeStarted))
  plan(sequential)
  
  compModelDf |>
    mutate(modelRef = map_chr(.x = modelRef,.f = class),
           modelComp = map_chr(.x = modelComp,.f = class),
           evaluatedDataPoints = map_dbl(.x = powerDf, .f = nrow)
    ) |> 
    filter(evaluatedDataPoints != 2) |> 
    bind_rows(
      refModelOrDf |> 
        mutate(modelComp = slice_head(compModelDf)$modelRef[[1]] |> class())
    ) |> 
    select(-c(modelRef,maxOR,evaluatedDataPoints), model = modelComp, prev = prevalance) |> 
    unnest(powerDf)
}



getPolygonDf <- function(powerDf, baseLineModel, sampleSize, verbose = TRUE){
  # browser()
  fn_getPowerQuantile <- function(model, df, q = 0.025){
    # browser()
    uniroot(f = \(or) predict(object = model,newdata = tibble(or = or)) - q, interval = c(1, max(df$or) + 20 ) )$root
  } 
  
  fn_getPowerQuantileU <- function(model, df, q = 0.975){
    # browser()
    if(max(df$power) < .975){
      max(df$or)
    }else {
      uniroot(f = \(or) predict(object = model,newdata = tibble(or = or)) - q, interval = c(1, max(df$or) + 20 ) )$root
    }
  } 
  
  suppressWarnings(suppressMessages({
  
    temp_dfForPolygon <- powerDf |> 
      nest(data = -c(mar,prev,model)) |> 
      arrange(prev,mar) |> 
      mutate(facet = sprintf("(m: %.2f)~vs~(pi: %.2f)", mar, prev),
             facet = as_factor(facet),
      ) |> 
      group_by(facet) |> 
      mutate(mps_model = map(.x = data,
                             .f = \(df) scam::scam(formula = power ~ s(or, bs="mpi"), data = df) ),
             or_025 = map2_dbl(.x = mps_model,.y = data,.f = ~ fn_getPowerQuantile(.x, .y, 0.025) ),
             or_975 = map2_dbl(.x = mps_model,.y = data,.f = ~ fn_getPowerQuantileU(.x, .y, 0.975) ),
             or_min = min(or_025),
             or_max = if_else(model == baseLineModel, NA, or_975),
      ) |> 
      mutate(model = str_replace(model,"_"," "),
             model = str_to_title(model),
             model = str_replace(model," ","-")
      ) |> 
      select(facet,mar,prev,model,mps_model,or_min,or_max) |> 
      mutate(base_model = mps_model[is.na(or_max)]) |> 
      na.omit() |> 
      transmute(mar, prev, model,
                polygon = pmap(.l = list(mps_model,base_model,or_min,or_max),
                               .f = \(comp_model,base_model,or_min,or_max){
                                 or_forward <- seq(or_min,or_max,length.out = 100)
                                 or_backward <- rev(or_forward)
                                 power_forward <- predict(object = comp_model,newdata = tibble(or = c(or_forward)))
                                 power_backward <- predict(object = base_model,newdata = tibble(or = c(or_backward)))
                                 bind_rows(
                                   tibble(or = or_forward,power = power_forward),
                                   tibble(or = or_backward,power = power_backward)
                                 ) |> 
                                   mutate(power = case_when(power < 0 ~ 0,
                                                            1 < power ~ 1,
                                                            TRUE ~ power))
                               }),
                power = map_dbl(.x = polygon,
                                .f = \(df){
                                  df |> 
                                    mutate(type = rep(c("forward","backward"),each = n()/2)) |> 
                                    pivot_wider(names_from = type, values_from = power) |> 
                                    summarise(LoPT = mean(abs(backward - forward)) * diff(range(or)) ) |> 
                                    pull()
                                }),
                power_std = map_dbl(.x = polygon,
                                    .f = \(df){
                                      df |> 
                                        mutate(type = rep(c("forward","backward"),each = n()/2)) |> 
                                        pivot_wider(names_from = type, values_from = power) |> 
                                        summarise(LoPT = mean(abs(backward - forward))) |> 
                                        mutate(LoPT = LoPT * .95) |> 
                                        pull()
                                    }),
                power_diffAtOne = map(.x = polygon,
                                      .f = \(df){
                                        df |> 
                                          slice_max(or) |> 
                                          pull(power) |> 
                                          diff()
                                      })
      )
    
    if(verbose){
      temp_dfForPolygon |> 
        ungroup() |> 
        arrange(mar,prev,power) |> 
        select(model,mar,prev,power,power_std) |> 
        print()
    }
    
    baseLineModel <- 
      tibble(model = baseLineModel) |> 
      mutate(model = str_replace(model,"_"," "),
             model = str_to_title(model),
             model = str_replace(model," ","-")
      ) |> 
      pull()
    
  }))
  
  invisible(return(list(baseLineModel = baseLineModel, dfForPolygon = temp_dfForPolygon, sampleSize = sampleSize)))
}

drawPolygon <- function(temp_dfForPolygon, filterMaxPowerDiffAtOne = .2){
  suppressWarnings(suppressMessages({
    baseLineModel <- temp_dfForPolygon$baseLineModel
    sampleSize <- temp_dfForPolygon$sampleSize
    temp_dfForPolygon <- temp_dfForPolygon$dfForPolygon |> 
      filter(power_diffAtOne <= filterMaxPowerDiffAtOne) |> 
      select(-power_diffAtOne)
    
    ggplot() +
      geom_polygon(data = temp_dfForPolygon |>  
                     # arrange(facet,-power_std) |> 
                     arrange(facet, -power/ min(power)) |> 
                     mutate(group = sprintf("%s %s",facet,model),
                            group = as_factor(group),
                            model = as.factor(model),
                            inc = (row_number()-1) * .25) |> 
                     unnest() |> 
                     group_by(model,facet) |> 
                     mutate(or = (or - min(or))/(max(or)-min(or)),
                            # or = or + inc,
                            power = power - inc
                     ),
                   mapping = aes(or, power, group = group, fill = model), 
                   col = "gray25",
                   alpha = .95,
                   # show.legend = F
      ) +
      geom_text(data = temp_dfForPolygon |> 
                  # arrange(power_std) |>
                  arrange(-power/ min(power)) |>
                  mutate(inc = (row_number()-1) * .25,
                         x = 1.05,
                         y = - .5 + inc * .5,
                         power_std = round(power_std,4)
                  ) |> 
                  select(model,power_std,x,y),
                mapping = aes(x,y,label = power_std),
                size = 3,hjust = 0,
                show.legend = FALSE) +
      geom_point(data = temp_dfForPolygon |> 
                   # arrange(power_std) |>
                   arrange(-power/ min(power)) |>
                   mutate(inc = (row_number()-1) * .25,
                          x = 1,
                          y = -.5 + inc * .5,
                          power_std = round(power_std,4)
                   ) |> 
                   select(model,power_std,x,y),
                 mapping = aes(x,y,fill = model),
                 shape = 21,
                 # size = 2.5,
                 show.legend = FALSE) +
      facet_wrap(~facet, scales = "free_x",ncol = length(unique(temp_dfForPolygon$mar)),labeller = label_parsed) +
      scale_color_viridis_d() +
      scale_fill_viridis_d() +
      lims(x = c(0,1.3)) +
      scale_y_continuous(breaks = NULL) +
      labs(title = glue::glue("Simulation results under a {baseLineModel} underlying genetic model for {sampleSize} cases and {sampleSize} control."),
           fill = "Model:",
           x = expression("odds ratio"~(rho))) +
      theme(legend.position = "bottom")
  
  }))
}






