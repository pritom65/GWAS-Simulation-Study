source(here::here("scripts","scripts.R"))





model <- ModelDominantCC()
modelComp <- list(ModelGenotypicCC(), ModelCoDominantCC(), ModelAdditiveCC(), ModelRecessiveCC())
parmsGrid <- expand_grid(mar = c(.1,.25,.4), prevalance = seq(.1,.5,.1))

nCases <- 500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 2500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 5000
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))





model <- ModelCoDominantCC()
modelComp <- list(ModelGenotypicCC(), ModelDominantCC(), ModelAdditiveCC(), ModelRecessiveCC())
parmsGrid <- expand_grid(mar = c(.1,.25,.4), prevalance = seq(.1,.5,.1))

nCases <- 500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 2500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 5000
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))







model <- ModelRecessiveCC()
modelComp <- list(ModelGenotypicCC(), ModelCoDominantCC(), ModelAdditiveCC(), ModelDominantCC())
parmsGrid <- expand_grid(mar = c(.1,.25,.4), prevalance = seq(.1,.5,.1)) 

nCases <- 500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 2500
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))

nCases <- 5000
model$setParms(n_case = nCases, n_control = nCases)
powerDf <- runFullExperiment(modelRef = model,modelComp = modelComp,parmsGrid = parmsGrid,leaveEmptyCore = 2)
write_rds(x = powerDf,file = glue::glue("rds_files/model_{class(model)}_sample_{2*nCases}.rds"))







