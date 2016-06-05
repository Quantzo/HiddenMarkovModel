#############Instaling and testing libraries#############
testPackage <- function(name) {
    if (!(name %in% rownames(installed.packages()))) {
        print(paste("Installing ", name, " package"))
        install.packages(c(name))
    }
}

testPackage("HMM")
testPackage("jsonlite")
testPackage("foreach")
testPackage("doParallel")

library(jsonlite)
library(HMM)
library(foreach)
library(doParallel)

#############Cluster initialization#############
initializeCluster <- function(numberOfCores) {
    cl <- makeCluster(numberOfCores)
    clusterSetRNGStream(cl, 9956)
    registerDoParallel(cl)
    return(cl)
}


#############Definitions#############
#Definition of symbols used in HMM
AminoAcidSymbols <- c("G", "P", "D", "E", "K", "R", "H", "S", "T", "N", "Q", "A", "M", "Y", "W", "V", "I", "L", "F", "C")


#############Data loading#############
#Loading Data to data frame
dataFrame <- fromJSON("Data/jsonProteins.json", flatten = TRUE)

#############Preparing training and test data sets#############
testRows <- sample(1:nrow(dataFrame), nrow(dataFrame)/3);

trainingSet <- dataFrame[testRows,]
testSet <- dataFrame[ - testRows,]

#############Preparing observation vector for learning#############
observationVecTraining <- unlist(strsplit(paste(trainingSet[, c("firstStructure")]), ''))

#############Baum Welch training#############
train <- function(model, observation, iter) {
    return(tmp <- baumWelch(hmm = model, observation = observation, maxIterations = iter))
}
#############Predicting results form model#############
#Funtion that cleans results vector from digits
cleanResultsVector <- function(resultsVector) {
    return(gsub("[[:digit:]]", "", resultsVector))
}

#Fuction that returns clean results vector
predictHmm <- function(model, observation) {
    return(cleanResultsVector(viterbi(model, observation)))
}

#############Evaluation functions#############
correctPositions <- function(predictedVector, labeledVector) {
    return(predictedVector[labeledVector == predictedVector])
}

numberOfStates <- function(vector, stateSymbol) {
    return(length(which(vector == stateSymbol)))
}



evaluateEntry <- function(dataSetEntry, model) {
    firstStructureVector <- unlist(strsplit(dataSetEntry$firstStructure, ''))
    secondaryStructureVector <- unlist(strsplit(dataSetEntry$secondaryStructure, ''))

    totalNumber <- length(firstStructureVector)
    totalHelixesNumber <- numberOfStates(secondaryStructureVector, "H")
    totalStrandsNumber <- numberOfStates(secondaryStructureVector, "B")
    totalCoilsNumber <- numberOfStates(secondaryStructureVector, "C")

    predictedVector <- predictHmm(model, firstStructureVector)

    correctVector <- correctPositions(predictedVector, secondaryStructureVector)

    correctNumber <- length(correctVector)

    correctHelixesNumber <- numberOfStates(correctVector, "H")
    correctStrandsNumber <- numberOfStates(correctVector, "B")
    correctCoilsNumber <- numberOfStates(correctVector, "C")

    result <- list()

    result$TotalNumber <- totalNumber
    result$CorrectNumber <- correctNumber

    result$TotalHelixesNumber <- totalHelixesNumber
    result$CorrectHelixesNumber <- correctHelixesNumber

    result$TotalStrandsNumber <- totalStrandsNumber
    result$CorrectStrandsNumber <- correctStrandsNumber

    result$TotalCoilsNumber <- totalCoilsNumber
    result$CorrectCoilsNumber <- correctCoilsNumber


    return(result)
}
combineEntriesResult <- function(ob1, ob2) {
    result <- list()

    result$TotalNumber <- ob1$TotalNumber + ob2$TotalNumber
    result$CorrectNumber <- ob1$CorrectNumber + ob2$CorrectNumber

    result$TotalHelixesNumber <- ob1$TotalHelixesNumber + ob2$TotalHelixesNumber
    result$CorrectHelixesNumber <- ob1$CorrectHelixesNumber + ob2$CorrectHelixesNumber

    result$TotalStrandsNumber <- ob1$TotalStrandsNumber + ob2$TotalStrandsNumber
    result$CorrectStrandsNumber <- ob1$CorrectStrandsNumber + ob2$CorrectStrandsNumber
     
    result$TotalCoilsNumber <- ob1$TotalCoilsNumber + ob2$TotalCoilsNumber
    result$CorrectCoilsNumber <- ob1$CorrectCoilsNumber + ob2$CorrectCoilsNumber

    return(result)
}

evaluateDataSet <- function(dataSet, model) {
    results <- foreach(i = 1:nrow(dataSet), .combine = 'combineEntriesResult', .export = c('evaluateEntry', 'numberOfStates', 'correctPositions', 'predictHmm', 'cleanResultsVector', 'combineEntriesResult'), .packages = c("HMM")) %dopar% {
        evaluateEntry(dataSet[i,], model)
    }
    return(results)
}
#############Helper functions#############
q3Score <- function(resultsOfTest) {
    testResults <- foreach(i = 1:ncol(resultsOfTest), .combine = 'cbind', .export = c('calculateQ3')) %do% {
        calculateQ3(resultsOfTest[,i])
        }
    return(testResults)
}

calculateQ3 <- function(resultOfTestIteration) {
    result <- list()

    result$TotalQ3 <- resultOfTestIteration$CorrectNumber / resultOfTestIteration$TotalNumber
    result$HelixesQ3 <- resultOfTestIteration$CorrectHelixesNumber / resultOfTestIteration$TotalHelixesNumber
    result$StrandsQ3 <- resultOfTestIteration$CorrectStrandsNumber / resultOfTestIteration$TotalStrandsNumber
    result$CoilsQ3 <- resultOfTestIteration$CorrectCoilsNumber / resultOfTestIteration$TotalCoilsNumber
    result$TestType <- resultOfTestIteration$TestType
    result$NumberOfStates <- resultOfTestIteration$NumberOfStates
    return(result)
}

saveTestResultToFile <- function(result, fileName) {
    write.table(t(result), file = fileName, row.names = FALSE, na = "", col.names = TRUE, sep = ",")
}

generateAllStatesVector <- function(numberOfStates = 0, helixNumber = numberOfStates, strandNumber = numberOfStates, coilNumber = numberOfStates) {
    return(c(generateStateVector("H", helixNumber), generateStateVector("B", strandNumber), generateStateVector("C", coilNumber)))
}

generateStateVector <- function(symbol, numberOfStates) {
    return(paste(symbol, seq(1:numberOfStates), sep=""))
}

#############Testing functions#############
performUniformGrowthTest <- function(testSet, symbols, learningObservations, maxStatesNumber, minStatesNumber = 1, learningIterations = 100, numberOfCores = 1, modelInitialStartProbs = TRUE, modelInitialTransProbs = TRUE, modelInitialEmissionProbs = TRUE) {
    cl <- initializeCluster(numberOfCores)
     testResults <- foreach(i = minStatesNumber:maxStatesNumber, .combine = 'cbind', .export = c('uniformGrowthTestEntry', 'generateAllStatesVector', 'generateStateVector', 'train', 'evaluateDataSet', 'evaluateEntry', 'numberOfStates', 'correctPositions', 'predictHmm', 'cleanResultsVector', 'combineEntriesResult'), .packages = c("HMM", "foreach", "doParallel")) %dopar% {
        uniformGrowthTestEntry(i, symbols, learningObservations, testSet, modelInitialStartProbs, modelInitialTransProbs, modelInitialEmissionProbs, learningIterations)
     }

    stopCluster(cl)
    return(testResults)    
}

uniformGrowthTestEntry <- function(numberOfStates, symbols, learningObservations, testSet, modelInitialStartProbs, modelInitialTransProbs , modelInitialEmissionProbs , learningIterations) {
    prob <- function(x) {
        x / sum(x)
        }

    startProbs <- NULL
    transProbsMatrix <- NULL
    emissionProbsMatrix <- NULL

    states <- generateAllStatesVector(numberOfStates = numberOfStates)

    if (modelInitialStartProbs) {
        startProbs <- prob(runif(length(states)))
    }
    if (modelInitialEmissionProbs) {
        emissionProbsMatrix <- apply(matrix(runif(length(states) * length(symbols)), nrow = length(states), ncol = length(symbols)), 1, prob)
    }
    if (modelInitialTransProbs) {
        transProbsMatrix <- apply(matrix(runif(length(states) * length(states)), nrow = length(states), ncol = length(states)), 1, prob)
    }


    initializedModel <- initHMM(States = states, Symbols = symbols, startProbs = startProbs, transProbs = transProbsMatrix, emissionProbs = emissionProbsMatrix)
    trainedModel <- train(initializedModel, learningObservations, learningIterations)
    result <- evaluateDataSet(testSet, trainedModel$hmm)
    result$TestType <- "Uniform Growth"
    result$NumberOfStates <- numberOfStates
    return(result)
}

performUniformGrowthWithConstHelixes <- function(testSet, symbols, learningObservations, maxStatesNumber, helixesNumber, minStatesNumber = 1, learningIterations = 100, numberOfCores = 1, modelInitialStartProbs = TRUE, modelInitialTransProbs = TRUE, modelInitialEmissionProbs = TRUE){
    cl <- initializeCluster(numberOfCores)
     testResults <- foreach(i = minStatesNumber:maxStatesNumber, .combine = 'cbind', .export = c('uniformGrowthTestEntryWithConstHelixes', 'generateAllStatesVector', 'generateStateVector', 'train', 'evaluateDataSet', 'evaluateEntry', 'numberOfStates', 'correctPositions', 'predictHmm', 'cleanResultsVector', 'combineEntriesResult'), .packages = c("HMM", "foreach", "doParallel")) %dopar% {
        uniformGrowthTestEntryWithConstHelixes(i, symbols, learningObservations, helixesNumber, testSet, modelInitialStartProbs, modelInitialTransProbs, modelInitialEmissionProbs, learningIterations)
     }

    stopCluster(cl)
    return(testResults)   
}

uniformGrowthTestEntryWithConstHelixes <- function(numberOfStates, symbols, learningObservations, helixesNumber, testSet, modelInitialStartProbs, modelInitialTransProbs, modelInitialEmissionProbs, learningIterations) {
    prob <- function(x) {
        x / sum(x)
        }

    startProbs <- NULL
    transProbsMatrix <- NULL
    emissionProbsMatrix <- NULL

    states <- c(generateStateVector("H", helixesNumber), generateStateVector("C", numberOfStates), generateStateVector("B", numberOfStates))

    if (modelInitialStartProbs) {
        startProbs <- prob(runif(length(states)))
    }
    if (modelInitialEmissionProbs) {
        emissionProbsMatrix <- apply(matrix(runif(length(states) * length(symbols)), nrow = length(states), ncol = length(symbols)), 1, prob)
    }
    if (modelInitialTransProbs) {
        transProbsMatrix <- apply(matrix(runif(length(states) * length(states)), nrow = length(states), ncol = length(states)), 1, prob)
    }


    initializedModel <- initHMM(States = states, Symbols = symbols, startProbs = startProbs, transProbs = transProbsMatrix, emissionProbs = emissionProbsMatrix)
    trainedModel <- train(initializedModel, learningObservations, learningIterations)
    result <- evaluateDataSet(testSet, trainedModel$hmm)
    result$TestType <- "Uniform Growth with Constant Helixes Number"
	result$HelixesNumber <- helixesNumber
    result$NumberOfStates <- numberOfStates
    return(result)
}

performUniformGrowthWithConstHelixesAndStrands <- function(testSet, symbols, learningObservations, maxStatesNumber, helixesNumber, strandsNumber, minStatesNumber = 1, learningIterations = 100, numberOfCores = 1, modelInitialStartProbs = TRUE, modelInitialTransProbs = TRUE, modelInitialEmissionProbs = TRUE){
    cl <- initializeCluster(numberOfCores)
     testResults <- foreach(i = minStatesNumber:maxStatesNumber, .combine = 'cbind', .export = c('uniformGrowthTestEntryWithConstHelixesAndStrands', 'generateAllStatesVector', 'generateStateVector', 'train', 'evaluateDataSet', 'evaluateEntry', 'numberOfStates', 'correctPositions', 'predictHmm', 'cleanResultsVector', 'combineEntriesResult'), .packages = c("HMM", "foreach", "doParallel")) %dopar% {
        uniformGrowthTestEntryWithConstHelixesAndStrands(i, symbols, learningObservations, helixesNumber, strandsNumber, testSet, modelInitialStartProbs, modelInitialTransProbs, modelInitialEmissionProbs, learningIterations)
     }

    stopCluster(cl)
    return(testResults)   
}

uniformGrowthTestEntryWithConstHelixesAndStrands <- function(numberOfStates, symbols, learningObservations, helixesNumber, strandsNumber, testSet, modelInitialStartProbs, modelInitialTransProbs, modelInitialEmissionProbs, learningIterations) {
    prob <- function(x) {
        x / sum(x)
        }

    startProbs <- NULL
    transProbsMatrix <- NULL
    emissionProbsMatrix <- NULL

    states <- c(generateStateVector("H", helixesNumber), generateStateVector("C", numberOfStates), generateStateVector("B", strandsNumber))

    if (modelInitialStartProbs) {
        startProbs <- prob(runif(length(states)))
    }
    if (modelInitialEmissionProbs) {
        emissionProbsMatrix <- apply(matrix(runif(length(states) * length(symbols)), nrow = length(states), ncol = length(symbols)), 1, prob)
    }
    if (modelInitialTransProbs) {
        transProbsMatrix <- apply(matrix(runif(length(states) * length(states)), nrow = length(states), ncol = length(states)), 1, prob)
    }


    initializedModel <- initHMM(States = states, Symbols = symbols, startProbs = startProbs, transProbs = transProbsMatrix, emissionProbs = emissionProbsMatrix)
    trainedModel <- train(initializedModel, learningObservations, learningIterations)
    result <- evaluateDataSet(testSet, trainedModel$hmm)
    result$TestType <- "Uniform Growth with Constant Helixes and Strands Number"
    result$HelixesNumber <- helixesNumber
    result$StrandsNumber <- strandsNumber
    result$NumberOfStates <- numberOfStates
    return(result)
}

#############Sample#############
test2 <- function() {
     x <- performUniformGrowthWithConstHelixesAndStrands(testSet, AminoAcidSymbols, observationVecTraining[1:7000], 18, 13, 14, minStatesNumber = 10, learningIterations = 50, numberOfCores = 8)
    saveTestResultToFile(x, "uniformGrowthConstHelixesAndStrands.csv")
    saveTestResultToFile(q3Score(x), "uniformGrowthConstHelixesAndStrandsQ3.csv")
    return(x)
}
