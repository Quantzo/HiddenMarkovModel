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

#Definitions of states used in HMM from figure 1
HelixStates <- c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15")
CoilStates <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12")
StrandStates <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")
AllStates <- c(HelixStates, CoilStates, StrandStates)

#Transition Matrix from figure 1
TransitionMatrix <- t(matrix(c(0, 0.25, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0.25, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16, 0, 0, 0, 0, 0.18, 0.16, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.25, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0.25, 0.25, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.375, 0, 0, 0, 0.25, 0, 0.375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0.375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.25, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0.25, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.375, 0, 0.25, 0, 0.375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.75, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.25), 36))

#Hidden Markov Model initialization
Markov <- initHMM(States = AllStates, Symbols = AminoAcidSymbols, transProbs = NULL, startProbs = NULL, emissionProbs = NULL)

#############Data loading#############
#Loading Data to data frame
dataFrame <- fromJSON("Data/jsonProteins.json", flatten = TRUE)

#############Preparing training and test data sets#############
testRows <- sample(1:nrow(dataFrame), nrow(dataFrame)/3);

trainingSet <- dataFrame[testRows,]
testSet <- dataFrame[ - testRows,]

#############Preparing observation vector for learning#############
observationTraining <- unlist(strsplit(paste(trainingSet[, c("firstStructure")]), ''))

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


performUniformGrowthTest <- function(testSet, symbols, learningObservations, maxStatesNumber, modelStartProbs = NULL, learningIterations = 100, numberOfCores = 1) {
    cl <- initializeCluster(numberOfCores)
     testResults <- foreach(i = 1:maxStatesNumber, .combine = 'cbind', .export = c('uniformGrowthTestEntry', 'generateAllStatesVector', 'generateStateVector', 'train', 'evaluateDataSet', 'evaluateEntry', 'numberOfStates', 'correctPositions', 'predictHmm', 'cleanResultsVector', 'combineEntriesResult' ), .packages = c("HMM", "foreach", "doParallel")) %dopar% {
        uniformGrowthTestEntry(i, symbols, learningObservations, testSet, modelStartProbs, learningIterations)
     }

    stopCluster(cl)
    return(testResults)    
}

generateAllStatesVector <- function(numberOfStates = 0, helixNumber = numberOfStates, strandNumber = numberOfStates, coilNumber = numberOfStates) {
    return(c(generateStateVector("H", helixNumber), generateStateVector("B", strandNumber), generateStateVector("C", coilNumber)))
}

generateStateVector <- function(symbol, numberOfStates) {
    return(paste(symbol, seq(1:numberOfStates), sep=""))
}


uniformGrowthTestEntry <- function(numberOfStates, symbols, learningObservations, testSet, modelStartProbs, learningIterations) {
    states <- generateAllStatesVector(numberOfStates = numberOfStates)
    initializedModel <- initHMM(States = states, Symbols = symbols, startProbs = modelStartProbs)
    trainedModel <- train(initializedModel, learningObservations, learningIterations)
    result <- evaluateDataSet(testSet, trainedModel$hmm)
    result$TestType = "Uniform Growth"
    result$NumberOfStates = numberOfStates
    return(result)
}