#Finalised InputOutput analysis

library("BoolNet")
library("ggplot2")
library("scales") #needed for "alpha" function for transparency
library("reshape2") #needed for melt() of dataframe

#inputs are gene names
#inputActivities is a length(inputs) x 2 matrix, each row specifies the start and end activity level of an input to be screened, in the same order as written in "inputs"
#The intervals provided by inputActivities are split in a total of nrActivityValues levels, e.g. 0->100 and 101 activity values will check every percentage from 0%,1%,...,100% active.
#Bitstrings are generated for the input nodes which have a given percentage of 1s across the entire trajectory. Random start states are chosen to conduct repeated state transitions
#where the state of the input nodes is overwritten every time according to the specified bitstrings. The activity of all nodes across the last nrAveragedStates states in the trajectory 
#is averaged and taken as the pseudo-continuous output value.
#Returns a list with [[1]] and [[2]] containing the mean + std.devs of the output activities in matrices of size nrActivityLevels x length(network$genes).

analyseInputOutputRelations <- function(network, inputs, inputActivities, filename = paste0("InputOutputAnalysisLinePlot.RData"),
                                        nrStartStates = 1000, nrTransitions = 1000, nrAveragedStates = 100, #Should be good default values, discovered by trial & error
                                        nrActivityValues = 101){
  
  ###CHECKING INPUT FOR CONSISTENCY
  stopifnot(inherits(network, "BooleanNetwork") || inherits(network, "SymbolicBooleanNetwork"))
  
  stopifnot(nrow(inputActivities) == length(inputs))
  if (is.numeric(inputs)){
    stopifnot(max(inputs) <= length(network$genes) && min(inputs) >= 1)
  }
  stopifnot(nrow(inputActivities) == length(inputs) && ncol(inputActivities) == 2)
  stopifnot(nrTransitions > nrAveragedStates)
  
  #If nodes are given as gene names, convert to indices
  if (typeof(inputs) == "character"){
    inputnames <- inputs
    for (i in 1:length(inputs)){
      geneIndex <- which(network$genes == inputs[i])
      if (is.na(geneIndex)){stop(paste0("There is no node named ", inputs[i], " in the network."))}
      inputs[i] <- geneIndex
    }
    inputs <- as.numeric(inputs)
  }
  
  ###GENERATE VALUES OF ACTIVITY PERCENTAGES FOR ALL INPUTS
  inputActivityList <- matrix(NA, nrow = length(inputs), ncol = nrActivityValues)
  for (i in 1:length(inputs)){
    if (inputActivities[i,1] == inputActivities[i,2]){#input is fixed at a given activity
      inputActivityList[i,] <- rep(inputActivities[i, 1], nrActivityValues)
    } else {#activity varies
      possibleActivityLevels <- seq(from = inputActivities[i,1],
                                    to = inputActivities[i,2], 
                                    by = (inputActivities[i,2]-inputActivities[i,1])/(nrActivityValues-1))
      inputActivityList[i,] <- possibleActivityLevels
    }
  }
  
  yaxis_means <- yaxis_stds <- matrix(0, nrow = nrActivityValues, ncol = length(network$genes))
  
  
  ###RUN SIMULATION
  #Loop over all activity levels
  pb <- txtProgressBar(min = 0, max = nrActivityValues, style = 3)
  for (a in 1:nrActivityValues){
    setTxtProgressBar(pb, a)
    
    #Generate length nrTransitions+1 input vectors for every input and given activity (inputActivityList[i,a])
    inputBitVectors <- rep(list(NA), length(inputs))
    # Binary vectors of length nrTransitions, having inputActivityList[i,a] percent 1s
    
    #Every network component is an output
    outputsByState <- matrix(NA, nrow = nrStartStates, ncol=length(network$genes))
    for (run in 1:nrStartStates){ #All start states get a different bitstring for a given activity%
      # --> More variation introduced, but likely makes simulation take longer -> Compare run times
      #called mode = random before: 1s randomly distributed along entire bitstring -> mode=block option removed!
      for (i in 1:length(inputs)){
        input_a <- inputActivityList[i,a]
        inputBitVectors[[i]] <- rep(0, nrTransitions+1)
        nrSample1s <- round((nrTransitions+1)*input_a/100)
        index1s <- sample(seq(1,nrTransitions+1), size = nrSample1s, replace = FALSE)
        inputBitVectors[[i]][index1s] <- 1
      }
      
      #Generate random initial state
      state <- sample(c(0,1), size=length(network$genes), replace=TRUE)
      outputActivities <- rep(0, length(network$genes))
      
      #Make initial state conform to activity specifications of predetermined inputBitVectors
      for (i in 1:length(inputs)){
        state[inputs[i]] <- inputBitVectors[[i]][1]
      }
      
      #Loop, do repeated state transitions for the chosen start state at this activity%
      for (s in 1:nrTransitions){
        oldstate <- state
        state <- stateTransition(network, state = state)
        
        #Re-fix input nodes in newly obtained state
        for (i in 1:length(inputs)){
          state[inputs[i]] <- inputBitVectors[[i]][s+1]
        }
        
        #Final nrAveragedStates are used for averaging output node activity
        if (s > nrTransitions-nrAveragedStates){
          outputActivities <- outputActivities + state
        }
        
      }#all state transitions done
      
      outputActivities <- outputActivities/nrAveragedStates
      outputsByState[run,] <- outputActivities #output activities given this particular start state
    }#all runs over random start states done
    
    #For every gene (column), get mean+std across all runs (rows) -> save to [a,g] of yaxis_means and yaxis_stds matrices
    for (g in 1:length(network$genes)){
      yaxis_means[a,g] <- mean(outputsByState[,g])
      yaxis_stds[a,g] <- sd(outputsByState[,g])
    }
  }#all runs over activity values done
  colnames(yaxis_means) <- colnames(yaxis_stds) <- network$genes
  rownames(inputActivityList) <- inputnames
  
  return(list(yaxis_means, yaxis_stds, inputActivityList))
}#[[1]] is matrix of nrActivities  x NrGenes with mean output activities, [[2]] has corresponding std.devs, [[3]] has levels of input activity percentages for all inputs


#Plot the pseudo-continuous output response of the plottedOutputGene vs the activity of plottedInputGene.
#Red line is mean activity, lightblue ribbon indicates +/- one std.dev around the mean.
#Works for an arbitrary number of varied inputs, although only one input can be shown for the x-axis.
#If x-axis i.e. [[3]] is inverted 100->0, x axis reads still 0->100 from left to right but plot becomes flip automatically to match the result
#Input-Output relations will be analysed if the optional argument simResult is not provided. This argument should have the structure of the return given by
#the function analyseInputOutputRelations(). If this is given, the function directly skips to generating a line plot of the output activity, showing mean+StdDev of the output
#as a function of the varying input activities. Any node in the network can be chosen as the output for the plot.
InputOutputLinePlot <- function(network, inputs, inputActivities, plottedInputGene, plottedOutputGene, nrStartStates = 1000, nrTransitions = 1000, nrAveragedStates = 100, 
                                nrActivityValues = 101, filename = "InputOutputAnalysisLinePlot.RData", meanLineColor="red", StdDevColor="lightblue", 
                                StdDevAlpha = 0.5, nr_xticks = 10, xlab = "Input activity", ylab = "Output activity", title = "Average activity response of output node", 
                                cex.lab = 1.3, cex.main=1.5, simResult){
  
  if(missing(simResult)){
    print(paste0("Analysing input-output relations, result will be saved under ", filename, "."))
    simResult <- analyseInputOutputRelations(network=network, inputs=inputs, inputActivities = inputActivities, filename = filename,
                                nrStartStates = nrStartStates, nrTransitions = nrTransitions, nrAveragedStates = nrAveragedStates,
                                nrActivityValues = nrActivityValues)
    save(simResult, file=filename)
  }
  
  #Convert gene names to indices
  plottedInputGeneIndex <- which(network$genes == plottedInputGene)
  plottedOutputGeneIndex <- which(network$genes == plottedOutputGene)
  
  x <- simResult[[3]][which(rownames(simResult[[3]]) == plottedInputGene),]
  means <- simResult[[1]][,plottedOutputGeneIndex]
  lowerribbon <- means - simResult[[2]][,plottedOutputGeneIndex]
  upperribbon <- means + simResult[[2]][,plottedOutputGeneIndex]
  plot(x, means, type = "l",
       xlim = c(min(x), max(x)), ylim = c(0, 1), xaxs="i",yaxs="i", cex.lab=cex.lab, cex.main=cex.main, yaxt="n", xaxt="n",
       xlab = xlab, ylab = ylab, main=title)
  axis(2, at = seq(0, 1, by = 0.1), las=2) 
  axis(1, at = seq(min(x), max(x), by = (max(x)-min(x))/nr_xticks), las=1)
  polygon(c(x, rev(x)), c(lowerribbon, rev(upperribbon)),
          col=alpha(StdDevColor, StdDevAlpha))
  lines(x, means, type = "l", col = meanLineColor)
  return(simResult)
}


#Plot a heatmap showing the pseudo-continuous response of the output gene. Only works if exactly two inputs have been varied.
#simResultList is an optional argument. If this is not given, a simulation will be performed varying both inputs over their respective ranges independently.
#The return of this function can be given as input for simResultList in which case the simulation will be skipped and only the plot will be modified.
InputOutputHeatmap <- function(network, inputs, plottedOutputGene, inputActivities, nrStartStates = 1000, nrTransitions = 1000, nrAveragedStates = 100, nrActivityValues = 101, 
                               colorpalette = "RdBu", filename = "InputOutputAnalysisHeatmap.RData",
                               xlab = "Input 1 activity", ylab = "Input 2 activity",
                               title = "Average activity response of output node", simResultList){
  
  input1 <- inputs[1]
  input2 <- inputs[2]
  
  stopifnot(length(inputs)==2)
  stopifnot(input1 %in% network$genes)
  stopifnot(input2 %in% network$genes)
  stopifnot(plottedOutputGene %in% network$genes)
  
  if (missing(simResultList)){#perform simulation if no result given as input
    simResultList <- rep(list(NA), nrActivityValues+1)
    
    #Generate input activity % levels, i.e. what would be saved in [[3]] in final result
    inputActivityList <- matrix(NA, nrow = length(inputs), ncol = nrActivityValues)
    for (i in 1:length(inputs)){
      if (inputActivities[i,1] == inputActivities[i,2]){#input is fixed at a given activity
        inputActivityList[i,] <- rep(inputActivities[i, 1], nrActivityValues)
      } else {#activity varies
        possibleActivityLevels <- seq(from = inputActivities[i,1],
                                      to = inputActivities[i,2], 
                                      by = (inputActivities[i,2]-inputActivities[i,1])/(nrActivityValues-1))
        inputActivityList[i,] <- possibleActivityLevels
      }
    }
    rownames(inputActivityList) <- inputs
    
    simResultList[[nrActivityValues+1]] <- inputActivityList #Access last list element for structure of axes in plot
    
    for (i in 1:nrActivityValues){
      inputActivities_firstInputFixed <- inputActivities
      inputActivities_firstInputFixed[1,] <- rep(inputActivityList[1,i],2)
      print(inputActivities_firstInputFixed)
      simResultList[[i]] <- analyseInputOutputRelations(network = network, inputs = inputs, inputActivities = inputActivities_firstInputFixed, 
                                               filename = paste0(as.character(network), ".RData"), nrStartStates = nrStartStates, nrTransitions = nrTransitions, 
                                               nrAveragedStates = nrAveragedStates, nrActivityValues = nrActivityValues)[[1]] #take only means from this result
    }
    save(simResultList, file=filename)
  }#end simulation
  
  
  #Generate matrix from means over all the different runs for the specified plottedOutputGene argument
  plottedOutputGeneIndex <- which(network$genes == plottedOutputGene)
  outputValues2plot <- matrix(NA, nrow = nrActivityValues, ncol = nrActivityValues)
  for (c in 1:nrActivityValues){
    outputValues2plot[c,] <- simResultList[[c]][,plottedOutputGeneIndex]
  }
  
  #Generate plot
  xactivities <- simResultList[[nrActivityValues+1]][1,]
  yactivities <- simResultList[[nrActivityValues+1]][2,]
  df <- melt(outputValues2plot)
  print(outputValues2plot)
  print(df)
  df$Var1 <- xactivities[df$Var1]
  df$Var2 <- yactivities[df$Var2]
  print(df)
  colnames(df) <- c("Input1activity", "Input2activity", "Output")
  p <- ggplot(df, aes(x = Input1activity, y = Input2activity, fill = Output)) + geom_tile() + coord_equal() + theme_bw() + 
    scale_fill_distiller(palette = colorpalette, direction = 1) + ggtitle(title) + xlab(xlab) + ylab(ylab)
  plot(p)
  
  return(simResultList)
}

#required functions for attractor fusion
#This function takes a specific attractor from an attractor object and calculates its averaged activities over all its states weighted for the fraction of its basin of attraction.
AttractorActivityWeightedBasin<- function(attractor, nr){
  totalattrs <- length(attractor$attractors)
  length <- length(attractor$stateInfo$genes)
  #calculateTotalSizeof Basin of attraction
  Basins<- c()
  for (a in 1:length(attractor$attractors)){
    BasinAttr <- attractor$attractors[[a]]$basinSize
    Basins<- c(Basins, BasinAttr)
    TotalSize<- sum(Basins)
  }
  activities <- apply(attractor[[2]][[nr]]$involvedStates, MARGIN = 2, BoolNet:::dec2bin, len = length)
  ProbofOne<- apply(activities, MARGIN = 1, function(x){(sum(x)/ncol(activities))*(((attractor$attractors[[nr]]$basinSize)/TotalSize)*100)})
  ProbofOne
  FinalMatrix <- cbind(ProbofOne)
  rownames(FinalMatrix) <- attractor$stateInfo$genes
  return(FinalMatrix)
  
}

#This fucntion does the same as above but sums over all existing attractors
OverallActivities<- function(attractor){
  AllIn<- matrix(0,nrow=length(attractor$stateInfo$genes) , ncol=1)
  for (a in 1:length(attractor$attractors)){
    NewCalculation<- AttractorActivityWeightedBasin(attractor, a)
    AllIn<- cbind(AllIn, NewCalculation[,1])
  }
  AllIn<-AllIn[,-1]
  #return(AllIn)
  return(as.data.frame(rowSums(AllIn)))
}

#perturbation/stability test
modTestNetworkProperties <- function (network, numRandomNets = 100, testFunction = "testIndegree",
                                      testFunctionParams = list(), accumulation = c("characteristic",
                                                                                    "kullback_leibler"), alternative = c("greater", "less"),
                                      sign.level = 0.05, drawSignificanceLevel = TRUE, klBins,
                                      klMinVal = 1e-05, linkage = c("uniform", "lattice"), functionGeneration = c("uniform",
                                                                                                                  "biased"), validationFunction, failureIterations = 10000,
                                      simplify = FALSE, noIrrelevantGenes = TRUE, d_lattice = 1,
                                      zeroBias = 0.5, title = "", xlab, xlim, breaks = 30, ...)
{
  stopifnot(inherits(network, "BooleanNetwork") || inherits(network,
                                                            "SymbolicBooleanNetwork"))
  if (is.character(testFunction))
    testFunctionName <- testFunction
  else testFunctionName <- ""
  testFunction <- match.fun(testFunction)
  accumulate <- (match.arg(accumulation) == "characteristic")
  origResult <- testFunction(network, accumulate, testFunctionParams)
  numGenes <- length(network$interactions)
  if (inherits(network, "SymbolicBooleanNetwork"))
    inputGenes <- sapply(network$interactions, function(interaction) length(getInputs(interaction)))
  else inputGenes <- sapply(network$interactions, function(interaction) length(interaction$input))
  if (missing(validationFunction))
    validationFunction <- NULL
  randomResults <- lapply(seq_len(numRandomNets), function(i) {
    randomNet <- generateRandomNKNetwork(n = numGenes, k = inputGenes,
                                         topology = "fixed", linkage = linkage, functionGeneration = functionGeneration,
                                         validationFunction = validationFunction, failureIterations = failureIterations,
                                         simplify = simplify, noIrrelevantGenes = noIrrelevantGenes,
                                         d_lattice = d_lattice, zeroBias = zeroBias)
    randomRes <- testFunction(randomNet, accumulate, testFunctionParams)
    return(randomRes)
  })
  if (accumulate)
    randomResults <- unlist(randomResults)
  args <- list(...)
  res <- switch(match.arg(accumulation, c("characteristic",
                                          "kullback_leibler")), characteristic = {
                                            if (missing(xlab)) {
                                              xlab <- switch(testFunctionName, testIndegree = "Gini index of state in-degrees",
                                                             testAttractorRobustness = "% of identical attractors",
                                                             testTransitionRobustness = "Normalized Hamming distance",
                                                             "accumulated results")
                                            }
                                            if (missing(xlim)) {
                                              xlim <- range(c(origResult, randomResults))
                                            }
                                            alternative <- match.arg(alternative, c("greater", "less"))
                                            if (alternative == "greater") pval <- sum(randomResults <
                                                                                        origResult)/length(randomResults) else pval <- sum(randomResults >
                                                                                                                                             origResult)/length(randomResults)
                                            if (testFunctionName == "testIndegree" | testFunctionName ==
                                                "testAttractorRobustness") {
                                              r <- hist(randomResults, xlim = xlim, xlab = xlab,
                                                        main = title, xaxt = "n", ...)
                                              axis(side = 1, at = seq(xlim[1], xlim[2], length.out = 11))
                                            } else {
                                              r <- hist(randomResults, xlim = xlim, xlab = xlab,
                                                        main = title, ...)
                                            }
                                            print("orig Result :")
                                            print(origResult)
                                            abline(v = origResult, col = "red")
                                            if (alternative == "greater") text(x = origResult, pos = 2,
                                                                               y = max(r$counts) * 0.75, labels = paste("> ", round(pval *
                                                                                                                                      100), "%\nof random results", sep = ""), col = "red",
                                                                               cex = 0.75) else text(x = origResult, pos = 4, y = max(r$counts) *
                                                                                                       0.75, labels = paste("< ", round(pval * 100), "%\nof random results",
                                                                                                                            sep = ""), col = "red", cex = 0.75)
                                            if (drawSignificanceLevel) {
                                              if (alternative == "greater") {
                                                quant <- quantile(randomResults, 1 - sign.level)
                                                abline(v = quant, col = "blue")
                                                text(x = quant, pos = 2, y = max(r$counts) *
                                                       0.85, labels = paste((1 - sign.level) * 100,
                                                                            "% quantile", sep = ""), col = "blue", cex = 0.75)
                                              } else {
                                                quant <- quantile(randomResults, sign.level)
                                                print("random Result :")
                                                print(quant)
                                                abline(v = quant, col = "blue")
                                                text(x = quant, pos = 4, y = max(r$counts) *
                                                       0.85, labels = paste(sign.level * 100, "% quantile",
                                                                            sep = ""), col = "blue", cex = 0.75)
                                              }
                                            }
                                            list(hist = r, pval = 1 - pval, significant = (1 - pval <=
                                                                                             sign.level))
                                          }, kullback_leibler = {
                                            if (missing(xlab)) xlab <- "Kullback-Leibler distance"
                                            if (missing(klBins)) {
                                              bins <- unique(c(origResult, unlist(randomResults)))
                                              bins <- c(bins, max(bins) + 1)
                                            } else {
                                              bins <- unique(c(origResult, unlist(randomResults)))
                                              if (klBins < length(bins)) bins <- seq(min(bins),
                                                                                     max(bins), length.out = klBins + 1) else bins <- c(bins,
                                                                                                                                        max(bins) + 1)
                                            }
                                            vals <- sapply(randomResults, function(results) kullbackLeiblerDistance(origResult,
                                                                                                                    results, bins = bins, minVal = klMinVal))
                                            r <- hist(vals, xlab = xlab, main = title, breaks = breaks,
                                                      ...)
                                            list(hist = r)
                                          }, stop("'accumulation' must be one of \"characteristic\",\"kullback_leibler\""))
  return(res)
}



