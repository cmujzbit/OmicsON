
# NEW PUBLIC API:
clusterUsingOntology <- function(chebiIdsDataFrame, rootColumnName, ontologyRepresentatnion) {
    ontologyDataFrame <- ontologyRepresentatnion(baseData = chebiIdsDataFrame, rootColumnName = rootColumnName)
    ontologyDataFrame
}


# NEW PUBLIC API:
mapReactomePathwaysUnderOrganism <- function(chebiOntologyIds, organismTaxonomyId='9606', idsColumnName = 'ontologyId', rootColumnName = 'root') {
    # x <- c("supp", "dose")
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(idsColumnName, rootColumnName)
    }
    chebiIdsToEnsembleIds <- ddply(.data = chebiOntologyIds, columnsUseInIteration, .fun = function(vectorElement) {

        # print(as.character(vectorElement[1, c(idsColumnName)]))
        idToCheck <- as.character(strsplit(as.character(vectorElement[1, c(idsColumnName)]), ":")[[1]][2])
        # print("idToCheck")
        # print(idToCheck)
        if (is.na(idToCheck)) {
            idToCheck <- as.character("0");
            # print(idToCheck)
        }
        pathwayIds <- getPathwaysIdsForChebiUnderOrganism(idToCheck, taxonIdToReactomeCodes[[organismTaxonomyId]]$speciesCode)
        ensembleIds <- getEnsemblIdsForPathwayIds(pathwayIds)
        uniProtIds <- getUniProtIdsForPathwayIds(pathwayIds)
        # print(pathwayIds)
        # print(ensembleIds)
        # print(uniProtIds)
        genesSymbolsFromEnsemble <- character(0)
        if (0 != length(ensembleIds)) {
            genesSymbolsFromEnsemble <- getSymbolsBaseOnEnsemblGensIdsUsingMyGenePackage(ensembleIds, organismTaxonomyId = organismTaxonomyId)
            # print(genesSymbolsFromEnsemble)
            # str(genesSymbolsFromEnsemble)
        }
        genesSymbolsFromUniProt <- character(0)
        if (0 != length(uniProtIds)) {
            genesSymbolsFromUniProt <- getSymbolsBaseOnUniProtIdsUsingMyGenePackage(uniProtIds, organismTaxonomyId = organismTaxonomyId)
        }

        chebiIdToEnsembleIds <- data.frame('ensembleIds' = I(list(ensembleIds)),
                                           'uniProtIds' = I(list(uniProtIds)),
                                           'reactomeIds' = I(list(pathwayIds)),
                                           'genesSymbolsFromEnsemble' = I(list(genesSymbolsFromEnsemble)),
                                           'genesSymbolsFromUniProt' = I(list(genesSymbolsFromUniProt)))
        chebiIdToEnsembleIds
    })
    chebiIdsToEnsembleIds[,c(2, 1, 3:7)]
}


#Static hashmap taxonId <-> (reactome$speciesName, reactome$spaciesCode) Code is used in resouce files.
taxonIdToReactomeCodes <- new.env()
taxonIdToReactomeCodes[['9606']] <- list(speciesName='Homo sapiens', speciesCode='HSA')
taxonIdToReactomeCodes[['3702']] <- list(speciesName='Arabidopsis thaliana', speciesCode='ATH')
taxonIdToReactomeCodes[['9913']] <- list(speciesName='Bos taurus', speciesCode='BTA')
taxonIdToReactomeCodes[['6239']] <- list(speciesName='Caenorhabditis elegans', speciesCode='CEL')
taxonIdToReactomeCodes[['9615']] <- list(speciesName='Canis familiaris', speciesCode='CFA')
taxonIdToReactomeCodes[['7955']] <- list(speciesName='Danio rerio', speciesCode='DRE')
taxonIdToReactomeCodes[['44689']] <- list(speciesName='Dictyostelium discoideum', speciesCode='DDI')
taxonIdToReactomeCodes[['7227']] <- list(speciesName='Drosophila melanogaster', speciesCode='DME')
taxonIdToReactomeCodes[['9031']] <- list(speciesName='Gallus gallus', speciesCode='GGA')
taxonIdToReactomeCodes[['10090']] <- list(speciesName='Mus musculus', speciesCode='MMU')
taxonIdToReactomeCodes[['1773']] <- list(speciesName='Mycobacterium tuberculosis', speciesCode='MTU')
taxonIdToReactomeCodes[['4530']] <- list(speciesName='Oryza sativa', speciesCode='OSA')
taxonIdToReactomeCodes[['5833']] <- list(speciesName='Plasmodium falciparum', speciesCode='PFA')
taxonIdToReactomeCodes[['10116']] <- list(speciesName='Rattus norvegicus', speciesCode='RNO')
taxonIdToReactomeCodes[['4932']] <- list(speciesName='Saccharomyces cerevisiae', speciesCode='SCE')
taxonIdToReactomeCodes[['4896']] <- list(speciesName='Schizosaccharomyces pombe', speciesCode='SPO')
taxonIdToReactomeCodes[['8364']] <- list(speciesName='Xenopus tropicalis', speciesCode='XTR')
taxonIdToReactomeCodes[['59729']] <- list(speciesName='Taeniopygia guttata', speciesCode='TGU')
taxonIdToReactomeCodes[['9823']] <- list(speciesName='Sus scrofa', speciesCode='SSC')


# NEW PUBLIC API
getStringNeighbours <- function(chebiIdsToReactomePathways, stringOrganismId = 9606, stringDbVersion = "10", idsColumnName = 'ontologyId', rootColumnName = 'root', listOfEnsembleIdColumnName = 'ensembleIds') {
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(idsColumnName, rootColumnName)
    }
    string_db <- STRINGdb$new(
        version = stringDbVersion,
        species = stringOrganismId,
        input_directory = path.expand("~"))
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[!chebiIdsToReactomePathways[idsColumnName] == 0, ]
    dfWithString <- ddply(.data = chebiIdsToRealReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
        returnNeighbourVector <- character(length = 0)
        stringGenesSymbols <- character(length = 0)
        if (0 == length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
        } else {
            proteinIds <- getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(
                dfElement[1,listOfEnsembleIdColumnName], organismTaxonomyId = stringOrganismId
            )
            stringId1 <- string_db$mp(proteinIds)
            returnNeighbourVector <- string_db$get_neighbors(stringId1)
            ensembleIdsFromStringDb <- mapFromStringIdsToEnsembleIds(returnNeighbourVector)
            stringGenesSymbols <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDb, organismTaxonomyId = stringOrganismId)
        }
        dffff <- data.frame('ensembleIds' = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(returnNeighbourVector))),
                            'stringGenesSymbols' = I(list(unique(stringGenesSymbols))) )
        dffff
    })
    dfWithString
}


# NEW API
mapFromStringIdsToEnsembleIds <- function(vactofOfStringIds) {
    ensembleIds <- laply(vactofOfStringIds, .fun = function(vectorElement){
        ensemblePeptideId <- strsplit(vectorElement, "[.]")[[1]][2]
        ensemblePeptideId
    })
    ensembleIds
}


# NEW API.
getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {

    additionalInformationBaseOnEnsemblGenId <- invisible(
        mygene::queryMany(qterms = gensIdsVector[[1]], scopes = c("ensembl.gene"),
                  fields = c("symbol","ensembl.protein")))
    equivalentEnsemlProteinsIds <- unlist(additionalInformationBaseOnEnsemblGenId$ensembl)

    equivalentEnsemlProteinsIdsVector <- lapply(equivalentEnsemlProteinsIds, function(characterListElement){
        characterListElement
    });
    equivalentEnsemlProteinsIdsVector <- as.character(unlist(equivalentEnsemlProteinsIdsVector))
    equivalentEnsemlProteinsIdsVector
}


# NEW API.
getSymbolsBaseOnEnsemblGensIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    # genes$symbol
    # genes$ensembl.protein
    # gensIdsVectorGlobal <<- gensIdsVector
    # print(gensIdsVector)
    additionalInformationBaseOnEnsemblGenId <- invisible(queryMany(gensIdsVector, scopes = 'ensembl.gene', fields = c("symbol","ensembl.protein"),
                                                         species = organismTaxonomyId))
    # print(additionalInformationBaseOnEnsemblGenId)
    equivalentEnsemlProteinsIdsVector <- unlist(
        additionalInformationBaseOnEnsemblGenId$symbol[!is.na(additionalInformationBaseOnEnsemblGenId$symbol)]
    )
    equivalentEnsemlProteinsIdsVector <- as.character(equivalentEnsemlProteinsIdsVector)
    equivalentEnsemlProteinsIdsVector
}

# NEW API.
getSymbolsBaseOnUniProtIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    # genes$symbol
    # genes$ensembl.protein
    additionalInformationBaseOnEnsemblGenId <- invisible(queryMany(gensIdsVector, scopes = 'uniprot', fields = c("symbol"),
                                                         species = organismTaxonomyId))
    equivalentEnsemlProteinsIdsVector <- unlist(
        additionalInformationBaseOnEnsemblGenId$symbol[!is.na(additionalInformationBaseOnEnsemblGenId$symbol)]
    )
    equivalentEnsemlProteinsIdsVector <- as.character(equivalentEnsemlProteinsIdsVector)
    equivalentEnsemlProteinsIdsVector
}


# NEW API.
getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    equivalentEnsemlProteinsIdsVector <- character(0)
    if (0 != length(gensIdsVector)){
        additionalInformationBaseOnEnsemblPeptidId <- invisible(queryMany(
            gensIdsVector, scopes = 'ensemblprotein'
            , fields = c("symbol","ensembl.protein"), species = organismTaxonomyId
        ))
        equivalentEnsemlProteinsIdsVector <- unlist(
            additionalInformationBaseOnEnsemblPeptidId$symbol[!is.na(additionalInformationBaseOnEnsemblPeptidId$symbol)]
        )
        equivalentEnsemlProteinsIdsVector <- as.character(equivalentEnsemlProteinsIdsVector)
    }
    equivalentEnsemlProteinsIdsVector
}

# NEW PUBLIC API
readGroupsAsDf <- function(pathToFileWithGroupDefinition ) {
    maxColLength <- max(count.fields(pathToFileWithGroupDefinition, sep = '\t'))
    model <- read.table(file = pathToFileWithGroupDefinition, header = TRUE, fill = TRUE,
                        stringsAsFactors = FALSE, sep = "\t", strip.white = TRUE)
}

# PRIVATE
reduceCorrelation <- function(matrix, cutoff = 1) {
    corMatrix <- cor(x = matrix)
    corMatrix[upper.tri(corMatrix)] <- 0
    diag(corMatrix) <- 0
    matrix <- matrix[,!apply(corMatrix,2,function(x) any(x > cutoff))]
    matrix
}

# NEW API CCA
makeCanonicalCorrelationAnalysis <- function(xNamesVector, yNamesVector, XDataFrame, YDataFrame,
                                             xCutoff = 1, yCutoff = 1, scalingFactor = 1) {

    XDataFrame[seq(2,length(XDataFrame))] <- XDataFrame[seq(2,length(XDataFrame))] * scalingFactor
    YDataFrame[seq(2,length(YDataFrame))] <- YDataFrame[seq(2,length(YDataFrame))] * scalingFactor

    commonColNames <- intersect(colnames(XDataFrame), colnames(YDataFrame))
    XDataFrame <- XDataFrame[
        ,c(colnames(XDataFrame)[1],commonColNames)]
    YDataFrame <- YDataFrame[
        ,c(colnames(YDataFrame)[1],commonColNames)]

    # Where XData = transcriptomicsData and YData = lipidomicsData.
    XData <- data.frame(XDataFrame[!duplicated(XDataFrame[1]), ], row.names = 1)
    YData <- data.frame(YDataFrame[!duplicated(YDataFrame[1]), ], row.names = 1)

    transposedXData <- as.data.frame(t(XData))
    transposedYData <- as.data.frame(t(YData))

    interX <- intersect(colnames(transposedXData), xNamesVector)
    interY <- intersect(colnames(transposedYData), yNamesVector)

    X <- transposedXData[as.character(interX)]
    Y <- transposedYData[as.character(interY)]

    X <- as.matrix(X)
    Y <- as.matrix(Y)

    X <- OmicsON:::reduceCorrelation(X, cutoff = xCutoff)
    Y <- OmicsON:::reduceCorrelation(Y, cutoff = yCutoff)

    cca.fit <- NULL

    if (!length(X)) {
        print("CCA is not possible - XDataFrame do not have specified names in xNamesVector")
    } else if (!length(Y)) {
        print("CCA is not possible - YDataFrame do not have specified names in yNamesVector")
    } else {
        cca.fit <- tryCatch(
            {
                yacca::cca(X, Y)
            },
            error = function(cond) {
                message("OmicsON - Original message from yacca:")
                message(cond$call, appendLF = TRUE)
                message(cond$message, appendLF = TRUE)
                message("OmicsON - CCA with used inputs can not solve task, proposed improvement:",
                        appendLF = TRUE)
                if (cond$call == "cor(x, o$canvarx, use = use)") {
                    message("Increase xCutoff or yCutoff value.", appendLF = TRUE)
                }
                # Choose a return value in case of error
                return(NA)
            },
            warning = function(cond) {
                message("OmicsON - Included CCA present warning.")
                message("Original message (yacca):")
                message(cond)
                # Choose a return value in case of warning
                return(NULL)
            },
            finally = {
                message("OmicsON - CCA (yacca) finished.")
            }
        )
    }
    cca.fit$xCutoff <- xCutoff
    cca.fit$yCutoff <- yCutoff
    cca.fit$functionalGrouping$xData <- interX
    cca.fit$functionalGrouping$yData <- interY
    cca.fit$correlationCutOff$xData <- colnames(X)
    cca.fit$correlationCutOff$yData <- colnames(Y)
    cca.fit
}


# NEW PUBLIC API
makePermutationTestOnCCA <- function(XDataFrame, YDataFrame, numberOfRowsForTestOnX, numberOfRowsForTestOnY, numberOfIterations = 100, countedCCA) {
    # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    # print(numberOfRowsForTestOnX)
    # print(numberOfRowsForTestOnY)
    # print()
    # print()
    vectorOfXrd <- as.numeric();
    vectorOfYrd <- as.numeric();
    for (i in 1:numberOfIterations) {
        xNV <- as.character(XDataFrame[sample(nrow(XDataFrame), numberOfRowsForTestOnX), ][,1])
        yNV <- as.character(YDataFrame[sample(nrow(YDataFrame), numberOfRowsForTestOnY), ][,1])
        ccaResult <- OmicsON::makeCanonicalCorrelationAnalysis(xNamesVector = xNV, yNamesVector = yNV, XDataFrame = XDataFrame, YDataFrame = YDataFrame)
        # print("++++++++++++++++++++++++++++++++++")
        # print.default(ccaResult)
        if (is.na(ccaResult) || is.null(ccaResult)) {

        } else {
            vectorOfXrd <- c(vectorOfXrd, ccaResult$xrd)
            vectorOfYrd <- c(vectorOfYrd, ccaResult$yrd)
        }
    }

    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print.default(countedCCA)
    # print("*****************************")
    # print(countedCCA$xrd)
    # print(countedCCA$yrd)
    # print(vectorOfXrd)
    # print(vectorOfYrd)
    # print(mean(vectorOfXrd))
    # print(mean(vectorOfYrd))
    meanOfXrd <- NA;
    meanOfYrd <- NA;
    if (0 != length(vectorOfXrd)) {
        meanOfXrd <- mean(vectorOfXrd);
    }
    if (0 != length(vectorOfYrd)) {
        meanOfYrd <- mean(vectorOfYrd);
    }
    testResult <- list("countedCCAOnX" = countedCCA$xrd,
                       "countedCCAOnY" = countedCCA$yrd,
                       "meanOnX" = meanOfXrd,
                       "meanOnY" = meanOfYrd)
    testResult
}


# NEW PUBLIC API
makeCCAOnGroups <- function(groupsDefinitionDF, mappingDF, leftMappingColumnName = 'root', rightMappingColumnName = 'genesSymbolsFromEnsemble', groupsDataDF, mappingDataDF){
    ddply(.data = groupsDefinitionDF['Molecules'], .(Molecules), .fun = function(dfElement) {
        # print("???????????????????")
        # print(dfElement)
        rightSideIdsToAnalys <- unlist(strsplit(as.character(dfElement), split = " "));
        # print(rightSideIdsToAnalys)
        leftSideIdsToAnalys <- mappingDF[mappingDF[[leftMappingColumnName]] %in% rightSideIdsToAnalys,][[rightMappingColumnName]]
        leftSideIdsToAnalys <- unique(unlist(leftSideIdsToAnalys))

        ccaResults <- OmicsON::makeCanonicalCorrelationAnalysis(
            xNamesVector = leftSideIdsToAnalys,
            yNamesVector = rightSideIdsToAnalys,
            XDataFrame = mappingDataDF,
            YDataFrame = groupsDataDF)

        # print("######################################")
        # print(ccaResults)
        # print("---------------------------")
        # print.default(ccaResults)
        #TODO : Use user defined column name instead of symbol. ChEBI column too.
        numberOfRowsForTestOnX <- nrow(mappingDataDF[mappingDataDF$symbol %in% leftSideIdsToAnalys, ])
        numberOfRowsForTestOnY <- nrow(groupsDataDF[groupsDataDF$ChEBI %in% rightSideIdsToAnalys, ])
        parmutationTestResult <- makePermutationTestOnCCA(XDataFrame = mappingDataDF, YDataFrame = groupsDataDF,
                                                          numberOfRowsForTestOnX = numberOfRowsForTestOnX,
                                                          numberOfRowsForTestOnY = numberOfRowsForTestOnY,
                                                          numberOfIterations = 50, countedCCA = ccaResults);

        dfWithCca <- data.frame('right' = I(list(rightSideIdsToAnalys)),
                                'left' = I(list(leftSideIdsToAnalys)),
                                'ccaResults' = I(list(ccaResults)),
                                'ccaPermutationTestResults' = I(list(parmutationTestResult)))
        dfWithCca
    })
}


# NEW PUBLIC API
plotCanonicalCorrelationAnalysisResults <- function(ccaResults,
                                                    x.name = "Transcriptomics",
                                                    y.name = "Lipidomics",
                                                    cv = 1, thirdLineText = "",
                                                    ...) {

    helio.plot(ccaResults, x.name = x.name, y.name = y.name,
               sub = paste("Canonical Variate", cv, "\n",
                           "xCutoff = ", ccaResults$xCutoff, ",",
                           "yCutoff = ", ccaResults$yCutoff, "\n",
                           thirdLineText, sep = " "), ...)
}


# NEW PUBLIC API
makePartialLeastSquaresRegression <- function(xNamesVector, yNamesVector,
                                              XDataFrame, YDataFrame,
                                              formula = Y ~ X,
                                              ncompValue = 10,
                                              xCutoff = 1, yCutoff = 1, ...) {

    commonColNames <- intersect(colnames(XDataFrame), colnames(YDataFrame))
    XDataFrame <- XDataFrame[
        ,c(colnames(XDataFrame)[1],commonColNames)]
    YDataFrame <- YDataFrame[
        ,c(colnames(YDataFrame)[1],commonColNames)]

    # Where XData = transcriptomicsData and YData = lipidomicsData.
    XData <- data.frame(XDataFrame[!duplicated(XDataFrame[1]), ], row.names = 1)
    YData <- data.frame(YDataFrame[!duplicated(YDataFrame[1]), ], row.names = 1)

    transposedXData <- as.data.frame(t(XData))
    transposedYData <- as.data.frame(t(YData))

    interX <- intersect(colnames(transposedXData), xNamesVector)
    interY <- intersect(colnames(transposedYData), yNamesVector)

    X <- transposedXData[as.character(interX)]
    Y <- transposedYData[as.character(interY)]

    X <- as.matrix(X)
    Y <- as.matrix(Y)

    X <- OmicsON:::reduceCorrelation(X, cutoff = xCutoff)
    Y <- OmicsON:::reduceCorrelation(Y, cutoff = yCutoff)

    Xmelt <- I(X)
    Ymelt <- I(Y)

    combined <- data.frame(X = I(Xmelt), Y = I(Ymelt))

    PLSResults <- NULL

    PLSResults <- tryCatch(
        {
            pls::plsr(formula = formula, data = combined, ncomp = ncompValue, ...)
        },
        error = function(cond) {
            message("OmicsON - Included PLS can not solve task.")
            message("Original message (pls):")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning = function(cond) {
            message("OmicsON - Included PLS present warning.")
            message("Original message (pls):")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally = {
            message("OmicsON - PLS (pls) finished.")
        }
    )

    PLSResults$xCutoff <- xCutoff
    PLSResults$yCutoff <- yCutoff
    PLSResults$functionalGrouping$xData <- interX
    PLSResults$functionalGrouping$yData <- interY
    PLSResults$correlationCutOff$xData <- colnames(X)
    PLSResults$correlationCutOff$yData <- colnames(Y)

    PLSResults
}

# NEW PUBLIC API
makePermutationTestOnPLS <- function(XDataFrame, YDataFrame, numberOfRowsForTestOnX, numberOfRowsForTestOnY, numberOfIterations = 100, countedPLS) {
    # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    # print(numberOfRowsForTestOnX)
    # print(numberOfRowsForTestOnY)
    # print()
    # print()
    vectorOfVarianceExplained <- as.numeric();
    for (i in 1:numberOfIterations) {
        xNV <- as.character(XDataFrame[sample(nrow(XDataFrame), numberOfRowsForTestOnX), ][,1])
        yNV <- as.character(YDataFrame[sample(nrow(YDataFrame), numberOfRowsForTestOnY), ][,1])
        plsResult <- OmicsON::makePartialLeastSquaresRegression(xNamesVector = xNV, yNamesVector = yNV, XDataFrame = XDataFrame, YDataFrame = YDataFrame)
        # print("++++++++++++++++++++++++++++++++++")
        # print.default(ccaResult)
        if (is.na(plsResult) || is.null(plsResult)) {

        } else {
            vectorOfVarianceExplained <- c(vectorOfVarianceExplained, as.numeric(plsResult$varianceExplained[1]))
        }
    }

    # print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print.default(countedCCA)
    # print("*****************************")
    # print(countedCCA$xrd)
    # print(countedCCA$yrd)
    # print(vectorOfVarianceExplained)
    # print(mean(vectorOfVarianceExplained))
    # print(countedPLS)

    meanOfVarianceExplained <- NA;
    if (0 != length(vectorOfVarianceExplained)) {
        meanOfVarianceExplained <- mean(vectorOfVarianceExplained);
    }

    countedPLSRecord <- NA;
    if (is.na(countedPLS) || is.null(countedPLS)) {

    } else {
        countedPLSRecord <- countedPLS$varianceExplained[1];
    }

    testResult <- list("countedPLSOnX" = countedPLSRecord,
                       "meanOnVarianceExplained" = meanOfVarianceExplained)
    testResult
}

# NEW PUBLIC API
makePLSOnGroups <- function(groupsDefinitionDF, mappingDF, leftMappingColumnName = 'root', rightMappingColumnName = 'genesSymbolsFromEnsemble', groupsDataDF, mappingDataDF){
    ddply(.data = groupsDefinitionDF['Molecules'], .(Molecules), .fun = function(dfElement) {
        # print("???????????????????")
        # print(dfElement)
        rightSideIdsToAnalys <- unlist(strsplit(as.character(dfElement), split = " "));
        # print(rightSideIdsToAnalys)
        leftSideIdsToAnalys <- mappingDF[mappingDF[[leftMappingColumnName]] %in% rightSideIdsToAnalys,][[rightMappingColumnName]]
        leftSideIdsToAnalys <- unique(unlist(leftSideIdsToAnalys))

        plsResults <- makePartialLeastSquaresRegression(
            xNamesVector = leftSideIdsToAnalys,
            yNamesVector = rightSideIdsToAnalys,
            XDataFrame = mappingDataDF,
            YDataFrame = groupsDataDF)

        # print("######################################")
        # print(ccaResults)
        # print("---------------------------")
        # print.default(ccaResults)
        #TODO : Use user defined column name instead of symbol. ChEBI column too.
        numberOfRowsForTestOnX <- nrow(mappingDataDF[mappingDataDF$symbol %in% leftSideIdsToAnalys, ])
        numberOfRowsForTestOnY <- nrow(groupsDataDF[groupsDataDF$ChEBI %in% rightSideIdsToAnalys, ])
        parmutationTestResult <- makePermutationTestOnPLS(XDataFrame = mappingDataDF, YDataFrame = groupsDataDF,
                                                          numberOfRowsForTestOnX = numberOfRowsForTestOnX,
                                                          numberOfRowsForTestOnY = numberOfRowsForTestOnY,
                                                          numberOfIterations = 50, countedPLS = plsResults);

        dfWithCca <- data.frame('right' = I(list(rightSideIdsToAnalys)),
                                'left' = I(list(leftSideIdsToAnalys)),
                                'plsResults' = I(list(plsResults)),
                                'plsPermutationTestResults' = I(list(parmutationTestResult)))
        dfWithCca
    })
}

# NEW PUBLIC API
plotRmsepForPLS <- function(PLSResult, resetToActualMfrow = TRUE,
                            thirdLineText = "",
                            selectionVector = colnames(PLSResult$coefficients), ...) {

    actualMfrow <- par("mfrow")

    rmsep <- pls::RMSEP(PLSResult)
    selectionIndexes <- which(names(rmsep$val[1,,1]) %in% selectionVector)
    rmsep$val <- rmsep$val[, selectionIndexes, , drop = FALSE]
    # rmsep$comps <- rmsep$comps[selectionIndexes]

    if ("" == thirdLineText) {
        sub = paste("number of components = ", length(rmsep$comps), ", ",
                    "xCutoff = ", PLSResult$xCutoff, ", ",
                    "yCutoff = ", PLSResult$yCutoff, sep = "")
    } else {
        sub = paste("number of components = ", length(rmsep$comps), ", ",
                    "xCutoff = ", PLSResult$xCutoff, ", ",
                    "yCutoff = ", PLSResult$yCutoff, ", ",
                    thirdLineText, sep = "")
    }

    plot(rmsep, legendpos = "topright", xlab = sub, ...)

    if (resetToActualMfrow) {
        par(mfrow = actualMfrow)
    }
}




# NEW PUBLIC API
plotRegression <- function(PLSResult, ncompValue = NULL) {
    if (is.null(ncompValue)) {
        plot(PLSResult, asp = 1, line = TRUE)
    } else {
        if (ncompValue > PLSResult$ncomp) {
            ncompValue <- PLSResult$ncomp
        }
        plot(PLSResult, ncomp = ncompValue, asp = 1, line = TRUE)
    }
}

# TODO: Analiza różnicowa, differencial analysis.

# TODO: Check exceptions hadling in methods.

# TODO: Check plots. :)
makePLSCharts <- function(PLS) {
    png("PLS_loadings.png", width = 640, height = 480)
        par(mfrow = c(2,2))
        biplot(PLS, which = "x") # Default
        biplot(PLS, which = "y")
        biplot(PLS, which = "scores")
        biplot(PLS, which = "loadings")
    dev.off()
}


# OLD PUBLIC API
createFunctionalInteractionsDataFrameOld <- function(chebiToReactomeDataFrame, singleIdColumnName = 'ontologyId', idsListColumnName = 'ensembleIds') {
    functionalInteractionsDataFrame <- ddply(.data = chebiToReactomeDataFrame, c(singleIdColumnName), .fun = function(dfElement) {
        functionalInteractionsRows <- adply(.data = dfElement[1,c(idsListColumnName)][[1]], .margins = 1, dfff = dfff, .fun = function(listElement, dfff) {
            functionalInteractionsRow <- data.frame("Gene1" = dfElement[1, c(singleIdColumnName)],
                               "Gene2" = listElement,
                               "Annotation" = "reactome",
                               "Direction" = "-",
                               "Score" = 1.00,
                               stringsAsFactors = FALSE)
            functionalInteractionsRow
        })
        functionalInteractionsRows
    })
    functionalInteractionsDataFrame[,c("Gene1", "Gene2", "Annotation", "Direction", "Score")]
}


# NEW PUBLIC API
createFunctionalInteractionsDataFrame <- function(chebiToReactomeDataFrame, singleIdColumnName = 'ontologyId', idsListColumnName = 'ensembleIds',
                                                  attachRootColumn = TRUE, rootIdCollumnName = 'root') {
    functionalInteractionsDataFrame <- ddply(.data = chebiToReactomeDataFrame, c(singleIdColumnName), .fun = function(dfElement) {
        functionalInteractionsRows <- adply(.data = dfElement[1,c(idsListColumnName)][[1]], .margins = 1, dfff = dfff, .fun = function(listElement, dfff) {
            if (attachRootColumn) {
                functionalInteractionsRow <- data.frame("root" = dfElement[1, c(rootIdCollumnName)],
                                                        "singleIdColumnName" = dfElement[1, c(singleIdColumnName)],
                                                        "idsListColumnName" = listElement,
                                                        "Annotation" = "-",
                                                        "Direction" = "-",
                                                        "Score" = 1.00,
                                                        stringsAsFactors = FALSE)
            } else {
                functionalInteractionsRow <- data.frame("singleIdColumnName" = dfElement[1, c(singleIdColumnName)],
                                                        "idsListColumnName" = listElement,
                                                        "Annotation" = "-",
                                                        "Direction" = "-",
                                                        "Score" = 1.00,
                                                        stringsAsFactors = FALSE)
            }

            functionalInteractionsRow
        })
        functionalInteractionsRows
    })

    if (attachRootColumn) {
        functionalInteractionsDataFrame <- functionalInteractionsDataFrame[,c("root", "singleIdColumnName", "idsListColumnName", "Annotation", "Direction", "Score")]
        colnames(functionalInteractionsDataFrame)[2:3] <- c(singleIdColumnName, idsListColumnName)
    } else {
        functionalInteractionsDataFrame <- functionalInteractionsDataFrame[,c("singleIdColumnName", "idsListColumnName", "Annotation", "Direction", "Score")]
        colnames(functionalInteractionsDataFrame)[1:2] <- c(singleIdColumnName, idsListColumnName)
    }

    functionalInteractionsDataFrame
}
