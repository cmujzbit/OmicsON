
### Helpers
get_var_name <- function(variable) {
    deparse(substitute(variable))
}


### Set up package
library(OmicsON)
OmicsON::setUpReactomeMapping(ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt",
                              Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt",
                              UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")


### Read input data
pathToFileWithLipidomicsData <- system.file(package="OmicsON",
                                            "extdata", "nm-lipidomics.txt")
lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)

pathToFileWithTranscriptomicsData <- system.file(package="OmicsON",
                                                 "extdata", "nm-transcriptomics.txt")
transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)


### Convert data for functions inputs
xNamesVector <- as.character(transcriptomicsInputData$symbol)
yNamesVector <- as.character(lipidomicsInputData$ChEBI)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData


### Analysis of raw data
CcaResultsRawData <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsRawData,
    main = "Structural Correlations (Transcriptomics vs Lipidomics)",
    thirdLineText = "",
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(CcaResultsRawData), ".png"))


PlsResultsRawData <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVector,
    yNamesVector = yNamesVector,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotRmsepForPLS(
    PLSResult = PlsResultsRawData,
    nCols = 6, nRows = 3, lty = c(2), ask=FALSE,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(PlsResultsRawData), ".png"))


### Analysis of reactome deorated data

### Human
decoratedByReactomeHuman <- OmicsON::decorateByReactomeData(
    chebiMoleculesDf = lipidomicsInputData,
    chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')

ontology2EnsembleReactomeHuman <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactomeHuman,
    singleIdColumnName = 'ontologyId', idsListColumnName = 'genesSymbolsFromEnsemble')

xNamesVectorReactomeHuman <- as.character(ontology2EnsembleReactomeHuman$genesSymbolsFromEnsemble)
yNamesVectorReactomeHuman <- as.character(ontology2EnsembleReactomeHuman$root)


CcaResultsReactomeEnsembleExtentionHuman <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVectorReactomeHuman,
    yNamesVector = yNamesVectorReactomeHuman,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionHuman,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(CcaResultsReactomeEnsembleExtentionHuman), ".png"))

# TODO : Add ncompValue manipulation possibility to documentation, vignetts.
PlsResultsReactomeEnsembleExtentionHuman <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVectorReactomeHuman,
    yNamesVector = yNamesVectorReactomeHuman,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7, ncompValue = 7)

OmicsON::plotRmsepForPLS(
    PLSResult = PlsResultsReactomeEnsembleExtentionHuman,
    nCols = 6, nRows = 3, lty = c(2), ask=FALSE,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(PlsResultsReactomeEnsembleExtentionHuman), ".png"))


### Mouse
decoratedByReactomeMouse <- OmicsON::decorateByReactomeData(
    chebiMoleculesDf = lipidomicsInputData,
    chebiIdsColumnName = "ChEBI", organismTaxonomyId = '10090')

ontology2EnsembleReactomeMouse <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByReactomeMouse,
    singleIdColumnName = 'ontologyId', idsListColumnName = 'genesSymbolsFromEnsemble')

xNamesVectorReactomeMouse <- as.character(ontology2EnsembleReactomeMouse$genesSymbolsFromEnsemble)
yNamesVectorReactomeMouse <- as.character(ontology2EnsembleReactomeMouse$root)


CcaResultsReactomeEnsembleExtentionMouse <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVectorReactomeMouse,
    yNamesVector = yNamesVectorReactomeMouse,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsReactomeEnsembleExtentionMouse,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(CcaResultsReactomeEnsembleExtentionMouse), ".png"))

PlsResultsReactomeEnsembleExtentionMouse <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVectorReactomeMouse,
    yNamesVector = yNamesVectorReactomeMouse,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotRmsepForPLS(
    PLSResult = PlsResultsReactomeEnsembleExtentionMouse,
    nCols = 6, nRows = 3, lty = c(2), ask=FALSE,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(PlsResultsReactomeEnsembleExtentionMouse), ".png"))


### Analysis of string deorated data

### Human
decoratedByStringHuman <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactomeHuman,
    listOfEnsembleIdColumnName = 'ensembleIds', stringOrganismId = '9606')

ontology2EnsembleStringHuman <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByStringHuman,
    singleIdColumnName = 'ontologyId', idsListColumnName = 'stringGenesSymbolsExpand')

xNamesVectorStringHuman <- as.character(ontology2EnsembleStringHuman$stringGenesSymbolsExpand)
yNamesVectorStringHuman <- as.character(ontology2EnsembleStringHuman$root)


CcaResultsStringEnsembleExtentionHuman <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVectorStringHuman,
    yNamesVector = yNamesVectorStringHuman,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionHuman,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(CcaResultsStringEnsembleExtentionHuman), ".png"))


PlsResultsStringEnsembleExtentionHuman <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVectorStringHuman,
    yNamesVector = yNamesVectorStringHuman,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotRmsepForPLS(
    PLSResult = PlsResultsStringEnsembleExtentionHuman,
    nCols = 6, nRows = 3, lty = c(2), ask=FALSE,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(PlsResultsStringEnsembleExtentionHuman), ".png"))


### Mouse
decoratedByStringMouse <- OmicsON::decorateByStringDbData(
    chebiIdsToReactomePathways = decoratedByReactomeMouse,
    listOfEnsembleIdColumnName = 'ensembleIds', stringOrganismId = '10090')

ontology2EnsembleStringMouse <- OmicsON::createFunctionalInteractionsDataFrame(
    chebiToReactomeDataFrame = decoratedByStringMouse,
    singleIdColumnName = 'ontologyId', idsListColumnName = 'stringGenesSymbolsExpand')

xNamesVectorStringMouse <- as.character(ontology2EnsembleStringMouse$stringGenesSymbolsExpand)
yNamesVectorStringMouse <- as.character(ontology2EnsembleStringMouse$root)


CcaResultsStringEnsembleExtentionMouse <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = xNamesVectorStringMouse,
    yNamesVector = yNamesVectorStringMouse,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotCanonicalCorrelationAnalysisResults(
    ccaResults = CcaResultsStringEnsembleExtentionMouse,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(CcaResultsStringEnsembleExtentionMouse), ".png"))


PlsResultsStringEnsembleExtentionMouse <- OmicsON::makePartialLeastSquaresRegression(
    xNamesVector = xNamesVectorStringMouse,
    yNamesVector = yNamesVectorStringMouse,
    XDataFrame = XDataFrame,
    YDataFrame = YDataFrame,
    xCutoff = 0.6, yCutoff = 0.7)

OmicsON::plotRmsepForPLS(
    PLSResult = PlsResultsStringEnsembleExtentionMouse,
    nCols = 6, nRows = 3, lty = c(2), ask=FALSE,
    image_file_path = paste0("D:/projects/science/tmp/plot_",
                             get_var_name(PlsResultsStringEnsembleExtentionMouse), ".png"))



