Prefered/supported version of R is 3.5.3.

Use this command to install OmicsON package:

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
install.packages("BiocManager")
install.packages("devtools")
BiocManager::install(pkgs = c("STRINGdb", "BiocCheck", "mygene"))
devtools::install_github(repo = "cmujzbit/OmicsON", dependencies = TRUE)


Short analysis example:

OmicsON::setUpReactomeMapping(
  ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt",
  Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt",
  UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")

pathToFileWithLipidomicsData <- system.file(package = "OmicsON", "extdata", "nm-lipidomics.txt")
lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)

pathToFileWithTranscriptomicsData <- system.file(package = "OmicsON", "extdata", "nm-transcriptomics.txt")
transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)

XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData
xNamesVector <- as.character(transcriptomicsInputData$symbol)
yNamesVector <- as.character(lipidomicsInputData$ChEBI)

CcaResultsNoExtention <- OmicsON::makeCanonicalCorrelationAnalysis(
  xNamesVector = xNamesVector,
  yNamesVector = yNamesVector,
  XDataFrame = XDataFrame,
  YDataFrame = YDataFrame,
  xCutoff = 0.6, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
  ccaResults = CcaResultsNoExtention,
  main = "Structural Correlations (Transcriptomics vs Lipidomics)",
  thirdLineText = "")

PlsResultsNoExtention <- OmicsON::makePartialLeastSquaresRegression(
  xNamesVector = xNamesVector,
  yNamesVector = yNamesVector,
  XDataFrame = XDataFrame,
  YDataFrame = YDataFrame)
OmicsON::plotRmsepForPLS(
  PLSResult = PlsResultsNoExtention,
  nCols = 3, nRows = 2, lty = c(2))

lipidomicsInputDataDecoratedByReactome <- OmicsON::decorateByReactomeData(
  chebiMoleculesDf = lipidomicsInputData,
  chebiIdsColumnName = "ChEBI",
  organismTaxonomyId = '9606')

ontology2GenesSymboleFromEnsembleFunctionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(
  chebiToReactomeDataFrame = lipidomicsInputDataDecoratedByReactome,
  singleIdColumnName = 'ontologyId',
  idsListColumnName = 'genesSymbolsFromEnsemble')

xNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$genesSymbolsFromEnsemble)
yNamesVector <- as.character(ontology2GenesSymboleFromEnsembleFunctionalInteractions$root)
XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

CcaResultsReactomeEnsembleExtentionOldData <- OmicsON::makeCanonicalCorrelationAnalysis(
  xNamesVector = xNamesVector,
  yNamesVector = yNamesVector,
  XDataFrame = XDataFrame,
  YDataFrame = YDataFrame,
  xCutoff = 0.6, yCutoff = 0.7)
OmicsON::plotCanonicalCorrelationAnalysisResults(
  ccaResults = CcaResultsReactomeEnsembleExtentionOldData,
  main = "Structural Correlations (Transcriptomics vs Lipidomics)",
  thirdLineText = "Reactome")

PlsResultsReactomeEnsembleExtentionOldData <- OmicsON::makePartialLeastSquaresRegression(
  xNamesVector = xNamesVector,
  yNamesVector = yNamesVector,
  XDataFrame = XDataFrame,
  YDataFrame = YDataFrame
)
OmicsON::plotRmsepForPLS(
  PLSResult = PlsResultsReactomeEnsembleExtentionOldData,
  nCols = 3, nRows = 2, lty = c(2))
