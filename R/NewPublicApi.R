
#' @title Set up space of mappings in package.
#' @description Function is used to set up initial space of mappings in package.
#' It is necessary to set up this space before any analysis. It is crutial for the process
#' of mapping ChEBI, Ensemble and UniProt ids to Reactome pthways ids.
#' @param ChEBI2ReactomeFileURL A string which represents
#' full URL to mapping file on Reactome page. This file include ChEBI ids
#' mapped to Reactome pathway ids.
#' @param Ensembl2ReactomeFileURL A string which represents
#' full URL to mapping file on Reactome page. This file include Ensemble ids
#' mapped to Reactome pathway ids.
#' @param UniProt2ReactomeFileURL A string which represents
#' full URL to mapping file on Reactome page. This file include UniProt ids
#' mapped to Reactome pathway ids.
#' @examples
#' OmicsON::setUpReactomeMapping(
#' ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt",
#' Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt",
#' UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")
#'
setUpReactomeMapping <- function(ChEBI2ReactomeFileURL, Ensembl2ReactomeFileURL, UniProt2ReactomeFileURL) {

    ChEBI2ReactomeFileURL <- RCurl::getURL(ChEBI2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    chEBIToReactomeLowestLevel <<- read.table(
        textConnection(ChEBI2ReactomeFileURL),
        header = FALSE
    )

    Ensembl2ReactomeFileURL <- RCurl::getURL(Ensembl2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    Ensembl2ReactomeLowestLevel <<- read.table(
        textConnection(Ensembl2ReactomeFileURL),
        header = FALSE
    )

    UniProt2ReactomeFileURL <- RCurl::getURL(UniProt2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    UniProt2ReactomeLowestLevel <<- read.table(
        textConnection(UniProt2ReactomeFileURL),
        header = FALSE
    )
}

#' @title Decorate data by data presented in Reactome database.
#' @description Use this function on you data.frame with ChEBI ids. It decorates data
#' by data presented in Reactome database. It is done by searching of ontologically related
#' molecules presented in Reactome's pathways. If ChEBI id has representation in
#' Reactome pathways then we relay on it.
#' If ChEBI doesn't occure in Reactome mapping, then we are try to find representat
#' of that group by ontology in ChEBI database. There are also ids witch are not presented and do not have any no representants in ChEBI ontology,
#' so mapping and decoration is empty for that id.
#' @param chebiMoleculesDf Data frame with column which represents ChEBI ids and set of data.
#' @param chebiIdsColumnName A string which represents a name of column with ChEBI ids.
#' @param organismTaxonomyId A string which represents an organism taxonomy id, for example human is '9606'.
#' @return Data frame with columns:
#' - root - ChEBIs ids given by user,
#' - ontologyId - ChEBI ids used in the calculation, it is taken from ChEBI ontology base on root,
#' - ensembleIds - List including vector of Ensemble's Ids,
#' - uniProtIds - List including vector of UniProt's Ids,
#' - reactomeIds - List including vector of pathway's ids from Reactome DB,
#' - genesSymbolsFromEnsemble - List including vector of gen's symbols from Reactome DB base on pathway and Ensemble's Ids,
#' - genesSymbolsFromUniProt - List including vector of gen's symbols from Reactome DB base on pathway and UniProt's Ids,
#' @examples
#' pathToFileWithLipidomicsData <- system.file(package = "OmicsON", "extdata", "nm-lipidomics.txt")
#' lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
#' lipidomicsInputDataDecoratedByReactome <- OmicsON::decorateByReactomeData(
#' chebiMoleculesDf = lipidomicsInputData,
#' chebiIdsColumnName = "ChEBI",
#' organismTaxonomyId = '9606')
#'
decorateByReactomeData <- function(chebiMoleculesDf, chebiIdsColumnName, organismTaxonomyId = '9606') {

    clusteredSmallMolecules <- OmicsON::clusterUsingOntology(
        chebiIdsDataFrame = chebiMoleculesDf,
        rootColumnName = chebiIdsColumnName,
        ontologyRepresentatnion = OmicsON::firstExistsInReactomeChebiOntology)

    mergedSmallMolecules <- OmicsON::mergeChEBIOntologyWithChildFavoring(
        clusteredSmallMolecules,
        rootColumnName = 'root')

    OmicsON::mapReactomePathwaysUnderOrganism(
        chebiOntologyIds = mergedSmallMolecules[, c("ontologyId", "root"), drop = FALSE],
        organismTaxonomyId = organismTaxonomyId,
        idsColumnName = "ontologyId",
        rootColumnName = "root")
}

#' @title Decorate data from previous step (Reactome decoration) by data presented in STRING DB database.
#' @description In this part you search for any extra interactions of gens which you
#' find in Reactome. STRING calls them neighbours. To do it just put results achieved
#' from Reactome decoration step as input to this function.
#' @param chebiIdsToReactomePathways Data frame returned from OmicsON::decorateByReactomeData().
#' @param listOfEnsembleIdColumnName A string which represents a name of column with Ensemble or UniProt ids.
#' @param stringOrganismId A string which represents an organism taxonomy id, for example human is '9606'.
#' @param stringDbVersion String with STRING DB version used to decoration.
#' @param idsColumnName Column name of ChEBI ids from ontology mapping.
#' @param rootColumnName Column name of original ChEBI ids.
#' @return Input data frame with presented extra columns:
#' - stringIds - List including vector of all STRING's ids used in computations.
#' - stringGenesSymbolsExpand - List including vector of all neighbours find in STRING database.
#' - stringGenesSymbolsNarrow - List including vector of intersection of all neighbours per id from set of ids used in search.
#' @examples
#' pathToFileWithLipidomicsData <- system.file(package = "OmicsON", "extdata", "nm-lipidomics.txt")
#' lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
#' lipidomicsInputDataDecoratedByReactome <- OmicsON::decorateByReactomeData(
#' chebiMoleculesDf = lipidomicsInputData,
#' chebiIdsColumnName = "ChEBI",
#' organismTaxonomyId = '9606')
#' decoratedByStringBaseOnEnsembleIds <- OmicsON::decorateByStringDbData(
#' chebiIdsToReactomePathways = lipidomicsInputDataDecoratedByReactome,
#' listOfEnsembleIdColumnName = 'ensembleIds')
decorateByStringDbData <- function(chebiIdsToReactomePathways, listOfEnsembleIdColumnName,
                                   stringOrganismId = '9606', stringDbVersion = "10",
                                   idsColumnName = 'ontologyId', rootColumnName = 'root') {
    stringOrganismId <- as.numeric(stringOrganismId)

    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(rootColumnName, idsColumnName)
    }

    string_db <- STRINGdb$new(
        version = stringDbVersion,
        species = stringOrganismId,
        input_directory = path.expand("~"))
    dfWithString <- ddply(.data = chebiIdsToReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
        extendedByStringAsVector <- character(length = 0)
        stringGenesSymbolsExpand <- character(length = 0)
        stringGenesSymbolsNarrow <- character(length = 0)
        if (0 != length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
            toTranslate <- data.frame("translate" = dfElement[1,listOfEnsembleIdColumnName][[1]])
            translated <- string_db$map( toTranslate, "translate", removeUnmappedRows = TRUE )
            stringGraph <- string_db$get_graph()

            extendedByString <- plyr::ddply(.data = translated, .variables = c("STRING_id"), .fun = function(r) {
                reksultDataFrame <- data.frame("res" = I(list(c(""))))
                try((function() {
                    reksultDataFrame <<- data.frame("res" = I(list(igraph::neighbors(stringGraph, r[,"STRING_id"])$name)))
                })())
                reksultDataFrame
            })

            extendedByStringAsVector <- unique(unlist(extendedByString[,"res"]))
            ensembleIdsFromStringDbExpand <- mapFromStringIdsToEnsembleIds(extendedByStringAsVector)
            stringGenesSymbolsExpand <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDbExpand,
                                                                                           organismTaxonomyId = stringOrganismId)

            interSect <- character(length = 0)
            if (1 >= length(extendedByString[,"res"])) {
                interSect <- character(length = 0)
            } else {
                interSect <- unique(extendedByString[,"res"][[1]])
                for(i in length(extendedByString[,"res"])) {
                    interSect <- intersect(interSect, extendedByString[,"res"][[i]])
                }
            }
            ensembleIdsFromStringDbNarrow <- mapFromStringIdsToEnsembleIds(interSect)
            stringGenesSymbolsNarrow <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDbNarrow,
                                                                                           organismTaxonomyId = stringOrganismId)
        }

        dffff <- data.frame(listOfEnsembleIdColumnName = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(extendedByStringAsVector))),
                            'stringGenesSymbolsExpand' = I(list(unique(stringGenesSymbolsExpand))),
                            'stringGenesSymbolsNarrow' = I(list(unique(stringGenesSymbolsNarrow))))
        dffff
    })
    names(dfWithString)[3] <- listOfEnsembleIdColumnName
    dfWithString
}
