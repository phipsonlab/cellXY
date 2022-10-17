#' Find male-female doublet cells
#'
#' This function will identify male-female doublet cells in scRNA-seq data. The
#' classifier is based on random forest models that have been trained
#' on mouse and human single cell RNA-seq data.
#'
#' @aliases findMfDoublet
#' @param x counts matrix, rows correspond to genes and columns correspond to
#' cells. Row names must be gene symbols.
#' @param genome the genome the data arises from. Current options are
#' human: genome = "Hs" or mouse: genome = "Mm".
#' @param qc logical, indicates whether to perform quality control or not.
#' qc = TRUE will predict cells that pass quality control only and the filtered
#' cells will not be classified. qc = FALSE will predict every cell except the
#' cells with zero counts on *XIST/Xist* and the sum of the Y genes.
#' Default is TRUE.
#'
#' @return a dataframe with predicted labels for each cell
#'
#' @importFrom stats predict
#' @export findMfDoublet
#'
#' @author Xinyi Jin
#'
#' @examples
#'
#' library(speckle)
#' library(SingleCellExperiment)
#' library(CellBench)
#' library(org.Hs.eg.db)
#'
#' sc_data <- load_sc_data()
#' sc_10x <- sc_data$sc_10x
#'
#' counts <- counts(sc_10x)
#' ann <- select(org.Hs.eg.db, keys=rownames(sc_10x),
#'              columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
#' m <- match(rownames(counts), ann$ENSEMBL)
#' rownames(counts) <- ann$SYMBOL[m]
#'
#' sex <- findMfDoublet(counts, genome="Hs")
#'
#' table(sex$prediction)
#' boxplot(counts["XIST",]~sex$prediction)
#'
findMfDoublet<-function(x, genome=NULL, qc = FALSE)
#    Classify cells as doublets or single-cell
#    Xinyi Jin
#    9 September 2022
{
    # Perform some checks on the data
    if (is.null(x)) stop("Counts matrix missing")
    x <- as.matrix(x)

    if (is.null(genome)){
        message("Genome not specified. Human genome used. Options are 'Hs' for
        human and 'Mm' for mouse. We currently don't support other genomes.")
    }
    # Default is Hs
    genome <- match.arg(genome,c("Hs","Mm"))

    # pre-process
    processed.data<-preprocessDb(x, genome = genome, qc = FALSE)

    # the processed transposed count matrix
    tcm <-processed.data$tcm.final

    # the normalised, transposed count matrix
    data.df <- processed.data$data.df

    # cells that filtered by QC
    # discarded.cells <- processed.data$discarded.cells

    # cells with zero count on superX and superY
    # zero.cells <- processed.data$zero.cells

    # store the final predictions
    final.pred<-data.frame(prediction=rep("NA", ncol(x)))
    row.names(final.pred)<- colnames(x)

    # load trained models
    if(genome == "Mm"){
        model <- Mm_db_model
    }
    else{
        model <- Hs_db_model
    }

    preds <- predict(model, newdata = data.df)
    final.pred[row.names(data.df), "prediction"]<- as.character(preds)

    final.pred
}

