utils::globalVariables(c("cor", "pvalue", "feature_id"))

#' shinyParse_peakGeneCor
#'
#' @import shiny
#' @importFrom FindIT2 plot_peakGeneCor
#' @importFrom stats cor
#'
#' @param mmAnnoCor the annotated GRange object from peakGeneCor
#'
#' @return shiny mode
#' @export
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     library("FindIT2")
#'     data("RNA_normCount")
#'     data("ATAC_normCount")
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'     peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)[1:100]
#'     mmAnno <- mm_geneScan(peak_GR, Txdb)
#'
#'     ATAC_colData <- data.frame(
#'         row.names = colnames(ATAC_normCount),
#'         type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
#'     )
#'
#'     ATAC_normCount_merge <- integrate_replicates(ATAC_normCount, ATAC_colData)
#'     RNA_colData <- data.frame(
#'         row.names = colnames(RNA_normCount),
#'         type = gsub("_R[0-9]", "", colnames(RNA_normCount))
#'     )
#'     RNA_normCount_merge <- integrate_replicates(RNA_normCount, RNA_colData)
#'     mmAnnoCor <- peakGeneCor(
#'         mmAnno = mmAnno,
#'         peakScoreMt = ATAC_normCount_merge,
#'         geneScoreMt = RNA_normCount_merge,
#'         parallel = FALSE
#'     )
#'
#'     # shinyParse_peakGeneCor(mmAnnoCor)
#'
#' }
shinyParse_peakGeneCor <- function(mmAnnoCor) {
    all_genes <- unique(mmAnnoCor$gene_id)

    mmAnnoCor$distanceToTSS <- as.integer(mmAnnoCor$distanceToTSS)

    ui <- shiny::fluidPage(
        fluidRow(
            column(
                2,
                suppressWarnings(selectizeInput("gene",
                                                "Gene",
                                                choices = all_genes
                )),
            ),
            column(
                2,
                numericInput("cor", "Cor_filter",
                             value = 0, min = 0, max = 1,
                             step = 0.01
                )
            ),
            column(
                2,
                numericInput("pvalue", "pvalue_filter",
                             value = 1, min = 0, max = 1,
                             step = 1e-50
                )
            )
        ),
        fluidRow(
            column(
                4,
                tableOutput("GRange"),
                tableOutput("geneScore"),
                tableOutput("peakScore")
            )
        ),
        fluidRow(
            column(
                12,
                actionButton("click", "plot"),
                actionButton("reset", "Clear"),
                plotOutput("plot_peakGene")
            )
        )
    )


    peakScoreMt <- metadata(mmAnnoCor)$peakScoreMt %>%
        tibble::as_tibble(rownames = "feature_id")
    geneScoreMt <- metadata(mmAnnoCor)$geneScoreMt %>%
        tibble::as_tibble(rownames = "gene_id")

    server <- function(input, output, session) {
        data <- reactive(subset(
            mmAnnoCor,
            gene_id %in% input$gene &
                abs(cor) > input$cor &
                pvalue < input$pvalue
        ))

        output$GRange <- renderTable(data.frame(data()))
        output$peakScore <- renderTable(
            dplyr::filter(peakScoreMt, feature_id %in% data()$feature_id)
        )
        output$geneScore <- renderTable(geneScoreMt %>%
                                            dplyr::filter(gene_id == input$gene))

        # https://shiny.rstudio.com/articles/action-buttons.html
        v <- reactiveValues(data = NULL)
        observeEvent(input$click, {
            v$data <- plot_peakGeneCor(
                mmAnnoCor = data(),
                select_gene = input$gene
            )
        })

        observeEvent(input$reset, {
            v$data <- NULL
        })


        output$plot_peakGene <- renderPlot(
            {
                if (is.null(v$data)) {
                    return()
                }
                v$data
            },
            res = 96
        )
    }

    shinyApp(ui, server)
}
