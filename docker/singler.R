#suppressMessages(suppressWarnings(library("singleCellTK")))
#suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("optparse")))

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
  make_option(
    c('-f','--input_file'),
    help='Path to the count matrix input.'
  ),
  make_option(
    c('-o','--output_file_prefix'),
    help='The prefix for the output file'
  ),
  make_option(
    c('--level'),
    help='A string that specifies the granularity of cell typing. Choose from main or fine.'
  ),
  make_option(
    c('--reference'),
    help='A string that specifies a reference provided by SingleR. Choose from hpca, bpe, mp, dice, immgen, mouse, zeisel.'
  ),
  make_option(
    c('--featureType'),
    help='A string for whether to use gene symbols or Ensembl IDs when using a SingleR built-in reference. Choose from symbol, ensembl.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}

# Import counts as a data.frame
cnts <- read.table(
    file = opt$input_file,
    sep = "\t",
    row.names = 1,
    header=T
)

# change to a sparse matrix representation, necessary for SCE
cnts <- as(as.matrix(cnts), "sparseMatrix")

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=cnts)
)

# SingleR requires log normalized counts to match its database
sce <- runNormalization(
    inSCE = sce, 
    normalizationMethod = "logNormCounts", 
    outAssayName = "logcounts", 
    useAssay = "counts"
)

# Run SingleR
sce <- runSingleR(
    inSCE = sce,
    useAssay = "logcounts",
    useBltinRef = opt$reference,
    level = opt$level,
    featureType = opt$featureType
)

# Build varname
varName <- paste(
    "SingleR", opt$reference, paste(opt$level, "labels", sep="."),
    sep="_"
)

df.final <- data.frame(
    cell_barcodes = as.vector(colnames(sce)),
    cell_types = as.vector(sce[[varName]])
)

# Write results to file
output_filename <- paste(opt$output_file_prefix, contrast_str, 'tsv', sep='.')
write.table(
    df.final,
    output_filename,
    sep='\t', 
    quote=F, 
    row.names = TRUE
)