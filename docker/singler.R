suppressMessages(suppressWarnings(library("singleCellTK")))
suppressMessages(suppressWarnings(library("Matrix")))
suppressMessages(suppressWarnings(library("tidyverse")))
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

# in case of other capitalization conventions, etc.
featureType = tolower(opt$featureType)
level = tolower(opt$level)

# change the human-readable name to the 'encoded' one:
if (opt$reference == 'Human_Primary_Cell_Atlas'){
  reference = 'hpca'
} else if(opt$reference == 'ENCODE_Blueprint'){
  reference = 'bpe'
} else if(opt$reference == 'Muraro_Pancreas'){
   reference = 'mp'
} else if(opt$reference == 'Database_of_Immune_Cell_Expression'){
  reference = 'dice'
} else if(opt$reference == 'Immunological_Genome_Project'){
  reference = 'immgen'
} else if(opt$reference == 'Mouse'){
  reference = 'mouse'
} else if(opt$reference == 'Zeisel_Mouse_Brain'){
  reference = 'zeisel'
}

# Import counts as a data.frame
cnts <- read.table(
    file = opt$input_file,
    sep = "\t",
    row.names = 1,
    header=T,
    check.names=FALSE
)

# above, we set check.names=F to prevent the mangling of the sample names.
# Now, we stash those original sample names and run make.names, so that any downstream
# functions, etc. don't run into trouble. In the end, we convert back to the original names
orig_cols = colnames(cnts)
new_colnames = make.names(orig_cols)
colnames(cnts) = new_colnames

colname_mapping = data.frame(
    orig_names = orig_cols,
    row.names=new_colnames,
    stringsAsFactors=F
)

# change to a sparse matrix representation, necessary for SCE
cnts <- Matrix(as.matrix(cnts), sparse=T)

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
    useBltinRef = reference,
    level = level,
    featureType = featureType
)

# Build varname which will be used to select the annotations
varName <- sprintf('SingleR_%s_%s_labels', reference, level)

df.final <- data.frame(
    cell_barcodes = as.vector(colnames(sce)),
    cell_types = as.vector(sce[[varName]])
)

m <- merge(df.final,colname_mapping, by.y=0, by.x='cell_barcodes')
m <- column_to_rownames(m, var='orig_names')

# Write results to file
output_filename <- paste(opt$output_file_prefix, 
  'cell_types',
  opt$reference,
  opt$level, 'tsv', sep='.')
write.table(
    m['cell_types'],
    output_filename,
    sep='\t', 
    quote=F, 
    row.names = TRUE,
)