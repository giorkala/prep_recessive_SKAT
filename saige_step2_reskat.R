library(SAIGE)
library(optparse)

option_list = list(
    make_option(c("-o", "--out"), type="character", default="saige_gene", help="prefix for any output"),
    make_option("--bfile", type="character", default="chr11", help="prefix for the bed/bim/fam"),
    make_option("--vcf", type="character", default="chr11", help="filepath for the VCF file"),
    make_option("--chrom", type="character", default="chr11", help="tag for chromosome"),
    make_option("--bgen", type="character", default="chr11", help="filepath for the BGEN file"),
    make_option("--sample", type="character", help="full path to a sample file"),
    make_option("--model", type="character", help="prefix for rda and varRatio files"),
    make_option("--grm", type="character", help="full path to the GRM"),
    make_option("--annot", "-a", type="character", help="full path to file with annotation")
)

# make_option("-wdir", type="character", default="./", help="working directory"),

opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)

# print(c(opts$wdir, opts$model))

SPAGMMATtest(
    bedFile=paste(opts$bfile, '.bed', sep=''),
    bimFile=paste(opts$bfile, '.bim', sep=''),
    famFile=paste(opts$bfile, '.fam', sep='')
    bgenFile = opts$bgen,
    bgenFileIndex = paste(opts$bgen, '.bgi', sep=''),
    subSampleFile = opts$sample,
    min_MAC = 0.5,
    min_MAF = 0.0,
    GMMATmodelFile = paste(opts$model, '.rda', sep=''),
    varianceRatioFile = paste(opts$model, '.varianceRatio.txt', sep=''),
    sparseGRMFile = opts$grm,
    sparseGRMSampleIDFile = paste(opts$grm, '.sampleIDs.txt', sep=''),
    groupFile = opts$annot,
    maxMAF_in_groupTest = c(0.05),
    MACCutoff_to_CollapseUltraRare = 5,
    annotation_in_groupTest = c('synonymous', 'pLoF|pLoF|pLoF|pLoF;damaging_missense;damaging_missense|damaging_missense;pLoF|damaging_missense'),
    is_Firth_beta = TRUE,
    pCutoffforSKAT = 0.01,
    is_output_moreDetails = TRUE,
    is_output_NA_in_groupTest = TRUE,
    is_single_in_groupTest = TRUE,
    LOCO = FALSE,
    is_fastTest = TRUE,
    SAIGEOutputFile_rvTest = TRUE,
    is_write_outcome = TRUE
)
