library(optparse)

option_list = list(
    make_option("--vcfIDs", type="character", help="filepath for the list of IDs in the VCF"),
    make_option("--model", type="character", help="the rda file with the model"),
    make_option("--out", type="character", default="samples.effective.txt", help="filepath for output")
)

opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)

# load(paste(opts$model, ".rda", sep=''))
load(opts$model)
vcf = read.csv(opts$vcfIDs, header=FALSE)

final = intersect(modglm$sampleID, vcf$V1)
cat(paste("Will write IDs for", length(final), "samples to", opts$out, "\n"))
write.table(final, file=opts$out, quote=FALSE, row.names=FALSE, col.names=FALSE)
