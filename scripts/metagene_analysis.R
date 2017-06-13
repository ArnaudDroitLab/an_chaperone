library(metagene)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

design = read.table("input/metagene_design.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
bam.files = design$Samples

# Get coordinates for all genes.
coordinates = AnnotationDbi::select(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,
                                    keys(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, "GENEID"),
                                    c("TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND"),
                                    "GENEID")

# Get gene bodies
gene.bodies = GRanges(coordinates)                                    
                                    
# Get region flanking TSS                                    
tss.regions.pre30 = flank(GRanges(coordinates), width=30, start=TRUE, both=FALSE)
tss.regions.post300 = flank(tss.regions.pre30, width=300, start=FALSE, both=FALSE)

tss.regions.start = ifelse(strand(tss.regions.pre30)=="+", start(tss.regions.pre30), start(tss.regions.post300))
tss.regions.end = ifelse(strand(tss.regions.pre30)=="+", end(tss.regions.post300), end(tss.regions.pre30))
tss.regions = GRanges(data.frame(seqnames=seqnames(tss.regions.pre30),
                                 start=tss.regions.start,
                                 end=tss.regions.end,
                                 strand=strand(tss.regions.pre30)))

region.list=list(TSS=tss.regions, GeneBody=gene.bodies)                                 
                                 
out.path = "output/metagenes"
dir.create("output/metagenes", recursive=TRUE, showWarnings=FALSE)
metagenes = list()
for(region.name in loop.values) {
    loaded.cache.filename = file.path(out.path, paste0(region.name, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = metagene$new(regions=GRangesList(region.list[[region.name]]), bam_files=design$Samples)
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(out.path, paste0(region.name, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        metagene.obj$produce_matrices(design = design, normalization="RPM", flip_regions=TRUE, bin_count=330)
        save(metagene.obj, file=matrix.cache.filename)
    } else {
        load(matrix.cache.filename)
    }
    
    df.cache.filename = file.path(out.path, paste0(region.name, " data-frame.RData"))
    if(!file.exists(df.cache.filename)) {
        metagene.obj$produce_data_frame(stat = "bootstrap")
        save(metagene.obj, file=df.cache.filename)
    } else {
        load(df.cache.filename)
    }
    
    metagenes[[region.name]] = metagene.obj
}

                                 