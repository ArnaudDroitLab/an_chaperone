require(csaw)
require(ChIPseeker)
require(edgeR)
require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
require(org.Sc.sgd.db)

# Load data
bam.files = Sys.glob("output/pipeline/alignment/*/*.sorted.dup.bam")
condition = gsub(".*\\/(.*).sorted.dup.bam", "\\1", bam.files)
design.file = read.table("input/csaw_design_2.txt", header=TRUE)
design.file = design.file[match(condition, design.file$Sample),]

param <- readParam(minq=50)
data <- windowCounts(bam.files, ext=200, width=10, param=param)

# 2. Filtering out uninteresting regions.
keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]

# 3. Calculating normalization factors.
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
normfacs <- normOffsets(binned)

# 4. Identifying DB windows.
y <- asDGEList(data, norm.factors=normfacs)
txdb = makeTxDbFromGFF("input/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf")

compare.subsets <- function(design.file, dge.object, data, indices, label, threshold=0.05) {
    excludeNAs = indices
    excludeNAs[is.na(excludeNAs)] = FALSE
    
    design.subset = design.file[excludeNAs,]
    dge.subset = y[,excludeNAs]
    
    design = c()
    if(length(unique(design.subset$Temp))==2) {
        raw.factors = factor(design.subset$Temp, c("Ctl", "39C"))
        design = model.matrix(~raw.factors)
    } else {
        other.factor = setdiff(design.subset$Strain, "WT")
        raw.factors = factor(design.subset$Strain, c("WT", other.factor))
        design = model.matrix(~raw.factors)
    }
    
    dge.disp <- estimateDisp(dge.subset, design)
    fit <- glmQLFit(dge.disp, design, robust=TRUE)
    results <- glmQLFTest(fit)
    
    merged <- mergeWindows(rowRanges(data[,excludeNAs]), tol=1000L)
    tabcom <- combineTests(merged$id, results$table)
    
    sig.regions = merged$region
    mcols(sig.regions) = tabcom
    sig.regions = sig.regions[sig.regions$FDR <= threshold]
    
    sig.regions = annotatePeak(sig.regions, tssRegion=c(-1, 1),
                               TxDb=txdb,
                               annoDb="org.Sc.sgd.db", level="gene")
    
    dir.create("output/csaw", recursive=TRUE, showWarnings=FALSE)
    write.table(as.data.frame(sig.regions), file.path("output/csaw", paste0(label, ".txt")),
                sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

for(target in c("Iws1", "Chd1")) {
    # dspt2 effect
    dspt2.subset = design.file$Target==target & design.file$Temp=="Ctl"
    compare.subsets(design.file, y, data, dspt2.subset, paste0(target, " - dspt2 vs WT"))
    
    # spt6 effect
    spt6.subset = design.file$Target==target & design.file$Temp=="39C"
    compare.subsets(design.file, y, data, spt6.subset, paste0(target, " - spt6 vs WT"))
    
    # Temp effect
    temp.subset = design.file$Target==target & design.file$Strain=="WT"
    compare.subsets(design.file, y, data, temp.subset, paste0(target, " - 39C vs Ctl in WT"))    
}

spt2.strain.subset = design.file$Target=="Spt2"
compare.subsets(design.file, y, data, spt2.strain.subset, paste0("Spt2 - spt6 vs WT"))   

spt6.strain.subset = design.file$Target=="Spt6"
compare.subsets(design.file, y, data, spt2.strain.subset, paste0("Spt6 - dspt2 vs WT"))   

# Downstream analysis
results = list()
for(contrast in list.files("output/csaw", pattern="*.txt"))  {
    results[[contrast]] = read.table(file.path("output/csaw", contrast), sep="\t", header=TRUE)
}
names(results) = gsub(".txt", "", names(results))

library(VennDiagram)
library(ef.utils)
for(contrast in c("dspt2 vs WT", "spt6 vs WT", "39C vs Ctl", "Chd1", "Iws1")) {
    # Subset the results to only keep those relevant to the comparison.
    relevant.results = results[grepl(contrast, names(results))]
    
    # Fix the names by removing the useless parts.
    fixed.names = names(relevant.results)
    fixed.names = gsub(contrast, "", fixed.names)
    fixed.names = gsub("\\s+-\\s+", "", fixed.names)
    names(relevant.results) = fixed.names
    
    # Get gene names/genomic ranges for relevant results
    gene.list=list()
    gr.list = list()
    for(i in names(relevant.results)) {
        gene.list[[i]] = relevant.results[[i]]$geneId
        gr.list[[i]] = GRanges(relevant.results[[i]])
    }
    
    # Plot venn diagram of gene names.
    venn.diagram(gene.list, filename=file.path("output/csaw", paste0("Affected gene comparison for ", contrast, ".tiff")))
    
    # Plot venn diagram of regions.
    intersect.obj = build_intersect(GRangesList(gr.list))
    intersect_venn_plot(intersect.obj, filename=file.path("output/csaw", paste0("Affected region comparison for ", contrast, ".tiff")))
}
