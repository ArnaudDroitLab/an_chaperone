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

