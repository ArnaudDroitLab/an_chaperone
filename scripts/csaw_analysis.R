require(csaw)
require(ChIPseeker)
require(edgeR)
require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
require(org.Sc.sgd.db)
require(rtracklayer)

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

# Get all genes.
coordinates = import("input/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf")
coordinates = coordinates[coordinates$type=="gene"]

# Generate intergenic gap, and add dummy annotation so we can concatenate it to the genes.
nostrands = coordinates
strand(nostrands)="+"

intergenic = gaps(nostrands)
for(i in names(mcols(coordinates))) {
    mcols(intergenic)[[i]] = NA
}

strands(intergenic) = "*"

# Concatenate genes and intergenic regions.
full.genome = c(coordinates, intergenic)
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

    # Genewise results
    olap <- findOverlaps(full.genome, rowRanges(data[,excludeNAs]))
    tab.by.gene <- combineOverlaps(olap, results$table)
    contrast.results = full.genome
    mcols(contrast.results) = cbind(mcols(contrast.results), tab.by.gene)
    
    sig.regions = contrast.results[contrast.results$FDR <= threshold & !is.na(contrast.results$FDR)]
    
    sig.overlap = findOverlaps(sig.regions, rowRanges(data[,excludeNAs]))
    fc.mean = rep(0, length(sig.regions))
    for(i in queryHits(sig.overlap)) {
        hitIndices = subjectHits(sig.overlap)[queryHits(sig.overlap)==i]
        fc.mean[i] = mean(results$coefficients[hitIndices,2])
    }
    
    sig.regions$Mean.FC = fc.mean
    
    dir.create("output/csaw", recursive=TRUE, showWarnings=FALSE)
    write.table(as.data.frame(sig.regions), file.path("output/csaw", paste0(label, ".txt")),
                sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

for(target in c("Iws1", "Chd1")) {
    # dspt2 effect
    dspt2.subset = design.file$Target==target & design.file$Temp=="Ctl"
    compare.subsets(design.file, y, data, dspt2.subset, paste0(target, "-dspt2_vs_WT"))
    
    # spt6 effect
    spt6.subset = design.file$Target==target & design.file$Temp=="39C"
    compare.subsets(design.file, y, data, spt6.subset, paste0(target, "-spt6_vs_WT"))
    
    # Temp effect
    temp.subset = design.file$Target==target & design.file$Strain=="WT"
    compare.subsets(design.file, y, data, temp.subset, paste0(target, "-39C_vs_Ctl_in_WT"))    
}

spt2.strain.subset = design.file$Target=="Spt2"
compare.subsets(design.file, y, data, spt2.strain.subset, paste0("Spt2-spt6_vs_WT"))   

spt6.strain.subset = design.file$Target=="Spt6"
compare.subsets(design.file, y, data, spt6.strain.subset, paste0("Spt6-dspt2_vs_WT"))   

# Downstream analysis
results = list()
for(contrast in list.files("output/csaw", pattern="*.txt$"))  {
    results[[contrast]] = read.table(file.path("output/csaw", contrast), sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
}
names(results) = gsub(".txt", "", names(results))

library(VennDiagram)
library(ef.utils)
for(contrast in c("dspt2_vs_WT", "spt6_vs_WT", "39C_vs_Ctl_in_WT", "Chd1", "Iws1")) {
    # Subset the results to only keep those relevant to the comparison.
    relevant.results = results[grepl(contrast, names(results))]
    
    # Fix the names by removing the useless parts.
    fixed.names = names(relevant.results)
    fixed.names = gsub(contrast, "", fixed.names)
    fixed.names = gsub("-", "", fixed.names)
    names(relevant.results) = fixed.names
    
    # Get gene names/genomic ranges for relevant results
    gene.list=list()
    gr.list = list()
    for(i in names(relevant.results)) {
        gene.list[[i]] = relevant.results[[i]]$gene_id[!is.na(relevant.results[[i]]$gene_id)]
        
        relevant.results[[i]]$start = relevant.results[[i]]$start + 1
        relevant.results[[i]]$end = relevant.results[[i]]$end
        
        gr.list[[i]] = GRanges(relevant.results[[i]])
    }
    
    # Plot venn diagram of gene names.
    venn.diagram(gene.list, filename=file.path("output/csaw", paste0("Affected gene comparison for ", contrast, ".tiff")))
    
    # Plot venn diagram of regions.
    intersect.obj = build_intersect(GRangesList(gr.list))
    intersect_venn_plot(intersect.obj, filename=file.path("output/csaw", paste0("Affected region comparison for ", contrast, ".tiff")))
}


target = "Chd1"
temp.subset = design.file$Target==target & design.file$Strain=="WT"
design.file = design.file
dge.object = y
data = data 
indices = temp.subset
label = paste0(target, "-39C_vs_Ctl_in_WT")
threshold=0.05
