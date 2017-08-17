library(metagene)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(ggplot2)

# Utility function to snap coordinates to chromosomes minimum/maximum
snapToChrStart <- function(regions) {
    start(regions) = pmax(start(regions), 1)
    return(regions)
}

snapToChrEnd <- function(regions, chr.lengths) {
    chr.end = chr.lengths[as.character(seqnames(regions))]
    end(regions) = pmin(end(regions), chr.end)
    return(regions)
}

snapToChrEdges <- function(regions, chr.lengths) {
    snapToChrEnd(snapToChrStart(regions), chr.lengths)
}

# Function to get region flanking TSS
tss.flank <- function(regions, before, after, chr.length, keep.smaller=FALSE) {
    # Only keep genes which have the minimum width
    kept.regions = regions
    kept.indices = TRUE
    if(!keep.smaller) {
        kept.indices = width(regions)>after
        kept.regions = regions[kept.indices] 
    }

    tss.regions.pre = flank(GRanges(kept.regions), width=before, start=TRUE, both=FALSE)
    tss.regions.post = flank(tss.regions.pre, width=after, start=FALSE, both=FALSE)
    
    tss.regions.start = ifelse(strand(tss.regions.pre)=="+", start(tss.regions.pre), start(tss.regions.post))
    tss.regions.end = ifelse(strand(tss.regions.pre)=="+", end(tss.regions.post), end(tss.regions.pre))
    tss.regions = GRanges(data.frame(seqnames=seqnames(tss.regions.pre),
                                     start=tss.regions.start,
                                     end=tss.regions.end,
                                     strand=strand(tss.regions.pre)))    
    tss.regions = snapToChrEdges(tss.regions, chr.length)
    mcols(tss.regions) = mcols(kept.regions)
    
    return(tss.regions)                                     
}
                           
# Read the design file, which contains the names of the BAM to import,
design = read.table("input/metagene_design.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Retrieve chromosome lengths so we won't go over.
chr.length = seqlengths(BSgenome.Scerevisiae.UCSC.sacCer3) - 1
names(chr.length) = gsub("^chr", "", names(chr.length))
names(chr.length) = gsub("^M", "Mito", names(chr.length))

# Get coordinates for all genes.
coordinates = import("input/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf")
coordinates = coordinates[coordinates$type=="gene"]

# Keep only tRNA
tRNA = coordinates[coordinates$gene_biotype=="tRNA"]
tRNA_30 = tRNA
length_tRNA = end(tRNA_30) - start(tRNA_30)
tRNA_offset = floor(length_tRNA * 0.3)
start(tRNA_30) = start(tRNA_30) - tRNA_offset
end(tRNA_30) = end(tRNA_30) + tRNA_offset
                           
# Extend gene bodies by 500
gene.bodies.500 = flank(coordinates, width=500, start=FALSE, both=TRUE)
gene.bodies.500 = snapToChrEdges(gene.bodies.500, chr.length)                     
                                 
# Get the length of sequences from the bam file.                                  
region.list=list(TSS_30  = tss.flank(coordinates, 30, 300, keep.smaller=FALSE, chr.length),
                 TSS_500 = tss.flank(coordinates, 500, 500, keep.smaller=TRUE, chr.length),
                 All_GeneBody=coordinates, 
                 All_GeneBody_500=gene.bodies.500)                                 

               
                 
# Read expression levels
expr.levels = read.table("input/Holstege.txt", header=TRUE, sep="\t")
expr.levels = expr.levels[!is.na(expr.levels$ExpressionLevel),]
expr.levels = expr.levels[order(expr.levels$ExpressionLevel),]

# Determine splitting points (Q1, Q2, ... Q5)
q.split = quantile(expr.levels$ExpressionLevel, seq(0,1,0.2))

# Associate expression values with their quantile.
expr.quintile = cut(expr.levels$ExpressionLevel, q.split, include.lowest=TRUE)
expr.levels$Q = expr.quintile
levels(expr.levels$Q) =  paste0("Q", 5:1)

split.by.quintile <- function(regions, expr.levels) {
    results = list()
    for(quintile in unique(expr.levels$Q)) {
        results[[quintile]] = regions[regions$gene_id %in% expr.levels$ORF[expr.levels$Q == quintile]]
    }
    
    return(GRangesList(results))
}
      
region.by.quintile = lapply(region.list, split.by.quintile, expr.levels=expr.levels)
names(region.by.quintile) = paste0(names(region.by.quintile), "_byQ")

region.list = c(region.list, tRNA=tRNA, tRNA_30=tRNA_30)
region.list = lapply(region.list, GRangesList)  
region.list = c(region.list, region.by.quintile)      
      
      
out.path = "output/metagenes"
dir.create("output/metagenes", recursive=TRUE, showWarnings=FALSE)
metagenes = list()
for(region.name in names(region.list)) {
    loaded.cache.filename = file.path(out.path, paste0(region.name, " Loaded.RData"))
    if(!file.exists(loaded.cache.filename)) {
        metagene.obj = metagene$new(regions=region.list[[region.name]], bam_files=design$Samples)
        save(metagene.obj, file=loaded.cache.filename)
    } else {
        load(loaded.cache.filename)
    }
    
    matrix.cache.filename = file.path(out.path, paste0(region.name, " matrix.RData"))
    if(!file.exists(matrix.cache.filename)) {
        bin_count = min(c(unlist(lapply(region.list[[region.name]], width)), 200))
        metagene.obj$produce_matrices(design = design, normalization="RPM", flip_regions=TRUE, bin_count=bin_count)
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

update.scale <- function(values, new.min, new.max) {
    # Make it start at 0.
    values = values - min(values)
    
    # Make it span the whole interval.
    values = values * (new.max - new.min) / (max(values) - min(values))
    
    # Yank it back to the requested minimum.
    values = values + new.min
}

# Plot
no.expr.split = names(metagenes)[!grepl("_byQ", names(metagenes))]
no.expr.split = setdiff(no.expr.split, "All_GeneBody_500")
noTAG.display = "REMOVE"
#noTAG.display = "SUBTRACT"
#noTAG.display = "SHOW"
for(region.name in no.expr.split) {
    plot.df = metagenes[[region.name]]$get_data_frame()
    
    # Fix group name
    plot.df$group = gsub("Chd1.dspt2", "Chd1_dspt2", plot.df$group)
    plot.df$group = gsub("Chd1_Spt6", "Chd1_spt6", plot.df$group)
    
    # Attach sample information.
    plot.df$Target = gsub("^(.*)_(.*)_(.*)_region_1", "\\1", plot.df$group)
    plot.df$Strain = gsub("^(.*)_(.*)_(.*)_region_1", "\\2", plot.df$group)
    plot.df$Temperature = gsub("^(.*)_(.*)_(.*)_region_1", "\\3", plot.df$group)
    plot.df$Condition = paste0(plot.df$Strain, "-", plot.df$Temperature)
    
    # Remove NoTAG
    noTAG = plot.df[plot.df$group=="NoTAG_region_1",]
    plot.df = plot.df[plot.df$group!="NoTAG_region_1",]
    
    # Either subtract NoTAG from all, or display it in all facets.
    if(noTAG.display == "SHOW") {
        noTAG$Condition="NoTAG"
        old.nrow = nrow(noTAG)
        noTAG = rbind(noTAG, noTAG, noTAG, noTAG)
        noTAG$Target = rep(unique(plot.df$Target), each=old.nrow)
        plot.df = rbind(plot.df, noTAG)
    } else if(noTAG.display == "SUBTRACT") {
        plot.df$value = plot.df$value - noTAG$value
        plot.df$qinf = plot.df$qinf - noTAG$value
        plot.df$qsup = plot.df$qsup - noTAG$value
    }
    
    # Fix positions
    if(region.name=="TSS_30") {
        # TSS should start at -30 and end at 300
        #plot.df$position = plot.df$position - min(plot.df$position) - 30
        plot.df$position = update.scale(plot.df$position, -30, 300)
        xlabel = "Distance from TSS (bp)"
    } else if(region.name=="TSS_500") {
        # TSS should start at -30 and end at 300
        plot.df$position = update.scale(plot.df$position, -500, 500)
        xlabel = "Distance from TSS (bp)"    
    } else if(region.name=="tRNA_30") {
        plot.df$position = update.scale(plot.df$position, -0.3, 1.3)
        xlabel = "Relative position between TSS (0) and TES (1)"  
    } else {
        # Gene body should start at 0 and end at 1
        plot.df$position = plot.df$position - min(plot.df$position)
        plot.df$position = plot.df$position / max(plot.df$position)
        xlabel = "Relative position between TSS (0) and TES (1)"
    }
    
    scale_values = c("dspt2-Cl"="#31a354", "spt6-39C"="#d95f0e", "WT-39C"="#3182bd", "WT-Cl"="#9ecae1", "NoTAG"="#000000")
    full.plot = ggplot(plot.df, aes(x=position, y=value, ymin=qinf, ymax=qsup, fill=Condition)) +
        geom_line(mapping=aes(color=Condition)) +
        geom_ribbon(alpha=0.6) +
        scale_fill_manual(values=scale_values) +
        scale_color_manual(values=scale_values) + 
        facet_grid(~Target) +
        labs(x=xlabel, y="RPM")
        
    if(region.name=="TSS_30" || region.name=="TSS_500") {
        # Add line at TSS.
        full.plot = full.plot + geom_vline(xintercept=0, color="black", linetype="dotted")
    }
    ggsave(filename=file.path("output/metagenes", paste0(region.name, ".pdf")), plot=full.plot, width=14)
}      

# Plot by expression quantile.
expr.split = names(metagenes)[grepl("_byQ", names(metagenes))]
#expr.split = setdiff(expr.split, "All_GeneBody_500")
noTAG.display = "REMOVE"
#noTAG.display = "SUBTRACT"
#noTAG.display = "SHOW"
for(region.name in expr.split) {
    plot.df = metagenes[[region.name]]$get_data_frame()
    
    # Fix group name
    plot.df$group = gsub("Chd1.dspt2", "Chd1_dspt2", plot.df$group)
    plot.df$group = gsub("Chd1_Spt6", "Chd1_spt6", plot.df$group)
    
    # Attach sample information.
    plot.df$Target = gsub("^(.*)_(.*)_(.*)_(Q.)", "\\1", plot.df$group)
    plot.df$Strain = gsub("^(.*)_(.*)_(.*)_(Q.)", "\\2", plot.df$group)
    plot.df$Temperature = gsub("^(.*)_(.*)_(.*)_(Q.)", "\\3", plot.df$group)
    plot.df$Quintile = gsub("^(.*)_(.*)_(.*)_(Q.)", "\\4", plot.df$group)
    plot.df$Condition = paste0(plot.df$Strain, "-", plot.df$Temperature)
    
    # Remove NoTAG
    noTAG = plot.df[grepl("NoTAG", plot.df$group),]
    plot.df = plot.df[!grepl("NoTAG", plot.df$group),]
    
    # Either subtract NoTAG from all, or display it in all facets.
    if(noTAG.display == "SHOW") {
        noTAG$Condition="NoTAG"
        old.nrow = nrow(noTAG)
        noTAG = rbind(noTAG, noTAG, noTAG, noTAG)
        noTAG$Target = rep(unique(plot.df$Target), each=old.nrow)
        plot.df = rbind(plot.df, noTAG)
    } else if(noTAG.display == "SUBTRACT") {
        plot.df$value = plot.df$value - noTAG$value
        plot.df$qinf = plot.df$qinf - noTAG$value
        plot.df$qsup = plot.df$qsup - noTAG$value
    }
    
    # Fix positions
    if(region.name=="TSS_30") {
        # TSS should start at -30 and end at 300
        #plot.df$position = plot.df$position - min(plot.df$position) - 30
        plot.df$position = update.scale(plot.df$position, -30, 300)
        xlabel = "Distance from TSS (bp)"
    } else if(region.name=="TSS_500") {
        # TSS should start at -30 and end at 300
        plot.df$position = update.scale(plot.df$position, -500, 500)
        xlabel = "Distance from TSS (bp)"    
    } else {
        # Gene body should start at 0 and end at 1
        plot.df$position = plot.df$position - min(plot.df$position)
        plot.df$position = plot.df$position / max(plot.df$position)
        xlabel = "Relative position between TSS (0) and TES (1)"
    }
    
    #scale_values = c("dspt2-Cl"="#31a354", "spt6-39C"="#d95f0e", "WT-39C"="#3182bd", "WT-Cl"="#9ecae1", "NoTAG"="#000000")
    full.plot = ggplot(plot.df, aes(x=position, y=value, ymin=qinf, ymax=qsup, fill=Quintile)) +
        geom_line(mapping=aes(color=Quintile)) +
        geom_ribbon(alpha=0.6) +
#        scale_fill_manual(values=scale_values) +
#        scale_color_manual(values=scale_values) + 
        facet_grid(Condition~Target) +
        labs(x=xlabel, y="RPM")
        
    if(region.name=="TSS_30" || region.name=="TSS_500") {
        # Add line at TSS.
        full.plot = full.plot + geom_vline(xintercept=0, color="black", linetype="dotted")
    }
    dir.create("output/metagenes/by_quantile", recursive=TRUE, showWarnings=FALSE)
    ggsave(filename=file.path("output/metagenes/by_quantile", paste0(region.name, ".pdf")), plot=full.plot, width=14)
}      