library(metagene)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(rtracklayer)

# Read the design file, which contains the names of the BAM to import,
design = read.table("input/metagene_design.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Get coordinates for all genes.
coordinates = import("input/Saccharomyces_cerevisiae.R64-1-1.Ensembl77.gtf")
coordinates = coordinates[coordinates$type=="gene"]
coordinates = coordinates[end(coordinates)-start(coordinates)>330]                  
                                    
# Get region flanking TSS                                    
tss.regions.pre30 = flank(GRanges(coordinates), width=30, start=TRUE, both=FALSE)
tss.regions.post300 = flank(tss.regions.pre30, width=300, start=FALSE, both=FALSE)

tss.regions.start = ifelse(strand(tss.regions.pre30)=="+", start(tss.regions.pre30), start(tss.regions.post300))
tss.regions.end = ifelse(strand(tss.regions.pre30)=="+", end(tss.regions.post300), end(tss.regions.pre30))
tss.regions = GRanges(data.frame(seqnames=seqnames(tss.regions.pre30),
                                 start=tss.regions.start,
                                 end=tss.regions.end,
                                 strand=strand(tss.regions.pre30)))

# Get the length of sequences from the bam file.                                  
                                 
                                 
region.list=list(TSS=tss.regions, GeneBody=coordinates)                                 
                                 
out.path = "output/metagenes"
dir.create("output/metagenes", recursive=TRUE, showWarnings=FALSE)
metagenes = list()
for(region.name in names(region.list)) {
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

# Plot
for(region.name in names(metagenes)) {
    plot.df = metagenes[[region.name]]$get_data_frame()
    
    # Fix group name
    plot.df$group = gsub("Chd1.dspt2", "Chd1_dspt2", plot.df$group)
    plot.df$group = gsub("Chd1_Spt6", "Chd1_spt6", plot.df$group)
    
    # Attach sample information.
    plot.df$Target = gsub("^(.*)_(.*)_(.*)_region_1", "\\1", plot.df$group)
    plot.df$Strain = gsub("^(.*)_(.*)_(.*)_region_1", "\\2", plot.df$group)
    plot.df$Temperature = gsub("^(.*)_(.*)_(.*)_region_1", "\\3", plot.df$group)
    plot.df$Condition = paste0(plot.df$Strain, "-", plot.df$Temperature)
    
    # Fix positions
    if(region.name=="TSS") {
        # TSS should start at -30 and end at 300
        plot.df$position = plot.df$position - min(plot.df$position) - 30
        xlabel = "Distance from TSS (bp)"
    } else {
        # Gene body should start at 0 and end at 1
        plot.df$position = plot.df$position - min(plot.df$position)
        plot.df$position = plot.df$position / max(plot.df$position)
        xlabel = "Relative position between TSS (0) and TES (1)"
    }
    
    scale_values = c("dspt2-Cl"="#31a354", "spt6-39C"="#d95f0e", "WT-39C"="#3182bd", "WT-Cl"="#9ecae1")
    full.plot = ggplot(plot.df, aes(x=position, y=value, ymin=qinf, ymax=qsup, fill=Condition)) +
        geom_line(mapping=aes(color=Condition)) +
        geom_ribbon(alpha=0.6) +
        scale_fill_manual(values=scale_values) +
        scale_color_manual(values=scale_values) + 
        facet_grid(~Target) +
        labs(x=xlabel, y="RPM")
        
    if(region.name=="TSS") {
        # Add line at TSS.
        full.plot = full.plot + geom_vline(xintercept=0, color="black", linetype="dotted")
    }
    ggsave(filename=file.path("output/metagenes", paste0(region.name, ".pdf")), plot=full.plot)
}      
