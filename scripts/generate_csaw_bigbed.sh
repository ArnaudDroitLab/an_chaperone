#module load mugqic/ucsc/v346
module load ucsc/20141124

fetchChromSizes sacCer3 > sacCer3.chrom.sizes

for i in output/csaw/*.txt
do
    sed -e 1d "$i" | sed -e 's/^/chr/' | sed -e 's/chrMito/chrM/' | cut -f 1-4 | sort -k1,1 -k2,2n  > "$i.sorted"
    bedToBigBed "$i.sorted" sacCer3.chrom.sizes "$i.bigbed"
    rm "$i.sorted"
done

