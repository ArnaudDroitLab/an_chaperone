for i in output/pipeline/alignment/*/*.sorted.dup.bam
do
    samplename=`basename $i .sorted.dup.bam`
    if [ ! -e output/pipeline/tracks/$samplename.bw ]
    then 
        mkdir -p output/pipeline/jobs
        script=output/pipeline/jobs/$samplename.make_bigwig.sh
        cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A eav-760-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.12
bamCoverage -e 200 --binSize 5 -p 16 --normalizeUsingRPKM \
    -b $i \
    -o output/pipeline/tracks/$samplename.bw
EOF
        workdir=`pwd`
        qsub $script -o $script.stdout -e $script.stderr -d $workdir
    fi
done
