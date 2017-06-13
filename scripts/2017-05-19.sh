#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-19T13:15:16
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_trimmomatic_stats: 0 job... skipping
#   bwa_mem_picard_sort_sam: 3 jobs
#   samtools_view_filter: 3 jobs
#   picard_merge_sam_files: 3 jobs
#   picard_mark_duplicates: 4 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 25 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 26 jobs
#   macs2_callpeak: 4 jobs
#   homer_annotate_peaks: 13 jobs
#   homer_find_motifs_genome: 13 jobs
#   annotation_graphs: 1 job
#   TOTAL: 98 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/project/eav-760-aa/Chaperone/output/pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS.cc67b5051abf415f58113da9de201db9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS.cc67b5051abf415f58113da9de201db9.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt2-Myc_spt6_39C_Cl2_RS	SM:Spt2-Myc_spt6_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS/Spt2-Myc_spt6_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS.cc67b5051abf415f58113da9de201db9.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS.1ce065e1db69320cdee0a0de98d4b6f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS.1ce065e1db69320cdee0a0de98d4b6f0.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_Wt_39C_Cl1_RS	SM:Chd1-Myc_Wt_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS/Chd1-Myc_Wt_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS.1ce065e1db69320cdee0a0de98d4b6f0.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.855ad85bcd07ad9fb618f370dc076b20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.855ad85bcd07ad9fb618f370dc076b20.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Saccharomyces_cerevisiae" \
  --variable assembly="R64-1-1" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.855ad85bcd07ad9fb618f370dc076b20.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS.c75d4f864f6a5c572baf09c882138589.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS.c75d4f864f6a5c572baf09c882138589.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS/Spt2-Myc_spt6_39C_Cl2_RS.sorted.bam \
  > alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS/Spt2-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS.c75d4f864f6a5c572baf09c882138589.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS.774806b0f198b4fa0cf2cdbc18ff248c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS.774806b0f198b4fa0cf2cdbc18ff248c.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS/Chd1-Myc_Wt_39C_Cl1_RS.sorted.bam \
  > alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS/Chd1-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS.774806b0f198b4fa0cf2cdbc18ff248c.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done'
mkdir -p alignment/Chd1-Myc_spt6_39C_Cl2 && \
ln -s -f Chd1-Myc_spt6_39C_Cl2_RS/Chd1-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.merged.bam
symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done'
mkdir -p alignment/Spt2-Myc_spt6_39C_Cl2 && \
ln -s -f Spt2-Myc_spt6_39C_Cl2_RS/Spt2-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.merged.bam
symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done'
mkdir -p alignment/Chd1-Myc_Wt_39C_Cl1 && \
ln -s -f Chd1-Myc_Wt_39C_Cl1_RS/Chd1-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.merged.bam
symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2.8a20eca04f89ba5e6ef1829bb3ebe434.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2.8a20eca04f89ba5e6ef1829bb3ebe434.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.merged.bam \
  OUTPUT=alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2.8a20eca04f89ba5e6ef1829bb3ebe434.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2.554ea16c20e446d4dfbc6ef2e45a367d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2.554ea16c20e446d4dfbc6ef2e45a367d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.merged.bam \
  OUTPUT=alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2.554ea16c20e446d4dfbc6ef2e45a367d.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1.7b0365eb377a9136f822d0788ab60ac5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1.7b0365eb377a9136f822d0788ab60ac5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.merged.bam \
  OUTPUT=alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1.7b0365eb377a9136f822d0788ab60ac5.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done'
mkdir -p report && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.7656bc386a3d9a3d5154eb93d8240d3d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.7656bc386a3d9a3d5154eb93d8240d3d.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam \
  > alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam \
  > alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam \
  > alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam \
  > alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam \
  > alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam \
  > alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  > alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  > alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  > alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  > alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam \
  > alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam \
  > alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  > alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  > alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  > alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  > alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam \
  > alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam \
  > alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam \
  > alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam \
  > alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam \
  > alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam \
  > alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam \
  > alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam \
  > alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  > alignment/No-TAG/No-TAG.sorted.dup.bam.flagstat
metrics.flagstat.7656bc386a3d9a3d5154eb93d8240d3d.mugqic.done
)
metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$metrics_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.07e87c51a237d303e4571e2d46c294cf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.07e87c51a237d303e4571e2d46c294cf.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in Chd1-Myc_dspt2Cl2 Chd1-Myc_dspt2Cl1 Iws1-Myc_dspt2Cl2 Iws1-Myc_dspt2Cl1 Spt6-Myc_dspt2Cl2 Spt6-Myc_dspt2Cl1 Chd1-Myc_spt6_39C_Cl2 Chd1-Myc_spt6_39C_Cl1 Iws1-Myc_spt6_39C_Cl2 Iws1-Myc_spt6_39C_Cl1 Spt2-Myc_spt6_39C_Cl2 Spt2-Myc_spt6_39C_Cl1 Chd1-Myc_Wt_39C_Cl2 Chd1-Myc_Wt_39C_Cl1 Iws1-Myc_Wt_39C_Cl2 Iws1-Myc_Wt_39C_Cl1 Spt2-Myc_Wt_39C_Cl2 Spt2-Myc_Wt_39C_Cl1 Chd1-Myc_WtCl2 Chd1-Myc_WtCl1 Iws1-Myc_WtCl2 Iws1-Myc_WtCl1 Spt6-Myc_WtCl2 Spt6-Myc_WtCl1 No-TAG
do
  flagstat_file=alignment/$sample/$sample.sorted.dup.bam.flagstat
  echo -e "$sample	`grep -P '^\d+ \+ \d+ mapped' $flagstat_file | grep -Po '^\d+'`	`grep -P '^\d+ \+ \d+ duplicate' $flagstat_file | grep -Po '^\d+'`"
done | \
awk -F"	" '{OFS="	"; print $0, $3 / $2 * 100}' | sed '1iSample	Aligned Filtered Reads	Duplicate Reads	Duplicate %' \
  > metrics/SampleMetrics.stats && \
mkdir -p report && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F "	" 'FNR==NR{trim_line[$1]=$0; surviving[$1]=$3; next}{OFS="	"; if ($1=="Sample") {print trim_line[$1], $2, "Aligned Filtered %", $3, $4} else {print trim_line[$1], $2, $2 / surviving[$1] * 100, $3, $4}}' metrics/trimSampleTable.tsv metrics/SampleMetrics.stats \
  > report/trimMemSampleTable.tsv
else
  cp metrics/SampleMetrics.stats report/trimMemSampleTable.tsv
fi && \
trim_mem_sample_table=`if [[ -f metrics/trimSampleTable.tsv ]] ; then LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8)}}' report/trimMemSampleTable.tsv ; else LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4)}}' report/trimMemSampleTable.tsv ; fi` && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.07e87c51a237d303e4571e2d46c294cf.mugqic.done
)
metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.Chd1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_dspt2Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_dspt2Cl2.2e99a73d76db07ad17ffa0897f140d43.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_dspt2Cl2.2e99a73d76db07ad17ffa0897f140d43.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_dspt2Cl2 \
  alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_dspt2Cl2.2e99a73d76db07ad17ffa0897f140d43.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.Chd1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_dspt2Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_dspt2Cl1.c77b483c9af6a84fe5a40e56904f7c9b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_dspt2Cl1.c77b483c9af6a84fe5a40e56904f7c9b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_dspt2Cl1 \
  alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_dspt2Cl1.c77b483c9af6a84fe5a40e56904f7c9b.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.Iws1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_dspt2Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_dspt2Cl2.c79c31a179c6a09456435ac0edf304c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_dspt2Cl2.c79c31a179c6a09456435ac0edf304c4.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_dspt2Cl2 \
  alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_dspt2Cl2.c79c31a179c6a09456435ac0edf304c4.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.Iws1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_dspt2Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_dspt2Cl1.86fd54d3f2c4bade43e0d45c352b17fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_dspt2Cl1.86fd54d3f2c4bade43e0d45c352b17fa.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_dspt2Cl1 \
  alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_dspt2Cl1.86fd54d3f2c4bade43e0d45c352b17fa.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.Spt6-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_dspt2Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt6-Myc_dspt2Cl2.44790af0c8301d016e1c167667c1152f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt6-Myc_dspt2Cl2.44790af0c8301d016e1c167667c1152f.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt6-Myc_dspt2Cl2 \
  alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt6-Myc_dspt2Cl2.44790af0c8301d016e1c167667c1152f.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.Spt6-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_dspt2Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt6-Myc_dspt2Cl1.9d22f633a6041cbc9e36358ef72d02db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt6-Myc_dspt2Cl1.9d22f633a6041cbc9e36358ef72d02db.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt6-Myc_dspt2Cl1 \
  alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt6-Myc_dspt2Cl1.9d22f633a6041cbc9e36358ef72d02db.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2.278e4d2973ecaf0200b222f6dd28f86e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2.278e4d2973ecaf0200b222f6dd28f86e.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_spt6_39C_Cl2 \
  alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2.278e4d2973ecaf0200b222f6dd28f86e.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl1.8279da5caedea2c61f03932e352d4deb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl1.8279da5caedea2c61f03932e352d4deb.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_spt6_39C_Cl1 \
  alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl1.8279da5caedea2c61f03932e352d4deb.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2.046b8d2d4081bb2625980e96ca0d67a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2.046b8d2d4081bb2625980e96ca0d67a9.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_spt6_39C_Cl2 \
  alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2.046b8d2d4081bb2625980e96ca0d67a9.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1.4de020fb2138147a9a457c50668d9d17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1.4de020fb2138147a9a457c50668d9d17.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_spt6_39C_Cl1 \
  alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1.4de020fb2138147a9a457c50668d9d17.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2.3dd70bfb20580c5ab85becb9583bd660.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2.3dd70bfb20580c5ab85becb9583bd660.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt2-Myc_spt6_39C_Cl2 \
  alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2.3dd70bfb20580c5ab85becb9583bd660.mugqic.done
)
homer_make_tag_directory_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl1.8af3df9ee391b8eda7d2919a9ece24a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl1.8af3df9ee391b8eda7d2919a9ece24a4.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt2-Myc_spt6_39C_Cl1 \
  alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl1.8af3df9ee391b8eda7d2919a9ece24a4.mugqic.done
)
homer_make_tag_directory_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_13_JOB_ID: homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2.2330a611835e51e6472efe3c87e6d568.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2.2330a611835e51e6472efe3c87e6d568.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_Wt_39C_Cl2 \
  alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2.2330a611835e51e6472efe3c87e6d568.mugqic.done
)
homer_make_tag_directory_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_14_JOB_ID: homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1.48b24de06236c4d6b3d7eb25b1e4d452.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1.48b24de06236c4d6b3d7eb25b1e4d452.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_Wt_39C_Cl1 \
  alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1.48b24de06236c4d6b3d7eb25b1e4d452.mugqic.done
)
homer_make_tag_directory_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_15_JOB_ID: homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl2.aa8eea7ef64f377062759c37c5fcde9a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl2.aa8eea7ef64f377062759c37c5fcde9a.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_Wt_39C_Cl2 \
  alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl2.aa8eea7ef64f377062759c37c5fcde9a.mugqic.done
)
homer_make_tag_directory_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_16_JOB_ID: homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1.a4c610c2f1a50657dbb2e61f3c06821e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1.a4c610c2f1a50657dbb2e61f3c06821e.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_Wt_39C_Cl1 \
  alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1.a4c610c2f1a50657dbb2e61f3c06821e.mugqic.done
)
homer_make_tag_directory_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_17_JOB_ID: homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2.99280c7dcfc75a867bd0593d690ebb5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2.99280c7dcfc75a867bd0593d690ebb5b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt2-Myc_Wt_39C_Cl2 \
  alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2.99280c7dcfc75a867bd0593d690ebb5b.mugqic.done
)
homer_make_tag_directory_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_18_JOB_ID: homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1.27e8ce39860ab3aa4b2b2d83416aa548.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1.27e8ce39860ab3aa4b2b2d83416aa548.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt2-Myc_Wt_39C_Cl1 \
  alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1.27e8ce39860ab3aa4b2b2d83416aa548.mugqic.done
)
homer_make_tag_directory_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_19_JOB_ID: homer_make_tag_directory.Chd1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_WtCl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_WtCl2.510de56cda48c48262ddb4080891c53c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_WtCl2.510de56cda48c48262ddb4080891c53c.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_WtCl2 \
  alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_WtCl2.510de56cda48c48262ddb4080891c53c.mugqic.done
)
homer_make_tag_directory_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_20_JOB_ID: homer_make_tag_directory.Chd1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_WtCl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Chd1-Myc_WtCl1.fa212d45b5cf69370087fd5a8d13f616.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Chd1-Myc_WtCl1.fa212d45b5cf69370087fd5a8d13f616.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Chd1-Myc_WtCl1 \
  alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Chd1-Myc_WtCl1.fa212d45b5cf69370087fd5a8d13f616.mugqic.done
)
homer_make_tag_directory_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_21_JOB_ID: homer_make_tag_directory.Iws1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_WtCl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_WtCl2.c98b8311ae074d4564f0250a89b166a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_WtCl2.c98b8311ae074d4564f0250a89b166a5.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_WtCl2 \
  alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_WtCl2.c98b8311ae074d4564f0250a89b166a5.mugqic.done
)
homer_make_tag_directory_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_22_JOB_ID: homer_make_tag_directory.Iws1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_WtCl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Iws1-Myc_WtCl1.c20dfd339f7768040e7619de783024e1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Iws1-Myc_WtCl1.c20dfd339f7768040e7619de783024e1.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Iws1-Myc_WtCl1 \
  alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Iws1-Myc_WtCl1.c20dfd339f7768040e7619de783024e1.mugqic.done
)
homer_make_tag_directory_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_23_JOB_ID: homer_make_tag_directory.Spt6-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_WtCl2
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt6-Myc_WtCl2.fb66390952497a22434af1e44168e964.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt6-Myc_WtCl2.fb66390952497a22434af1e44168e964.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt6-Myc_WtCl2 \
  alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt6-Myc_WtCl2.fb66390952497a22434af1e44168e964.mugqic.done
)
homer_make_tag_directory_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_24_JOB_ID: homer_make_tag_directory.Spt6-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_WtCl1
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.Spt6-Myc_WtCl1.c4e6a14663beca521bd43e5193ac8d9b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.Spt6-Myc_WtCl1.c4e6a14663beca521bd43e5193ac8d9b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/Spt6-Myc_WtCl1 \
  alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.Spt6-Myc_WtCl1.c4e6a14663beca521bd43e5193ac8d9b.mugqic.done
)
homer_make_tag_directory_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_25_JOB_ID: homer_make_tag_directory.No-TAG
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.No-TAG
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.No-TAG.446ceb9460267778154750ec6fd2c254.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.No-TAG.446ceb9460267778154750ec6fd2c254.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/No-TAG \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  -checkGC -genome R64-1-1
homer_make_tag_directory.No-TAG.446ceb9460267778154750ec6fd2c254.mugqic.done
)
homer_make_tag_directory_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID:$homer_make_tag_directory_11_JOB_ID:$homer_make_tag_directory_12_JOB_ID:$homer_make_tag_directory_13_JOB_ID:$homer_make_tag_directory_14_JOB_ID:$homer_make_tag_directory_15_JOB_ID:$homer_make_tag_directory_16_JOB_ID:$homer_make_tag_directory_17_JOB_ID:$homer_make_tag_directory_18_JOB_ID:$homer_make_tag_directory_19_JOB_ID:$homer_make_tag_directory_20_JOB_ID:$homer_make_tag_directory_21_JOB_ID:$homer_make_tag_directory_22_JOB_ID:$homer_make_tag_directory_23_JOB_ID:$homer_make_tag_directory_24_JOB_ID:$homer_make_tag_directory_25_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.cbfc6d840267524918ddc1141ec9591c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.cbfc6d840267524918ddc1141ec9591c.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/design.txt \
  /gs/project/eav-760-aa/Chaperone/output/pipeline && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in Chd1-Myc_dspt2Cl2 Chd1-Myc_dspt2Cl1 Iws1-Myc_dspt2Cl2 Iws1-Myc_dspt2Cl1 Spt6-Myc_dspt2Cl2 Spt6-Myc_dspt2Cl1 Chd1-Myc_spt6_39C_Cl2 Chd1-Myc_spt6_39C_Cl1 Iws1-Myc_spt6_39C_Cl2 Iws1-Myc_spt6_39C_Cl1 Spt2-Myc_spt6_39C_Cl2 Spt2-Myc_spt6_39C_Cl1 Chd1-Myc_Wt_39C_Cl2 Chd1-Myc_Wt_39C_Cl1 Iws1-Myc_Wt_39C_Cl2 Iws1-Myc_Wt_39C_Cl1 Spt2-Myc_Wt_39C_Cl2 Spt2-Myc_Wt_39C_Cl1 Chd1-Myc_WtCl2 Chd1-Myc_WtCl1 Iws1-Myc_WtCl2 Iws1-Myc_WtCl1 Spt6-Myc_WtCl2 Spt6-Myc_WtCl1 No-TAG
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.cbfc6d840267524918ddc1141ec9591c.mugqic.done
)
qc_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.Chd1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_dspt2Cl2.e7943cabf7cf0a24586c3b26644dd7f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_dspt2Cl2.e7943cabf7cf0a24586c3b26644dd7f9.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_dspt2Cl2 && \
makeUCSCfile \
  tags/Chd1-Myc_dspt2Cl2 | \
gzip -1 -c > tracks/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_dspt2Cl2.e7943cabf7cf0a24586c3b26644dd7f9.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.Chd1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_dspt2Cl1.4dba69180ffed21a2e6584de75f24933.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_dspt2Cl1.4dba69180ffed21a2e6584de75f24933.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_dspt2Cl1 && \
makeUCSCfile \
  tags/Chd1-Myc_dspt2Cl1 | \
gzip -1 -c > tracks/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_dspt2Cl1.4dba69180ffed21a2e6584de75f24933.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.Iws1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_dspt2Cl2.8c68c9a7deef3e919fefaf6e83df9203.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_dspt2Cl2.8c68c9a7deef3e919fefaf6e83df9203.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_dspt2Cl2 && \
makeUCSCfile \
  tags/Iws1-Myc_dspt2Cl2 | \
gzip -1 -c > tracks/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_dspt2Cl2.8c68c9a7deef3e919fefaf6e83df9203.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.Iws1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_dspt2Cl1.3fe13903d1c2311bacd4054011ff4b26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_dspt2Cl1.3fe13903d1c2311bacd4054011ff4b26.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_dspt2Cl1 && \
makeUCSCfile \
  tags/Iws1-Myc_dspt2Cl1 | \
gzip -1 -c > tracks/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_dspt2Cl1.3fe13903d1c2311bacd4054011ff4b26.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.Spt6-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt6-Myc_dspt2Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt6-Myc_dspt2Cl2.e64a4736445b00aa5d8b6524104a7ebc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt6-Myc_dspt2Cl2.e64a4736445b00aa5d8b6524104a7ebc.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt6-Myc_dspt2Cl2 && \
makeUCSCfile \
  tags/Spt6-Myc_dspt2Cl2 | \
gzip -1 -c > tracks/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt6-Myc_dspt2Cl2.e64a4736445b00aa5d8b6524104a7ebc.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.Spt6-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt6-Myc_dspt2Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt6-Myc_dspt2Cl1.cb5dcf53263d6f43003253c066c82037.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt6-Myc_dspt2Cl1.cb5dcf53263d6f43003253c066c82037.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt6-Myc_dspt2Cl1 && \
makeUCSCfile \
  tags/Spt6-Myc_dspt2Cl1 | \
gzip -1 -c > tracks/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt6-Myc_dspt2Cl1.cb5dcf53263d6f43003253c066c82037.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl2.da77cf7bf3829137861c3b85c22132d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl2.da77cf7bf3829137861c3b85c22132d1.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_spt6_39C_Cl2 && \
makeUCSCfile \
  tags/Chd1-Myc_spt6_39C_Cl2 | \
gzip -1 -c > tracks/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl2.da77cf7bf3829137861c3b85c22132d1.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl1.9f94e36d17afb7e258877a21beb5c0b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl1.9f94e36d17afb7e258877a21beb5c0b0.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_spt6_39C_Cl1 && \
makeUCSCfile \
  tags/Chd1-Myc_spt6_39C_Cl1 | \
gzip -1 -c > tracks/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_spt6_39C_Cl1.9f94e36d17afb7e258877a21beb5c0b0.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl2.5da2df694524e40f16aa72ac90797122.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl2.5da2df694524e40f16aa72ac90797122.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_spt6_39C_Cl2 && \
makeUCSCfile \
  tags/Iws1-Myc_spt6_39C_Cl2 | \
gzip -1 -c > tracks/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl2.5da2df694524e40f16aa72ac90797122.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl1.bd699333b8ad84020ef91b6852aef340.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl1.bd699333b8ad84020ef91b6852aef340.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_spt6_39C_Cl1 && \
makeUCSCfile \
  tags/Iws1-Myc_spt6_39C_Cl1 | \
gzip -1 -c > tracks/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_spt6_39C_Cl1.bd699333b8ad84020ef91b6852aef340.mugqic.done
)
homer_make_ucsc_file_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl2.c402b58a0ef957aa813a79281675af00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl2.c402b58a0ef957aa813a79281675af00.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt2-Myc_spt6_39C_Cl2 && \
makeUCSCfile \
  tags/Spt2-Myc_spt6_39C_Cl2 | \
gzip -1 -c > tracks/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl2.c402b58a0ef957aa813a79281675af00.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_12_JOB_ID: homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_12_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl1.19557041754dd42af62c8fd648af5b28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl1.19557041754dd42af62c8fd648af5b28.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt2-Myc_spt6_39C_Cl1 && \
makeUCSCfile \
  tags/Spt2-Myc_spt6_39C_Cl1 | \
gzip -1 -c > tracks/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt2-Myc_spt6_39C_Cl1.19557041754dd42af62c8fd648af5b28.mugqic.done
)
homer_make_ucsc_file_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_13_JOB_ID: homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_13_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl2.104785d296e16d32881b01200e83afa3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl2.104785d296e16d32881b01200e83afa3.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_Wt_39C_Cl2 && \
makeUCSCfile \
  tags/Chd1-Myc_Wt_39C_Cl2 | \
gzip -1 -c > tracks/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl2.104785d296e16d32881b01200e83afa3.mugqic.done
)
homer_make_ucsc_file_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_14_JOB_ID: homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_14_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl1.33ce008256e37b7b8d34aded252dce73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl1.33ce008256e37b7b8d34aded252dce73.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_Wt_39C_Cl1 && \
makeUCSCfile \
  tags/Chd1-Myc_Wt_39C_Cl1 | \
gzip -1 -c > tracks/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_Wt_39C_Cl1.33ce008256e37b7b8d34aded252dce73.mugqic.done
)
homer_make_ucsc_file_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_15_JOB_ID: homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_15_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl2.6b2d21e2ce954d9896a7163cbb628171.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl2.6b2d21e2ce954d9896a7163cbb628171.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_Wt_39C_Cl2 && \
makeUCSCfile \
  tags/Iws1-Myc_Wt_39C_Cl2 | \
gzip -1 -c > tracks/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl2.6b2d21e2ce954d9896a7163cbb628171.mugqic.done
)
homer_make_ucsc_file_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_16_JOB_ID: homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_16_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl1.abb924cb94c7570094f6f01be0bcf163.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl1.abb924cb94c7570094f6f01be0bcf163.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_Wt_39C_Cl1 && \
makeUCSCfile \
  tags/Iws1-Myc_Wt_39C_Cl1 | \
gzip -1 -c > tracks/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_Wt_39C_Cl1.abb924cb94c7570094f6f01be0bcf163.mugqic.done
)
homer_make_ucsc_file_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_17_JOB_ID: homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$homer_make_tag_directory_17_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl2.40f006a529a25be67dcf23ac11529de4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl2.40f006a529a25be67dcf23ac11529de4.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt2-Myc_Wt_39C_Cl2 && \
makeUCSCfile \
  tags/Spt2-Myc_Wt_39C_Cl2 | \
gzip -1 -c > tracks/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl2.40f006a529a25be67dcf23ac11529de4.mugqic.done
)
homer_make_ucsc_file_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_18_JOB_ID: homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$homer_make_tag_directory_18_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl1.721253e1068741616478c8f1dccda6d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl1.721253e1068741616478c8f1dccda6d8.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt2-Myc_Wt_39C_Cl1 && \
makeUCSCfile \
  tags/Spt2-Myc_Wt_39C_Cl1 | \
gzip -1 -c > tracks/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt2-Myc_Wt_39C_Cl1.721253e1068741616478c8f1dccda6d8.mugqic.done
)
homer_make_ucsc_file_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_19_JOB_ID: homer_make_ucsc_file.Chd1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_WtCl2
JOB_DEPENDENCIES=$homer_make_tag_directory_19_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_WtCl2.a9bccfac57e61acbb0d8ee3f457e8ec0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_WtCl2.a9bccfac57e61acbb0d8ee3f457e8ec0.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_WtCl2 && \
makeUCSCfile \
  tags/Chd1-Myc_WtCl2 | \
gzip -1 -c > tracks/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_WtCl2.a9bccfac57e61acbb0d8ee3f457e8ec0.mugqic.done
)
homer_make_ucsc_file_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_20_JOB_ID: homer_make_ucsc_file.Chd1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Chd1-Myc_WtCl1
JOB_DEPENDENCIES=$homer_make_tag_directory_20_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Chd1-Myc_WtCl1.585a113af61c7673f22675f0dc351e28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Chd1-Myc_WtCl1.585a113af61c7673f22675f0dc351e28.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Chd1-Myc_WtCl1 && \
makeUCSCfile \
  tags/Chd1-Myc_WtCl1 | \
gzip -1 -c > tracks/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Chd1-Myc_WtCl1.585a113af61c7673f22675f0dc351e28.mugqic.done
)
homer_make_ucsc_file_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_21_JOB_ID: homer_make_ucsc_file.Iws1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_WtCl2
JOB_DEPENDENCIES=$homer_make_tag_directory_21_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_WtCl2.2c2044e4bfd550971b4d92423c11dcd8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_WtCl2.2c2044e4bfd550971b4d92423c11dcd8.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_WtCl2 && \
makeUCSCfile \
  tags/Iws1-Myc_WtCl2 | \
gzip -1 -c > tracks/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_WtCl2.2c2044e4bfd550971b4d92423c11dcd8.mugqic.done
)
homer_make_ucsc_file_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_22_JOB_ID: homer_make_ucsc_file.Iws1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Iws1-Myc_WtCl1
JOB_DEPENDENCIES=$homer_make_tag_directory_22_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Iws1-Myc_WtCl1.2766dfbb9abb474a919fe077b632efd3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Iws1-Myc_WtCl1.2766dfbb9abb474a919fe077b632efd3.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Iws1-Myc_WtCl1 && \
makeUCSCfile \
  tags/Iws1-Myc_WtCl1 | \
gzip -1 -c > tracks/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Iws1-Myc_WtCl1.2766dfbb9abb474a919fe077b632efd3.mugqic.done
)
homer_make_ucsc_file_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_23_JOB_ID: homer_make_ucsc_file.Spt6-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt6-Myc_WtCl2
JOB_DEPENDENCIES=$homer_make_tag_directory_23_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt6-Myc_WtCl2.b88b82a908a7372b3205e030cf05c5ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt6-Myc_WtCl2.b88b82a908a7372b3205e030cf05c5ac.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt6-Myc_WtCl2 && \
makeUCSCfile \
  tags/Spt6-Myc_WtCl2 | \
gzip -1 -c > tracks/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt6-Myc_WtCl2.b88b82a908a7372b3205e030cf05c5ac.mugqic.done
)
homer_make_ucsc_file_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_24_JOB_ID: homer_make_ucsc_file.Spt6-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.Spt6-Myc_WtCl1
JOB_DEPENDENCIES=$homer_make_tag_directory_24_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.Spt6-Myc_WtCl1.5a2271ebad94eb11c6ec98723bff4ae6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.Spt6-Myc_WtCl1.5a2271ebad94eb11c6ec98723bff4ae6.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/Spt6-Myc_WtCl1 && \
makeUCSCfile \
  tags/Spt6-Myc_WtCl1 | \
gzip -1 -c > tracks/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.ucsc.bedGraph.gz
homer_make_ucsc_file.Spt6-Myc_WtCl1.5a2271ebad94eb11c6ec98723bff4ae6.mugqic.done
)
homer_make_ucsc_file_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_25_JOB_ID: homer_make_ucsc_file.No-TAG
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.No-TAG
JOB_DEPENDENCIES=$homer_make_tag_directory_25_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.No-TAG.a2c13dc03d8be6706050ad89f00d9059.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.No-TAG.a2c13dc03d8be6706050ad89f00d9059.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/No-TAG && \
makeUCSCfile \
  tags/No-TAG | \
gzip -1 -c > tracks/No-TAG/No-TAG.ucsc.bedGraph.gz
homer_make_ucsc_file.No-TAG.a2c13dc03d8be6706050ad89f00d9059.mugqic.done
)
homer_make_ucsc_file_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_26_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_25_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-spt6-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-spt6-39C.8c641a6ac0a8d8cf5d29050c89fe1f99.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-spt6-39C.8c641a6ac0a8d8cf5d29050c89fe1f99.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-spt6-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-spt6-39C/Chd1-spt6-39C \
  >& peak_call/Chd1-spt6-39C/Chd1-spt6-39C.diag.macs.out
macs2_callpeak.Chd1-spt6-39C.8c641a6ac0a8d8cf5d29050c89fe1f99.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt2-spt6-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt2-spt6-39C.9dc2cc4bcf20672e42f76ca464288cb5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt2-spt6-39C.9dc2cc4bcf20672e42f76ca464288cb5.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt2-spt6-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt2-spt6-39C/Spt2-spt6-39C \
  >& peak_call/Spt2-spt6-39C/Spt2-spt6-39C.diag.macs.out
macs2_callpeak.Spt2-spt6-39C.9dc2cc4bcf20672e42f76ca464288cb5.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.Chd1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-WT-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-WT-39C.b895d94e31253f62b5998903309db099.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-WT-39C.b895d94e31253f62b5998903309db099.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-WT-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-WT-39C/Chd1-WT-39C \
  >& peak_call/Chd1-WT-39C/Chd1-WT-39C.diag.macs.out
macs2_callpeak.Chd1-WT-39C.b895d94e31253f62b5998903309db099.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.dcbd0355ed87aef14b6caf32d9da8654.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.dcbd0355ed87aef14b6caf32d9da8654.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in Chd1-dspt2-Ctl Iws1-dspt2-Ctl Spt6-dspt2-Ctl Chd1-spt6-39C Iws1-spt6-39C Spt2-spt6-39C Chd1-WT-39C Iws1-WT-39C Spt2-WT-39C Chd1-WT-Ctl Iws1-WT-Ctl Spt6-WT-Ctl
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.dcbd0355ed87aef14b6caf32d9da8654.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.Chd1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Chd1-dspt2-Ctl.da6038c9474c56701367d556241956cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Chd1-dspt2-Ctl.da6038c9474c56701367d556241956cd.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl && \
annotatePeaks.pl \
  peak_call/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl \
  -genomeOntology annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl \
  > annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl.annotated.csv",
  "annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Chd1-dspt2-Ctl.da6038c9474c56701367d556241956cd.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.Iws1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Iws1-dspt2-Ctl.3e68bd5699de3276abf8f90b1758580a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Iws1-dspt2-Ctl.3e68bd5699de3276abf8f90b1758580a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl && \
annotatePeaks.pl \
  peak_call/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl \
  -genomeOntology annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl \
  > annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl.annotated.csv",
  "annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Iws1-dspt2-Ctl.3e68bd5699de3276abf8f90b1758580a.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.Spt6-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt6-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Spt6-dspt2-Ctl.4064c7922e04b342e6950c3c9e890307.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Spt6-dspt2-Ctl.4064c7922e04b342e6950c3c9e890307.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl && \
annotatePeaks.pl \
  peak_call/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl \
  -genomeOntology annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl \
  > annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl.annotated.csv",
  "annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Spt6-dspt2-Ctl.4064c7922e04b342e6950c3c9e890307.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Chd1-spt6-39C.f35289d95f9088f71af4e38a049a8d17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Chd1-spt6-39C.f35289d95f9088f71af4e38a049a8d17.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Chd1-spt6-39C/Chd1-spt6-39C && \
annotatePeaks.pl \
  peak_call/Chd1-spt6-39C/Chd1-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Chd1-spt6-39C/Chd1-spt6-39C \
  -genomeOntology annotation/Chd1-spt6-39C/Chd1-spt6-39C \
  > annotation/Chd1-spt6-39C/Chd1-spt6-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Chd1-spt6-39C/Chd1-spt6-39C.annotated.csv",
  "annotation/Chd1-spt6-39C/Chd1-spt6-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Chd1-spt6-39C.f35289d95f9088f71af4e38a049a8d17.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.Iws1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-spt6-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Iws1-spt6-39C.5712bac30da1af8bd65a81f049dcc26a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Iws1-spt6-39C.5712bac30da1af8bd65a81f049dcc26a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Iws1-spt6-39C/Iws1-spt6-39C && \
annotatePeaks.pl \
  peak_call/Iws1-spt6-39C/Iws1-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Iws1-spt6-39C/Iws1-spt6-39C \
  -genomeOntology annotation/Iws1-spt6-39C/Iws1-spt6-39C \
  > annotation/Iws1-spt6-39C/Iws1-spt6-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Iws1-spt6-39C/Iws1-spt6-39C.annotated.csv",
  "annotation/Iws1-spt6-39C/Iws1-spt6-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Iws1-spt6-39C.5712bac30da1af8bd65a81f049dcc26a.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt2-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Spt2-spt6-39C.915739db63fa65f36df7cb8968f58e1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Spt2-spt6-39C.915739db63fa65f36df7cb8968f58e1b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Spt2-spt6-39C/Spt2-spt6-39C && \
annotatePeaks.pl \
  peak_call/Spt2-spt6-39C/Spt2-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Spt2-spt6-39C/Spt2-spt6-39C \
  -genomeOntology annotation/Spt2-spt6-39C/Spt2-spt6-39C \
  > annotation/Spt2-spt6-39C/Spt2-spt6-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Spt2-spt6-39C/Spt2-spt6-39C.annotated.csv",
  "annotation/Spt2-spt6-39C/Spt2-spt6-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Spt2-spt6-39C.915739db63fa65f36df7cb8968f58e1b.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.Chd1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-WT-39C
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Chd1-WT-39C.3990e66b1971a03c35474d8843e5a896.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Chd1-WT-39C.3990e66b1971a03c35474d8843e5a896.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Chd1-WT-39C/Chd1-WT-39C && \
annotatePeaks.pl \
  peak_call/Chd1-WT-39C/Chd1-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Chd1-WT-39C/Chd1-WT-39C \
  -genomeOntology annotation/Chd1-WT-39C/Chd1-WT-39C \
  > annotation/Chd1-WT-39C/Chd1-WT-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Chd1-WT-39C/Chd1-WT-39C.annotated.csv",
  "annotation/Chd1-WT-39C/Chd1-WT-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Chd1-WT-39C.3990e66b1971a03c35474d8843e5a896.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.Iws1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-WT-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Iws1-WT-39C.1f0cbb810e1170f19286b4a2d44292aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Iws1-WT-39C.1f0cbb810e1170f19286b4a2d44292aa.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Iws1-WT-39C/Iws1-WT-39C && \
annotatePeaks.pl \
  peak_call/Iws1-WT-39C/Iws1-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Iws1-WT-39C/Iws1-WT-39C \
  -genomeOntology annotation/Iws1-WT-39C/Iws1-WT-39C \
  > annotation/Iws1-WT-39C/Iws1-WT-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Iws1-WT-39C/Iws1-WT-39C.annotated.csv",
  "annotation/Iws1-WT-39C/Iws1-WT-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Iws1-WT-39C.1f0cbb810e1170f19286b4a2d44292aa.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.Spt2-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt2-WT-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Spt2-WT-39C.622edac71864c149f87cf4f849894447.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Spt2-WT-39C.622edac71864c149f87cf4f849894447.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Spt2-WT-39C/Spt2-WT-39C && \
annotatePeaks.pl \
  peak_call/Spt2-WT-39C/Spt2-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Spt2-WT-39C/Spt2-WT-39C \
  -genomeOntology annotation/Spt2-WT-39C/Spt2-WT-39C \
  > annotation/Spt2-WT-39C/Spt2-WT-39C.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Spt2-WT-39C/Spt2-WT-39C.annotated.csv",
  "annotation/Spt2-WT-39C/Spt2-WT-39C",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Spt2-WT-39C.622edac71864c149f87cf4f849894447.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks.Chd1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Chd1-WT-Ctl.46e8af83812ce4d5069e7f0874ab4745.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Chd1-WT-Ctl.46e8af83812ce4d5069e7f0874ab4745.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Chd1-WT-Ctl/Chd1-WT-Ctl && \
annotatePeaks.pl \
  peak_call/Chd1-WT-Ctl/Chd1-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Chd1-WT-Ctl/Chd1-WT-Ctl \
  -genomeOntology annotation/Chd1-WT-Ctl/Chd1-WT-Ctl \
  > annotation/Chd1-WT-Ctl/Chd1-WT-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Chd1-WT-Ctl/Chd1-WT-Ctl.annotated.csv",
  "annotation/Chd1-WT-Ctl/Chd1-WT-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Chd1-WT-Ctl.46e8af83812ce4d5069e7f0874ab4745.mugqic.done
)
homer_annotate_peaks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_11_JOB_ID: homer_annotate_peaks.Iws1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Iws1-WT-Ctl.7b538806c277853d602e6fdf6a73dd9f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Iws1-WT-Ctl.7b538806c277853d602e6fdf6a73dd9f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Iws1-WT-Ctl/Iws1-WT-Ctl && \
annotatePeaks.pl \
  peak_call/Iws1-WT-Ctl/Iws1-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Iws1-WT-Ctl/Iws1-WT-Ctl \
  -genomeOntology annotation/Iws1-WT-Ctl/Iws1-WT-Ctl \
  > annotation/Iws1-WT-Ctl/Iws1-WT-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Iws1-WT-Ctl/Iws1-WT-Ctl.annotated.csv",
  "annotation/Iws1-WT-Ctl/Iws1-WT-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Iws1-WT-Ctl.7b538806c277853d602e6fdf6a73dd9f.mugqic.done
)
homer_annotate_peaks_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_12_JOB_ID: homer_annotate_peaks.Spt6-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt6-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.Spt6-WT-Ctl.8cef75ad391f363b5f50b24ba583aead.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.Spt6-WT-Ctl.8cef75ad391f363b5f50b24ba583aead.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/Spt6-WT-Ctl/Spt6-WT-Ctl && \
annotatePeaks.pl \
  peak_call/Spt6-WT-Ctl/Spt6-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  -gsize R64-1-1 \
  -cons -CpG \
  -go annotation/Spt6-WT-Ctl/Spt6-WT-Ctl \
  -genomeOntology annotation/Spt6-WT-Ctl/Spt6-WT-Ctl \
  > annotation/Spt6-WT-Ctl/Spt6-WT-Ctl.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/Spt6-WT-Ctl/Spt6-WT-Ctl.annotated.csv",
  "annotation/Spt6-WT-Ctl/Spt6-WT-Ctl",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.Spt6-WT-Ctl.8cef75ad391f363b5f50b24ba583aead.mugqic.done
)
homer_annotate_peaks_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_13_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID:$homer_annotate_peaks_11_JOB_ID:$homer_annotate_peaks_12_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.069bb177345a0c054f1d04748bde6eb0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.069bb177345a0c054f1d04748bde6eb0.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in Chd1-dspt2-Ctl Iws1-dspt2-Ctl Spt6-dspt2-Ctl Chd1-spt6-39C Iws1-spt6-39C Spt2-spt6-39C Chd1-WT-39C Iws1-WT-39C Spt2-WT-39C Chd1-WT-Ctl Iws1-WT-Ctl Spt6-WT-Ctl
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.069bb177345a0c054f1d04748bde6eb0.mugqic.done
)
homer_annotate_peaks_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.Chd1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Chd1-dspt2-Ctl.acd78bd80f2b58ee2982605d2a5de855.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Chd1-dspt2-Ctl.acd78bd80f2b58ee2982605d2a5de855.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl && \
findMotifsGenome.pl \
  peak_call/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl \
  -preparsedDir annotation/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Chd1-dspt2-Ctl.acd78bd80f2b58ee2982605d2a5de855.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.Iws1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Iws1-dspt2-Ctl.3d8f2f4aebebc5eed6fc8dbaa6043831.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Iws1-dspt2-Ctl.3d8f2f4aebebc5eed6fc8dbaa6043831.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl && \
findMotifsGenome.pl \
  peak_call/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl \
  -preparsedDir annotation/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Iws1-dspt2-Ctl.3d8f2f4aebebc5eed6fc8dbaa6043831.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.Spt6-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt6-dspt2-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Spt6-dspt2-Ctl.0f45cb1832965936f6d55efa237984b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Spt6-dspt2-Ctl.0f45cb1832965936f6d55efa237984b8.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl && \
findMotifsGenome.pl \
  peak_call/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl \
  -preparsedDir annotation/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Spt6-dspt2-Ctl.0f45cb1832965936f6d55efa237984b8.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Chd1-spt6-39C.7f7bf74d51529e8e16cd9d51017cd133.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Chd1-spt6-39C.7f7bf74d51529e8e16cd9d51017cd133.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Chd1-spt6-39C/Chd1-spt6-39C && \
findMotifsGenome.pl \
  peak_call/Chd1-spt6-39C/Chd1-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Chd1-spt6-39C/Chd1-spt6-39C \
  -preparsedDir annotation/Chd1-spt6-39C/Chd1-spt6-39C/preparsed \
  -p 4
homer_find_motifs_genome.Chd1-spt6-39C.7f7bf74d51529e8e16cd9d51017cd133.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.Iws1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-spt6-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Iws1-spt6-39C.33b8ffbe0b661d55146c087551cc62ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Iws1-spt6-39C.33b8ffbe0b661d55146c087551cc62ef.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Iws1-spt6-39C/Iws1-spt6-39C && \
findMotifsGenome.pl \
  peak_call/Iws1-spt6-39C/Iws1-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Iws1-spt6-39C/Iws1-spt6-39C \
  -preparsedDir annotation/Iws1-spt6-39C/Iws1-spt6-39C/preparsed \
  -p 4
homer_find_motifs_genome.Iws1-spt6-39C.33b8ffbe0b661d55146c087551cc62ef.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt2-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Spt2-spt6-39C.2ddca48ec093b10b66635ff0e953171a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Spt2-spt6-39C.2ddca48ec093b10b66635ff0e953171a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Spt2-spt6-39C/Spt2-spt6-39C && \
findMotifsGenome.pl \
  peak_call/Spt2-spt6-39C/Spt2-spt6-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Spt2-spt6-39C/Spt2-spt6-39C \
  -preparsedDir annotation/Spt2-spt6-39C/Spt2-spt6-39C/preparsed \
  -p 4
homer_find_motifs_genome.Spt2-spt6-39C.2ddca48ec093b10b66635ff0e953171a.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.Chd1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-WT-39C
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Chd1-WT-39C.973efeb15bcf46bff895d0ccaed5f40e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Chd1-WT-39C.973efeb15bcf46bff895d0ccaed5f40e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Chd1-WT-39C/Chd1-WT-39C && \
findMotifsGenome.pl \
  peak_call/Chd1-WT-39C/Chd1-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Chd1-WT-39C/Chd1-WT-39C \
  -preparsedDir annotation/Chd1-WT-39C/Chd1-WT-39C/preparsed \
  -p 4
homer_find_motifs_genome.Chd1-WT-39C.973efeb15bcf46bff895d0ccaed5f40e.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.Iws1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-WT-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Iws1-WT-39C.cf3c9f6701cba96b40c9c60db624b525.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Iws1-WT-39C.cf3c9f6701cba96b40c9c60db624b525.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Iws1-WT-39C/Iws1-WT-39C && \
findMotifsGenome.pl \
  peak_call/Iws1-WT-39C/Iws1-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Iws1-WT-39C/Iws1-WT-39C \
  -preparsedDir annotation/Iws1-WT-39C/Iws1-WT-39C/preparsed \
  -p 4
homer_find_motifs_genome.Iws1-WT-39C.cf3c9f6701cba96b40c9c60db624b525.mugqic.done
)
homer_find_motifs_genome_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome.Spt2-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt2-WT-39C
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Spt2-WT-39C.b5f512c1e4223585e54f96ff547cf5bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Spt2-WT-39C.b5f512c1e4223585e54f96ff547cf5bb.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Spt2-WT-39C/Spt2-WT-39C && \
findMotifsGenome.pl \
  peak_call/Spt2-WT-39C/Spt2-WT-39C_peaks.narrowPeak \
  R64-1-1 \
  annotation/Spt2-WT-39C/Spt2-WT-39C \
  -preparsedDir annotation/Spt2-WT-39C/Spt2-WT-39C/preparsed \
  -p 4
homer_find_motifs_genome.Spt2-WT-39C.b5f512c1e4223585e54f96ff547cf5bb.mugqic.done
)
homer_find_motifs_genome_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_10_JOB_ID: homer_find_motifs_genome.Chd1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Chd1-WT-Ctl.1a665079569f9760bb888bf2eebb5951.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Chd1-WT-Ctl.1a665079569f9760bb888bf2eebb5951.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Chd1-WT-Ctl/Chd1-WT-Ctl && \
findMotifsGenome.pl \
  peak_call/Chd1-WT-Ctl/Chd1-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Chd1-WT-Ctl/Chd1-WT-Ctl \
  -preparsedDir annotation/Chd1-WT-Ctl/Chd1-WT-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Chd1-WT-Ctl.1a665079569f9760bb888bf2eebb5951.mugqic.done
)
homer_find_motifs_genome_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_11_JOB_ID: homer_find_motifs_genome.Iws1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Iws1-WT-Ctl.48bf32f20aacf14839231043a4da8b7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Iws1-WT-Ctl.48bf32f20aacf14839231043a4da8b7e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Iws1-WT-Ctl/Iws1-WT-Ctl && \
findMotifsGenome.pl \
  peak_call/Iws1-WT-Ctl/Iws1-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Iws1-WT-Ctl/Iws1-WT-Ctl \
  -preparsedDir annotation/Iws1-WT-Ctl/Iws1-WT-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Iws1-WT-Ctl.48bf32f20aacf14839231043a4da8b7e.mugqic.done
)
homer_find_motifs_genome_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_12_JOB_ID: homer_find_motifs_genome.Spt6-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt6-WT-Ctl
JOB_DEPENDENCIES=
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.Spt6-WT-Ctl.b865c0d1321337126e091e45aafa98fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.Spt6-WT-Ctl.b865c0d1321337126e091e45aafa98fc.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/Spt6-WT-Ctl/Spt6-WT-Ctl && \
findMotifsGenome.pl \
  peak_call/Spt6-WT-Ctl/Spt6-WT-Ctl_peaks.narrowPeak \
  R64-1-1 \
  annotation/Spt6-WT-Ctl/Spt6-WT-Ctl \
  -preparsedDir annotation/Spt6-WT-Ctl/Spt6-WT-Ctl/preparsed \
  -p 4
homer_find_motifs_genome.Spt6-WT-Ctl.b865c0d1321337126e091e45aafa98fc.mugqic.done
)
homer_find_motifs_genome_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_13_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID:$homer_find_motifs_genome_6_JOB_ID:$homer_find_motifs_genome_7_JOB_ID:$homer_find_motifs_genome_8_JOB_ID:$homer_find_motifs_genome_9_JOB_ID:$homer_find_motifs_genome_10_JOB_ID:$homer_find_motifs_genome_11_JOB_ID:$homer_find_motifs_genome_12_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.7adceae880cf72feefff8e7c04341786.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.7adceae880cf72feefff8e7c04341786.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in Chd1-dspt2-Ctl Iws1-dspt2-Ctl Spt6-dspt2-Ctl Chd1-spt6-39C Iws1-spt6-39C Spt2-spt6-39C Chd1-WT-39C Iws1-WT-39C Spt2-WT-39C Chd1-WT-Ctl Iws1-WT-Ctl Spt6-WT-Ctl
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.7adceae880cf72feefff8e7c04341786.mugqic.done
)
homer_find_motifs_genome_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs_1_JOB_ID: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID:$homer_annotate_peaks_11_JOB_ID:$homer_annotate_peaks_12_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.e0c9ff4f662eaec3932888081d00e73e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.e0c9ff4f662eaec3932888081d00e73e.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../raw/design.txt \
  /gs/project/eav-760-aa/Chaperone/output/pipeline && \
mkdir -p report/annotation/ && \
if [[ -f annotation/peak_stats.csv ]]
then
  cp annotation/peak_stats.csv report/annotation/
peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, $2,  sprintf("%\47d", $3), $4, sprintf("%\47.1f", $5), sprintf("%\47.1f", $6), sprintf("%\47.1f", $7), sprintf("%\47.1f", $8)}}' annotation/peak_stats.csv`
else
  peak_stats_table=""
fi
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.annotation_graphs.md \
  --variable peak_stats_table="$peak_stats_table" \
  --variable proximal_distance="2" \
  --variable distal_distance="10" \
  --variable distance5d_lower="10" \
  --variable distance5d_upper="100" \
  --variable gene_desert_size="100" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.annotation_graphs.md \
  > report/ChipSeq.annotation_graphs.md && \
for contrast in Chd1-dspt2-Ctl Iws1-dspt2-Ctl Spt6-dspt2-Ctl Chd1-spt6-39C Iws1-spt6-39C Spt2-spt6-39C Chd1-WT-39C Iws1-WT-39C Spt2-WT-39C Chd1-WT-Ctl Iws1-WT-Ctl Spt6-WT-Ctl
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.e0c9ff4f662eaec3932888081d00e73e.mugqic.done
)
annotation_graphs_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$annotation_graphs_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r14-n04&ip=10.241.129.4&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=25" --quiet --output-document=/dev/null

