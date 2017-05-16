#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-10T14:18:30
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 25 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 26 jobs
#   samtools_view_filter: 26 jobs
#   picard_merge_sam_files: 25 jobs
#   picard_mark_duplicates: 26 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 25 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 26 jobs
#   macs2_callpeak: 13 jobs
#   homer_annotate_peaks: 13 jobs
#   homer_find_motifs_genome: 13 jobs
#   annotation_graphs: 1 job
#   TOTAL: 223 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/project/eav-760-aa/Chaperone/output/pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.Chd1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_dspt2Cl2_RS.559d93e2e7f0e1de0b96c92802ee4029.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_dspt2Cl2_RS.559d93e2e7f0e1de0b96c92802ee4029.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_dspt2Cl2 && \
`cat > trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_8.Chd1-Myc_dspt2Cl2_R1.fastq.gz \
  trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.log
trimmomatic.Chd1-Myc_dspt2Cl2_RS.559d93e2e7f0e1de0b96c92802ee4029.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.Chd1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_dspt2Cl1_RS.18e963d78a39266a05ae680bf378ccfc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_dspt2Cl1_RS.18e963d78a39266a05ae680bf378ccfc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_dspt2Cl1 && \
`cat > trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_9.Chd1-Myc_dspt2Cl1_R1.fastq.gz \
  trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.log
trimmomatic.Chd1-Myc_dspt2Cl1_RS.18e963d78a39266a05ae680bf378ccfc.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.Iws1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_dspt2Cl2_RS.aed169afbf074a08fdae0f68bb58da20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_dspt2Cl2_RS.aed169afbf074a08fdae0f68bb58da20.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_dspt2Cl2 && \
`cat > trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_15.Iws1-Myc_dspt2Cl2_R1.fastq.gz \
  trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.log
trimmomatic.Iws1-Myc_dspt2Cl2_RS.aed169afbf074a08fdae0f68bb58da20.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.Iws1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_dspt2Cl1_RS.0b04f3413c8f9a2f9d4d02f3a1b760d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_dspt2Cl1_RS.0b04f3413c8f9a2f9d4d02f3a1b760d0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_dspt2Cl1 && \
`cat > trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_6.Iws1-Myc_dspt2Cl1_R1.fastq.gz \
  trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.log
trimmomatic.Iws1-Myc_dspt2Cl1_RS.0b04f3413c8f9a2f9d4d02f3a1b760d0.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.Spt6-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt6-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt6-Myc_dspt2Cl2_RS.8fc01831771c03b29b7eae9929860f08.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt6-Myc_dspt2Cl2_RS.8fc01831771c03b29b7eae9929860f08.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt6-Myc_dspt2Cl2 && \
`cat > trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_27.Spt6-Myc_dspt2Cl2_R1.fastq.gz \
  trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.log
trimmomatic.Spt6-Myc_dspt2Cl2_RS.8fc01831771c03b29b7eae9929860f08.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.Spt6-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt6-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt6-Myc_dspt2Cl1_RS.ca01436df9bfc50cc75be318b7dff65c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt6-Myc_dspt2Cl1_RS.ca01436df9bfc50cc75be318b7dff65c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt6-Myc_dspt2Cl1 && \
`cat > trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_23.Spt6-Myc_dspt2Cl1_R1.fastq.gz \
  trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.log
trimmomatic.Spt6-Myc_dspt2Cl1_RS.ca01436df9bfc50cc75be318b7dff65c.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.Chd1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_spt6_39C_Cl2_RS.f3a58b89e71627ec96025aadce7e5d2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_spt6_39C_Cl2_RS.f3a58b89e71627ec96025aadce7e5d2b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_spt6_39C_Cl2 && \
`cat > trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_22.Chd1-Myc_spt6_39C_Cl2_R1.fastq.gz \
  trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.log
trimmomatic.Chd1-Myc_spt6_39C_Cl2_RS.f3a58b89e71627ec96025aadce7e5d2b.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.Chd1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_spt6_39C_Cl1_RS.6eece0bffa1d97c739c931009f5fbe7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_spt6_39C_Cl1_RS.6eece0bffa1d97c739c931009f5fbe7f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_spt6_39C_Cl1 && \
`cat > trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_20.Chd1-Myc_spt6_39C_Cl1_R1.fastq.gz \
  trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.log
trimmomatic.Chd1-Myc_spt6_39C_Cl1_RS.6eece0bffa1d97c739c931009f5fbe7f.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.Iws1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_spt6_39C_Cl2_RS.8ef81a8009d5e979d8b901187070695b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_spt6_39C_Cl2_RS.8ef81a8009d5e979d8b901187070695b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_spt6_39C_Cl2 && \
`cat > trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_16.Iws1-Myc_spt6_39C_Cl2_R1.fastq.gz \
  trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.log
trimmomatic.Iws1-Myc_spt6_39C_Cl2_RS.8ef81a8009d5e979d8b901187070695b.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.Iws1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_spt6_39C_Cl1_RS.4819ceb03ad5e40bbd2a3ea07c457e1d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_spt6_39C_Cl1_RS.4819ceb03ad5e40bbd2a3ea07c457e1d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_spt6_39C_Cl1 && \
`cat > trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_14.Iws1-Myc_spt6_39C_Cl1_R1.fastq.gz \
  trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.log
trimmomatic.Iws1-Myc_spt6_39C_Cl1_RS.4819ceb03ad5e40bbd2a3ea07c457e1d.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.Spt2-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt2-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt2-Myc_spt6_39C_Cl2_RS.cbf41ee2209c97554b73ba0934e9d2d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt2-Myc_spt6_39C_Cl2_RS.cbf41ee2209c97554b73ba0934e9d2d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt2-Myc_spt6_39C_Cl2 && \
`cat > trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_19.Spt2-Myc_spt6_39C_Cl2_R1.fastq.gz \
  trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.log
trimmomatic.Spt2-Myc_spt6_39C_Cl2_RS.cbf41ee2209c97554b73ba0934e9d2d8.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.Spt2-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt2-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt2-Myc_spt6_39C_Cl1_RS.c0558a3f95fc8cc809a98474867574e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt2-Myc_spt6_39C_Cl1_RS.c0558a3f95fc8cc809a98474867574e9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt2-Myc_spt6_39C_Cl1 && \
`cat > trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_12.Spt2-Myc_spt6_39C_Cl1_R1.fastq.gz \
  trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.log
trimmomatic.Spt2-Myc_spt6_39C_Cl1_RS.c0558a3f95fc8cc809a98474867574e9.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.Chd1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_Wt_39C_Cl2_RS.26766ce7542476841d330d343e36bdb0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_Wt_39C_Cl2_RS.26766ce7542476841d330d343e36bdb0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_Wt_39C_Cl2 && \
`cat > trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_11.Chd1-Myc_Wt_39C_Cl2_R1.fastq.gz \
  trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.log
trimmomatic.Chd1-Myc_Wt_39C_Cl2_RS.26766ce7542476841d330d343e36bdb0.mugqic.done
)
trimmomatic_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_14_JOB_ID: trimmomatic.Chd1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_Wt_39C_Cl1_RS.05624ec96c12e953083370923397fb47.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_Wt_39C_Cl1_RS.05624ec96c12e953083370923397fb47.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_Wt_39C_Cl1 && \
`cat > trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_10.Chd1-Myc_Wt_39C_Cl1_R1.fastq.gz \
  trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.log
trimmomatic.Chd1-Myc_Wt_39C_Cl1_RS.05624ec96c12e953083370923397fb47.mugqic.done
)
trimmomatic_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_15_JOB_ID: trimmomatic.Iws1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_Wt_39C_Cl2_RS.83c631c6e897bf2bdcf79211ac5a42fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_Wt_39C_Cl2_RS.83c631c6e897bf2bdcf79211ac5a42fc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_Wt_39C_Cl2 && \
`cat > trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_18.Iws1-Myc_Wt_39C_Cl2_R1.fastq.gz \
  trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.log
trimmomatic.Iws1-Myc_Wt_39C_Cl2_RS.83c631c6e897bf2bdcf79211ac5a42fc.mugqic.done
)
trimmomatic_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_16_JOB_ID: trimmomatic.Iws1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_Wt_39C_Cl1_RS.4fea594f04b6577633fb98def4d2e6a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_Wt_39C_Cl1_RS.4fea594f04b6577633fb98def4d2e6a5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_Wt_39C_Cl1 && \
`cat > trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_7.Iws1-Myc_Wt_39C_Cl1_R1.fastq.gz \
  trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.log
trimmomatic.Iws1-Myc_Wt_39C_Cl1_RS.4fea594f04b6577633fb98def4d2e6a5.mugqic.done
)
trimmomatic_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_17_JOB_ID: trimmomatic.Spt2-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt2-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt2-Myc_Wt_39C_Cl2_RS.2ce6e3234659b43705fc6e89f774fa2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt2-Myc_Wt_39C_Cl2_RS.2ce6e3234659b43705fc6e89f774fa2b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt2-Myc_Wt_39C_Cl2 && \
`cat > trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_5.Spt2-Myc_Wt_39C_Cl2_R1.fastq.gz \
  trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.log
trimmomatic.Spt2-Myc_Wt_39C_Cl2_RS.2ce6e3234659b43705fc6e89f774fa2b.mugqic.done
)
trimmomatic_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_18_JOB_ID: trimmomatic.Spt2-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt2-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt2-Myc_Wt_39C_Cl1_RS.5c85eb309f587c2840ff24dff259b5ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt2-Myc_Wt_39C_Cl1_RS.5c85eb309f587c2840ff24dff259b5ea.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt2-Myc_Wt_39C_Cl1 && \
`cat > trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_4.Spt2-Myc_Wt_39C_Cl1_R1.fastq.gz \
  trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.log
trimmomatic.Spt2-Myc_Wt_39C_Cl1_RS.5c85eb309f587c2840ff24dff259b5ea.mugqic.done
)
trimmomatic_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_19_JOB_ID: trimmomatic.Chd1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_WtCl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_WtCl2_RS.74dca7189ffd5f2498228d1a453b9d24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_WtCl2_RS.74dca7189ffd5f2498228d1a453b9d24.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_WtCl2 && \
`cat > trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_3.Chd1-Myc_WtCl2_R1.fastq.gz \
  trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.log
trimmomatic.Chd1-Myc_WtCl2_RS.74dca7189ffd5f2498228d1a453b9d24.mugqic.done
)
trimmomatic_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_20_JOB_ID: trimmomatic.Chd1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Chd1-Myc_WtCl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Chd1-Myc_WtCl1_RS.4c27cbf03a2fd47dbdd9da26ac715921.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Chd1-Myc_WtCl1_RS.4c27cbf03a2fd47dbdd9da26ac715921.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Chd1-Myc_WtCl1 && \
`cat > trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_1.Chd1-Myc_WtCl1_R1.fastq.gz \
  trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.log
trimmomatic.Chd1-Myc_WtCl1_RS.4c27cbf03a2fd47dbdd9da26ac715921.mugqic.done
)
trimmomatic_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_21_JOB_ID: trimmomatic.Iws1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_WtCl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_WtCl2_RS.6aa251769053cb5fc0fca9103de90561.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_WtCl2_RS.6aa251769053cb5fc0fca9103de90561.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_WtCl2 && \
`cat > trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_13.Iws1-Myc_WtCl2_R1.fastq.gz \
  trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.log
trimmomatic.Iws1-Myc_WtCl2_RS.6aa251769053cb5fc0fca9103de90561.mugqic.done
)
trimmomatic_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_22_JOB_ID: trimmomatic.Iws1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Iws1-Myc_WtCl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Iws1-Myc_WtCl1_RS.2d43dbdbb760376fbaf3578598f54859.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Iws1-Myc_WtCl1_RS.2d43dbdbb760376fbaf3578598f54859.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Iws1-Myc_WtCl1 && \
`cat > trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4097.008.Index_2.Iws1-Myc_WtCl1_R1.fastq.gz \
  trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.log
trimmomatic.Iws1-Myc_WtCl1_RS.2d43dbdbb760376fbaf3578598f54859.mugqic.done
)
trimmomatic_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_23_JOB_ID: trimmomatic.Spt6-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt6-Myc_WtCl2_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt6-Myc_WtCl2_RS.f68300c91b3bb5e5a423987e0cbfbcd6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt6-Myc_WtCl2_RS.f68300c91b3bb5e5a423987e0cbfbcd6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt6-Myc_WtCl2 && \
`cat > trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_21.Spt6-Myc_WtCl2_R1.fastq.gz \
  trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.log
trimmomatic.Spt6-Myc_WtCl2_RS.f68300c91b3bb5e5a423987e0cbfbcd6.mugqic.done
)
trimmomatic_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_24_JOB_ID: trimmomatic.Spt6-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.Spt6-Myc_WtCl1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.Spt6-Myc_WtCl1_RS.08d6fa8a31077a3d4b63207416ce1011.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.Spt6-Myc_WtCl1_RS.08d6fa8a31077a3d4b63207416ce1011.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/Spt6-Myc_WtCl1 && \
`cat > trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_25.Spt6-Myc_WtCl1_R1.fastq.gz \
  trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.log
trimmomatic.Spt6-Myc_WtCl1_RS.08d6fa8a31077a3d4b63207416ce1011.mugqic.done
)
trimmomatic_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_25_JOB_ID: trimmomatic.No-TAG_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.No-TAG_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.No-TAG_RS.f8fec815b2fe33a748d3a844b4ee516e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.No-TAG_RS.f8fec815b2fe33a748d3a844b4ee516e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/No-TAG && \
`cat > trim/No-TAG/No-TAG_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/eav-760-aa/Chaperone/raw/HI.4117.006.Index_2.No-TAG_R1.fastq.gz \
  trim/No-TAG/No-TAG_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/No-TAG/No-TAG_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/No-TAG/No-TAG_RS.trim.log
trimmomatic.No-TAG_RS.f8fec815b2fe33a748d3a844b4ee516e.mugqic.done
)
trimmomatic_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID:$trimmomatic_13_JOB_ID:$trimmomatic_14_JOB_ID:$trimmomatic_15_JOB_ID:$trimmomatic_16_JOB_ID:$trimmomatic_17_JOB_ID:$trimmomatic_18_JOB_ID:$trimmomatic_19_JOB_ID:$trimmomatic_20_JOB_ID:$trimmomatic_21_JOB_ID:$trimmomatic_22_JOB_ID:$trimmomatic_23_JOB_ID:$trimmomatic_24_JOB_ID:$trimmomatic_25_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.56ea0d277a51d60cb76a7ea1a1a5d964.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.56ea0d277a51d60cb76a7ea1a1a5d964.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_dspt2Cl2	Chd1-Myc_dspt2Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_dspt2Cl1	Chd1-Myc_dspt2Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_dspt2Cl2	Iws1-Myc_dspt2Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_dspt2Cl1	Iws1-Myc_dspt2Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt6-Myc_dspt2Cl2	Spt6-Myc_dspt2Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt6-Myc_dspt2Cl1	Spt6-Myc_dspt2Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_spt6_39C_Cl2	Chd1-Myc_spt6_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_spt6_39C_Cl1	Chd1-Myc_spt6_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_spt6_39C_Cl2	Iws1-Myc_spt6_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_spt6_39C_Cl1	Iws1-Myc_spt6_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt2-Myc_spt6_39C_Cl2	Spt2-Myc_spt6_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt2-Myc_spt6_39C_Cl1	Spt2-Myc_spt6_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_Wt_39C_Cl2	Chd1-Myc_Wt_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_Wt_39C_Cl1	Chd1-Myc_Wt_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_Wt_39C_Cl2	Iws1-Myc_Wt_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_Wt_39C_Cl1	Iws1-Myc_Wt_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt2-Myc_Wt_39C_Cl2	Spt2-Myc_Wt_39C_Cl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt2-Myc_Wt_39C_Cl1	Spt2-Myc_Wt_39C_Cl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_WtCl2	Chd1-Myc_WtCl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Chd1-Myc_WtCl1	Chd1-Myc_WtCl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_WtCl2	Iws1-Myc_WtCl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Iws1-Myc_WtCl1	Iws1-Myc_WtCl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt6-Myc_WtCl2	Spt6-Myc_WtCl2_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/Spt6-Myc_WtCl1	Spt6-Myc_WtCl1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/No-TAG/No-TAG_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/No-TAG	No-TAG_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.56ea0d277a51d60cb76a7ea1a1a5d964.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl2_RS.3c1389d3741189d8d1f0fa45244871b7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl2_RS.3c1389d3741189d8d1f0fa45244871b7.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_dspt2Cl2_RS	SM:Chd1-Myc_dspt2Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS/Chd1-Myc_dspt2Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl2_RS.3c1389d3741189d8d1f0fa45244871b7.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl1_RS.2dea08f04788472f056bd934fc65913a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl1_RS.2dea08f04788472f056bd934fc65913a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_dspt2Cl1_RS	SM:Chd1-Myc_dspt2Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS/Chd1-Myc_dspt2Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_dspt2Cl1_RS.2dea08f04788472f056bd934fc65913a.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl2_RS.8d954f42656be9f7799fe451982e0ef7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl2_RS.8d954f42656be9f7799fe451982e0ef7.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_dspt2Cl2_RS	SM:Iws1-Myc_dspt2Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS/Iws1-Myc_dspt2Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl2_RS.8d954f42656be9f7799fe451982e0ef7.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl1_RS.499d23482ab13edf7cfcdae8f6033305.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl1_RS.499d23482ab13edf7cfcdae8f6033305.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_dspt2Cl1_RS	SM:Iws1-Myc_dspt2Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS/Iws1-Myc_dspt2Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_dspt2Cl1_RS.499d23482ab13edf7cfcdae8f6033305.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl2_RS.f58855c6e30e0b36bdb33e4dde7ae0db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl2_RS.f58855c6e30e0b36bdb33e4dde7ae0db.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt6-Myc_dspt2Cl2_RS	SM:Spt6-Myc_dspt2Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS/Spt6-Myc_dspt2Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl2_RS.f58855c6e30e0b36bdb33e4dde7ae0db.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl1_RS.d6da991a80f9fe510f5dbd4f3039f210.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl1_RS.d6da991a80f9fe510f5dbd4f3039f210.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt6-Myc_dspt2Cl1_RS	SM:Spt6-Myc_dspt2Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS/Spt6-Myc_dspt2Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt6-Myc_dspt2Cl1_RS.d6da991a80f9fe510f5dbd4f3039f210.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl2_RS.06375309d6defd194c21ab3415e52d13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl2_RS.06375309d6defd194c21ab3415e52d13.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_spt6_39C_Cl2_RS	SM:Chd1-Myc_spt6_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS/Chd1-Myc_spt6_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl2_RS.06375309d6defd194c21ab3415e52d13.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl1_RS.b48dad198681acc1fa3c6d1da7fca14a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl1_RS.b48dad198681acc1fa3c6d1da7fca14a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_spt6_39C_Cl1_RS	SM:Chd1-Myc_spt6_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS/Chd1-Myc_spt6_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_spt6_39C_Cl1_RS.b48dad198681acc1fa3c6d1da7fca14a.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl2_RS.e655e73963c2c471fdc849ab2753f67e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl2_RS.e655e73963c2c471fdc849ab2753f67e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_spt6_39C_Cl2_RS	SM:Iws1-Myc_spt6_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS/Iws1-Myc_spt6_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl2_RS.e655e73963c2c471fdc849ab2753f67e.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl1_RS.9a6ffabcdbd34c1d5f1f0cdd8ddc94c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl1_RS.9a6ffabcdbd34c1d5f1f0cdd8ddc94c6.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_spt6_39C_Cl1_RS	SM:Iws1-Myc_spt6_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS/Iws1-Myc_spt6_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_spt6_39C_Cl1_RS.9a6ffabcdbd34c1d5f1f0cdd8ddc94c6.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
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
bwa_mem_picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl1_RS.7a3caac063e0a40fd08bcc7573e91faf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl1_RS.7a3caac063e0a40fd08bcc7573e91faf.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt2-Myc_spt6_39C_Cl1_RS	SM:Spt2-Myc_spt6_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS/Spt2-Myc_spt6_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt2-Myc_spt6_39C_Cl1_RS.7a3caac063e0a40fd08bcc7573e91faf.mugqic.done
)
bwa_mem_picard_sort_sam_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_13_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl2_RS.ee4d038d92baf8fef4a0b9394d2bf6e5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl2_RS.ee4d038d92baf8fef4a0b9394d2bf6e5.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_Wt_39C_Cl2_RS	SM:Chd1-Myc_Wt_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS/Chd1-Myc_Wt_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl2_RS.ee4d038d92baf8fef4a0b9394d2bf6e5.mugqic.done
)
bwa_mem_picard_sort_sam_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_14_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_14_JOB_ID
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
bwa_mem_picard_sort_sam_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_15_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_15_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl2_RS.2aa44c1cbf0a7077d17d9ca1bb6dda93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl2_RS.2aa44c1cbf0a7077d17d9ca1bb6dda93.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_Wt_39C_Cl2_RS	SM:Iws1-Myc_Wt_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS/Iws1-Myc_Wt_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl2_RS.2aa44c1cbf0a7077d17d9ca1bb6dda93.mugqic.done
)
bwa_mem_picard_sort_sam_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_16_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_16_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl1_RS.8c0f1a16e99adc8eabd3cfec92bff8d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl1_RS.8c0f1a16e99adc8eabd3cfec92bff8d6.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_Wt_39C_Cl1_RS	SM:Iws1-Myc_Wt_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS/Iws1-Myc_Wt_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_Wt_39C_Cl1_RS.8c0f1a16e99adc8eabd3cfec92bff8d6.mugqic.done
)
bwa_mem_picard_sort_sam_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_17_JOB_ID: bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$trimmomatic_17_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl2_RS.fd33fa7e0be16cf4948bddf51c7e09ec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl2_RS.fd33fa7e0be16cf4948bddf51c7e09ec.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt2-Myc_Wt_39C_Cl2_RS	SM:Spt2-Myc_Wt_39C_Cl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS/Spt2-Myc_Wt_39C_Cl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl2_RS.fd33fa7e0be16cf4948bddf51c7e09ec.mugqic.done
)
bwa_mem_picard_sort_sam_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_18_JOB_ID: bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$trimmomatic_18_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl1_RS.5d189ffd2722fdfdbfb5e3b6830bca78.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl1_RS.5d189ffd2722fdfdbfb5e3b6830bca78.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt2-Myc_Wt_39C_Cl1_RS	SM:Spt2-Myc_Wt_39C_Cl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS/Spt2-Myc_Wt_39C_Cl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt2-Myc_Wt_39C_Cl1_RS.5d189ffd2722fdfdbfb5e3b6830bca78.mugqic.done
)
bwa_mem_picard_sort_sam_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_19_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_WtCl2_RS
JOB_DEPENDENCIES=$trimmomatic_19_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_WtCl2_RS.53f12eabe1f50f159f20384485c19449.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_WtCl2_RS.53f12eabe1f50f159f20384485c19449.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_WtCl2_RS	SM:Chd1-Myc_WtCl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS/Chd1-Myc_WtCl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_WtCl2_RS.53f12eabe1f50f159f20384485c19449.mugqic.done
)
bwa_mem_picard_sort_sam_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_20_JOB_ID: bwa_mem_picard_sort_sam.Chd1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Chd1-Myc_WtCl1_RS
JOB_DEPENDENCIES=$trimmomatic_20_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Chd1-Myc_WtCl1_RS.0d46c14711942e32e505edf45292a37b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Chd1-Myc_WtCl1_RS.0d46c14711942e32e505edf45292a37b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Chd1-Myc_WtCl1_RS	SM:Chd1-Myc_WtCl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS/Chd1-Myc_WtCl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Chd1-Myc_WtCl1_RS.0d46c14711942e32e505edf45292a37b.mugqic.done
)
bwa_mem_picard_sort_sam_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_21_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_WtCl2_RS
JOB_DEPENDENCIES=$trimmomatic_21_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_WtCl2_RS.ea9ad0c69116d740a0a47faec4a3e3b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_WtCl2_RS.ea9ad0c69116d740a0a47faec4a3e3b1.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_WtCl2_RS	SM:Iws1-Myc_WtCl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS/Iws1-Myc_WtCl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_WtCl2_RS.ea9ad0c69116d740a0a47faec4a3e3b1.mugqic.done
)
bwa_mem_picard_sort_sam_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_22_JOB_ID: bwa_mem_picard_sort_sam.Iws1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Iws1-Myc_WtCl1_RS
JOB_DEPENDENCIES=$trimmomatic_22_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Iws1-Myc_WtCl1_RS.903fe14ac37f588d020910ae09363354.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Iws1-Myc_WtCl1_RS.903fe14ac37f588d020910ae09363354.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Iws1-Myc_WtCl1_RS	SM:Iws1-Myc_WtCl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS/Iws1-Myc_WtCl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Iws1-Myc_WtCl1_RS.903fe14ac37f588d020910ae09363354.mugqic.done
)
bwa_mem_picard_sort_sam_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_23_JOB_ID: bwa_mem_picard_sort_sam.Spt6-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt6-Myc_WtCl2_RS
JOB_DEPENDENCIES=$trimmomatic_23_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt6-Myc_WtCl2_RS.8b5ba785a1a119b9aa7ba44a8f45146d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt6-Myc_WtCl2_RS.8b5ba785a1a119b9aa7ba44a8f45146d.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt6-Myc_WtCl2_RS	SM:Spt6-Myc_WtCl2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS/Spt6-Myc_WtCl2_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt6-Myc_WtCl2_RS.8b5ba785a1a119b9aa7ba44a8f45146d.mugqic.done
)
bwa_mem_picard_sort_sam_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_24_JOB_ID: bwa_mem_picard_sort_sam.Spt6-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.Spt6-Myc_WtCl1_RS
JOB_DEPENDENCIES=$trimmomatic_24_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.Spt6-Myc_WtCl1_RS.52c5e1033b0b63672e821ccdfa91236e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.Spt6-Myc_WtCl1_RS.52c5e1033b0b63672e821ccdfa91236e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:Spt6-Myc_WtCl1_RS	SM:Spt6-Myc_WtCl1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS/Spt6-Myc_WtCl1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.Spt6-Myc_WtCl1_RS.52c5e1033b0b63672e821ccdfa91236e.mugqic.done
)
bwa_mem_picard_sort_sam_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_25_JOB_ID: bwa_mem_picard_sort_sam.No-TAG_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.No-TAG_RS
JOB_DEPENDENCIES=$trimmomatic_25_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.No-TAG_RS.a36e13fd5ab0949e14f1cf20bfd817b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.No-TAG_RS.a36e13fd5ab0949e14f1cf20bfd817b6.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/No-TAG/No-TAG_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:No-TAG_RS	SM:No-TAG	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Saccharomyces_cerevisiae.R64-1-1/genome/bwa_index/Saccharomyces_cerevisiae.R64-1-1.fa \
  trim/No-TAG/No-TAG_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/No-TAG/No-TAG_RS/No-TAG_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.No-TAG_RS.a36e13fd5ab0949e14f1cf20bfd817b6.mugqic.done
)
bwa_mem_picard_sort_sam_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_26_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID:$bwa_mem_picard_sort_sam_11_JOB_ID:$bwa_mem_picard_sort_sam_12_JOB_ID:$bwa_mem_picard_sort_sam_13_JOB_ID:$bwa_mem_picard_sort_sam_14_JOB_ID:$bwa_mem_picard_sort_sam_15_JOB_ID:$bwa_mem_picard_sort_sam_16_JOB_ID:$bwa_mem_picard_sort_sam_17_JOB_ID:$bwa_mem_picard_sort_sam_18_JOB_ID:$bwa_mem_picard_sort_sam_19_JOB_ID:$bwa_mem_picard_sort_sam_20_JOB_ID:$bwa_mem_picard_sort_sam_21_JOB_ID:$bwa_mem_picard_sort_sam_22_JOB_ID:$bwa_mem_picard_sort_sam_23_JOB_ID:$bwa_mem_picard_sort_sam_24_JOB_ID:$bwa_mem_picard_sort_sam_25_JOB_ID
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
bwa_mem_picard_sort_sam_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.Chd1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_dspt2Cl2_RS.28daba30178545bc32e124aaa9499c5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_dspt2Cl2_RS.28daba30178545bc32e124aaa9499c5b.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS/Chd1-Myc_dspt2Cl2_RS.sorted.bam \
  > alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2_RS/Chd1-Myc_dspt2Cl2_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_dspt2Cl2_RS.28daba30178545bc32e124aaa9499c5b.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.Chd1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_dspt2Cl1_RS.e672fa8fc106f81185ff2b37f775d3b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_dspt2Cl1_RS.e672fa8fc106f81185ff2b37f775d3b1.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS/Chd1-Myc_dspt2Cl1_RS.sorted.bam \
  > alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1_RS/Chd1-Myc_dspt2Cl1_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_dspt2Cl1_RS.e672fa8fc106f81185ff2b37f775d3b1.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.Iws1-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_dspt2Cl2_RS.985a17adb4c5195e42394697949af3ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_dspt2Cl2_RS.985a17adb4c5195e42394697949af3ba.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS/Iws1-Myc_dspt2Cl2_RS.sorted.bam \
  > alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2_RS/Iws1-Myc_dspt2Cl2_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_dspt2Cl2_RS.985a17adb4c5195e42394697949af3ba.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.Iws1-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_dspt2Cl1_RS.2578a224069c8a7c0510f767696662df.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_dspt2Cl1_RS.2578a224069c8a7c0510f767696662df.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS/Iws1-Myc_dspt2Cl1_RS.sorted.bam \
  > alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1_RS/Iws1-Myc_dspt2Cl1_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_dspt2Cl1_RS.2578a224069c8a7c0510f767696662df.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.Spt6-Myc_dspt2Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt6-Myc_dspt2Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt6-Myc_dspt2Cl2_RS.f21c15731dcc893fe4f1f4c0db056286.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt6-Myc_dspt2Cl2_RS.f21c15731dcc893fe4f1f4c0db056286.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS/Spt6-Myc_dspt2Cl2_RS.sorted.bam \
  > alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2_RS/Spt6-Myc_dspt2Cl2_RS.sorted.filtered.bam
samtools_view_filter.Spt6-Myc_dspt2Cl2_RS.f21c15731dcc893fe4f1f4c0db056286.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.Spt6-Myc_dspt2Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt6-Myc_dspt2Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt6-Myc_dspt2Cl1_RS.8aef4428ed2dde18b83e3970523d8306.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt6-Myc_dspt2Cl1_RS.8aef4428ed2dde18b83e3970523d8306.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS/Spt6-Myc_dspt2Cl1_RS.sorted.bam \
  > alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1_RS/Spt6-Myc_dspt2Cl1_RS.sorted.filtered.bam
samtools_view_filter.Spt6-Myc_dspt2Cl1_RS.8aef4428ed2dde18b83e3970523d8306.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.Chd1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_spt6_39C_Cl2_RS.2fc7a2c58394736746f9a79d2c9a982a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_spt6_39C_Cl2_RS.2fc7a2c58394736746f9a79d2c9a982a.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS/Chd1-Myc_spt6_39C_Cl2_RS.sorted.bam \
  > alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2_RS/Chd1-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_spt6_39C_Cl2_RS.2fc7a2c58394736746f9a79d2c9a982a.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.Chd1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_spt6_39C_Cl1_RS.7788c0b344b2847aa6304b977807cdbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_spt6_39C_Cl1_RS.7788c0b344b2847aa6304b977807cdbd.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS/Chd1-Myc_spt6_39C_Cl1_RS.sorted.bam \
  > alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1_RS/Chd1-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_spt6_39C_Cl1_RS.7788c0b344b2847aa6304b977807cdbd.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.Iws1-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_spt6_39C_Cl2_RS.5aa64891ffba46d6eb2c703009fef9dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_spt6_39C_Cl2_RS.5aa64891ffba46d6eb2c703009fef9dd.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS/Iws1-Myc_spt6_39C_Cl2_RS.sorted.bam \
  > alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2_RS/Iws1-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_spt6_39C_Cl2_RS.5aa64891ffba46d6eb2c703009fef9dd.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.Iws1-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_spt6_39C_Cl1_RS.ad1de0f4dfe055110025c2191703b0d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_spt6_39C_Cl1_RS.ad1de0f4dfe055110025c2191703b0d9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS/Iws1-Myc_spt6_39C_Cl1_RS.sorted.bam \
  > alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1_RS/Iws1-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_spt6_39C_Cl1_RS.ad1de0f4dfe055110025c2191703b0d9.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt2-Myc_spt6_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
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
samtools_view_filter_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.Spt2-Myc_spt6_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt2-Myc_spt6_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt2-Myc_spt6_39C_Cl1_RS.01ed9ccb5da216f6e2ae520f714167bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt2-Myc_spt6_39C_Cl1_RS.01ed9ccb5da216f6e2ae520f714167bb.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS/Spt2-Myc_spt6_39C_Cl1_RS.sorted.bam \
  > alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1_RS/Spt2-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Spt2-Myc_spt6_39C_Cl1_RS.01ed9ccb5da216f6e2ae520f714167bb.mugqic.done
)
samtools_view_filter_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_13_JOB_ID: samtools_view_filter.Chd1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_13_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_Wt_39C_Cl2_RS.fd8e399de6aaf9491a02fd778415203e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_Wt_39C_Cl2_RS.fd8e399de6aaf9491a02fd778415203e.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS/Chd1-Myc_Wt_39C_Cl2_RS.sorted.bam \
  > alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2_RS/Chd1-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_Wt_39C_Cl2_RS.fd8e399de6aaf9491a02fd778415203e.mugqic.done
)
samtools_view_filter_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_14_JOB_ID: samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_14_JOB_ID
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
samtools_view_filter_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_15_JOB_ID: samtools_view_filter.Iws1-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_15_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_Wt_39C_Cl2_RS.596ce8feb972bf5a1a744485a53f85c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_Wt_39C_Cl2_RS.596ce8feb972bf5a1a744485a53f85c3.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS/Iws1-Myc_Wt_39C_Cl2_RS.sorted.bam \
  > alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2_RS/Iws1-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_Wt_39C_Cl2_RS.596ce8feb972bf5a1a744485a53f85c3.mugqic.done
)
samtools_view_filter_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_16_JOB_ID: samtools_view_filter.Iws1-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_16_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_Wt_39C_Cl1_RS.6cafcc9b5c7ae0d1c07cd8b66b2bbcf9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_Wt_39C_Cl1_RS.6cafcc9b5c7ae0d1c07cd8b66b2bbcf9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS/Iws1-Myc_Wt_39C_Cl1_RS.sorted.bam \
  > alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1_RS/Iws1-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_Wt_39C_Cl1_RS.6cafcc9b5c7ae0d1c07cd8b66b2bbcf9.mugqic.done
)
samtools_view_filter_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_17_JOB_ID: samtools_view_filter.Spt2-Myc_Wt_39C_Cl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt2-Myc_Wt_39C_Cl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_17_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt2-Myc_Wt_39C_Cl2_RS.a5ab79a86cc54efb57adee6cee706af0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt2-Myc_Wt_39C_Cl2_RS.a5ab79a86cc54efb57adee6cee706af0.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS/Spt2-Myc_Wt_39C_Cl2_RS.sorted.bam \
  > alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2_RS/Spt2-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam
samtools_view_filter.Spt2-Myc_Wt_39C_Cl2_RS.a5ab79a86cc54efb57adee6cee706af0.mugqic.done
)
samtools_view_filter_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_18_JOB_ID: samtools_view_filter.Spt2-Myc_Wt_39C_Cl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt2-Myc_Wt_39C_Cl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_18_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt2-Myc_Wt_39C_Cl1_RS.c368e79350e555dd759254e0564ee1f5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt2-Myc_Wt_39C_Cl1_RS.c368e79350e555dd759254e0564ee1f5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS/Spt2-Myc_Wt_39C_Cl1_RS.sorted.bam \
  > alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1_RS/Spt2-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam
samtools_view_filter.Spt2-Myc_Wt_39C_Cl1_RS.c368e79350e555dd759254e0564ee1f5.mugqic.done
)
samtools_view_filter_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_19_JOB_ID: samtools_view_filter.Chd1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_WtCl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_19_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_WtCl2_RS.adf32e04467c75d6716c764db2d3947a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_WtCl2_RS.adf32e04467c75d6716c764db2d3947a.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS/Chd1-Myc_WtCl2_RS.sorted.bam \
  > alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2_RS/Chd1-Myc_WtCl2_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_WtCl2_RS.adf32e04467c75d6716c764db2d3947a.mugqic.done
)
samtools_view_filter_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_20_JOB_ID: samtools_view_filter.Chd1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Chd1-Myc_WtCl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_20_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Chd1-Myc_WtCl1_RS.8e7c21e10fb0641ebcd08cadeb5baa00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Chd1-Myc_WtCl1_RS.8e7c21e10fb0641ebcd08cadeb5baa00.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS/Chd1-Myc_WtCl1_RS.sorted.bam \
  > alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1_RS/Chd1-Myc_WtCl1_RS.sorted.filtered.bam
samtools_view_filter.Chd1-Myc_WtCl1_RS.8e7c21e10fb0641ebcd08cadeb5baa00.mugqic.done
)
samtools_view_filter_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_21_JOB_ID: samtools_view_filter.Iws1-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_WtCl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_21_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_WtCl2_RS.cdfebb243a01eb99090bee998530ecdd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_WtCl2_RS.cdfebb243a01eb99090bee998530ecdd.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS/Iws1-Myc_WtCl2_RS.sorted.bam \
  > alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2_RS/Iws1-Myc_WtCl2_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_WtCl2_RS.cdfebb243a01eb99090bee998530ecdd.mugqic.done
)
samtools_view_filter_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_22_JOB_ID: samtools_view_filter.Iws1-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Iws1-Myc_WtCl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_22_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Iws1-Myc_WtCl1_RS.72c10dcc7a77264e0d51c051c30c2edc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Iws1-Myc_WtCl1_RS.72c10dcc7a77264e0d51c051c30c2edc.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS/Iws1-Myc_WtCl1_RS.sorted.bam \
  > alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1_RS/Iws1-Myc_WtCl1_RS.sorted.filtered.bam
samtools_view_filter.Iws1-Myc_WtCl1_RS.72c10dcc7a77264e0d51c051c30c2edc.mugqic.done
)
samtools_view_filter_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_23_JOB_ID: samtools_view_filter.Spt6-Myc_WtCl2_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt6-Myc_WtCl2_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_23_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt6-Myc_WtCl2_RS.15e9894ca7c76e4175bb364bab67dc03.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt6-Myc_WtCl2_RS.15e9894ca7c76e4175bb364bab67dc03.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS/Spt6-Myc_WtCl2_RS.sorted.bam \
  > alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2_RS/Spt6-Myc_WtCl2_RS.sorted.filtered.bam
samtools_view_filter.Spt6-Myc_WtCl2_RS.15e9894ca7c76e4175bb364bab67dc03.mugqic.done
)
samtools_view_filter_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_24_JOB_ID: samtools_view_filter.Spt6-Myc_WtCl1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.Spt6-Myc_WtCl1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_24_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.Spt6-Myc_WtCl1_RS.a8a6ee7c3d5a902e1d9cdbbd49b7258c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.Spt6-Myc_WtCl1_RS.a8a6ee7c3d5a902e1d9cdbbd49b7258c.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS/Spt6-Myc_WtCl1_RS.sorted.bam \
  > alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1_RS/Spt6-Myc_WtCl1_RS.sorted.filtered.bam
samtools_view_filter.Spt6-Myc_WtCl1_RS.a8a6ee7c3d5a902e1d9cdbbd49b7258c.mugqic.done
)
samtools_view_filter_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_25_JOB_ID: samtools_view_filter.No-TAG_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.No-TAG_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_25_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.No-TAG_RS.086a7ddf22a8976f91a9e5985e891011.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.No-TAG_RS.086a7ddf22a8976f91a9e5985e891011.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/No-TAG/No-TAG_RS/No-TAG_RS.sorted.bam \
  > alignment/No-TAG/No-TAG_RS/No-TAG_RS.sorted.filtered.bam
samtools_view_filter.No-TAG_RS.086a7ddf22a8976f91a9e5985e891011.mugqic.done
)
samtools_view_filter_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_26_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID:$samtools_view_filter_11_JOB_ID:$samtools_view_filter_12_JOB_ID:$samtools_view_filter_13_JOB_ID:$samtools_view_filter_14_JOB_ID:$samtools_view_filter_15_JOB_ID:$samtools_view_filter_16_JOB_ID:$samtools_view_filter_17_JOB_ID:$samtools_view_filter_18_JOB_ID:$samtools_view_filter_19_JOB_ID:$samtools_view_filter_20_JOB_ID:$samtools_view_filter_21_JOB_ID:$samtools_view_filter_22_JOB_ID:$samtools_view_filter_23_JOB_ID:$samtools_view_filter_24_JOB_ID:$samtools_view_filter_25_JOB_ID
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
samtools_view_filter_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_dspt2Cl2.08695a95482a565cf662aa1051209819.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_dspt2Cl2.08695a95482a565cf662aa1051209819.mugqic.done'
mkdir -p alignment/Chd1-Myc_dspt2Cl2 && \
ln -s -f Chd1-Myc_dspt2Cl2_RS/Chd1-Myc_dspt2Cl2_RS.sorted.filtered.bam alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.merged.bam
symlink_readset_sample_bam.Chd1-Myc_dspt2Cl2.08695a95482a565cf662aa1051209819.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_dspt2Cl1.9e54c6f428008b780f092c083c45fdaf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_dspt2Cl1.9e54c6f428008b780f092c083c45fdaf.mugqic.done'
mkdir -p alignment/Chd1-Myc_dspt2Cl1 && \
ln -s -f Chd1-Myc_dspt2Cl1_RS/Chd1-Myc_dspt2Cl1_RS.sorted.filtered.bam alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.merged.bam
symlink_readset_sample_bam.Chd1-Myc_dspt2Cl1.9e54c6f428008b780f092c083c45fdaf.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_dspt2Cl2.b05db6778d851b8ee7a21a9e03dd9ccc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_dspt2Cl2.b05db6778d851b8ee7a21a9e03dd9ccc.mugqic.done'
mkdir -p alignment/Iws1-Myc_dspt2Cl2 && \
ln -s -f Iws1-Myc_dspt2Cl2_RS/Iws1-Myc_dspt2Cl2_RS.sorted.filtered.bam alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.merged.bam
symlink_readset_sample_bam.Iws1-Myc_dspt2Cl2.b05db6778d851b8ee7a21a9e03dd9ccc.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_dspt2Cl1.08d896e74c28874269c56f1265c372f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_dspt2Cl1.08d896e74c28874269c56f1265c372f9.mugqic.done'
mkdir -p alignment/Iws1-Myc_dspt2Cl1 && \
ln -s -f Iws1-Myc_dspt2Cl1_RS/Iws1-Myc_dspt2Cl1_RS.sorted.filtered.bam alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.merged.bam
symlink_readset_sample_bam.Iws1-Myc_dspt2Cl1.08d896e74c28874269c56f1265c372f9.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.Spt6-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt6-Myc_dspt2Cl2
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt6-Myc_dspt2Cl2.092145fad36710113d76940360665842.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt6-Myc_dspt2Cl2.092145fad36710113d76940360665842.mugqic.done'
mkdir -p alignment/Spt6-Myc_dspt2Cl2 && \
ln -s -f Spt6-Myc_dspt2Cl2_RS/Spt6-Myc_dspt2Cl2_RS.sorted.filtered.bam alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.merged.bam
symlink_readset_sample_bam.Spt6-Myc_dspt2Cl2.092145fad36710113d76940360665842.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.Spt6-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt6-Myc_dspt2Cl1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt6-Myc_dspt2Cl1.66532f6bf238d720701602a9a4e8e4b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt6-Myc_dspt2Cl1.66532f6bf238d720701602a9a4e8e4b5.mugqic.done'
mkdir -p alignment/Spt6-Myc_dspt2Cl1 && \
ln -s -f Spt6-Myc_dspt2Cl1_RS/Spt6-Myc_dspt2Cl1_RS.sorted.filtered.bam alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.merged.bam
symlink_readset_sample_bam.Spt6-Myc_dspt2Cl1.66532f6bf238d720701602a9a4e8e4b5.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done'
mkdir -p alignment/Chd1-Myc_spt6_39C_Cl2 && \
ln -s -f Chd1-Myc_spt6_39C_Cl2_RS/Chd1-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.merged.bam
symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl2.2d56b9777f6c5a88b9fdc5c07d956ee1.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl1.82e9a621f7e8a092c101fdbf541766a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl1.82e9a621f7e8a092c101fdbf541766a2.mugqic.done'
mkdir -p alignment/Chd1-Myc_spt6_39C_Cl1 && \
ln -s -f Chd1-Myc_spt6_39C_Cl1_RS/Chd1-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.merged.bam
symlink_readset_sample_bam.Chd1-Myc_spt6_39C_Cl1.82e9a621f7e8a092c101fdbf541766a2.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl2.a22249122153c2f387a13076ce2e830f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl2.a22249122153c2f387a13076ce2e830f.mugqic.done'
mkdir -p alignment/Iws1-Myc_spt6_39C_Cl2 && \
ln -s -f Iws1-Myc_spt6_39C_Cl2_RS/Iws1-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.merged.bam
symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl2.a22249122153c2f387a13076ce2e830f.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl1.c0b9b97cce97d82ee2e17ab541e91774.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl1.c0b9b97cce97d82ee2e17ab541e91774.mugqic.done'
mkdir -p alignment/Iws1-Myc_spt6_39C_Cl1 && \
ln -s -f Iws1-Myc_spt6_39C_Cl1_RS/Iws1-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.merged.bam
symlink_readset_sample_bam.Iws1-Myc_spt6_39C_Cl1.c0b9b97cce97d82ee2e17ab541e91774.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done'
mkdir -p alignment/Spt2-Myc_spt6_39C_Cl2 && \
ln -s -f Spt2-Myc_spt6_39C_Cl2_RS/Spt2-Myc_spt6_39C_Cl2_RS.sorted.filtered.bam alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.merged.bam
symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl2.43d5d22712411a73a141c5c48be205cc.mugqic.done
)
picard_merge_sam_files_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_12_JOB_ID: symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl1.398f1ce112e52dacf4c8b082c9684adc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl1.398f1ce112e52dacf4c8b082c9684adc.mugqic.done'
mkdir -p alignment/Spt2-Myc_spt6_39C_Cl1 && \
ln -s -f Spt2-Myc_spt6_39C_Cl1_RS/Spt2-Myc_spt6_39C_Cl1_RS.sorted.filtered.bam alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.merged.bam
symlink_readset_sample_bam.Spt2-Myc_spt6_39C_Cl1.398f1ce112e52dacf4c8b082c9684adc.mugqic.done
)
picard_merge_sam_files_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_13_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_13_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl2.4ee753ecb220148a512d122cdb277585.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl2.4ee753ecb220148a512d122cdb277585.mugqic.done'
mkdir -p alignment/Chd1-Myc_Wt_39C_Cl2 && \
ln -s -f Chd1-Myc_Wt_39C_Cl2_RS/Chd1-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.merged.bam
symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl2.4ee753ecb220148a512d122cdb277585.mugqic.done
)
picard_merge_sam_files_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_14_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_14_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done'
mkdir -p alignment/Chd1-Myc_Wt_39C_Cl1 && \
ln -s -f Chd1-Myc_Wt_39C_Cl1_RS/Chd1-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.merged.bam
symlink_readset_sample_bam.Chd1-Myc_Wt_39C_Cl1.b8300dd0add3aa6123d5e295b88caf38.mugqic.done
)
picard_merge_sam_files_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_15_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_15_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl2.e77e543a8c567f82c134feff79b8f0bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl2.e77e543a8c567f82c134feff79b8f0bd.mugqic.done'
mkdir -p alignment/Iws1-Myc_Wt_39C_Cl2 && \
ln -s -f Iws1-Myc_Wt_39C_Cl2_RS/Iws1-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.merged.bam
symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl2.e77e543a8c567f82c134feff79b8f0bd.mugqic.done
)
picard_merge_sam_files_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_16_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_16_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl1.f3c6fbbe6de6b974135d4dadaf4c7ec3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl1.f3c6fbbe6de6b974135d4dadaf4c7ec3.mugqic.done'
mkdir -p alignment/Iws1-Myc_Wt_39C_Cl1 && \
ln -s -f Iws1-Myc_Wt_39C_Cl1_RS/Iws1-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.merged.bam
symlink_readset_sample_bam.Iws1-Myc_Wt_39C_Cl1.f3c6fbbe6de6b974135d4dadaf4c7ec3.mugqic.done
)
picard_merge_sam_files_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_17_JOB_ID: symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$samtools_view_filter_17_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl2.946bf69493e641be083b3b44b0e70422.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl2.946bf69493e641be083b3b44b0e70422.mugqic.done'
mkdir -p alignment/Spt2-Myc_Wt_39C_Cl2 && \
ln -s -f Spt2-Myc_Wt_39C_Cl2_RS/Spt2-Myc_Wt_39C_Cl2_RS.sorted.filtered.bam alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.merged.bam
symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl2.946bf69493e641be083b3b44b0e70422.mugqic.done
)
picard_merge_sam_files_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_18_JOB_ID: symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$samtools_view_filter_18_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl1.3e42f92f122c9b2d192f7648db9aca16.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl1.3e42f92f122c9b2d192f7648db9aca16.mugqic.done'
mkdir -p alignment/Spt2-Myc_Wt_39C_Cl1 && \
ln -s -f Spt2-Myc_Wt_39C_Cl1_RS/Spt2-Myc_Wt_39C_Cl1_RS.sorted.filtered.bam alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.merged.bam
symlink_readset_sample_bam.Spt2-Myc_Wt_39C_Cl1.3e42f92f122c9b2d192f7648db9aca16.mugqic.done
)
picard_merge_sam_files_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_19_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_WtCl2
JOB_DEPENDENCIES=$samtools_view_filter_19_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_WtCl2.1ebf6bc21f19a010b021191f5f081a91.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_WtCl2.1ebf6bc21f19a010b021191f5f081a91.mugqic.done'
mkdir -p alignment/Chd1-Myc_WtCl2 && \
ln -s -f Chd1-Myc_WtCl2_RS/Chd1-Myc_WtCl2_RS.sorted.filtered.bam alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.merged.bam
symlink_readset_sample_bam.Chd1-Myc_WtCl2.1ebf6bc21f19a010b021191f5f081a91.mugqic.done
)
picard_merge_sam_files_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_20_JOB_ID: symlink_readset_sample_bam.Chd1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Chd1-Myc_WtCl1
JOB_DEPENDENCIES=$samtools_view_filter_20_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Chd1-Myc_WtCl1.8fee36755ee2cfa5b712d10b8ba5e23d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Chd1-Myc_WtCl1.8fee36755ee2cfa5b712d10b8ba5e23d.mugqic.done'
mkdir -p alignment/Chd1-Myc_WtCl1 && \
ln -s -f Chd1-Myc_WtCl1_RS/Chd1-Myc_WtCl1_RS.sorted.filtered.bam alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.merged.bam
symlink_readset_sample_bam.Chd1-Myc_WtCl1.8fee36755ee2cfa5b712d10b8ba5e23d.mugqic.done
)
picard_merge_sam_files_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_21_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_WtCl2
JOB_DEPENDENCIES=$samtools_view_filter_21_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_WtCl2.b67983cf979150bf2c604d62fde72a20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_WtCl2.b67983cf979150bf2c604d62fde72a20.mugqic.done'
mkdir -p alignment/Iws1-Myc_WtCl2 && \
ln -s -f Iws1-Myc_WtCl2_RS/Iws1-Myc_WtCl2_RS.sorted.filtered.bam alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.merged.bam
symlink_readset_sample_bam.Iws1-Myc_WtCl2.b67983cf979150bf2c604d62fde72a20.mugqic.done
)
picard_merge_sam_files_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_22_JOB_ID: symlink_readset_sample_bam.Iws1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Iws1-Myc_WtCl1
JOB_DEPENDENCIES=$samtools_view_filter_22_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Iws1-Myc_WtCl1.692d9e3f6086b46c1430ae0ccb8a38d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Iws1-Myc_WtCl1.692d9e3f6086b46c1430ae0ccb8a38d6.mugqic.done'
mkdir -p alignment/Iws1-Myc_WtCl1 && \
ln -s -f Iws1-Myc_WtCl1_RS/Iws1-Myc_WtCl1_RS.sorted.filtered.bam alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.merged.bam
symlink_readset_sample_bam.Iws1-Myc_WtCl1.692d9e3f6086b46c1430ae0ccb8a38d6.mugqic.done
)
picard_merge_sam_files_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_23_JOB_ID: symlink_readset_sample_bam.Spt6-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt6-Myc_WtCl2
JOB_DEPENDENCIES=$samtools_view_filter_23_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt6-Myc_WtCl2.296a03a111350a7499b776c1df0833cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt6-Myc_WtCl2.296a03a111350a7499b776c1df0833cc.mugqic.done'
mkdir -p alignment/Spt6-Myc_WtCl2 && \
ln -s -f Spt6-Myc_WtCl2_RS/Spt6-Myc_WtCl2_RS.sorted.filtered.bam alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.merged.bam
symlink_readset_sample_bam.Spt6-Myc_WtCl2.296a03a111350a7499b776c1df0833cc.mugqic.done
)
picard_merge_sam_files_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_24_JOB_ID: symlink_readset_sample_bam.Spt6-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.Spt6-Myc_WtCl1
JOB_DEPENDENCIES=$samtools_view_filter_24_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.Spt6-Myc_WtCl1.86c52104f445c463684d8f2662355982.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.Spt6-Myc_WtCl1.86c52104f445c463684d8f2662355982.mugqic.done'
mkdir -p alignment/Spt6-Myc_WtCl1 && \
ln -s -f Spt6-Myc_WtCl1_RS/Spt6-Myc_WtCl1_RS.sorted.filtered.bam alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.merged.bam
symlink_readset_sample_bam.Spt6-Myc_WtCl1.86c52104f445c463684d8f2662355982.mugqic.done
)
picard_merge_sam_files_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_25_JOB_ID: symlink_readset_sample_bam.No-TAG
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.No-TAG
JOB_DEPENDENCIES=$samtools_view_filter_25_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.No-TAG.1d11f6cab239a528a1556b9b3117d813.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.No-TAG.1d11f6cab239a528a1556b9b3117d813.mugqic.done'
mkdir -p alignment/No-TAG && \
ln -s -f No-TAG_RS/No-TAG_RS.sorted.filtered.bam alignment/No-TAG/No-TAG.merged.bam
symlink_readset_sample_bam.No-TAG.1d11f6cab239a528a1556b9b3117d813.mugqic.done
)
picard_merge_sam_files_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.Chd1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_dspt2Cl2.135d5d687d580ee6abe98a7f0b17e4af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_dspt2Cl2.135d5d687d580ee6abe98a7f0b17e4af.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.merged.bam \
  OUTPUT=alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_dspt2Cl2.135d5d687d580ee6abe98a7f0b17e4af.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.Chd1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_dspt2Cl1.b5e4f77b54eb544b323720dcb3081634.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_dspt2Cl1.b5e4f77b54eb544b323720dcb3081634.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.merged.bam \
  OUTPUT=alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_dspt2Cl1.b5e4f77b54eb544b323720dcb3081634.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.Iws1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_dspt2Cl2.82ce42a8205be61204765064719210fe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_dspt2Cl2.82ce42a8205be61204765064719210fe.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.merged.bam \
  OUTPUT=alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_dspt2Cl2.82ce42a8205be61204765064719210fe.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.Iws1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_dspt2Cl1.0a655ca0bbc17648359b545652e0e297.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_dspt2Cl1.0a655ca0bbc17648359b545652e0e297.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.merged.bam \
  OUTPUT=alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_dspt2Cl1.0a655ca0bbc17648359b545652e0e297.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.Spt6-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt6-Myc_dspt2Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt6-Myc_dspt2Cl2.e5292aa2ead8c79bba43b77abdc5db9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt6-Myc_dspt2Cl2.e5292aa2ead8c79bba43b77abdc5db9c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.merged.bam \
  OUTPUT=alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt6-Myc_dspt2Cl2.e5292aa2ead8c79bba43b77abdc5db9c.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.Spt6-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt6-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt6-Myc_dspt2Cl1.cc1a45e26be1107cde92c583c07fb77f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt6-Myc_dspt2Cl1.cc1a45e26be1107cde92c583c07fb77f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.merged.bam \
  OUTPUT=alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt6-Myc_dspt2Cl1.cc1a45e26be1107cde92c583c07fb77f.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
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
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl1.6286c53efa10518c8c7b6a7bd1b7f8ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl1.6286c53efa10518c8c7b6a7bd1b7f8ae.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.merged.bam \
  OUTPUT=alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_spt6_39C_Cl1.6286c53efa10518c8c7b6a7bd1b7f8ae.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl2.f40150da5b49f53da9d28120c9e6ad1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl2.f40150da5b49f53da9d28120c9e6ad1b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.merged.bam \
  OUTPUT=alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl2.f40150da5b49f53da9d28120c9e6ad1b.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl1.a63d8944afa35937376151329c011257.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl1.a63d8944afa35937376151329c011257.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.merged.bam \
  OUTPUT=alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_spt6_39C_Cl1.a63d8944afa35937376151329c011257.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
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
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_12_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl1.10f0604faba97a3d49483432dae9743c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl1.10f0604faba97a3d49483432dae9743c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.merged.bam \
  OUTPUT=alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt2-Myc_spt6_39C_Cl1.10f0604faba97a3d49483432dae9743c.mugqic.done
)
picard_mark_duplicates_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_13_JOB_ID: picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_13_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl2.8c7a11e5b18c72f3b7df6410a8ffc6eb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl2.8c7a11e5b18c72f3b7df6410a8ffc6eb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.merged.bam \
  OUTPUT=alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl2.8c7a11e5b18c72f3b7df6410a8ffc6eb.mugqic.done
)
picard_mark_duplicates_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_14_JOB_ID: picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_14_JOB_ID
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
picard_mark_duplicates_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_15_JOB_ID: picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_15_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl2.9ede7a427bd8871d189a5ec26e827abb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl2.9ede7a427bd8871d189a5ec26e827abb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.merged.bam \
  OUTPUT=alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl2.9ede7a427bd8871d189a5ec26e827abb.mugqic.done
)
picard_mark_duplicates_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_16_JOB_ID: picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_16_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl1.c9e7d9ede1965bf705afb095fadc5bfb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl1.c9e7d9ede1965bf705afb095fadc5bfb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.merged.bam \
  OUTPUT=alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_Wt_39C_Cl1.c9e7d9ede1965bf705afb095fadc5bfb.mugqic.done
)
picard_mark_duplicates_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_17_JOB_ID: picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$picard_merge_sam_files_17_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl2.2b27efaa3a75cd0535c21ab750e5b7ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl2.2b27efaa3a75cd0535c21ab750e5b7ed.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.merged.bam \
  OUTPUT=alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam \
  METRICS_FILE=alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl2.2b27efaa3a75cd0535c21ab750e5b7ed.mugqic.done
)
picard_mark_duplicates_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_18_JOB_ID: picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_merge_sam_files_18_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl1.fcf2765694188141c4f248777dbc5d48.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl1.fcf2765694188141c4f248777dbc5d48.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.merged.bam \
  OUTPUT=alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam \
  METRICS_FILE=alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt2-Myc_Wt_39C_Cl1.fcf2765694188141c4f248777dbc5d48.mugqic.done
)
picard_mark_duplicates_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_19_JOB_ID: picard_mark_duplicates.Chd1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_WtCl2
JOB_DEPENDENCIES=$picard_merge_sam_files_19_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_WtCl2.0fa3bbc2ddc63abea045321cd8b8152a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_WtCl2.0fa3bbc2ddc63abea045321cd8b8152a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.merged.bam \
  OUTPUT=alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_WtCl2.0fa3bbc2ddc63abea045321cd8b8152a.mugqic.done
)
picard_mark_duplicates_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_20_JOB_ID: picard_mark_duplicates.Chd1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Chd1-Myc_WtCl1
JOB_DEPENDENCIES=$picard_merge_sam_files_20_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Chd1-Myc_WtCl1.8a7f65bad53f6c4c4e27bd6e60c3f0c7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Chd1-Myc_WtCl1.8a7f65bad53f6c4c4e27bd6e60c3f0c7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.merged.bam \
  OUTPUT=alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam \
  METRICS_FILE=alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Chd1-Myc_WtCl1.8a7f65bad53f6c4c4e27bd6e60c3f0c7.mugqic.done
)
picard_mark_duplicates_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_21_JOB_ID: picard_mark_duplicates.Iws1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_WtCl2
JOB_DEPENDENCIES=$picard_merge_sam_files_21_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_WtCl2.13173625fd408c7990b0e25fb7291dfe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_WtCl2.13173625fd408c7990b0e25fb7291dfe.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.merged.bam \
  OUTPUT=alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_WtCl2.13173625fd408c7990b0e25fb7291dfe.mugqic.done
)
picard_mark_duplicates_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_22_JOB_ID: picard_mark_duplicates.Iws1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Iws1-Myc_WtCl1
JOB_DEPENDENCIES=$picard_merge_sam_files_22_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Iws1-Myc_WtCl1.5d43f3e4f8106e4a689b6193292a24c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Iws1-Myc_WtCl1.5d43f3e4f8106e4a689b6193292a24c6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.merged.bam \
  OUTPUT=alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam \
  METRICS_FILE=alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Iws1-Myc_WtCl1.5d43f3e4f8106e4a689b6193292a24c6.mugqic.done
)
picard_mark_duplicates_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_23_JOB_ID: picard_mark_duplicates.Spt6-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt6-Myc_WtCl2
JOB_DEPENDENCIES=$picard_merge_sam_files_23_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt6-Myc_WtCl2.c16f83ce98850d116df9f34f1570f28b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt6-Myc_WtCl2.c16f83ce98850d116df9f34f1570f28b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.merged.bam \
  OUTPUT=alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam \
  METRICS_FILE=alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt6-Myc_WtCl2.c16f83ce98850d116df9f34f1570f28b.mugqic.done
)
picard_mark_duplicates_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_24_JOB_ID: picard_mark_duplicates.Spt6-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.Spt6-Myc_WtCl1
JOB_DEPENDENCIES=$picard_merge_sam_files_24_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.Spt6-Myc_WtCl1.6a709bb88683dbb03dc4d55a55e02560.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.Spt6-Myc_WtCl1.6a709bb88683dbb03dc4d55a55e02560.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.merged.bam \
  OUTPUT=alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam \
  METRICS_FILE=alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.Spt6-Myc_WtCl1.6a709bb88683dbb03dc4d55a55e02560.mugqic.done
)
picard_mark_duplicates_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_25_JOB_ID: picard_mark_duplicates.No-TAG
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.No-TAG
JOB_DEPENDENCIES=$picard_merge_sam_files_25_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.No-TAG.6ac1d07054540cb1d84f513fc4caa6d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.No-TAG.6ac1d07054540cb1d84f513fc4caa6d9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/No-TAG/No-TAG.merged.bam \
  OUTPUT=alignment/No-TAG/No-TAG.sorted.dup.bam \
  METRICS_FILE=alignment/No-TAG/No-TAG.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.No-TAG.6ac1d07054540cb1d84f513fc4caa6d9.mugqic.done
)
picard_mark_duplicates_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_26_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID:$picard_mark_duplicates_15_JOB_ID:$picard_mark_duplicates_16_JOB_ID:$picard_mark_duplicates_17_JOB_ID:$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_19_JOB_ID:$picard_mark_duplicates_20_JOB_ID:$picard_mark_duplicates_21_JOB_ID:$picard_mark_duplicates_22_JOB_ID:$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID
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
picard_mark_duplicates_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID:$picard_mark_duplicates_15_JOB_ID:$picard_mark_duplicates_16_JOB_ID:$picard_mark_duplicates_17_JOB_ID:$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_19_JOB_ID:$picard_mark_duplicates_20_JOB_ID:$picard_mark_duplicates_21_JOB_ID:$picard_mark_duplicates_22_JOB_ID:$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID
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
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.Chd1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.Iws1-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_dspt2Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.Iws1-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.Spt6-Myc_dspt2Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_dspt2Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.Spt6-Myc_dspt2Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_dspt2Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
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
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_spt6_39C_Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_spt6_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
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
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_13_JOB_ID: homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_14_JOB_ID: homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_14_JOB_ID
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
JOB_DEPENDENCIES=$picard_mark_duplicates_15_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_16_JOB_ID: homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_16_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_17_JOB_ID: homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl2
JOB_DEPENDENCIES=$picard_mark_duplicates_17_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_18_JOB_ID: homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt2-Myc_Wt_39C_Cl1
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_19_JOB_ID: homer_make_tag_directory.Chd1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_WtCl2
JOB_DEPENDENCIES=$picard_mark_duplicates_19_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_20_JOB_ID: homer_make_tag_directory.Chd1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Chd1-Myc_WtCl1
JOB_DEPENDENCIES=$picard_mark_duplicates_20_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_21_JOB_ID: homer_make_tag_directory.Iws1-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_WtCl2
JOB_DEPENDENCIES=$picard_mark_duplicates_21_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_22_JOB_ID: homer_make_tag_directory.Iws1-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Iws1-Myc_WtCl1
JOB_DEPENDENCIES=$picard_mark_duplicates_22_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_23_JOB_ID: homer_make_tag_directory.Spt6-Myc_WtCl2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_WtCl2
JOB_DEPENDENCIES=$picard_mark_duplicates_23_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_24_JOB_ID: homer_make_tag_directory.Spt6-Myc_WtCl1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.Spt6-Myc_WtCl1
JOB_DEPENDENCIES=$picard_mark_duplicates_24_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_25_JOB_ID: homer_make_tag_directory.No-TAG
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.No-TAG
JOB_DEPENDENCIES=$picard_mark_duplicates_25_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.Chd1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-dspt2-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-dspt2-Ctl.31adb995c5ae6966221f177552a86681.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-dspt2-Ctl.31adb995c5ae6966221f177552a86681.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-dspt2-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl \
  >& peak_call/Chd1-dspt2-Ctl/Chd1-dspt2-Ctl.diag.macs.out
macs2_callpeak.Chd1-dspt2-Ctl.31adb995c5ae6966221f177552a86681.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.Iws1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-dspt2-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-dspt2-Ctl.2b908ce38f5b6ffc5961a3571f3d7dd6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-dspt2-Ctl.2b908ce38f5b6ffc5961a3571f3d7dd6.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-dspt2-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl \
  >& peak_call/Iws1-dspt2-Ctl/Iws1-dspt2-Ctl.diag.macs.out
macs2_callpeak.Iws1-dspt2-Ctl.2b908ce38f5b6ffc5961a3571f3d7dd6.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.Spt6-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt6-dspt2-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt6-dspt2-Ctl.e8ca49dc7b37836e3326c1832ea65062.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt6-dspt2-Ctl.e8ca49dc7b37836e3326c1832ea65062.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt6-dspt2-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl \
  >& peak_call/Spt6-dspt2-Ctl/Spt6-dspt2-Ctl.diag.macs.out
macs2_callpeak.Spt6-dspt2-Ctl.e8ca49dc7b37836e3326c1832ea65062.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-spt6-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_25_JOB_ID
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
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.Iws1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-spt6-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-spt6-39C.ed8e0397c1cde906c45f51f7d98c7884.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-spt6-39C.ed8e0397c1cde906c45f51f7d98c7884.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-spt6-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-spt6-39C/Iws1-spt6-39C \
  >& peak_call/Iws1-spt6-39C/Iws1-spt6-39C.diag.macs.out
macs2_callpeak.Iws1-spt6-39C.ed8e0397c1cde906c45f51f7d98c7884.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt2-spt6-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_25_JOB_ID
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
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.Chd1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-WT-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID:$picard_mark_duplicates_25_JOB_ID
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
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.Iws1-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-WT-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_15_JOB_ID:$picard_mark_duplicates_16_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-WT-39C.359dc26f024ecbc1f5be4c9768ec1578.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-WT-39C.359dc26f024ecbc1f5be4c9768ec1578.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-WT-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-WT-39C/Iws1-WT-39C \
  >& peak_call/Iws1-WT-39C/Iws1-WT-39C.diag.macs.out
macs2_callpeak.Iws1-WT-39C.359dc26f024ecbc1f5be4c9768ec1578.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.Spt2-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt2-WT-39C
JOB_DEPENDENCIES=$picard_mark_duplicates_17_JOB_ID:$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt2-WT-39C.41c725639e6d09c3f189062be55cdff0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt2-WT-39C.41c725639e6d09c3f189062be55cdff0.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt2-WT-39C && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt2-WT-39C/Spt2-WT-39C \
  >& peak_call/Spt2-WT-39C/Spt2-WT-39C.diag.macs.out
macs2_callpeak.Spt2-WT-39C.41c725639e6d09c3f189062be55cdff0.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.Chd1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-WT-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_19_JOB_ID:$picard_mark_duplicates_20_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-WT-Ctl.f01a33bdea0bbf703590d67f4a3e3ba1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-WT-Ctl.f01a33bdea0bbf703590d67f4a3e3ba1.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-WT-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam \
  alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-WT-Ctl/Chd1-WT-Ctl \
  >& peak_call/Chd1-WT-Ctl/Chd1-WT-Ctl.diag.macs.out
macs2_callpeak.Chd1-WT-Ctl.f01a33bdea0bbf703590d67f4a3e3ba1.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.Iws1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-WT-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_21_JOB_ID:$picard_mark_duplicates_22_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-WT-Ctl.f5ed9226a17cbb7af561125a00aa384c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-WT-Ctl.f5ed9226a17cbb7af561125a00aa384c.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-WT-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam \
  alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-WT-Ctl/Iws1-WT-Ctl \
  >& peak_call/Iws1-WT-Ctl/Iws1-WT-Ctl.diag.macs.out
macs2_callpeak.Iws1-WT-Ctl.f5ed9226a17cbb7af561125a00aa384c.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak.Spt6-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt6-WT-Ctl
JOB_DEPENDENCIES=$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt6-WT-Ctl.cda25e2a11efea5db5facf81478d0c4d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt6-WT-Ctl.cda25e2a11efea5db5facf81478d0c4d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt6-WT-Ctl && \
macs2 callpeak --format BAM --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam \
  alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt6-WT-Ctl/Spt6-WT-Ctl \
  >& peak_call/Spt6-WT-Ctl/Spt6-WT-Ctl.diag.macs.out
macs2_callpeak.Spt6-WT-Ctl.cda25e2a11efea5db5facf81478d0c4d.mugqic.done
)
macs2_callpeak_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_10_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_12_JOB_ID
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
macs2_callpeak_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.Chd1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-dspt2-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.Iws1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-dspt2-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.Spt6-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt6-dspt2-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt2-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.Spt2-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt2-WT-39C
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks.Chd1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Chd1-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_11_JOB_ID: homer_annotate_peaks.Iws1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Iws1-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_12_JOB_ID: homer_annotate_peaks.Spt6-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.Spt6-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_12_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.Iws1-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-dspt2-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.Spt6-dspt2-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt6-dspt2-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.Chd1-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.Spt2-spt6-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt2-spt6-39C
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
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
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome.Spt2-WT-39C
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt2-WT-39C
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_10_JOB_ID: homer_find_motifs_genome.Chd1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Chd1-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_11_JOB_ID: homer_find_motifs_genome.Iws1-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Iws1-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_12_JOB_ID: homer_find_motifs_genome.Spt6-WT-Ctl
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.Spt6-WT-Ctl
JOB_DEPENDENCIES=$macs2_callpeak_12_JOB_ID
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n04&ip=10.241.129.14&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=25" --quiet --output-document=/dev/null

