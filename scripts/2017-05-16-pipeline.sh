#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-26T21:02:16
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_trimmomatic_stats: 0 job... skipping
#   bwa_mem_picard_sort_sam: 0 job... skipping
#   samtools_view_filter: 0 job... skipping
#   picard_merge_sam_files: 0 job... skipping
#   picard_mark_duplicates: 0 job... skipping
#   metrics: 0 job... skipping
#   macs2_callpeak: 13 jobs
#   TOTAL: 13 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/project/eav-760-aa/Chaperone/output/pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.Chd1-dspt2-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-dspt2-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-dspt2-Ctl_B.5efadc841a6cc6934b7f44c1e0f0783b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-dspt2-Ctl_B.5efadc841a6cc6934b7f44c1e0f0783b.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-dspt2-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_dspt2Cl2/Chd1-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_dspt2Cl1/Chd1-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-dspt2-Ctl_B/Chd1-dspt2-Ctl_B \
  >& peak_call/Chd1-dspt2-Ctl_B/Chd1-dspt2-Ctl_B.diag.macs.out
macs2_callpeak.Chd1-dspt2-Ctl_B.5efadc841a6cc6934b7f44c1e0f0783b.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.Iws1-dspt2-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-dspt2-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-dspt2-Ctl_B.4a116067b63366a6068a9f5346ece2af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-dspt2-Ctl_B.4a116067b63366a6068a9f5346ece2af.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-dspt2-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_dspt2Cl2/Iws1-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_dspt2Cl1/Iws1-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-dspt2-Ctl_B/Iws1-dspt2-Ctl_B \
  >& peak_call/Iws1-dspt2-Ctl_B/Iws1-dspt2-Ctl_B.diag.macs.out
macs2_callpeak.Iws1-dspt2-Ctl_B.4a116067b63366a6068a9f5346ece2af.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.Spt6-dspt2-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt6-dspt2-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt6-dspt2-Ctl_B.d756f885d5bb4fcd95d17f255820b0fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt6-dspt2-Ctl_B.d756f885d5bb4fcd95d17f255820b0fc.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt6-dspt2-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt6-Myc_dspt2Cl2/Spt6-Myc_dspt2Cl2.sorted.dup.bam \
  alignment/Spt6-Myc_dspt2Cl1/Spt6-Myc_dspt2Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt6-dspt2-Ctl_B/Spt6-dspt2-Ctl_B \
  >& peak_call/Spt6-dspt2-Ctl_B/Spt6-dspt2-Ctl_B.diag.macs.out
macs2_callpeak.Spt6-dspt2-Ctl_B.d756f885d5bb4fcd95d17f255820b0fc.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.Chd1-spt6-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-spt6-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-spt6-39C_B.badc0e5e9af2c8f23d9544d1156becb3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-spt6-39C_B.badc0e5e9af2c8f23d9544d1156becb3.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-spt6-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_spt6_39C_Cl2/Chd1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_spt6_39C_Cl1/Chd1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-spt6-39C_B/Chd1-spt6-39C_B \
  >& peak_call/Chd1-spt6-39C_B/Chd1-spt6-39C_B.diag.macs.out
macs2_callpeak.Chd1-spt6-39C_B.badc0e5e9af2c8f23d9544d1156becb3.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.Iws1-spt6-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-spt6-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-spt6-39C_B.3889aa34145e5d547909364549115433.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-spt6-39C_B.3889aa34145e5d547909364549115433.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-spt6-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_spt6_39C_Cl2/Iws1-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_spt6_39C_Cl1/Iws1-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-spt6-39C_B/Iws1-spt6-39C_B \
  >& peak_call/Iws1-spt6-39C_B/Iws1-spt6-39C_B.diag.macs.out
macs2_callpeak.Iws1-spt6-39C_B.3889aa34145e5d547909364549115433.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.Spt2-spt6-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt2-spt6-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt2-spt6-39C_B.9e6eec97952f54f8f0643622fba3ceb1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt2-spt6-39C_B.9e6eec97952f54f8f0643622fba3ceb1.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt2-spt6-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt2-Myc_spt6_39C_Cl2/Spt2-Myc_spt6_39C_Cl2.sorted.dup.bam \
  alignment/Spt2-Myc_spt6_39C_Cl1/Spt2-Myc_spt6_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt2-spt6-39C_B/Spt2-spt6-39C_B \
  >& peak_call/Spt2-spt6-39C_B/Spt2-spt6-39C_B.diag.macs.out
macs2_callpeak.Spt2-spt6-39C_B.9e6eec97952f54f8f0643622fba3ceb1.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.Chd1-WT-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-WT-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-WT-39C_B.0e44c618b332cc4fe0839a64f6501a64.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-WT-39C_B.0e44c618b332cc4fe0839a64f6501a64.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-WT-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_Wt_39C_Cl2/Chd1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Chd1-Myc_Wt_39C_Cl1/Chd1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-WT-39C_B/Chd1-WT-39C_B \
  >& peak_call/Chd1-WT-39C_B/Chd1-WT-39C_B.diag.macs.out
macs2_callpeak.Chd1-WT-39C_B.0e44c618b332cc4fe0839a64f6501a64.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.Iws1-WT-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-WT-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-WT-39C_B.8169458d4496173629b6083c8359acb0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-WT-39C_B.8169458d4496173629b6083c8359acb0.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-WT-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_Wt_39C_Cl2/Iws1-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Iws1-Myc_Wt_39C_Cl1/Iws1-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-WT-39C_B/Iws1-WT-39C_B \
  >& peak_call/Iws1-WT-39C_B/Iws1-WT-39C_B.diag.macs.out
macs2_callpeak.Iws1-WT-39C_B.8169458d4496173629b6083c8359acb0.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.Spt2-WT-39C_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt2-WT-39C_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt2-WT-39C_B.e4ca9726325b2cc31bb70cc31e94c270.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt2-WT-39C_B.e4ca9726325b2cc31bb70cc31e94c270.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt2-WT-39C_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt2-Myc_Wt_39C_Cl2/Spt2-Myc_Wt_39C_Cl2.sorted.dup.bam \
  alignment/Spt2-Myc_Wt_39C_Cl1/Spt2-Myc_Wt_39C_Cl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt2-WT-39C_B/Spt2-WT-39C_B \
  >& peak_call/Spt2-WT-39C_B/Spt2-WT-39C_B.diag.macs.out
macs2_callpeak.Spt2-WT-39C_B.e4ca9726325b2cc31bb70cc31e94c270.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.Chd1-WT-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Chd1-WT-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Chd1-WT-Ctl_B.581147a41b1d77f98b8a2740ad349d32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Chd1-WT-Ctl_B.581147a41b1d77f98b8a2740ad349d32.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Chd1-WT-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Chd1-Myc_WtCl2/Chd1-Myc_WtCl2.sorted.dup.bam \
  alignment/Chd1-Myc_WtCl1/Chd1-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Chd1-WT-Ctl_B/Chd1-WT-Ctl_B \
  >& peak_call/Chd1-WT-Ctl_B/Chd1-WT-Ctl_B.diag.macs.out
macs2_callpeak.Chd1-WT-Ctl_B.581147a41b1d77f98b8a2740ad349d32.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.Iws1-WT-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Iws1-WT-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Iws1-WT-Ctl_B.67a362dc3136b7d3585df9b9c303e347.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Iws1-WT-Ctl_B.67a362dc3136b7d3585df9b9c303e347.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Iws1-WT-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Iws1-Myc_WtCl2/Iws1-Myc_WtCl2.sorted.dup.bam \
  alignment/Iws1-Myc_WtCl1/Iws1-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Iws1-WT-Ctl_B/Iws1-WT-Ctl_B \
  >& peak_call/Iws1-WT-Ctl_B/Iws1-WT-Ctl_B.diag.macs.out
macs2_callpeak.Iws1-WT-Ctl_B.67a362dc3136b7d3585df9b9c303e347.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak.Spt6-WT-Ctl_B
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.Spt6-WT-Ctl_B
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.Spt6-WT-Ctl_B.6212afb668b0d23dcec0667248c8fa6e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.Spt6-WT-Ctl_B.6212afb668b0d23dcec0667248c8fa6e.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/Spt6-WT-Ctl_B && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 9725684.0 \
  --treatment \
  alignment/Spt6-Myc_WtCl2/Spt6-Myc_WtCl2.sorted.dup.bam \
  alignment/Spt6-Myc_WtCl1/Spt6-Myc_WtCl1.sorted.dup.bam \
  --control \
  alignment/No-TAG/No-TAG.sorted.dup.bam \
  --name peak_call/Spt6-WT-Ctl_B/Spt6-WT-Ctl_B \
  >& peak_call/Spt6-WT-Ctl_B/Spt6-WT-Ctl_B.diag.macs.out
macs2_callpeak.Spt6-WT-Ctl_B.6212afb668b0d23dcec0667248c8fa6e.mugqic.done
)
macs2_callpeak_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_10_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.d08a7691a3dd12a8cd3ae35f04fafc4c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.d08a7691a3dd12a8cd3ae35f04fafc4c.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in Chd1-dspt2-Ctl_B Iws1-dspt2-Ctl_B Spt6-dspt2-Ctl_B Chd1-spt6-39C_B Iws1-spt6-39C_B Spt2-spt6-39C_B Chd1-WT-39C_B Iws1-WT-39C_B Spt2-WT-39C_B Chd1-WT-Ctl_B Iws1-WT-Ctl_B Spt6-WT-Ctl_B
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.d08a7691a3dd12a8cd3ae35f04fafc4c.mugqic.done
)
macs2_callpeak_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n03&ip=10.241.129.13&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,macs2_callpeak&samples=25" --quiet --output-document=/dev/null

