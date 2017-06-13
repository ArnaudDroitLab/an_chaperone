export RAP_ID="eav-760-aa"

mkdir -p output/pipeline

$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.py -s '1-8,12' \
    -l debug \
    -r raw/readset.txt \
    -d raw/design_B.txt \
    -o output/pipeline \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.guillimin.ini \
        $MUGQIC_PIPELINES_HOME/resources/genomes/config/Saccharomyces_cerevisiae.R64-1-1.ini \
        input/chipseq.numpy.bug.ini        