#!/usr/bin/env nextflow

/*
 * Generate alignment statistics for a BAM file
 */
process SAMTOOLS_FLAGSTAT {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    tuple path(input_bam), path(input_bam_index)

    output:
    path "${input_bam}.flagstat.txt"

    script:
    """
    samtools flagstat '$input_bam' > '${input_bam}.flagstat.txt'
    """
}
