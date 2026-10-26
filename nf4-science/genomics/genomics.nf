#!/usr/bin/env nextflow

// Module INCLUDE statements
include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
include { GATK_VARIANTFILTRATION } from './modules/gatk_variantfiltration.nf'

/*
 * Pipeline parameters
 */
params {
    // Primary input (file of input files, one per line)
    input: Path

    // Accessory files
    reference: Path
    reference_index: Path
    reference_dict: Path
    intervals: Path

    // Base name for final output file
    cohort_name: String
}

workflow {

    main:
    // Create input channel from a CSV file listing input file paths
    reads_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> file(row.reads_bam) }

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )

    // Generate alignment statistics for each indexed BAM file
    SAMTOOLS_FLAGSTAT(SAMTOOLS_INDEX.out)

    // Collect variant calling outputs across samples
    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

    // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )

    // Apply hard filters to the joint-called variants
    GATK_VARIANTFILTRATION(
        GATK_JOINTGENOTYPING.out.vcf,
        GATK_JOINTGENOTYPING.out.idx,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )

    publish:
    indexed_bam = SAMTOOLS_INDEX.out
    gvcf = GATK_HAPLOTYPECALLER.out.vcf
    gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    filtered_vcf = GATK_VARIANTFILTRATION.out.vcf
    filtered_vcf_idx = GATK_VARIANTFILTRATION.out.idx
    flagstat = SAMTOOLS_FLAGSTAT.out
}

output {
    indexed_bam {
        path 'indexed_bam'
    }
    gvcf {
        path 'gvcf'
    }
    gvcf_idx {
        path 'gvcf'
    }
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
    filtered_vcf {
        path '.'
    }
    filtered_vcf_idx {
        path '.'
    }
    flagstat {
        path 'stats'
    }
}
