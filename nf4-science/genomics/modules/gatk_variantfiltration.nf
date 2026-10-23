#!/usr/bin/env nextflow

/*
 * Apply hard filters to cohort-level variant calls
 */
process GATK_VARIANTFILTRATION {

    container 'community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867'

    input:
    path input_vcf
    path input_vcf_idx
    val cohort_name
    path ref_fasta
    path ref_index
    path ref_dict

    output:
    path "${cohort_name}.filtered.vcf"     , emit: vcf
    path "${cohort_name}.filtered.vcf.idx" , emit: idx

    script:
    """
    gatk VariantFiltration \
        -R ${ref_fasta} \
        -V ${input_vcf} \
        -O ${cohort_name}.filtered.vcf \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "hard_filter"
    """
}
