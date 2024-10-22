#!/usr/bin/env nextflow

/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 'nf-after-sarek' - A Nextflow to get gVCF filtered per sample instead of joint calls (only variants)
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Judith Ballesteros Villascán
 GitHub: https://github.com/jbv2/after-sarek
 ----------------------------------------------------------------------------------------
 */

/* 
 Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */ 

params.input_tsv = null  // This will hold the TSV file path from the command line

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_FILTER   } from './modules/local/bcftools/filter/main'
include { BCFTOOLS_FILTER_XX   } from './modules/local/bcftools/filter_xx/main'
include { LIFTOVER   } from './modules/local/liftover/main'
include { BCFTOOLS_CONVERT     } from './modules/local/bcftools/convert/main'
include { BCFTOOLS_VIEW     } from './modules/local/bcftools/view/main'
include { BCFTOOLS_STATS     } from './modules/local/bcftools/stats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//


workflow {

    if (!params.input_tsv) {
        error "Please provide a TSV file using the '--input_tsv' parameter."
    }

    // Read the TSV file and extract columns
    ch_samples = Channel
        .fromPath(params.input_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def meta = [ id: row.sample_id ]
            def vcf  = row.inputvcf
            def tbi  = row.tbi
            def sex  = row.sex
            return [meta, vcf, tbi, sex]
        }
    
    // Filter samples by sex
    ch_xx_samples = ch_samples.filter { meta, vcf, tbi, sex -> sex == "xx" }
    ch_xy_samples = ch_samples.filter { meta, vcf, tbi, sex -> sex == "xy" }

    // Set params
    ch_ref_fasta19 = Channel.from(params.ref_fasta19)
    ch_fai19 = Channel.from(params.ref_fai19)
    ch_ref_fasta37 = Channel.from(params.ref_fasta37)
    ch_fai37 = Channel.from(params.ref_fai37)
    ch_chain = Channel.from(params.chain)
    ch_accesibility_bed = Channel.from(params.accesibility_bed)
    ch_chromosomes = Channel.from(params.chromosomes)

    // Filter

    ch_input_filter = ch_xy_samples
    .map{ meta, vcf, tbi, sex -> [meta, vcf]}

    ch_input_filter_xx = ch_xx_samples
    .map{ meta, vcf, tbi, sex -> [meta, vcf]}

    ch_filtered = BCFTOOLS_FILTER(ch_input_filter).filtered_vcf
    .mix(BCFTOOLS_FILTER_XX(ch_input_filter_xx).filtered_vcf)

    // LiftOver 

    ch_lift_input = ch_filtered
    .combine(ch_ref_fasta19)
    .combine(ch_fai19)
    .combine(ch_chain)
    .multiMap{ meta, vcf, fasta19, fai, chain ->
    vcf: [meta, vcf]
    fasta: [fasta19, fai]
    chain: chain
    }

    LIFTOVER(ch_lift_input)

    ch_convert_input = LIFTOVER.out.tmp_vcf
    .combine(LIFTOVER.out.tmp_tbi, by: 0)
    .combine(ch_chromosomes)
    .combine(ch_ref_fasta37)
    .multiMap{
        meta, vcf, tbi, chrom_file, fasta ->
        vcf: [meta, vcf, tbi]
        chrom_file: chrom_file
        fasta: fasta
    }

    // Convert
    BCFTOOLS_CONVERT(ch_convert_input)
    ch_view_input = BCFTOOLS_CONVERT.out.converted_vcf
    .combine(BCFTOOLS_CONVERT.out.converted_tbi, by: 0)
    .combine(ch_accesibility_bed)
    .multiMap{
        meta, vcf, tbi, bed ->
        vcf: [meta, vcf, tbi]
        bed: bed
    }

    // Apply accesiblity filter 1KGP
    BCFTOOLS_VIEW(ch_view_input)
    ch_stats_input = BCFTOOLS_VIEW.out.final_vcf
    .combine(BCFTOOLS_VIEW.out.final_tbi, by: 0)
    //.multiMap{meta, vcf, tbi ->
    //vcf: [meta, vcf, tbi]}
    
    BCFTOOLS_STATS(ch_stats_input)

}