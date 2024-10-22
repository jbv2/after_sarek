process BCFTOOLS_VIEW {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(bed)

    output:
    tuple val(meta), path("${meta.id}.haplotypecaller.g_hs37d5_gq20.accesibility.vcf.gz"), emit: final_vcf
    tuple val(meta), path("${meta.id}.haplotypecaller.g_hs37d5_gq20.accesibility.vcf.gz.tbi"), emit: final_tbi

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools view \\
    --regions-file $bed \\
    $vcf \\
    -Oz \\
    -o ${prefix}.haplotypecaller.g_hs37d5_gq20.accesibility.vcf.gz \\
    && bcftools index --tbi ${prefix}.haplotypecaller.g_hs37d5_gq20.accesibility.vcf.gz --threads ${task.cpus}
    """
}