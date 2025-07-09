#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    if (params.build_controls) {
        TARGETS(
            params.coverage_bed,
            params.genome_fasta,
        )

        ch_input_controls = Channel.fromPath(params.ctrls_csv_file, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample = row.sample
                def sample_sex = row.Chr_sex
                def bam_file = file(row.bam_file)

                return [sample, sample_sex, bam_file]
            }

        COVERAGE(
            TARGETS.out.targets,
            TARGETS.out.antitargets,
            params.genome_fasta,
            ch_input_controls,
        )

        ch_coverage_with_sex = ch_input_controls
            .map { sample, sex, _bam -> [sample, sex] }
            .combine(COVERAGE.out.coverage.map { cov_file ->
                def sample = cov_file.name.replaceAll(/coverage_(.+)\.cnn/, '$1')
                [sample, cov_file]
            }, by: 0)
            .map { _sample, sex, cov_file -> [sex, cov_file] }
            .groupTuple(by: 0)

        REFERENCE(
            ch_coverage_with_sex,
            params.genome_fasta,
        )
    }

    ch_input_samples = Channel.fromPath(params.samples_csv_file, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def bam_file = file(row.bam_file)
            def sample = row.sample
            def sample_sex = row.sex
            def ref_cnn_file = row.ref_cnn_file
            def probes = row.probes

            return [bam_file, sample, sample_sex, ref_cnn_file, probes]
        }
}

process TARGETS {
    publishDir "reference/", mode: 'copy'

    input:
    path coverage_bed
    path genome_fasta

    output:
    path ("targets.bed"), emit: targets
    path ("access.bed"), emit: access
    path ("antitargets.bed"), emit: antitargets

    script:
    """
    cnvkit.py target ${coverage_bed} --short-names --split -o targets.bed
    cnvkit.py access ${genome_fasta} -o access.bed
    cnvkit.py antitarget targets.bed -g access.bed -o antitargets.bed
    """
}

process COVERAGE {
    publishDir "reference/coverages", mode: 'copy'

    input:
    path targets
    path antitargets
    path genome_fasta
    tuple val(sample), val(sample_sex), path(bam_file)

    output:
    path "coverage_${sample}.cnn", emit: coverage

    script:
    """
    cnvkit.py coverage ${bam_file} ${targets} -p ${task.cpus} -o coverage_${sample}.cnn
    """
}

process REFERENCE {
    publishDir "reference", mode: 'copy'

    input:
    tuple val(sex), path(coverages)
    path genome_fasta

    output:
    path "reference_${sex == 'm' ? 'male' : 'female'}.cnn", emit: reference

    script:
    def sex_param = sex == 'm' ? 'male' : 'female'
    def output_name = "reference_${sex_param}.cnn"
    """
    cnvkit.py reference ${coverages} --sample-sex ${sex_param} --fasta ${genome_fasta} --output ${output_name}
    """
}
