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

        ch_all_coverages = COVERAGE.out.target_coverage
            .mix(COVERAGE.out.antitarget_coverage)
            .map { cov_file ->
                def sample = cov_file.name.replaceAll(/(target|antitarget)coverage_(.+)\.cnn/, '$2')
                [sample, cov_file]
            }
            .combine(ch_input_controls.map { sample, sex, _bam -> [sample, sex] }, by: 0)
            .map { _sample, cov_file, sex -> [sex, cov_file] }
            .groupTuple(by: 0)

        REFERENCE(
            ch_all_coverages,
            params.genome_fasta,
        )
    }

    if (params.call_cnvs) {
        ch_input_samples = Channel.fromPath(params.samples_csv_file, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sample = row.sample
                def sample_sex = row.Chr_sex
                def bam_file = file(row.bam_file)

                return [sample, sample_sex, bam_file]
            }

        ch_ref = Channel.value(
            [
                file(params.targets_file),
                file(params.antitargets_file),
                file(params.reference_m_file),
                file(params.reference_f_file),
            ]
        )

        CNV_CALLS(
            ch_input_samples,
            ch_ref,
        )
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
    tag "${sample}"
    publishDir "reference/coverages", mode: 'copy'

    input:
    path targets
    path antitargets
    path genome_fasta
    tuple val(sample), val(sample_sex), path(bam_file)

    output:
    path "targetcoverage_${sample}.cnn", emit: target_coverage
    path "antitargetcoverage_${sample}.cnn", emit: antitarget_coverage

    script:
    """
    cnvkit.py coverage ${bam_file} ${targets} -p ${task.cpus} -o targetcoverage_${sample}.cnn
    cnvkit.py coverage ${bam_file} ${antitargets} -p ${task.cpus} -o antitargetcoverage_${sample}.cnn
    """
}

process REFERENCE {
    tag "${sex}"
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

process CNV_CALLS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/", mode: 'copy'

    input:
    tuple val(sample), val(sample_sex), path(bam_file)
    tuple path(targets), path(antitargets), path(male_ref), path(female_ref)

    output:
    path ("${sample}.targetcoverage.cnn"), emit: target_cov
    path ("${sample}.antitargetcoverage.cnn"), emit: antitarget_cov
    path ("${sample}.ratios.cnr"), emit: ratios
    path ("${sample}.segments.cns"), emit: segments
    path ("${sample}.calls.cns"), emit: calls
    path ("${sample}.scatter.pdf"), emit: scatterplot
    path ("${sample}.diagram.pdf"), emit: diagram

    script:
    def ref = sample_sex == 'm' ? male_ref : female_ref
    """
    cnvkit.py coverage ${bam_file} ${targets} -o ${sample}.targetcoverage.cnn -p ${task.cpus}
    cnvkit.py coverage ${bam_file} ${antitargets} -o ${sample}.antitargetcoverage.cnn -p ${task.cpus}

    cnvkit.py fix ${sample}.targetcoverage.cnn ${sample}.antitargetcoverage.cnn ${ref} -o ${sample}.ratios.cnr

    cnvkit.py segment ${sample}.ratios.cnr -o ${sample}.segments.cns -p ${task.cpus}
    cnvkit.py call ${sample}.segments.cns --sample-sex ${sample_sex} -o ${sample}.calls.cns

    cnvkit.py scatter ${sample}.ratios.cnr -s ${sample}.segments.cns --title ${sample} -o ${sample}.scatter.pdf
    cnvkit.py diagram -s ${sample}.segments.cns -t 99 --sample-sex ${sample_sex} --title ${sample} -o ${sample}.diagram.pdf
    """
}
