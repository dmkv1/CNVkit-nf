# NextFlow wrapper for CNVkit

This pipeline was tested on the local Linux machine.

# Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [docker](https://docs.docker.com/desktop/setup/install/linux/)
- `etal/cnvkit:0.9.11` docker image

# Inputs
## Controls

To build CNVkit reference from whole exome sequenced set of controls, set `build_controls = true` in `nextflow.config` and prepare controls sample sheet with these columns:

| sample | Chr_sex | bam_file |
|--------|---------|----------|
| ctrl1  | f       | /path/to/ctrl1.bam |
| ctrl2  | m       | /path/to/ctrl2.bam |
| ctrl3  | f       | /path/to/ctrl3.bam |
| ...  | ...       | ... |

Column `Chr_sex` can be set to either `m` of `f`, estimation is not tested.

You will also need to specify baits file (`coverage_bed` in the config), usually provided by the exome library preparation kit manufacturer. It should correspond to genome version used for mapping.
Fasta file for the latter have to be specified in the `genome_fasta` field of the config.

## Samples

To call CNVs in .bam files you have to specify the reference and provide a sample sheet.

First, set `call_cnvs` in the config to `true`. Provide targets, antitargets, male, and/or female .cnn reference files in the corresponding fields.

Create `samples.csv` file with the following columns:

|  sample  | Chr_sex |       bam_file       |
|----------|---------|----------------------|
| sample1  | f       | /path/to/sample1.bam |
| sample2  | m       | /path/to/sample2.bam |
|...|...|...|

# Running the pipeline

To run the pipeline in either reference build or sample processing mode use:

```bash
nextflow run main.nf -resume
```

You can add `-bg` flag to run it detached from the current terminal.

# Outputs

## Reference

The result of reference building will be stored in the `reference` directory. Rename it after the succesfull build if you need to compile several references (e.g. for different sets of probes). The directory will contain `targets.bed`, `antitargets.bed`, `reference_male.cnn`, `reference_female.cnn`, `access.bed` and, for each sample: target and antitarget coverage file. Only first four files are needed for sample analysis.

## Sample analysis

The results of sample analysis would be stored in the `results` directory and would include, for each sample:
- `targetcoverage.cnn` and `antitargetcoverage.cnn`
- `ratios.cnr` - copy number rations
- `segments.cns` - copy number segments
- `calls.cns` - copy number segments
- `diagram.pdf` and `scatter.pdf` - default CNVkit plots. Scatter is usually very heavy and should be abringed or rasterized.

# Acknowledgement

[CNVkit](https://cnvkit.readthedocs.io/en/stable/) is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data.

Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014). [CNVkit: Genome-wide copy number detection and visualization from targeted sequencing.](http://dx.doi.org/10.1371/journal.pcbi.1004873) PLOS Computational Biology 12(4):e1004873
