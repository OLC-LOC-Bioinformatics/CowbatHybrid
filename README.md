# COWBAT-hybrid

This pipeline is designed to perform hybrid assembly on a set of Illumina (paired-end) and MinION (or other long reads)
for bacterial genomes, and then perform typing on the resulting assemblies.

### Usage

`cowbat-hybrid-assembly.py -i template_csv -o output_dir -r database_folder`

As input, this pipeline needs a CSV file that tells it where to find Illumina and MinION files.
The `template.csv` file in this repository has the appropriate headers: `MinION`, `Illumina_R1`, `Illumina_R2`, and `OutName`.
For the MinION column and the Illumina columns, put the absolute path to the raw fastq files you want to use.
`OutName` can be anything - it's what the resulting assembly will be called.

You'll also need to provide a path to the databases that the typing portion of the pipeline needs to run - 
these can be downloaded and set up following the instructions at `ADD DOCS FOR THIS TO COWBAT REPO`

### Installation

Conda recipe to be built... doing this any other way promises to be excruciatingly difficult.

### Output Files

To be finalized...