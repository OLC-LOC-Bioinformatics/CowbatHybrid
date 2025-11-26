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
these can be downloaded (see `Installation` to get info on getting them.)

### Installation

##### Cowbat-hybrid Pipeline and Dependencies

This is a pain to install - here's how it seems to work best:

1) Create a new conda environment (must by python 3.5!): `conda create -n cowbat_hybrid python=3.5`
2) Activate your new conda env: `conda activate cowbat_hybrid`
3) sudo apt install -y build-essential pkg-config libcurl4-openssl-dev libssl-dev zlib1g-dev
4) sudo apt install mummer
5) clone the git environment like this git clone https://github.com/OLC-LOC-Bioinformatics/CowbatHybrid.git
6) cd CowbatHybrid
7) Install some requirements via pip: `pip install -r requirements.txt`
8) Install other things via conda/mamba (you should have the `conda-forge` and `bioconda` channels set up): `mamba install racon bbmap blast mob_suite clark nanoplot prodigal filtlong porechop sistr_cmd mash unicycler=0.4.4`
9) One final thing to do before you're good to go: by default the bioconda version of pilon allows
a max of 1 GB RAM. To rectify, find the pilon executable on your system (`which pilon`) and then open
the resulting file in your text editor of choice. Change the -Xmx1g in `default_jvm_mem_opts = ['-Xms512m', '-Xmx1g']` (should be around line 16) 
to a higher amount of max RAM (12 GB or so should work `default_jvm_mem_opts = ['-Xms512m', '-Xmx12g']`)

You should now be good to go! Typing `cowbat-hybrid-assembly.py` will let you run the pipeline.

You might get an error saying that `famap` and `fahash` can't be found. The binaries for those should be in the
python olctools package - you can add `/miniconda_install_location/envs/your_cowbathybrid_env/lib/python3.5/site-packages/spadespipeline/ePCR/`
to your `$PATH` to fix this.

##### Databases

Databases can be downloaded by clicking the following link: `https://scist01.blob.core.windows.net/olc/0.3.4.tar.gz?sp=r&st=2018-11-23T18:43:10Z&se=2021-04-01T01:43:10Z&spr=https&sv=2017-11-09&sig=aIqmlsFgnz1VCpyocJGHzJqz0t5seHdgYguTmwE6bcs%3D&sr=b`

You'll then need to unzip and untar the files, but after that you're ready to go.

### Output Files

Your output folder will have an assembly for each input set of sequences, and a `reports` folder
with lots of information - the most important of these is `combinedMetadata.csv`, which summarizes
most of the analyses done.
