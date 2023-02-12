# FLAQ-SC2-ClearLabs (Florida Assembly Quality - SARS-CoV-2 - Clear Labs)
FL BPHL's SARS-CoV-2 (SC2) analysis pipeline for whole-genome tiled-amplicon data generated on the [Clear Labs Clear Dx instrument](https://www.clearlabs.com/cleardx/#fully-automated-next-generation-sequencing-platform). 

## About
FLAQ-SC2-ClearLabs was developed to analyze whole-genome tiled-amplicon data (i.e., [ARTIC protocol](https://artic.network/ncov-2019)) from SC2-positive clinical specimens sequenced on the [Clear Labs Clear Dx instrument](https://www.clearlabs.com/cleardx/#fully-automated-next-generation-sequencing-platform). This pipeline provides a secondary analysis using the outputs from Clear Dx BIP pipeline to generate a report including read/mapping/assembly quality metrics, Pango lineage, and quality flags to support screening samples prior to submissions to public repositories (e.g., GISAID, NCBI Genbank). The current version will run only on [HiPerGator](https://www.rc.ufl.edu/about/hipergator/)(HPG) using local Singularity containers for each pipeline process.

Stay tuned for FLAQ-SC2-ClearLabs's upgrade to [Daytona](https://github.com/BPHL-Molecular/Daytona), a platform agnostic [Nextflow](https://www.nextflow.io/) workflow. Daytona is currently under active development.

## Dependencies
- Python3
- Singularity/Apptainer
- Git

To load python3 into your current environment on HiPerGator, either use `module load python` to get the lastest version of python or activate your base conda environment. For more information on how to set up your base conda environment on HPG, see the [HiPerGator Analysis Reference Guide](https://github.com/StaPH-B/southeast-region/tree/master/hipergator)).

Singularity/Apptainer will be loaded as a module during your job execution on HPG using the sbatch job script in this repository. 

Git is already installed in your HPG environment upon login.

## Usage

For first time use, clone this repository to a directory in blue on HPG, such as in /blue/bphl-\<state\>/\<user\>/repos/bphl-molecular/.
```
cd /blue/bphl-<state>/<user>/repos/bphl-molecular/
git clone https://github.com/BPHL-Molecular/flaq_sc2_clearlabs.git
```
For future use, update any changes to your local repository on HPG by navigating to the flaq_sc2_clearlabs repository and pulling any changes.
```
cd flaq_sc2_clearlabs/
git pull
```
To run the FLAQ-SC2-ClearLabs pipeline, copy all files from the flaq_sc2_clearlabs local repository to your analysis folder. Make an input directory and copy your fastqs.
```
mkdir <analysis_dir>
cd <analysis_dir>
cp /blue/bphl-<state>/<user>/repos/bphl-molecular/flaq_sc2_clearlabs/* .
mkdir fastqs/ bams/ assemblies/
cp /path/to/fastqs/*.fastq fastqs/
cp /path/to/bams/*.bam bams/
cp /path/to/assemblies/*.fasta assemblies/
```
Rename your fastq files to the following format: sample.fastq. See below for a helpful command to rename your fastq files.
```
cd fastqs/
for i in *.fastq; do mv -- "$i" "${i%[PATTERN to REMOVE]}.fastq"; done
```
Edit your sbatch job submission script to include your email to receive an email notification upon job END or FAIL. Replace ENTER EMAIL in `#SBATCH --mail-user=ENTER EMAIL` with your email address. Make sure there is no space between = and your email address. Edit additional sbatch parameters as needed to run your job succesfully, such as the length of time the job will run.

Submit your job.
```
sbatch sbatch_clearlabs_report.sh
```

## Main processes
- [Samtools](https://github.com/samtools/samtools)
- [Pangolin](https://github.com/cov-lineages/pangolin)
- [VADR](https://github.com/ncbi/vadr)

## Primary outputs

Outputs from each process for each individual sample can be found in a sample-specific subdirectory within the FLAQ-SC2-ClearLabs analysis directory. Report.txt contains the main summary report with read/mapping/assembly quality metrics, Pango lineage, and quality flags to support screening samples prior to submissions to public repositories (e.g., GISAID, NCBI Genbank). Additional details can be found in the report outputs from each process. Final passing  assemblies (.fasta) and annotation alerts are copied into the run directory for easier access for use in downstream analyses or for samples that require manual review.

```
analysis_dir/
|__ <date>_flaq_run/
     |__ report.txt
     |__ sample1/
     |__ sample2/
|__ assemblies/
|__ vadr_error_reports/
```

## Developed by:
[@SESchmedes](https://www.github.com/SESchmedes)<br />

Please email bphl16BioInformatics@flhealth.gov for questions and support.
