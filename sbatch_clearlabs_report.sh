#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=clearlabs_report
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ENTER EMAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=3-00
#SBATCH --output=clearlabs_report.%j.out
#SBATCH --error=clearlabs_report.%j.err

module load apptainer

#Run script/command and use $SLURM_CPUS_ON_NODE
python flaq_sc2_clearlabs.py --fastqs fastqs/ --assemblies assemblies/ --bams bams/ --threads $SLURM_CPUS_ON_NODE --sotc S:L452R,S:E484K --pango_path /apps/staphb-toolkit/containers/pangolin_4.1.3-pdata-1.17.sif --pangolin v4.1.3 --pangolin_data v1.17
