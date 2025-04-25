#!/usr/bin/env python

#Author: Sarah Schmedes
#Email: sarah.schmedes@flhealth.gov

'''
This program takes in ClearLabs fastqs using the ARTIC primer schemes for SARS-CoV-2, bams,
and fastas and generates a read and assembly report (with VADR flags) along with pangolin lineage.
'''

import os
import sys
import subprocess
import argparse
import datetime
import pandas as pd
import re
import os.path
from Bio import SeqIO

#Parse arguments, get path for fastqs, primer version
parser = argparse.ArgumentParser(usage='flaq_sc2_clearlabs.py [options]')
parser.add_argument('--fastqs', help='path to clearlabs fastqs dir')
parser.add_argument('--bams', help='path to clearlabs bams dir')
parser.add_argument('--assemblies', help='path to clearlabs assemblies dir') 
parser.add_argument('--threads', default=8, dest='threads', help='specify number of threads, (default: %(default)s)')
parser.add_argument('--sotc', help='comma separated list of SOTCs to screen (e.g., S:L452R,S:E484K')
parser.add_argument('--pango_path', help='path to pangolin container', required=False)
parser.add_argument('--pangolin', help='pangolin version (e.g., v2.3)', required=False)
parser.add_argument('--pangolin_data', help='pangolin-data version (e.g., v1.3)', required=False)
parser.add_argument('--version', action='version', version='This is flaq_sc2_clearlabs: Version 1.2', help='print version')

if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

fastqs = os.path.abspath(args.fastqs) + '/'
bams = os.path.abspath(args.bams) + '/'
assemblies = os.path.abspath(args.assemblies) +'/'
threads = str(args.threads)
sotc_v = args.sotc
sotc_v = sotc_v.split(',')
pango_path = args.pango_path
pango_v = args.pangolin
pdata = args.pangolin_data
cwd = os.getcwd() + '/'

output_dir = cwd + datetime.date.today().strftime('%Y-%m-%d') + '_clearlabs_lineage_report'
subprocess.run('mkdir -p ' + output_dir, shell=True, check=True) #make output directory
subprocess.run('mkdir -p vadr_error_reports_clearlabs', shell=True, check=True) #folder for vadr error report for easier review
subprocess.run('mkdir -p assemblies_pass', shell=True, check=True) #folder for passing clearlabs assemblies
if os.path.exists('nextclade_latest.sif'):
    subprocess.run('rm nextclade_latest.sif', shell=True, check=True)
subprocess.run('singularity pull nextclade_latest.sif docker://nextstrain/nextclade:latest', shell=True, check=True) #pull latest nextclade container
subprocess.run('mkdir -p data/nextclade', shell=True, check=True) #create nextclade dataset folder
subprocess.run('singularity exec nextclade_latest.sif nextclade dataset get --name sars-cov-2 --output-dir data/nextclade', shell=True, check=True) #get latest nextclade dataset

# Pull the latest pangolin container if path not provided
if not args.pango_path:
    if os.path.exists('pangolin_latest.sif'):
        subprocess.run('rm pangolin_latest.sif', shell=True, check=True)
    subprocess.run('singularity pull pangolin_latest.sif docker://staphb/pangolin:latest', shell=True, check=True)
    # Get the container version info
    try:
        proc = subprocess.run('singularity exec pangolin_latest.sif pangolin --all-versions', shell=True, capture_output=True, text=True, check=True)
        version_info = proc.stdout.strip()
        # Extract pangolin version
        pango_match = re.search(r'pangolin: ([\d.]+)', version_info)
        # Extract pangolin-data version
        pdata_match = re.search(r'pangolin-data: ([\d.]+)', version_info)
        if pango_match and pdata_match:
            pango_v = "v" + pango_match.group(1)
            pdata = pdata_match.group(1)
            pangolin = pango_v + '_pdata-' + pdata
        else:
            pango_v = "latest"
            pdata = "latest"
            pangolin = "latest_auto-pulled"
    except Exception as e:
        print(f"Warning: Error extracting pangolin version: {e}")
        pango_v = "latest"
        pdata = "latest"
        pangolin = "latest_auto-pulled"
    pango_path = "pangolin_latest.sif"
else:
    if pango_v is not None and pdata is not None:
        pangolin = pango_v + '_pdata-' + pdata
    else:
        pangolin = "user-specified-container"

samples = []

#Look at some code examples to get fastq names R1, _1 or R1_001 (make work for more sample types later)
for f in os.listdir(fastqs):
    if f.endswith('.fastq'):
        sn = f.split(".")
        sn = sn[0]
        samples.append(sn)
unique = set(samples)
samples = list(unique)
samples.sort()

#Create output file
report = open(output_dir + '/report.txt', 'w')
header = ['sampleID', 'reference', 'start', 'end', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'numN', 'percent_ref_genome_cov', 'VADR_flag', 'QC_flag', 'pangolin_version', 'lineage', 'SOTC']
report.write('\t'.join(map(str,header)) + '\n')

#Run pipeline for each sample
for s in samples:
    sample_dir = output_dir + '/' + s + '/'
    subprocess.run('mkdir -p ' + sample_dir, shell=True, check=True) #mkdir for each sample name

    out_log = open(sample_dir + s + '.out', 'w')
    err_log = open(sample_dir + s + '.err', 'w')

    #Get number of raw reads
    proc_r = subprocess.run('cat ' + fastqs + s + '.fastq | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_r.stdout.rstrip()
    clean_reads = int(wc_out_1) / 4

    #Get bam file name
    proc_b = subprocess.run('ls ' + bams + s + '*.bam', shell=True, capture_output=True, text=True, check=True)
    bam_file = proc_b.stdout.rstrip()

    #Generate bam with only mapped reads
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools view -F 4 -u -h ' + bam_file + ' | singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools sort > ' + sample_dir + s + '.mapped.sorted.bam', shell=True, check=True)
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools index ' + sample_dir + s + '.mapped.sorted.bam', shell=True, check=True)    

    #Run samtools coverage to get map stats
    subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/samtools_1.12.sif samtools coverage ' + sample_dir + s + '.mapped.sorted.bam -o ' + sample_dir + s + '.coverage.txt', shell=True, stdout=out_log, stderr=err_log, check=True)

    #Get map stats
    with open(sample_dir + s + '.coverage.txt', 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        mapq = stats[8]

    #Get percentage of mapped reads/reads
    percent_map = "%0.4f"%(((int(reads_mapped)/int(clean_reads)))*100)
    
    #Get fasta file name
    proc_f = subprocess.run('ls ' + assemblies + s + '*.fasta', shell=True, capture_output=True, text=True, check=True)
    fasta_file = proc_f.stdout.rstrip()
    
    #Gather QC metrics for consensus assembly
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        num_bases = len(seq)
        ns = seq.lower().count('n')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)
    
    #Copy assembly to sample dir in order to rename and run VADR
    subprocess.run('cp ' + assemblies + s + '*.fasta ' + sample_dir + s + '.consensus.fa', shell=True, check=True)
    #Rename header in fasta to just sample name
    subprocess.run('sed -i \'s/^>.*/>' + s + '/\' ' + sample_dir + s + '.consensus.fa', shell=True, check=True)
        
    #QC flag
    pg_flag = ''
    dp_flag = ''
    qc_flag = ''
    if float(pg) < 79.5:
        pg_flag = 'FAIL: Percent genome < 80%'
        qc_flag = qc_flag + pg_flag
    else:
        if float(depth) < 100:
            dp_flag = 'FAIL: Mean read depth < 100x'
            qc_flag = qc_flag + dp_flag
        if qc_flag == '':
            qc_flag = qc_flag + 'PASS'
   
    #Copy passing assemblies to assemblies folder and vcf files to variants folder
    if qc_flag == 'PASS':
        subprocess.run('cp ' + sample_dir + s + '.consensus.fa assemblies_pass/', shell=True, check=True)
        
        #Run VADR
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/vadr_1.3.sif /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 ' + sample_dir + s + '.consensus.fa > ' + sample_dir + s + '.trimmed.fasta', shell=True, stdout=out_log, stderr=err_log, check=True)
        subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/vadr_1.3.sif v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir /opt/vadr/vadr-models/ ' + sample_dir + s + '.trimmed.fasta -f ' + sample_dir + 'vadr_results', shell=True, stdout=out_log, stderr=err_log, check=True)
        #Run VADR
        #subprocess.run('singularity exec -B $(pwd):/data /apps/staphb-toolkit/containers/vadr_1.1.3.sif v-annotate.pl --mxsize 64000 -s -r --nomisc --mkey NC_045512 --lowsim5term 2 --lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf,insertnn,deletinn --mdir /opt/vadr/vadr-models/ ' + sample_dir + s + '.consensus.fa ' + sample_dir + 'vadr_results', shell=True, stdout=out_log, stderr=err_log, check=True)

        #Parse through VADR outputs to get PASS or REVIEW flag
        vadr_flag = ''

        with open(sample_dir + 'vadr_results/vadr_results.vadr.pass.list', 'r') as p_list:
            result = p_list.readline()
            result = result.rstrip()
            if result == s:
                vadr_flag = 'PASS'

        with open(sample_dir + 'vadr_results/vadr_results.vadr.fail.list', 'r') as f_list:
            result = f_list.readline()
            result = result.rstrip()
            if result == s:
                vadr_flag = 'REVIEW'

        #Copy VADR error report to main analysis folder for easier review
        if vadr_flag == 'REVIEW':
            subprocess.run('cp ' + sample_dir + 'vadr_results/vadr_results.vadr.alt.list vadr_error_reports_clearlabs/', shell=True, check=True)
            subprocess.run('mv vadr_error_reports_clearlabs/vadr_results.vadr.alt.list vadr_error_reports_clearlabs/' + s + '.vadr.alt.list', shell=True, check=True)

        #Run pangolin
        subprocess.run('singularity exec -B $(pwd):/data ' + pango_path + ' pangolin ' + sample_dir + s + '.consensus.fa', shell=True, check=True)
        #Get lineage
        proc = subprocess.run('tail -n 1 lineage_report.csv | cut -d \',\' -f 2', shell=True, check=True, capture_output=True, text=True)
        lineage = proc.stdout.rstrip()
        subprocess.run('mv lineage_report.csv ' + sample_dir, shell=True, check=True)

        #Run nextclade
        subprocess.run('singularity exec -B $(pwd):/data nextclade_latest.sif nextclade run --output-csv ' + sample_dir + 'nextclade_report.csv --jobs ' + threads + ' --input-dataset data/nextclade ' + sample_dir + s + '.consensus.fa', shell=True, check=True)

        #Parse nextclade output and screen for sotc
        with open(sample_dir + 'nextclade_report.csv', 'r') as nc:
            header = nc.readline()
            c_results = nc.readline()
            c_results = c_results.rstrip()
            data = c_results.split(',')
            sotc = []
            for v in sotc_v:
                if v in data:
                    sotc.append(v)
            sotc_out = (',').join(sotc)
    else:
        vadr_flag = 'NA'
        lineage = 'NA'
        sotc_out = 'NA'

    #Write to output file
    results = [s, ref_name, start, end, int(clean_reads), reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg, vadr_flag, qc_flag, pangolin, lineage, sotc_out]
    report.write('\t'.join(map(str,results)) + '\n')
    out_log.close()
    err_log.close()

report.close()

#Create multi-fasta file of only passing assemblise
subprocess.run('cat assemblies_pass/*.fa > assemblies_pass.fasta', shell=True, check=True)

# Get nextclade version
try:
    proc = subprocess.run('singularity exec nextclade_latest.sif nextclade --version', shell=True, capture_output=True, text=True, check=True)
    nextclade_version_output = proc.stdout.strip()
    version_match = re.search(r'nextclade\s+([\d.]+)', nextclade_version_output)
    if version_match:
        nextclade_version = version_match.group(1)
    else:
        nextclade_version = "latest"
except:
    nextclade_version = "latest"

#Run nextclade
subprocess.run('singularity exec -B $(pwd):/data nextclade_latest.sif nextclade run --output-tsv nextclade_report_clearlabs.tsv --jobs ' + threads + ' --input-dataset data/nextclade assemblies_pass.fasta', shell=True, check=True)

# Add nextclade version to the report
if os.path.exists('nextclade_report_clearlabs.tsv'):
    df = pd.read_csv('nextclade_report_clearlabs.tsv', sep='\t')
    df.insert(df.columns.get_loc("clade"), "nextclade_version", nextclade_version)
    df.to_csv('nextclade_report_clearlabs.tsv', sep='\t', index=False)
