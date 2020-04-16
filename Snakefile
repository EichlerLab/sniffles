"""
Run ngmlr and sniffles.
"""

############
### Init ###
############

#
# Import
#

import os
import numpy as np
import pandas as pd


#
# Init
#

# Read config
configfile: 'config.json'


# Read cell table
cell_table_list = list()

with open(config['fofn']) as in_file:
    for line in in_file:
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        if not line.endswith('.subreads.bam'):
            raise RuntimeError('Expected subreads BAM in input FOFN (see config): Found: {}'.format(line))

        cell_table_list.append(pd.Series(
            [
                os.path.basename(line).rstrip('.subreads.bam'),
                line
            ],
            index=['CELL', 'FILE']
        ))

df_cell = pd.concat(cell_table_list, 1).T
df_cell = df_cell.set_index('CELL').squeeze()


# Temp directory

if 'TMPDIR' in os.environ:
    temp_dir = os.path.abspath(os.environ['TMPDIR'])
else:
    temp_dir = os.path.abspath('temp/rule_temp')


# Set LD_LIBRARY_PATH from config options (cluster wipes it out when distributing)

if 'ldpath' in config:
    os.environ['LD_LIBRARY_PATH'] = config['ldpath']



#############
### Rules ###
#############


#
# Call
#

# sniffles_call_compress
#
# Compress variant calls.
rule sniffles_call_compress:
    input:
        vcf='temp/variants.vcf'
    output:
        vcf='variants.vcf.gz',
        tbi='variants.vcf.gz.tbi'
    params:
        temp_dir=os.path.join(temp_dir, 'sniffles_call_compress')
    shell:
        """mkdir -p {params.temp_dir}; """
        """bcftools sort -T {params.temp_dir} -O z -o {output.vcf} {input.vcf}; """
        """tabix {output.vcf}"""

# sniffles_call
#
# Call variants.
rule sniffles_call:
    input:
        bam='mapped_reads.bam'
    output:
        vcf=temp('temp/variants.vcf')
    params:
        min_svlen=config.get('min_svlen', '30')
    shell:
        """sniffles -t 12 -m {input.bam} -l {params.min_svlen} -v {output.vcf} --genotype --cluster"""


#
# Merge BAM
#

# sniffles_merge_bam
#
# Merge aligned BAM files.
rule sniffles_merge_bam:
    input:
        bam=expand('temp/align/bam/sorted/{cell}.bam', cell=list(df_cell.index))
    output:
        bam='mapped_reads.bam',
        bai='mapped_reads.bam.bai'
    shell:
        """samtools merge {output.bam} {input.bam}; """
        """samtools index {output.bam}"""


#
# Map
#

# sniffles_map_sort
#
# Sort mapped reads.
rule sniffles_map_sort:
    input:
        bam='temp/align/bam/unsorted/{cell}.bam'
    output:
        bam='temp/align/bam/sorted/{cell}.bam'
    params:
        threads='4',
        mem='1.5G',
        temp_dir=lambda wildcards: os.path.join(temp_dir, 'sniffles_map_sort_{}'.format(wildcards.cell))
    shell:
        """echo "Hostname: $(hostname)"; """
        """echo "Temp directory: {params.temp_dir}"; """
        """mkdir -p {params.temp_dir}; """
        """samtools sort -@ {params.threads} -T {params.temp_dir} -o {output.bam} {input.bam}"""

# sniffles_map
#
# Map one cell.
rule sniffles_map:
    input:
        fastq='temp/align/fastq/{cell}.fastq.gz'
    output:
        bam=temp('temp/align/bam/unsorted/{cell}.bam')
    params:
        ref=config['reference'],
        threads='18',
        mem='8G',
        preset='pacbio'
    shell:
        """ngmlr --bam-fix -x {params.preset} -t {params.threads} -r {params.ref} -q {input.fastq} -o {output.bam}"""

# sniffles_bam_to_fq
#
# Convert a PacBio cell BAM to a FASTQ.
rule sniffles_bam_to_fq:
    input:
        bam=lambda wildcards: df_cell[wildcards.cell]
    output:
        fastq=temp('temp/align/fastq/{cell}.fastq.gz')
    shell:
        """bam2fastq -o temp/align/fastq/{wildcards.cell} {input.bam}"""
