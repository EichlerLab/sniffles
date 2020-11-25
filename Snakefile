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
def get_cell_table(sample=None):
    """
    Get a Pandas Series of subread BAMs. If `sample` is `None`, one record for each cell from all samples is returned
    with key (sample, cell). If `sample` is not `None`, records for that sample are returned keyed only by the cell.

    :param sample: Subset to this sample if not `None`.

    :return: Pandas series of subread BAMs.
    """

    # Read FOFN table
    df_fofn = pd.read_csv('samples.tsv', sep='\t', index_col=['SAMPLE'], usecols=('SAMPLE', 'FOFN'), squeeze=True)

    if sample is None:
        sample_list = list(df_fofn.keys())
    else:
        sample_list = [sample]

    # Bulid cell table
    cell_table_list = list()

    for sample_name in sample_list:
        with open(df_fofn[sample_name]) as in_file:
            for line in in_file:

                # Process FOFN line
                line = line.strip()

                if not line or line.startswith('#'):
                    continue

                if not line.endswith('.subreads.bam'):
                    raise RuntimeError(
                        'Expected subreads BAM in input FOFN (see config): Found: "{}": {}'.format(
                            line, df_fofn[sample_name]
                        )
                    )

                # Add to cell table
                cell_table_list.append(pd.Series(
                    [
                        sample_name,
                        os.path.basename(line).rstrip('.subreads.bam'),
                        line
                    ],
                    index=['SAMPLE', 'CELL', 'FILE']
                ))

    df_cell = pd.concat(cell_table_list, 1).T

    if sample is not None:
        del(df_cell['SAMPLE'])
        df_cell = df_cell.set_index('CELL').squeeze()
    else:
        df_cell = df_cell.set_index(['SAMPLE', 'CELL']).squeeze()

    return df_cell


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

localrules: sniffles_all

#
# Call
#

# sniffles_all
#
# Run all
rule sniffles_all:
    input:
        vcf=expand('{sample}.sniffles.vcf.gz', sample=set(get_cell_table().index.get_level_values('SAMPLE')))

# sniffles_call_compress
#
# Compress variant calls.
rule sniffles_call_compress:
    input:
        vcf='temp/{sample}/variants.vcf'
    output:
        vcf='{sample}.sniffles.vcf.gz',
        tbi='{sample}.sniffles.vcf.gz.tbi'
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
        bam='align/{sample}/mapped_reads.bam'
    output:
        vcf=temp('temp/{sample}/variants.vcf')
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
        bam=lambda wildcards: [
            'temp/{sample}/align/bam/sorted/{cell}.bam'.format(
                sample=wildcards.sample, cell=cell
            ) for cell in list(get_cell_table(wildcards.sample).keys())
         ]
    output:
        bam='align/{sample}/mapped_reads.bam',
        bai='align/{sample}/mapped_reads.bam.bai'
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
        bam='temp/{sample}/align/bam/unsorted/{cell}.bam'
    output:
        bam='temp/{sample}/align/bam/sorted/{cell}.bam'
    params:
        threads=config.get('sort_cores', '4'),
        mem=config.get('sort_mem', '1.5G'),
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
        fastq='temp/{sample}/align/fastq/{cell}.fastq.gz'
    output:
        bam=temp('temp/{sample}/align/bam/unsorted/{cell}.bam')
    params:
        ref=config['reference'],
        threads=config.get('map_cores', '18'),
        mem=config.get('map_mem', '2G'),
        preset=config.get('map_preset', 'pacbio')
    shell:
        """ngmlr --bam-fix -x {params.preset} -t {params.threads} -r {params.ref} -q {input.fastq} -o {output.bam}"""

# sniffles_bam_to_fq
#
# Convert a PacBio cell BAM to a FASTQ.
rule sniffles_bam_to_fq:
    input:
        bam=lambda wildcards: get_cell_table(wildcards.sample)[wildcards.cell]
    output:
        fastq=temp('temp/{sample}/align/fastq/{cell}.fastq.gz')
    shell:
        """bam2fastq -o temp/{wildcards.sample}/align/fastq/{wildcards.cell} {input.bam}"""
