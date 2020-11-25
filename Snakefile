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
import re


#
# Init
#

# Read config
configfile: 'config.json'


# Read cell table
def get_cell_table():
    """
    Get a DataFrame of input data with (sample, cell) index.

    :return: Pandas DataFrame of input data.
    """

    # Read FOFN table
    df_fofn = pd.read_csv('samples.tsv', sep='\t', index_col=['SAMPLE'], usecols=('SAMPLE', 'FOFN'), squeeze=True)

    sample_list = list(df_fofn.keys())

    # Build cell table
    cell_table_list = list()

    for sample_name in sample_list:
        with open(df_fofn[sample_name]) as in_file:
            for line in in_file:

                # Process FOFN line
                line = line.strip()

                if not line or line.startswith('#'):
                    continue

                if line.endswith('.bam'):
                    file_format = 'BAM'

                elif line.endswith('.fastq') or line.endswith('.fastq.gz'):
                    file_format = 'FASTQ'

                else:
                    raise RuntimeError(
                        'Expected BAM or FASTQ in input FOFN (see config): Found: "{}": {}'.format(
                            line, df_fofn[sample_name]
                        )
                    )

                # Add to cell table
                cell_table_list.append(pd.Series(
                    [
                        sample_name,
                        re.sub(r'\.(bam|fastq(\.gz)?)$', '', os.path.basename(line), flags=re.IGNORECASE),
                        file_format,
                        line
                    ],
                    index=['SAMPLE', 'CELL', 'FORMAT', 'FILE']
                ))

    df_cell = pd.concat(
        cell_table_list, axis=1
    ).T.set_index(
        ['SAMPLE', 'CELL']
    )

    return df_cell

def list_cells(sample):
    """
    Get a list of cells for a sample. Used by expand rules to find all input for one sample.

    :param sample: sample.

    :return: List of cells.
    """
    df_cell = get_cell_table().reset_index()

    return list(df_cell.loc[df_cell['SAMPLE'] == sample, 'CELL'])


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
            ) for cell in list_cells(wildcards.sample)
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
    run:

        cell_record = get_cell_table().loc[(wildcards.sample, wildcards.cell)]

        if cell_record['FORMAT'] == 'BAM':
            # Get converted FASTQ
            input_fastq = f'temp/{wildcards.sample}/align/fastq/{wildcards.cell}.fastq.gz'

        elif cell_record['FORMAT'] == 'FASTQ':
            # Get original FASTQ (no BAM > FASTQ conversion was needed)
            input_fastq = cell_record['FILE']

        else:
            # Unknown format. FORMAT field is out of sync with get_cell_table()
            raise RuntimeError(
                'Unknown file format: {} (Program bug: should have been caught by get_cell_table())'.format(
                    cell_record['FORMAT']
                )
            )

        # Run alignment
        shell(
            """ngmlr --bam-fix -x {params.preset} -t {params.threads} -r {params.ref} -q {input_fastq} -o {output.bam}"""
        )

# sniffles_bam_to_fq
#
# Convert a PacBio cell BAM to a FASTQ.
rule sniffles_bam_to_fq:
    input:
        bam=lambda wildcards: get_cell_table().loc[(wildcards.sample, wildcards.cell), 'FILE']
    output:
        fastq=temp('temp/{sample}/align/fastq/{cell}.fastq.gz')
    run:

        cell_record = get_cell_table().loc[(wildcards.sample, wildcards.cell)]

        if cell_record['FORMAT'] == 'BAM':
            # Generate FASTQ for ngmlr
            shell(
                """bam2fastq -o temp/{wildcards.sample}/align/fastq/{wildcards.cell} {input.bam}"""
            )

        elif cell_record['FORMAT'] == 'FASTQ':
            # Write an empty file: File is a flag when data is already FASTQ from the source
            with open(output.fastq, 'w') as out_file:
                pass

        else:
            # Unknown format. FORMAT field is out of sync with get_cell_table()
            raise RuntimeError(
                'Unknown file format: {} (Program bug: should have been caught by get_cell_table())'.format(
                    cell_record['FORMAT']
                )
            )

