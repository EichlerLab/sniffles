# Sniffles Pipeline #

Runs `ngmlr` and `sniffles` on a set of input PacBio reads.

## Modules ##



## Run ##

### Write config ###

Write a configuration file to `config.json` pointing to the reference sequnece and the sample sequence reads.

Example:
```
{
  "reference": "/net/eichler/vol27/projects/autism_genome_assembly/nobackups/sv/reference/hg38.no_alt.fa",
}
```

### Execute ###

Define a variable that gives the full path to the pbsv pipeline code, which is the directory that contains `Snakefile`
and this `README.md` file. The pipeline itself does not use the variable, but commands in this README will.

Example:
`PIPELINE_DIR=/net/eichler/vol27/projects/structural_variation/nobackups/pipelines/sniffles/202006`

Load required modules (may work with later versions of these modules):
```
module load ngmlr/0.2.7
module load sniffles/202006
module load pbconda/202004
module load samtools/1.10
module load htslib/1.9
module load bcftools/1.9
module load miniconda/4.5.12
```

Run distributed:
`mkdir -p log; snakemake -s ${PIPELINE_DIR}/Snakefile --ri -j 30 -k --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -e ./log -o ./log -pe serial {cluster.cpu} -l mfree={cluster.mem} -l h_rt={cluster.rt} -j y -w n -S /bin/bash" -w 60 -u ${PIPELINE_DIR}/config/cluster.json --config ldpath=$LD_LIBRARY_PATH`


### Notes ###

Samtools 1.10 does not to work with ngmlr alignments, which is probably due to trimming alignments down to only soft-
clipped bases (See https://github.com/philres/ngmlr/issues/51). Seems to work with Samtools 1.9, and pipeline would need
to be modified to remove these alignment records before moving to 1.10.
