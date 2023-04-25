# Computer Setup

## Create the environment
The conda envoronment is created by running

`$ conda env create -f environment.yml`

The tools listed in environment.yml are installed

***Note: We are going to create a docker container to run our analyses to prevent issues with OS differences. For 
example, an added bonus is that I am running a Mac M1 with an ARM Silicone chip (not intel) so there are limited 
distributions of some of the software we need. This is where docker can help because there are Ubuntu distributions
available that will run on my chip and can be installed using the Ubuntu inbuilt package mamager, apt-get***

### Additional dependencies
See requirements.txt

### Files

#### Genome fasta file
GRCh38 fasta file [downloaded from here](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
The file should be placed in the is located in the `data/genome` directory

#### CRAM files
VCF files for Quality Score Recalibartion were downloaded manually from this [bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
and the required files are documented [here](https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json)
Files are downloaded into the `data/qsrecal_vcf` directory

#### Gene Panel file & Human Phenotype Ontology file
The gene panel used in Exomiser configuration is compiled from Genomics England Panel App
[Parkinson Disease and Complex Parkinsonism (Version 1.110)](https://panelapp.genomicsengland.co.uk/panels/39/)

The HPO terms used in Exomiser configuration were compiled from phenotypes provided by the [MIM](https://www.omim.org/entry/607060) entry 
for LRRK2 related Parkinsons. This would potentially require user involvement to define the exact terms they want in the 
analysis. The phenotypic descriptions from MIM mim were mapped to [HPO Identifiers](https://hpo.jax.org/app/about)

Place these files in the `data/exomiser` directory

