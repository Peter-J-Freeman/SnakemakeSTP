# Computer Setup

## Create the environment
The conda envoronment is created by running

`$ conda env create -f environment.yml`

The tools listed in environment.yml are installed

***Note: I am running a mac M1 silicone chip and it seems there is no conda distribution for arm chips currently for 
samtools. I created by environment by hashing out samtools and installing using homebrew `brew install samtools
`***

### Additional dependencies
See requirements.txt

### Files

#### CRAM files
Patient CRAM and index files were downloaded using wget e.g.

```
$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/MXL/NA19750/alignment/NA19750.alt_bwamem_GRCh38DH.20150718.MXL.low_coverage.cram.crai
```

The web addresses are sourced from the `Sano_-_Bioinformatician_-_Take-Home_Task (1)` pdf. Files are downloaded into the `cram` directory

VCF files for Base Quality Score Recalibartion (BQSR) were downloaded manually from this [bucket](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
and the required files are documented [here](https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.json)
Files are downloaded into the `bqsr_vcf` directory

#### fasta file
GRCh38 fasta file [downloaded from here](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false)
The fasta file is located in the `fasta` directory

#### Gene Panel file
The gene panel used in Exomiser configuration is compiled from Genomics England Panel App
[Parkinson Disease and Complex Parkinsonism (Version 1.110)](https://panelapp.genomicsengland.co.uk/panels/39/)

#### Human Phenotype Ontology file
The HPO terms used in Exomiser configuration were compiled from phenotypes provided by the [MIM](https://www.omim.org/entry/607060) entry 
for LRRK2 related Parkinsons. This would potentially require user involvement to define the exact terms they want in the 
analysis. The phenotypic descriptions from MIM mim were mapped to [HPO Identifiers](https://hpo.jax.org/app/about)

# Analysis workflow
***Commands were run in detached mode using tmux***

This allows a degree of manual parrallelisation and speeds up analysis. For deployment, automated parallelisation would 
be developed, depending on HPC resources

See the [tmux cheat sheet](https://tmuxcheatsheet.com/) for multiple sessions

## Merge CRAM files
Merging the low coverage WGS file into the high coverage WES file makes it more likely that phasing can be assigned 
during variant calling. Phasing would allow us to determine if 2 closely (or more) associated variants on the same 
chromosome need to be considered together in order to calculate the clinical effect. Merging also provides additional 
data to improve the quality of variant calling and may improve calling at poorly covered variants in the wgs crams. 

Merging is run using the following commands

```
$ samtools merge -r -c -p -X \
 ./cram/NA19749_merged.cram \
 ./cram/NA19749.alt_bwamem_GRCh38DH.20150826.MXL.exome.cram \
 ./cram/NA19749.alt_bwamem_GRCh38DH.20150718.MXL.low_coverage.cram \
 ./cram/NA19749.alt_bwamem_GRCh38DH.20150826.MXL.exome.cram.crai \
 ./cram/NA19749.alt_bwamem_GRCh38DH.20150718.MXL.low_coverage.cram.crai
```

```
$ samtools merge -r -c -p -X \
 ./cram/NA19750_merged.cram \
 ./cram/NA19750.alt_bwamem_GRCh38DH.20150826.MXL.exome.cram \
 ./cram/NA19750.alt_bwamem_GRCh38DH.20150718.MXL.low_coverage.cram \
 ./cram/NA19750.alt_bwamem_GRCh38DH.20150826.MXL.exome.cram.crai \
 ./cram/NA19750.alt_bwamem_GRCh38DH.20150718.MXL.low_coverage.cram.crai
```

***Recommended option for deployment: We would look at adding threading using the `--threads` option to increase performance***

Check files using 

```
$ samtools quickcheck ./cram/NA19749_merged.cram
$ samtools quickcheck ./cram/NA19750_merged.cram
```

## Base Quality Score Recalibration (BQSR)

***BQSR is potentially obsolete for later versions of GATK tools and improved sequencing techniques. However, 
1000 genomes project data are very old now and may not have been sequenced using the most up-to-date techniques (so are more likely to have degradation in the read quality over time). I 
decided to include this step because the reads are likely to be shorter than modern NGS reads and the sequencing 
chemistry used to be less robust, so the data are more prone to error***

BQSR is run using the following commands

```
$ gatk BaseRecalibrator \
   --add-output-sam-program-record true \
   -I ./cram/NA19749_merged.cram \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   --known-sites ./bqsr_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   -O ./bqsr_vcf/NA19749_bqsr.table
```

```
$ gatk BaseRecalibrator \
   --add-output-sam-program-record true \
   -I ./cram/NA19750_merged.cram \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   --known-sites ./bqsr_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   -O ./bqsr_vcf/NA19750_bqsr.table
```


***Note: vcf files for indels were downloaded from the recommended source, but when run had errors in the VCF headers
With more time available I would look for alternates. The following were removed from the command above. However, I 
later managed to bug fix the issue and followed a different `download protocol` which was more
reliable for these files. I would put these commands back into the script for future iterations. The issue
turned out to be related to the URL provided to wget. The same website provided different URLs. In future, use the
`google public URL`***
```
   --known-sites ./bqsr_vcf/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites ./bqsr_vcf/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf\
```


Indexes created using 
```
$ samtools index -b -M ./cram/NA19749_merged.cram ./cram/NA19749_merged.crai
$ samtools index -b -M ./cram/NA19750_merged.cram ./cram/NA19750_merged.crai
```

## Mark duplicates
Illumina data frequently contain read duplicates. There are multiple reasons, but it is usually a result of an 
amplification bias for some molecules during the PCR stages of sample preparation. Again, 1000 genomes data is pretty 
old so may be more prone to such biases, so I chose to include this step

```
$ picard MarkDuplicates \
      I=./cram/NA19749_merged.cram \
      O=./cram/NA19749_merged.cram \
      M=./cram/NA19749_marked_metrics.txt
```

```
$ picard MarkDuplicates \
      I=./cram/NA19750_merged.cram \
      O=./cram/NA19750_merged.cram \
      M=./cram/NA19750_marked_metrics.txt
```

## Variant Calling
Variant calling utilises the [gatk HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
 workflow

```
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -I ./cram/NA19749_merged.cram \
   -O ./vcf/NA19749_merged.vcf.gz \
   -bamout ./vcf/NA19749_merged_bamout.bam
```

```
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -I ./cram/NA19750_merged.cram \
   -O ./vcf/NA19750_merged.vcf.gz \
   -bamout ./vcf/NA19750_merged_bamout.bam
```

https://gatk.broadinstitute.org/hc/en-us/articles/360035890551?id=9622


## Variant Quality Recalibration

Phase 1: Recalibrate SNPs
```
$ gatk VariantRecalibrator \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19749_merged.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./bqsr_vcf/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ./bqsr_vcf/1000G_omni2.5.hg38.vcf.gz   \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ./bqsr_vcf/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./bqsr_vcf/Homo_sapiens_assembly38.dbsnp138.vcf  \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -mode SNP \
   -O ./vcf/NA19749_snp.recal \
   --tranches-file ./vcf/NA19749_snp.tranches \
   --rscript-file ./vcf/NA19749_snp.plots.R
```

```
$ gatk ApplyVQSR \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19749_merged.vcf.gz \
   -O ./vcf/NA19749_snp.recal.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ./vcf/NA19749_snp.tranches \
   --recal-file ./vcf/NA19749_snp.recal \
   -mode SNP
```

```
$ gatk VariantRecalibrator \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19750_merged.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./bqsr_vcf/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ./bqsr_vcf/1000G_omni2.5.hg38.vcf.gz   \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ./bqsr_vcf/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./bqsr_vcf/Homo_sapiens_assembly38.dbsnp138.vcf  \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   -mode SNP \
   -O ./vcf/NA19750_snp.recal \
   --tranches-file ./vcf/NA19750_snp.tranches \
   --rscript-file ./vcf/NA19750_snp.plots.R
```

```
$ gatk ApplyVQSR \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19750_merged.vcf.gz \
   -O ./vcf/NA19750_snp.recal.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ./vcf/NA19750_snp.tranches \
   --recal-file ./vcf/NA19750_snp.recal \
   -mode SNP
```


Phase 2: Recalibrate INDELs

```
$ gatk VariantRecalibrator \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19749_snp.recal.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./bqsr_vcf/Homo_sapiens_assembly38.known_indels.vcf.gz  \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 ./bqsr_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  \
   -an QD -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode INDEL \
   -O ./vcf/NA19749_snp_indel.recal \
   --tranches-file ./vcf/NA19749_snp_indel.tranches \
   --rscript-file ./vcf/NA19749_snp_indel.plots.R
```

```
$ gatk ApplyVQSR \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19749_snp.recal.vcf.gz \
   -O ./vcf/NA19749_snp_indel.recal.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ./vcf/NA19749_snp_indel.tranches \
   --recal-file ./vcf/NA19749_snp_indel.recal \
   -mode INDEL
```


```
$ gatk VariantRecalibrator \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19750_snp.recal.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./bqsr_vcf/Homo_sapiens_assembly38.known_indels.vcf.gz  \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 ./bqsr_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  \
   -an QD -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode INDEL \
   -O ./vcf/NA19750_snp_indel.recal \
   --tranches-file ./vcf/NA19750_snp_indel.tranches \
   --rscript-file ./vcf/NA19750_snp_indel.plots.R
```

```
$ gatk ApplyVQSR \
   -R ./fasta/hg38_v0_Homo_sapiens_assembly38.fasta \
   -V ./vcf/NA19750_snp.recal.vcf.gz \
   -O ./vcf/NA19750_snp_indel.recal.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ./vcf/NA19750_snp_indel.tranches \
   --recal-file ./vcf/NA19750_snp_indel.recal \
   -mode INDEL
```

## Exomiser variant prioritisation

```
$ java -Xms4g -Xmx8g -jar ./exomizer/exomiser-cli-13.1.0/exomiser-cli-13.1.0.jar --analysis ./exomiser_analysis_configs/NA19749.yml 
$ java -Xms4g -Xmx8g -jar ./exomizer/exomiser-cli-13.1.0/exomiser-cli-13.1.0.jar --analysis ./exomiser_analysis_configs/NA19750.yml 
```

Exomiser configurations for each sample are located in the directory `exomiser_analysis_configs`
Exomiser results are written to `exomiser_results`

The setting are also written to the results and can be re-loaded for other samples that we may wish to 
analyse