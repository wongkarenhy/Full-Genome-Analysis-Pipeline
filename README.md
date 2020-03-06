# README 
# Clinical natural language processing and automated interpretation
## By Karen Wong and Michelle Verghese

**Required software:**<br>
1. jq-linux64 (for parsing json files)
2. clinphen (for NLP)
4. python2.7
3. python3.6 or above
4. annovar 
5. intervar (to annotate small variants)
6. bcftools 

**Required python packages:** allel,argparse,collections,datetime,io,itertools,networkx,numpy,obonet,os,subprocess,pandas,pyranges,re,sys

**Descriptions:** <br>
This tool parses SNPs, indels, and structural variations (SVs) from 10x Genomics linked-read and Bionano optical mapping data based on trio sequencing (singleton is allowed). SNPs/indels analysis can be done alone or in combination with SV analysis. In general, this tool parses a patient's electronic health record in JSON format and outputs a clinically relevant gene list. This gene list is then used to inform how genetic variants are prioritized. Genetic variants (SNPs, indels, and SVs) are vetted against a set of controls and parents. For SNPs and indels, variants are filtered based on allele frequencies reported by gnomad(?). Small variants reported as likely benign or benign by either Clinvar or Intervar are discarded from the pipeline. For SVs, the prevalent of these variants are compared against a set of 1KGP + CIAPM control sequenced previously by the Kwok lab. See below for more details. <br>

A pre-processing step is required in order to run this software. This pre-processing step takes the 10xG GATK output and applies filters based on GQ, DP, and PASS. This step removes the bulk of the variants that are likely to be artifacts. Remaining variants are annotated using Intervar, which is a wrapper for Annovar and it assigns ACMG pathogeneicity to all variants. Variants are additionally filtered for frameshift, nonframeshift, nonsynonymous, stopgain, stoploss, and splicing. Filtered variants are overlapped with the ranked gene list generated previously. Remaining variatns are ready for manual evaluation. <br>

For SVs, insertions, deletions, and duplications are annotated with known exons (exon-level not gene-level). Duplications and deletions are additionally used to search for known microdeletion and microduplication syndromes. Inversions and translocations are annotated with known genes (gene-level) and every call in these two categories are always reported. <br>

SV scripts have been designed to analyze BioNano Optical Mapping Data (.smap file format) and 10x Linked Reads Data (.vcf file format). All SV scripts read a proband file, mother file, father file, and a reference file which consists of SV calls from the 1000 Genome Project cohort as well as SV calls from all other parents in the study other than the parents of the proband being analyzed. The scripts output filtered proband calls with additional descriptor columns as a tab-delimited txt file. <br>

BioNanoDeletions, BioNanoInsertions, BioNanoDuplications selects calls of the SV type and eliminates calls below the inputted confidence threshold (default: 0.5). It performs a 50% reciprocal overlap with the reference file and removes calls that overlap. It performs a 50% reciprocal overlap with the inputted mother and father file separately and appends columns (Found_in_Mother, Found_in_Father) to describe overlap (True/False). It overlaps with exons and phenotypes and appends columns (Gene, Phenotype) with gene name and phenotype if found. <br>

BioNanoInversions, BioNanoTranslocations selects calls of the SV type and does not filter for confidence. It creates 20kb intervals around the start point and end point of the call. It overlaps the start and end intervals with reference file and removes calls that overlap. It overlaps the start and end intervals with the mother and father file separately and appends columns (Found_in_Mother, Found_in_Father) to describe overlap (True/False). It overlaps start and end intervals with exons and phenotypes and appends columns (Gene, Phenotype for start point; Gene2, Phenotype2 for end point) with gene name and phenotype if found. <br>

tenxDeletions reads the 10x Deletion calls (“dels.vcf”) and performs a 50% reciprocal overlap with reference file and removes calls that overlap. It performs a 50% reciprocal overlap with the inputted mother and father file separately and appends columns (Found_in_Mother, Found_in_Father) to describe overlap (True/False). It overlaps with exons and phenotypes and appends columns (Gene, Phenotype) with gene name and phenotype if found. <br>

tenxLargeSvDeletions, tenxLargeSvDuplications reads the 10x Large SV calls (“large_svs.vcf”) and selects calls of the SV type. It performs a 50% reciprocal overlap with reference file and removes calls that overlap. It performs a 50% reciprocal overlap with the inputted mother and father file separately and appends columns (Found_in_Mother, Found_in_Father) to describe overlap (True/False). It overlaps with exons and phenotypes and appends columns (Gene, Phenotype) with gene name and phenotype if found. <br>

tenxLargeSVInversions, tenxLargeSvUnknown, tenxLargeSvBreakends reads the 10x Large SV calls (“large_svs.vcf”) and selects calls of the SV type. It creates 10kb intervals around the start point and end point of the call. It overlaps the start and end intervals with reference file and removes calls that overlap. It overlaps the start and end intervals with the mother and father file separately and appends columns (Found_in_Mother, Found_in_Father) to describe overlap (True/False). It overlaps start and end intervals with exons and phenotypes and appends columns (Gene, Phenotype for start point; Gene2, Phenotype2 for end point) with gene name and phenotype if one is found. For unknwon and breakends types, only variants with quality score > 1 standard deviation above the mean are reported. <br>

All coordinates are based on hg38.<br>


**Command (must run pre-processing before this):**<br>
```
bash run_clinical_interpretation_pipeline.sh [path_to_json] [work_dir] [sample_id] [hpo_database_path][path_to_intervar] [bionano:true or false] [DLE/BspQI/None] [linkedReadSV:true or false] [trio or singleton]
```

**General assumptions about this program:**
1. All samples must be named BC0XX01, BC0XX02, and BC0XX03, where XX is the family ID shared acorss a trio<br>
    The last two digits (01/02/03) indicates father, mother, and proband respectively. Please use 04 and onward if there are more than one probands.<br>  
2. **[path_to_json]**<br>
    This is the path to the EHR file in JSON format. 
3. **[work_dir]**<br>
    By default, all output are written to this work directory. <br>
4. **[hpo_database_path]** <br>
    This is the HPO database path. Please refer to "Database files from HPO" section  <br>
5. **[path_to_intervar]** <br>
    This is the Intervar directory. There should be a folder named 'example' in this directory. <br>
6. **[bionano:true or false]** <br>
    Indicate true to analyze bionano data for SVs <br>
7. **[DLE/BspQI/None]** <br>
    Must use this option if bionano is true. BssSI is not supported as we don't have a large control database using this label. <br>
8. **[linkedReadSV:true or false]** <br>
    Indicate true to analyze linked-read SVs <br>
9. **[trio or singleton]** <br>
    Indicate whether this is a trio or a singleton case. If this is a singleton case, name the file as BC0XX03. For quad cases, submit each proband as a separate job.<br>
    
**Database files from HPO**<br>
http://purl.obolibrary.org/obo/hp.obo
http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab (build 1270)<br>
http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/genes_to_phenotype.txt (build 1270)<br>
http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/util/annotation/phenotype_to_genes.txt (build 1270)<br>

Please keep these four files in the hpo_database_path. <br>
I've noticed issues with columns and headers when you change the file version... <br>

**Getting started** <br>
To get started, pull the github repo and create two additional directories (bionano_sv and linkedRead_sv). Your work directory should look like this:<br>
```
.<br>
├── annotatedExon.bed
├── annotatedGene.bed
├── bionano_sv
├── cytoband.txt
├── human_pheno_ontology_b1270 
├── linkedRead_sv
├── README.md
└── scripts
```
**bionano_sv folder should have the following:** <br>
```
.
├── cases
│   ├── BspQI
│   │   ├── BC00101
│   │   │   └── exp_refineFinal1_merged_filter.smap 
│   │   ├── BC00102
│   │   │   └── exp_refineFinal1_merged_filter.smap 
│   │   └── BC00103
│   │       └── exp_refineFinal1_merged_filter.smap 
│   └── DLE<br>
│       ├── BC02101
│       │   └── exp_refineFinal1_merged_filter.smap 
│       ├── BC02102
│       │   └── exp_refineFinal1_merged_filter.smap 
│       └── BC02103
│           └── exp_refineFinal1_merged_filter.smap 
└── controls
    ├── BspQI
    │   └── HG00251
    │       └── exp_refineFinal1_merged_filter.smap 
    └── DLE
        └── HG00250 
            └── exp_refineFinal1_merged_filter.smap 
```            
            
            
**linkedRead_sv folder should have the following** <br>
```
.
├── cases
│   ├── BC00101
│   │   ├── dels.vcf.gz 
│   │   ├── dels.vcf.gz.tbi 
│   │   ├── large_svs.vcf.gz 
│   │   ├── large_svs.vcf.gz.tbi 
│   │   ├── phased_variants.vcf.gz 
│   │   └── phased_variants.vcf.gz.tbi 
│   └── BC00102
│       ├── dels.vcf.gz 
│       ├── dels.vcf.gz.tbi 
│       ├── large_svs.vcf.gz 
│       ├── large_svs.vcf.gz.tbi 
│       ├── phased_variants.vcf.gz 
│       └── phased_variants.vcf.gz.tbi 
└── controls<br>
    └── GM06986<br>
        ├── dels.vcf.gz <br>
        ├── dels.vcf.gz.tbi <br>
        ├── large_svs.vcf.gz <br>
        ├── large_svs.vcf.gz.tbi <br>
        ├── phased_variants.vcf.gz <br>
        └── phased_variants.vcf.gz.tbi <br>
 ```
    
    
**Output files explanations:**<br>
Score and normalized scores are two important indicators that show how relevant the variants are based on the proband's clinical phenome. The raw scores are normalized between 0-100 among all variants. 100 is assigned to highest ranking variant found in the proband. An html report combining all the result files in confident_set is generated for quick review. <br>

**Files of highest priority:**<br>
$SAMPLEID_confident_deletion_exons.txt<br>
$SAMPLEID_confident_deletion_syndrome.txt<br>
$SAMPLEID_confident_duplication_exons.txt<br>
$SAMPLEID_confident_duplication_syndrome.txt<br>
$SAMPLEID_smallVariants_compoundhet_candidates.txt<br>
$SAMPLEID_smallVariants_denovo_candidates.txt<br>
$SAMPLEID_smallVariants_recessive_candidates.txt<br>
$SAMPLEID_Bionano_translocations.txt<br>
$SAMPLEID_Bionano_inversions.txt/$SAMPLEID_10x_inversions_largeSV.txt<br>

```
.
├── BC05303_10x_deletion_largeSV_syndrome.txt
├── BC05303_10x_deletions_exons.txt
├── BC05303_10x_deletions_largeSV_exons.txt
├── BC05303_10x_deletion_syndrome.txt
├── BC05303_10x_duplication_largeSV_syndrome.txt
├── BC05303_10x_duplications_largeSV_exons.txt
├── BC05303_Bionano_deletions_exons.txt
├── BC05303_Bionano_deletion_syndrome.txt
├── BC05303_Bionano_duplications_exons.txt
├── BC05303_Bionano_duplication_syndrome.txt
├── BC05303_Bionano_insertions.txt
├── confident_set
│   ├── BC05303_10x_breakends_largeSV.txt
│   ├── BC05303_10x_inversions_largeSV.txt
│   ├── BC05303_10x_SV_SNPsIndels_compounthet_candidates.txt
│   ├── BC05303_10x_unknown_largeSV.txt
│   ├── BC05303_Bionano_inversions.txt
│   ├── BC05303_Bionano_SV_SNPsIndels_compounthet_candidates.txt
│   ├── BC05303_Bionano_translocations.txt
│   ├── BC05303_confident_deletion_exons.txt
│   ├── BC05303_confident_deletion_syndrome.txt
│   ├── BC05303_confident_duplication_exons.txt
│   ├── BC05303_confident_duplication_syndrome.txt
│   ├── BC05303_dominant_inherited_smallVariants_candidates.txt
│   ├── BC05303_smallVariants_compoundhet_candidates.txt
│   ├── BC05303_smallVariants_denovo_candidates.txt
│   ├── BC05303_smallVariants_recessive_candidates.txt 
│   └── report.html
└── misc
```

**Things to watch out for:**<br>
Small indels in the phased_variants.vcf.gz (specifically insertions larger than 5bp) are likely to be false positives. These are usually reported as de novo variants. Similar clipped read signature can be found in the parents usually. <br>
Large number of false positives in 10x SV calls. <br>

## Pre-processing step for SNPs and indels (must run this separately before running the clinical interpretation pipeline):<br>
**Step 1.** Change directory to the InterVar installation location<br>
```
cd /media/KwokRaid02/karen/software/InterVar<br>

bcftools view -i 'MIN(FMT/DP)>10 & MIN(FMT/GQ)>30' \
    -f PASS /media/KwokRaid04/CIAPM/CIAPM_longranger/BC00103_longranger/outs/phased_variants.vcf.gz | \
    bgzip -c > /media/KwokRaid02/karen/software/InterVar/input_vcf/BC00103_filtered.vcf.gz" 
```
**Step 2.** Generate the config files according to instructions provided by Intervar<br>
Run Intervar<br>
```
python2.7 ./Intervar.py -c ./configFiles/BC00103_config.ini
```
**Step 3.** Subset the Intervar output <br>
Grep functional variants and keep variants with maf <=0.05<br>
```
grep -w -E 'frameshift|nonframeshift|nonsynonymous|stopgain|stoploss|splicing' \
    ./example/BC00103.hg38_multianno.txt.intervar | \
    grep -i -v benign | awk -F'\t' '$15<=0.05' > ./example/BC00103.hg38_multianno.txt.intervar.FINAL

awk '{print $6}' ./example/BC00103.hg38_multianno.txt.intervar.FINAL | \
    sort -u > ./example/BC00103_smallVariant_geneList.txt
```
## Actual analysis:<br>
**Example command** <br>
```
bash /media/KwokRaid05/karen/ciapm/FGA/scripts/run_clinical_interpretation_pipeline.sh \
    /media/KwokRaid05/karen/ciapm/jsons/BC00103.json \
    /media/KwokRaid05/karen/ciapm/FGA \
    BC00103\
    /media/KwokRaid05/karen/ciapm/FGA/human_pheno_ontology_b1270/ \
    /media/KwokRaid02/karen/software/InterVar/ \
    true DLE true trio 
```








