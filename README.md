# 1011-barcodes


Retrieved table of genes and gene predictions for S. cerevisiae from UCSC Table Browser

Search settings for gene coordinates:

| Field     | Value                                 |
| -----     | -----                                 |
| clade     | other                                 |
| genome    | S. cerevisiae                         |
| assembly  | Apr. 2011 (SacCer_Apr2011/sacCer3)    |
| group     | Genes and Gene Predictions            |
| track     | SGD Genes                             |
| table     | sgdGene                               |
| region    | genome                                |


Search settings for expression coordinates:

| Field     | Value                                 |
| -----     | -----                                 |
| clade    | other                                  |
| genome    | S. cerevisiae                         |
| assembly  | Apr. 2011 (SacCer_Apr2011/sacCer3)    |
| group     | Expression and Regulation             |
| track     | Regulatory Code                       |
| table     | transRegCode                          |
| region    | genome                                |

# bgzip and tabix-index chromosome vcf files
```
(cd genomes && for file in *.vcf; do singularity exec ../src/singularity.sif bgzip $file; done)
(cd genomes && for file in *.vcf.gz; do singularity exec ../src/singularity.sif tabix -p vcf $file; done)

```

# extract regions from vcf
singularity exec src/singularity.sif tabix genomes/chromosome1.vcf.gz chromosome1:1-100

# test region

singularity exec src/singularity.sif Rscript stuff.R 1 500 1000