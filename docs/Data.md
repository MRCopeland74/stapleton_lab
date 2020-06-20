---
title: "Data"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---
<style>
th {
  font-size: 8px
}
td{
  font-size: 10px
}
</style>

# Data Collected

## Curve fitting on qPCR
[Hierarchical](https://osf.io/wqrsd/) and the [Sampling Plan](https://osf.io/mu5qs/) Data are available at this link. Here is a sample of the sampling plan dataset.

## vQTL Data
The [hyrbid data](https://github.com/AustinGratton/vQTL/blob/master/qPCR/FullHybvqtlinput.csv) and the [inbred data](https://github.com/AustinGratton/vQTL/blob/master/qPCR/FullInbvqtlinput.csv) we used for the vQTL data was all in these two csv file known as the FullHybvqtlinput and FullInbvqtlinput Data. This file contained a phenotype which is the height of the corn crop. There are 8 columns of different environment combinations of low water, low nitrogen, or presence of a pathogen. There is an environment column numbering the different combinations from 1-8. These combinations are either 1 or 0.  Then there is 3235 columns of different gene names. These are either A, B or NA. The for the rows we have a row indicating the chromosomes that the genes are on. There are 10 different chromosomes. There is another row that indicates the distance the gene is on the chromosome. The next 6672 rows are different tests with varying gene combinations and varying environmental combinations. The data is very similar however not the same for each environmental combination, and for some environmental combinations more tests have been done.

This is the hybrid sample dataset:

|stress     |Barcode   |month|BreedType|Genotype   |gpm27|tub1|gpm113b|gpm705a|gpm325a|dmt103b|gpm699d|gpm319|IDP1447|
|-----------|----------|-----|---------|-----------|-----|----|-------|-------|-------|-------|-------|------|-------|
|           |          |     |         |           |1    |1   |1      |1      |1      |1      |1      |1     |1      |
|           |          |     |         |           |0.01 |0.81|0.82   |0.84   |4.97   |6.89   |8.67   |14.82 |19.92  |
|1647.259062|201_B_1270|june |Hybrid   |PH207xMo010|B    |B   |B      |A      |B      |B      |B      |B     |B      |
|2138.410467|201_D_1272|june |Hybrid   |PH207xMo010|B    |B   |B      |A      |B      |B      |B      |B     |B      |
|1989.199039|206_A_1329|june |Hybrid   |PH207xMo017|A    |A   |-      |B      |A      |A      |A      |A     |A      |
|3174.316121|209_C_1367|june |Hybrid   |PH207xMo025|A    |A   |-      |A      |A      |A      |A      |A     |B      |
|387.3855291|212_A_1401|june |Hybrid   |PH207xMo029|A    |A   |A      |B      |A      |A      |A      |A     |A      |
|792.9015938|212_C_1403|june |Hybrid   |PH207xMo029|A    |A   |A      |B      |A      |A      |A      |A     |A      |
|1670.068523|221_B_1510|june |Hybrid   |PH207xMo055|B    |B   |B      |A      |B      |B      |B      |A     |A      |
|4797.379163|221_C_1511|june |Hybrid   |PH207xMo055|B    |B   |B      |A      |B      |B      |B      |A     |A      |
|1520.648232|231_A_1629|june |Hybrid   |PH207xMo276|A    |A   |A      |A      |A      |A      |A      |A     |A      |
|6226.595893|233_B_1654|june |Hybrid   |PH207xMo288|A    |A   |A      |B      |A      |A      |A      |A     |A      |
|6824.756723|233_C_1655|june |Hybrid   |PH207xMo288|A    |A   |A      |B      |A      |A      |A      |A     |A      |
|643.1691161|236_D_1692|june |Hybrid   |PH207xMo311|A    |-   |A      |A      |A      |-      |A      |A     |B      |
|1490.150142|242_A_1761|june |Hybrid   |PH207xMo354|A    |A   |A      |A      |A      |A      |A      |A     |A      |

This is inbred sample dataset:

|stress     |Barcode  |month|BreedType|Genotype|gpm27|tub1|gpm113b|gpm705a|gpm325a|dmt103b|gpm699d|gpm319|IDP1447|
|-----------|---------|-----|---------|--------|-----|----|-------|-------|-------|-------|-------|------|-------|
|           |         |     |         |        |1    |1   |1      |1      |1      |1      |1      |1     |1      |
|           |         |     |         |        |0.01 |0.81|0.82   |0.84   |4.97   |6.89   |8.67   |14.82 |19.92  |
|1852.195947|10_A_65  |june |Inbred   |Mo010   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|1767.175217|10_C_67  |june |Inbred   |Mo010   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|1479.159326|10_D_68  |june |Inbred   |Mo010   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|3865.705955|162_A_941|june |Inbred   |Mo355   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|4661.55073 |162_B_942|june |Inbred   |Mo355   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|2387.588833|162_D_944|june |Inbred   |Mo355   |B    |B   |B      |A      |B      |B      |B      |B     |B      |
|1662.058585|166_C_967|june |Inbred   |Mo360   |A    |A   |A      |B      |A      |-      |A      |B     |B      |
|234.0880247|170_A_989|june |Inbred   |Mo365   |-    |A   |-      |A      |A      |A      |A      |A     |A      |
|4751.784865|170_B_990|june |Inbred   |Mo365   |-    |A   |-      |A      |A      |A      |A      |A     |A      |
|775.0447083|170_C_991|june |Inbred   |Mo365   |-    |A   |-      |A      |A      |A      |A      |A     |A      |
|1341.596692|170_D_992|june |Inbred   |Mo365   |-    |A   |-      |A      |A      |A      |A      |A     |A      |
|120.6342255|19_A_125 |june |Inbred   |Mo017   |A    |A   |-      |B      |A      |A      |A      |A     |A      |
|2568.814877|19_B_126 |june |Inbred   |Mo017   |A    |A   |-      |B      |A      |A      |A      |A     |A      |
