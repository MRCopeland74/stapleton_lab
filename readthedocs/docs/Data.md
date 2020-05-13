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

## qPCR Data
Not sure what this looks like

## vQTL Data
The [data](https://github.com/MRCopeland74/stapleton_lab/blob/master/vQTL/ManchingStressData_Covar.csv) we used for the vQTL data was all in one csv file known as the Manching Stress Product Data. This file conatained a phenotype which is the height of the corn crop. There are 8 columns of different environment combinations of low water, low nitrogen, or presence of a pathogen. There is an environment column numbering the different combinations from 1-8. These combinations are either 1 or 0.  Then there is 3235 columns of different gene names. These are either A, B or NA. The for the rows we have a row indicating the chromosones that the genes are on. There are 10 different chromosones. There is another row that indicates the distance the gene is on the chromosone. The next 6672 rows are different tests with varying gene combinations and varying envronmental combinations. The data is very similar however not the same for each environmental combination, and for some environmental combinations more tests have been done.

|Height|Low.Water|Low.Nitrogen|Pathogen|Low_W_N|Low_W_P|Low_N_P|All|None|Env|gpm27|tub1|gpm113b|gpm705a|gpm325a|dmt103b|gpm699d|gpm319|IDP1447|bnl5.62a|
|------|---------|------------|--------|-------|-------|-------|---|----|---|-----|----|-------|-------|-------|-------|-------|------|-------|--------|
|      |         |            |        |       |       |       |   |    |   |1    |1   |1      |1      |1      |1      |1      |1     |1      |1       |
|      |         |            |        |       |       |       |   |    |   |0.01 |0.81|0.82   |0.84   |4.97   |6.89   |8.67   |14.82 |19.92  |21.12   |
|71    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|79    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|48    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|85    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|36    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|72    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|58    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|76    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|68    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|76    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|81    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|41    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |A      |-      |B      |B      |B      |-     |A      |A       |
|90    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|75    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|82    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|79    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|80    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|83    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|85    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|65    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|83    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|61    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|86    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|53    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|40    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|52    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
|84    |0        |0           |0       |0      |0      |0      |0  |1   |1  |B    |B   |B      |B      |B      |B      |B      |B     |B      |-       |
