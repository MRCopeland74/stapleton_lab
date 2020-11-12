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
[Hierarchical](https://osf.io/wqrsd/) and the [Sampling Plan](https://osf.io/mu5qs/) Data are available at this link. Here is a sample of the sampling plan dataset. The Hierarchical_exp_data_stress#3 is the most recent heirarchical entry.

## vQTL Data
The hyrbid and inbred data used for our vQTL analysis originates from `[FullHybvqtlinput.csv]`(https://github.com/AustinGratton/vQTL/blob/master/qPCR/FullHybvqtlinput.csv) and `[FullInbvqtlinput.csv]`(https://github.com/AustinGratton/vQTL/blob/master/qPCR/FullInbvqtlinput.csv) respectively. The list below details the factors included in these datasets.
 * Our phenotype of interest, stress gene amount (`stress`)
 * A unique identifier for each observation (`Barcode`)
 * The month that the sample was taken (`month`: *June, August, November*)
 * The breed type of the maize corresponding to the sample (`BreedType`: *Hybrid, Inbred*)
 * The genotype of the sample (`Genotype`)
The remaining 362 columns describe different gene names. These are either A, B or NA. The `Genotype` variable corresponds to a unique arrangement of these 362 genes.
For more information regarding `stress`, please see our [qPCR Analysis](https://stapleton-lab.readthedocs.io/en/latest/qPCR%20Analysis/).


Below is a truncated sample from the combined inbred and hybrid dataset. The first row after the headers is the chromosome number and the second is the distance on each chromosome. Some columns have been omitted below as they are irrelevant to the analysis process.

| stress    | BreedType | gpm27 | tub1 | gpm113b | gpm705a | gpm325a | dmt103b | gpm699d | gpm319 |
|-----------|-----------|-------|------|---------|---------|---------|---------|---------|--------|
|           |           | 1     | 1    | 1       | 1       | 1       | 1       | 1       | 1      |
|           |           | 0.01  | 0.81 | 0.82    | 0.84    | 4.97    | 6.89    | 8.67    | 14.82  |
| 2915.9323 | Inbred    | B     | B    | A       | -       | B       | B       | B       | -      |
| 5088.4966 | Inbred    | B     | B    | A       | -       | B       | B       | B       | -      |
| 5006.9534 | Inbred    | B     | B    | A       | -       | B       | B       | B       | -      |
| 4047.5168 | Inbred    | B     | B    | A       | -       | B       | B       | B       | -      |
| 4941.6801 | Inbred    | B     | B    | B       | B       | B       | B       | B       | B      |
| 3515.6499 | Inbred    | B     | B    | B       | B       | B       | B       | B       | B      |
| 2500.0000 | Inbred    | B     | B    | B       | B       | B       | B       | B       | B      |
| 3996.7554 | Inbred    | B     | B    | B       | B       | B       | B       | B       | B      |
