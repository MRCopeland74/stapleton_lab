---
title: "qPCR Analysis"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# qPCR Analysis
## Introduction
In order to conduct our research, two factors must be assessed for each of our samples of maize: the breed type and the amount of stress-product present. Fortunately, the former was recorded upon collection. Determining the latter, however, is a much more involved process. To this end we use Quantitative Polymerase Chain Reaction (qPCR), a technique that utilizes spectroscopy to measure gene expression. Given a genetic sample, qPCR determines the amount of a specified RNA sequence by assessing the fluorescence in the DNA. The fluorescence is measured by qPCR in 40 cycles and the cycle in which the amount of fluorescence passes a specified threshold is known as the cycle threshold (Ct) value.

A brief outline of qPCR is to follow.

## Cycle Threshold
There are 40 measured cycles in the raw data. Through these 40 cycles the threshhold is supposed to undergo an s curve in which the threshhold increases quickly then platos. qPCR analysis determines which cycle has the maximum first derivative of the threshhold which gives us the point in which the threshhold increases the most. However, the data we have does not always accurately follow this pattern. First of all some cases have very large initial thresholds. This drastically changes what the curve and maximum threshold increase is. Therefore, in order to prevent this analysis from being corrupted, the first 10 cycles of all data is excluded. This should only remove outliers and should not corrupt output as the maximum first derivative will still be the same. The second problem is that it seems the data is right censored. This means that by the end of the 40 cycles we still have not reached the maximum first derivative -- which would instead be in cycle number >40. This poses a difficult question into how to analyse this and we question whether we should throw out all of the samples that have a ct values of 40.

## Calibrated Data
The calibrated data does not have a sample ID however it does have a starting quantity. This is the data that we are going to use to determine the unknown starting quantity and ultimately the stress product of the experimental data. The calibrated data becomes a control for the project.
The calibrated data had the starting quantity and then the all-products and test1 ct values which were calculated during the qPCR analysis. We also have data from three experiments: June, August, November. The all-products and test1 are paired based on starting quantity. We are also able to determine the ratio between them. A boxplot based on the starting quantity shows the averages decreasing as the starting quantity increases. We have run two sets of analysis: one is analysis by month and the other is to put all the data together. The analysis entails a catagorical regressional analysis.

## Experimental Data
The experimental data all has unique sample ID that pairs an all-products and test1 ct value from the qPCR analysis. From the calibrated data we can calculate the starting quantity by predicting the probability of the experimental data having the same starting quantity compared to the calibrated data. We can then choose the highest possibility and alter the answer based on how probable it is. This will give us an adjusted ct threshold value. This value is added to the test1 value (which is always less than the all-products value) and then subtracted from the all products value. We multiply by 1000 to put the result into femtograms and we are left with the respective stress product for each sample ID which we can run vQTL/dglm on.


## Outputs
## Plots
