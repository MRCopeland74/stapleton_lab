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
Heirarchical output:

|    | sampleID.exp | test1.exp | allP.exp | ztest1             | zallP               | month | june | aug | exp.adjust       | exp.adjustTest1  | stress           |
|----|--------------|-----------|----------|--------------------|---------------------|-------|------|-----|------------------|------------------|------------------|
| 1  | 10_A_65      | 32090     | 39090    | -0.605598971136247 | 0.108396497246005   | june  | 1    | 0   | 1599.26301534119 | 33689.2630153412 | 5400.7369846588  |
| 2  | 10_B_66      | 35330     | 40000    | 0.754359904826838  | 0.360371525404632   | june  | 1    | 0   | 1798.6137576376  | 37128.6137576376 | 2871.3862423624  |
| 3  | 10_C_67      | 31410     | 38550    | -0.891022438930971 | -0.0411271458371387 | june  | 1    | 0   | 1562.788114847   | 32972.788114847  | 5577.211885153   |
| 4  | 10_D_68      | 32030     | 38510    | -0.630783394765194 | -0.0522029712507045 | june  | 1    | 0   | 1570.78464848193 | 33600.7846484819 | 4909.21535151807 |
| 5  | 162_A_941    | 30870     | 38220    | -1.11768225159149  | -0.132502705499058  | june  | 1    | 0   | 1502.52768283218 | 32372.5276828322 | 5847.47231716782 |
| 6  | 162_B_942    | 30780     | 39040    | -1.1554588870349   | 0.0945517154790463  | june  | 1    | 0   | 1552.52587247149 | 32332.5258724715 | 6707.47412752851 |
| 7  | 162_C_943    | 30790     | 38390    | -1.15126148309675  | -0.0854304474914021 | june  | 1    | 0   | 1507.15631873733 | 32297.1563187373 | 6092.84368126268 |
| 8  | 162_D_944    | 30780     | 38570    | -1.1554588870349   | -0.0355892331303548 | june  | 1    | 0   | 1520.52530965986 | 32300.5253096599 | 6269.47469034014 |
| 9  | 166_A_965    | 33190     | 39450    | -0.143884537938905 | 0.2080789259681     | june  | 1    | 0   | 1762.11799353673 | 34952.1179935367 | 4497.88200646327 |
| 10 | 166_C_967    | 33160     | 39390    | -0.156476749753379 | 0.19146518784775    | june  | 1    | 0   | 1761.35250936969 | 34921.3525093697 | 4468.64749063031 |

The most important value is the stress that has been calculated. This will allow us to attach it the inbred/hybrid dataframes by sampleId and run analysis on it.

## Plots
