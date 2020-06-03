---
title: "vQTL Analysis"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# vQTL Analysis and Scanonevar function
## Introduction
The goal of the vQTL analysis is to run vQTL using the scanonevar functionon on the set of data of corn crops. This process will give the LODs and p-values of each respective gene for the mean, variance and mean-variance. The high LODs or low p-values are the genes that are most significant and most useful to the project.
## Scanonevar Function
The Scanonevar function was developed by Robert Corty and is am extension of the Scanone function. It has the main goal of cross examining large data bases to to determine LODs and P-Values. The scanone function can easily be run on a local computer while the scanonevar function can take hours to run.
### Interactive model
For our project we need the interactive model. The additive model, y = b0+b1x1+b2x2...bnxn is not accuarate enough for our model as we need the combinations on how each environmental factor interacts with each gene. So instead we are using the interactive model which is y = b0+b1x1+b2x2+b3x1x2.
## Scanonevar.perm Function
This function runs numerous permutations of the scanonevar fucntion to identify the accuracy of the scanonevar function. This function, however, takes a lot longer to run and is not doable on a home computer and requires a super computer. We send these to TACC which is the HPC in Texas.
## Effectsizes function
## Outputs
## Plots
