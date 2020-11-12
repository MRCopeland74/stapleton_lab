---
title: "DGLM Analysis"
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

# DGLM Analysis

The data that we have accumulated is based on a double generalised linear model. Furthermore, the data is right censored as the maximum CP value obtainable is 40. This means that if a value of greater than 40 is obtained, it is recorded as 40. The stress is then calculated based off that CP value. It is double generalised because the CP value is based off a calculation and then the stress value is based off that calculation. The DGLM is therefore dependent on the varience of the first calculation.
