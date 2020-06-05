---
title: "Using R"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# Using R
R is a very convenient program for running cross analysis for large data pools. It is even more effective at running regression models with large matrices and tables. R is also very effective at reading in csv files and csv cross files. It also has many specific and quick table functions like `tapply()` and `sapply()` which are much more efficient than using for loops. Finally, R is also the language in which the vQTL package is created in.
### R Studio
R Studio is the IDE we use to run R in. It provides a much more user friendly interface for R. For example, one can have several tabs open to work between or hold different files and one can create and link new projects to Github for easy updating. It also remembers functions and variable names to increase calculation speed. It can create plots effectively which can easily be saved onto a computer. Its help window on the bottom right hand corner is also very useful and can link to other options that one may not know exist.
### Tidyverse package
Using the `tidyverse` is not a requirement, but for data manipulation they provide standards that can be used across any type of R project and Hadley Wickham (the tidyverse core contributor) diligently researches best software development practices. Install via ```install.packages("tidyverse")``` and check out the following free ebook by Hadley Wickham here [R for Data Science](https://r4ds.had.co.nz). Running the ```library(tidyverse)``` command will load several packages used commonly like tidyr, dplyr, readr, and purrr.
