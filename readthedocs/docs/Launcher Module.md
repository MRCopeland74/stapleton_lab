---
title: "Launcher Module"
author: "Michael Copeland"
date: "23/04/2020"
output: html_document
---
# Launcher Tutorial

Launcher is a utility for performing simple, data parallel, high throughput computing (HTC) workflows on clusters, massively parallel processor (MPP) systems, workgroups of computers, and personal machines. `Wilson:2014:LSF:2616498.2616533`

This tutorial is assuming the use of the launcher module found in the TACC Lmod's and some familiarity with TACC HPC systems.

## Purpose
There are many types of jobs for running with the **launcher** module but for the scope of this tutorial we are going to practice with R scripts. You are going to generate some random datasets with R and print them to csv files. Then you will create a launcher file to run 100 jobs on 10 nodes.  Finally you will use the **launcher** module to run all 100 jobs with very little typing involved.

## Loading Launcher

Login to one of the TACC systems such as Stampede2.

`ssh username@stampede2.tacc.utexas.edu`

This logs you into your `$HOME` directory.  We typically run analysis in the `$WORK` directory.

Just type following to move to your `$WORK` directory `cdw`.

To see the modules you currently have loaded:

`module list`

To load the launcher module if not listed in the output:

`module load launcher`

You can optionally save the launcher module as a default every time you login to the system:

`module list` confirm launcher is loaded then: `module save`

We also need to use the Rstats module.

`module load Rstats`

## Neccessary Files

You need to copy 4 files to use this tutorial:
* *randcsvGen.R* - This is an Rscript that randomly generates 100 csv files using the mtcars dataset.
* *launcherFile.sh* - This bash script will parse all the csv filenames into an array and return them to *launcherFile*.  *launcherFile* is a list of  Rscript command line arguments to run each job.
* *launcher.slurm* - This is the file which actually batches the jobs to different HPC nodes.
* *headDF.R* - A simple R script which prints the top of the dataset read in from each csv file.

All these files should be in the same directory in your `$WORK` directory.

## Modify launcher.slurm

Set `LAUNCHER_JOB_FILE` to point to your job file. We are going to use Vim, a built-in editor.

Type `vim launcher.slurm`

Type `i` and the prompt at the bottom will say `Insert`.
Use your arrow keys and change `<path_to_directory>` to the directory you are currently in, which is the `$WORK` directory if you did not place the files somewhere else. Change `<allocation>` to your group's allocation for the system.
m
Then hit `Esc` and type `:wq` which tells Vim to save and quit.

## Generate Some Data

The first file you will run is the *launcherFile.sh*.  This will run the R script to generate the csv files and created your *launcherFile*.

Type `chmod +x launcherFile.sh` to give your user permission to execute.

Then type `./launcherFile.sh`

This might take a few minutes to install the required R packages generating a lot of output in the terminal.

You will have 3 new directorys created:
* *randcsv/* - This is where all the random csv files are stored.
* *output/* - For the results from the analysis.
* *error/* - There will be some general warnings that are read as error by R.

The *launcherFile* should also be in your current directory now as well.
## Submit the Jobs

Your ready to submit your batch of 100 jobs.

Type `sbatch launcher.slurm`

A new file will be created for the job:  *slurm-job##.out*

You can check the status of your jobs.

`squeue -u <username>`

When the jobs are done you should have 100 files in the *output/* and *error/* directories.  Remember you can use Vim if you want to open the file or `cat <filename>` to peak inside.


## Wrap-up

Again you can you **launcher** to run jobs with python or Agave apps among many other uses.  Hopefully you now have some understanding how to use the **launcher** module to batch many jobs at once.


#### References
https://github.com/TACC/launcher#referencing-launcher

https://portal.tacc.utexas.edu/user-guides/stampede2#using-modules-to-manage-your-environment
