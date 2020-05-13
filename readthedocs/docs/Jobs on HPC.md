---
title: "Premium Manual"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

## The Batch Script and Jobs

This section of the tutorial will assume some basic knowledge of how to log into the HPC via `ssh` and to navigate directories via Command Line or Terminal. Feel free to go to the last section of the GitHub tutorial for a brief refresher of the `ssh` log-in process. Likewise, we will briefly discuss how to navigate directories in this tutorial, but more information is given in the "Reminders and Useful Commands" tutorial.

### Batch Script

Unlike a personal computer the HPC will not simply run an R or Python script. The high performance computer must be "told" what to do. This is because endusers have the ability to request how much and what type of resources are utilized when a job is submitted through the system. For example, you will have one option to choose how many cores the HPC uses to process your request and yet another option is to have the system email you when a job is submitted or completed. The batch script is what communicates your requests to the HPC. You may find an example batch script below:

`#!/bin/bash` (Never change this line) <br/>
`#` <br/>
`#` <br/>
`#SBATCH -J job_name` <br/>
`#SBATCH -o job_name.%j.out` <br/>
`#SBATCH -N number_of_nodes` <br/>
`#SBATCH -n total_number_of_MPI_tasks` <br/>
`#SBATCH -p process_queues` <br/>
`#SBATCH -t max_time` <br/>
`#SBATCH -A Envriotype` <br/>
`#SBATCH --mail-user=<email@host.com>` <br/>
`#SBATCH --mail-type=all` <br/>
`#------------------------------------------------------` <br/>
`mkdir -p output` <br/>
`Rscript --verbose ./script_to_be_run.R > ./output.Rout` <br/>

Notice, however, you must change a few generic identifiers with your own specifics. For example, "job_name" should be replaced with a name which is less than 8 characters long that will help you identify it in the queue. The expressions "number_of_nodes" (N) and "total_number_of_MPI_tasks" (n) must both be replaced with at least 1. The value "max_time" must be in the format hh:mm:ss, where hh is the number of hours, mm is the number of minutes, and ss is the number of seconds you think are needed for your job to complete. The expression "<email@host.com>" should be replaced with your preferred email address (without the arrows on the ends), and "script_to_be_run.R" should be replaced with the name of a file in the same folder as the batch script you wish the HPC to run. The only line which requires much thought is the one where you choose how to replace "process_queues". The typical selection is to replace it with "normal"; but you may also choose to select queues like "large", "long", or "development". There are several other options in addition to these, but know that the process queue you choose will determine when your job begins, how long it may run before being automatically terminated, and the number of compute resources available to you.

If you wish to know more about batch scripts, there are plenty of other tutorials online. One great tutorial is provided by Cornell University's Center for Advanced Computing. This [tutorial](https://www.cac.cornell.edu/education/training/StampedeJan2017/Envi.pdf) provides a good introduction to how Stampede2 operates and gives several great examples of batch scripts and other useful SLURM codes.

### Push and Pull Requests

When you wish to edit a batch script, open RStudio and go to the EnviroTyping Project created earlier. You can edit batch scripts in RStudio just like any other file, so it will be assumed you are working in this IDE. You may choose to edit an existing batch script or create a new one; ensure the file name is in the form of `file_name.sh` because the `.sh` extension is one way the HPC knows how to read the script.

Once you finish editing the batch script (or any other file), you must commit the change to your GitHub repository and submit a push request. First, open RStudio and look around the top-middle of the RStudio window for a smaller icon with a green "+", red "-", and grey "G" symbol vertically stacked.

<img src="https://github.com/TACC/EnviroTyping/blob/Tutorial_Additions/doc_files/docs/img/GitLogo.png" width="600"/>

Click on the symbol, then select "Commit". This starts the staging process. A new window should open where you select the files you want to stage, select all the files you want to merge into your repository and write a short message that describes the changes made in each file. Once all changes have been staged, click the green arrow that says "Push" to add the new changes to your GitHub repository. If there are no complications, your repository should be updated in less than a minute.

Next, you must make sure the repository cloned to your Stampede2 drive is updated as well. Log-in to Stampede2 as before and enter the command `git pull`. If there are no complications and few files to update, your command line editor will give descriptions of the files added, changed, or deleted in a matter of a few seconds. Make sure to update each drive ($home, $work, etc.) in Stampede2 as necessary.

This process is the same for any file -- not just batch scripts. Once your changes have been pushed and pulled without issue, you may continue to submit your batch script to the HPC for processing.

### Job Submittal Process

Assuming you are already logged into the system, there are three steps to submitting a job request to Stampede2. First, load the modules needed to run the job. Second, navigate to the directory of the batch script you wish to submit. Third, submit the batch script.

Loading modules requires the use of a single command: `module load <module>` or `ml <module>`. EnviroTyping requires the extensive use of R, so you will more than likely need to use `ml Rstats`. If need be, you can also load Python or other modules in similar manner.

Navigating directories via command line editors is simple. If you're unsure of the layout of the folders you can see all the folders and files on the current level of the directory by entering the command `ls`. You can then go down a level by typing the command `cd <folder_name>`. If you accidentally go into the wrong folder, the command `cd ..` will take you back up one level of your directory. You may repeat this process until you find the location of the batch script.

Finally, you are ready to submit the job request! Type `sbatch <file_name.sh>`. Stampede should output information about the job request, such as the request number and some other identifying information. From here you may do a number of things, such as check on the status of the request or even cancel the job. More information about your options is found in the "Reminders and Useful Commands" tutorial.

## Stampede Rules and Practices

One of the most important rules for using Stampede2 resources is to practice good citizenship. Once you create your TACC account, you become one of the thousands of users sharing the same compute power so it's best to start learning to "do unto others as you'd have them do unto you." The following are some basic rules to follow.

First, TACC recommends submitting job requests in $WORK. This is accomplished by entering the command `cdw` at the beginning of your session before you follow the steps to submit a job. Use $HOME (or `cdh`) to manage files. Don't submit an R (or other program) script without first submitting a batch script because you will inadvertently drain compute resources from the log-in nodes. Doing so will result in a nicely worded email from TACC.

Second, when requesting time in your batch script, make sure you request an appropriate amount of time for your job and the appropriate process queue. For example, if a job should complete in 2 hours, don't request 10 or ask for the "long" process queue. It would be best to ask for, maybe, 5 hours and the "normal" queue. This increases the chances your job will be accepted sooner and allows other resources to be properly allocated. Stampede2 uses a formula that "bills" users for the amount of resources they need for a job; however, for most users, the workings of the computation is irrelevant, so long as they know it exists. In other words, the more time and nodes you request the longer you will wait for the job to start processing. You can read more about the various queues and their properties in the [User Guide](https://portal.tacc.utexas.edu/user-guides/stampede2#running-queues).

There are many aspects to utilizing the HPC, so it never hurts to ask if you have any questions about the process or how to be a good citizen. Never hesitate to ask for help if you can't find what you need in the User Guide or another resource.

## RStudio Interface on HPC

If there comes a time when you'd like to work in the RStudio IDE while staying connected to Stampede and its compute resources, you can do that. Enter the command `sbatch /share/doc/slurm/job.Rstudio` after you are logged into Stampede, and give the HPC a moment to process the input. You will be given a URL to connect to the web interface of RStudio where you will be required to submit your TACC log-in credentials. After submitting them correctly you will then have access to the familiar IDE with a cleaner visualization of your directories and any output which may be produced. This is a helpful tool if you would like to keep an eye on any work as it is processed.
