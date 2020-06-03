---
title: "Reminders and Useful Commands"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

## Reminders

  First, recall that when utilizing the HPC, you are not actually processing information in your local machine. As such, be sure to not type commands only Command Line (on Windows) or Terminal (on Mac OS) recognizes. Stampede2 is a collection of Lustre servers that runs on Linux. The differences between languages may be small, but recognizing there could be a distinction will assist you should you ever need to look for extra help online.

  Second, it cannot be stressed enough that you must load modules anytime you want to run an analysis through Stampede2 *before* the batch script is submitted. As mentioned previously, you may do so by typing `module load <module>` or `ml <module>`. An example of this statement is `ml Rstats`. If you don't load the module, you may exit the HPC thinking your jobs are running, only to return at a later time to find the HPC did not know what to do with your batch script.

  Third, if your job requires the use of packages not already included in your base module, you must install the packages in the HPC first. It is easiest to do this by loading your preferred module, calling it, then installing the necessary packages. For example, if you need to install *tidyvese*, you would employ the following steps:

 `ml Rstats` <br/>
 `R # to open R` <br/>
 `install.packages("tidyverse")` <br/>
 `q() # to quit R` <br/>

  Once you finish installing the packages, you may continue to submit jobs as needed.

  Fourth, don't forget to push or pull updates from your Git respository. Many headaches are avoided if you get into the habit of typing `git pull` in the command line as soon as you push an updated file from your local machine and are logged into Stampede2.

 ## Other Useful Commands

  Navigating directories through a command line can be daunting the first couple of times you do it, but it becomes easier with time and experience. Here are some basic commands and their effects:

 `ls` returns the collection of folders and files in the current level of your directory, <br/>
 `cd <directory>` will take you to a lower level (subfolder) of the directory, <br/>
 `cd ..` will take you one level higher, <br/>
 `sbatch <batch script filename>` submits your job, <br/>
 `squeue -u <Stampede2 username>` tells you the details of *all* your current jobs, <br/>
 `scancel <job number>` cancels the requested job, <br/>
 `cat <job output filename>` reads the output file of a module (e.g., a .out or .Rout file), and <br/>
 `exit` disconnects your SSH session to Stampede2.

  There may be a time when you need to remove existing files or whole folders to free up space in your `$HOME` or `$WORK` directories. Doing so is straight-forward, but extreme caution should be practiced because data is not recoverable once removed using these methods.

 To remove a file: `rm <filename>`. <br/>
 To remove an empty folder: `rmdir <foldername>`. <br/>
 To remove a non-empty folder: `rmdir -R <foldername>`. <br/>
 To remove a non-empty folder but receive a prompt before removal: `rmdir -iR <foldername>`. <br/>

  It is easy to remove a file when you are in the proper directory already; but if you are not, simply specify the file's full address when using the above commands. Likewise, you may delete multiple files at once by listing the unwanted file names beside each other: `rm <filename1> <filename2>`.

  Command line programming is powerful. However, such power is hard to tame without the proper knowledge. If you ever need extra information for working with Stampede2, the Internet is a great resource; but the best help will be given by XSEDE's <a href=https://portal.xsede.org/web/xup/help-desk].>Help Desk</a>.
