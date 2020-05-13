---
title: "Using GitHub"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# Using Github Summary
Github becomes a very useful file sharing and storage system. In Github one can upload local files to the internet for other people to access. Github will override old files with new files such that one can always use the latest versioin. Github uses a push and pull system to upload and download the file.
To upload the files on Github, one can use the terminal on their computer. The first thing to do is to change the directory to the directory they want to upload on Github. Typically, one will have a folder on Github with the same name. Furthermore, one needs to clone that repository to Github. This can easily be done by copying the file path into Github where it says clone repository. Thereafter, they will type in their terminal *(git add .)*, then type *(git commit)*. After that, a screen will come up asking the user to name the upload. Type *(i)* to insert and then the user can type the name of the upload. To save and quit press *(esc)* then *(:w)* to save and *(:q)* to quit. Then type *(git push)*. this will update the git reposirtory and make all the updated work in that folder available to other people in the project.
To download the files from Github, the user needs again to have a destination folder cloned as a Git repository on their computer. Then change the directory to that directory and type *(git pull)* and all the files from that repository will be on copied to their computer. It is important to note that a specific github location/user must be cloned to a specific directory on the computer. It is not like once a directory is a github repository for one user, that it is a github repository for all.
Github is designed to work directly with members in group with several directories cloned for regular updating. R projects can also be cloned and directly linked to Github. On R Studio one can directly push project files to github without using the terminal.

## Intro

The files necessary for using EnviroTyping are located in a repository; therefore, it will be increasingly helpful over time if you create your own account on [GitHub](https://github.com). This allows you to easily access not only the EnviroTyping datasets available at the time of this writing, but also any changes made over time. If you have used GitHub previously, you can skim (or simply skip) this tutorial and jump to the information you need. However, if you are new to the whole experience following the below outline will provide you with a strong foundation to understand the tools GitHub provides.

## Creating the Account and Forking

This first step is straight-forward. Once you arrive on [GitHub](https://github.com), simply create an account by inputting a few identifying details, such as your email address. After your account is created you will be taken to your personal page. Here you can see your repositories, other projects you follow, and personal details. Feel free to customize it or add a picture because this account is yours, and becoming active on GitHub is a great networking opportunity as it connects you to the global research community.

Now, direct your attention to the search bar on the top left-hand side of the portal. In order to find the EnviroTyping repository, search for "TACC/EnviroTyping". There should be only one result, so click on its name. A quick glance at the directory may be overwhelming at first, but do not worry about navigating the GitHub just yet. Instead, look near the top of the portal for an icon that says "Fork". Click on it to have control over your own copy of the EnviroTyping repository. Verify the repository was forked by clicking on your profile icon on the top right-hand side of the portal, going to "Your Repositories", and checking that EnviroTyping is an option in the list.

This step is important because "forking" allows you to read and write the files in the master repository without fear of losing the original data. It also grants you the opportunity to contribute to the project by submitting "Pull" requests when you create (or edit) a file that helps further the project's goals.

## Navigating Directories

As previously alluded to, there are a plethora of directories within the EnviroTyping GitHub. Needless to say, it can be confusing when trying to navigate the massive number of files. In order to alleviate some of the initial confusion when navigating the files, below are a couple of tables which briefly describe some important folders in the directory and imply whether or not they are commonly used. However, do keep in mind these directories may change and these tables are meant solely as a guide. The first table describes overall parent folders.

|Directory|Description|Commonly Used?|
|---------|:---------:|:------------:|
|.ipynb_checkpoints|Holds one file|No|
|Premium @ 6750f8a|Collection of files for PReMiuM package|No|
|Simulated-data-for-performance-testing|Collection of files for simulated datasets|No|
|data|Central hub for data pertaining to each year's hybrids|Yes|
|doc_files|Collection of Markdown files for ReadtheDocs|Yes|
|Notes|Various notes|No|
|References|Collection of various articles for research|Yes|
|sandbox|Collection of data and files in current use|Yes|
|tools|Some benchmark functons|No|

You will find yearly GxE hybrid data in the folder entitled "data". Immediately within this folder are two more: "external" and "interim". The former consists of files and raw data used to make the cleaned data in the "interim" folder. Often you will utilize the data found in the "interim" folder. Within it, you will see more folders subset by years. Within each year folder, you may see a small selection of datasets with different suffixes in their names. Each dataset is slightly different. For example, the files with the suffix "wth_nas" have missing observations; and those with "shifted" move observations from later months to earlier months to remove any such missing observations. Having so many datasets with such minute differences is necessary to our research as it helps us establish any inconsistencies across different statistical tests.

It should be noted, however, the most commonly used files (and, consequently, most work) are found within "sandbox". Within "sandbox" are several folders, but below are the most frequently used.

|Directory|Description|
|---------|:---------:|
|posthoc_group_analysis|Various years' files used to create the output in the "Post hoc Anaylsis" section of the ReadtheDocs|
|shifted_data_analysis|Various years' workflows with different workflows for minimum, mean/median, or maximum variables|
|working_with_plots|Collection of files utilized to make different figures to visualize EnviroTyping data|

Your specific tasks will determine where you find (and place) your work. Nonetheless, feel free to navigate the GitHub and utilize any file as needed. This part of the tutorial is meant solely to help you understand the basic structure of the EnviroTyping GitHub. Understanding the nature of the directories and where to find the most important files will help tremendously when you need to reference multiple datasets or simply find a figure.

## Setting EnviroTyping Project in RStudio

Most of the code for EnviroTyping is written in R therefore, you should consider creating a "Project" in RStudio that you utilize when working with EnviroTyping files. Projects allow you to consistently work in a specific working directory on your local machine, easily connect to GitHub repositories, and streamline your work. If you do not already have RStudio downloaded on your computer, go to RStudio's [website](https://www.rstudio.com/products/rstudio/download) and select the installer for your OS.

To first create a project you will need the URL of your forked EnviroTyping repository. In GitHub go to your EnviroTyping repository. Click on the green button that reads "Clone or Download" on the right-hand side of the page, then press the "Copy" button beside the URL that appears.

Next, in RStudio go to "Projects" on the top right-hand side of the screen. Select "New Project", then "Version Control", then "Git". From here, paste the URL of your repository into the "Repository URL" field and type "EnviroTyping" as the project directory name. You may choose where you place the subdirector; this process will download the GitHub repository to your local machine, you are encouraged to choose a parent directory where you may easily find the EnviroTyping files. Finally, you may also select whether or not to have RStudio open a new session when you work in the Project (most opt for the new session).

Now, the Project setup in RStudio is complete! When you need to use R for EnviroTyping tasks, simply open RStudio, go to "Projects", and select "EnviroTyping". From there, you can choose to navigate the EnviroTyping files in the directory on the bottom left-hand side of the screen (depending on your selected layout of RStudio), create new files, or do whatever else you need.

## Cloning Repository into Stampede

You are almost finished! Now you need to ensure your portal on the Stampede servers has cloned the data from your personal GitHub fork of the EnviroTyping repository. Assuming you have already created an account with TACC, cloning your repository is straightforward.

First, you must be connected to a secure Wi-Fi or ethernet network. If you are not, TACC will not allow you to connect to the Stampede servers. Next, open up Terminal (Mac OS/Linux) or Command Line (Windows), and log onto the servers by SSH. Type into your command line prompt `ssh <TACC username>@stampede2.tacc.utexas.edu`. After inputting your password the system will ask you to verify your login information by sending you a text or email, depending upon your settings, with a unique code. Inputting the code correctly will finalize the log-in process. Make sure you copy the URL from your forked repository, then type: `git clone <repository URL>`. The URL will be the same as the one used for the R Project previously. You will know your command was accepted if Terminal or Command Line tells you that files are being downloaded. Once your repository is cloned, you may now submit jobs to be processed on the HPC. If you are unfamiliar with the SLURM commands used to communicate with Stampede2, feel free to reference the "Reminders and Useful Commands" section of the ReadtheDocs.
