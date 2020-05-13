---
title: "HPC Setup and Summary"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# HPC (High Performance Computers) Summary
As discussed in the Scanonevar.perm section, some functions when running giant data sets cannot be run on a home computer. Therefore, they need to be sent to an HPC. We use the TACC (Texas Advanced Computing Center). We are able to login to the HPC remotely and send files via githup to the HPC. We setup a github in the HPC and then our own Github then pull and push

## Why HPC?

Sometimes data can get so complex or large that a personal computer is no longer a viable option to run the analysis. Luckily there are high performance computing clusters available which cater to different requirements, from raw CPU or GPU power to large memory needs. These clusters can be utilized to handle analyses using R, Python, or other languages. The National Science Foundation has partnered with universities and organizations around the country to create HPC clusters for students at universities to freely use (in moderation) via [XSEDE](https://www.xsede.org).

## Necessary Accounts

There are a couple platforms essential to the mission of the EnviroTyping research, namely Cyverse and XSEDE/TACC. You will be asked to create an account on both platforms.

### Cyverse

Cyverse is an initiative led by the University of Arizona to improve how data is stored and shared amongst researchers. You may read more about its mission [here](http://www.cyverse.org/about) and create an account on [this page](https://user.cyverse.org/register).

### XSEDE/TACC

XSEDE is a virtual organization that serves to connect researchers to resources necessary to navigate the shear amount of data in the modern age. More background information is found on their [About](https://www.xsede.org/about/what-we-do) page. Create an account [here](https://portal.xsede.org/?p_p_id=58&p_p_lifecycle=0&p_p_state=maximized&p_p_mode=view&saveLastPath=0&_58_struts_action=%2Flogin%2Fcreate_account). An XSEDE account will allow you to access any resource in the XSEDE network.

Likewise, one such resource is the Texas Advanced Computing Center (TACC). As with the previously mentioned programs, more information is found on their [About](https://www.tacc.utexas.edu/about/overview) page. TACC holds the tools used to complete jobs too large to run on personal or local machines (such as the PReMiuM profile regressions discussed in the "Workflow" section of this ReadtheDocs). Specifically, you will most likely utilize the new Stampede2 servers. However, you must first create an
[account](https://portal.tacc.utexas.edu/account-request?p_p_id=createaccount_WAR_createaccountportlet&p_p_lifecycle=1&p_p_state=normal&p_p_mode=view&p_p_col_id=column-1&p_p_col_count=1&_createaccount_WAR_createaccountportlet_action=continue). It is encouraged you carefully read each page in the creation process thoroughly because TACC tends to be strict in its allocation of compute resources, and the more you know in the beginning the more it will help in the long run. We will discuss how to properly use Stampede2 in more detail, but feel free to briefly peruse its [User Guide](https://portal.tacc.utexas.edu/user-guides/stampede2) once your account is created.
