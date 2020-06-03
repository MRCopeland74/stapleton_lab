---
title: "Updating the Docs"
author: "Michael Copeland"
date: "22/04/2020"
output: html_document
---

# Updating the Docs

Updating the documentation is easy and should be done as users discover useful tips and tricks along their own workflows. All documentation is stored on GitHub in plain-text at [EnviroTyping Docs](https://github.com/TACC/EnviroTyping/tree/master/md_doc_files/md_files).

## Accessing the source

Make a working copy of the documentation.

**Using command line via terminal**

From your working directory, download the project from GitHub::

`git clone https://github.com/TACC/EnviroTyping.git`

After a change has been made to the master repository, [ReadtheDocs](https://readthedocs.org) automatically builds fresh html documentation hosted on their servers.

**Using a desktop web browser**

Browse to [EnviroTyping Project on GitHub](https://github.com/TACC/EnviroTyping) and click "Clone or Download" at the right. You can use the [GitHub desktop app](https://desktop.github.com/) or work with the compressed folder using the text editor of your choice. For comprehensive edits you may wish to use the [Atom Editor](https://atom.io) with the markdown preview package enabled with the documentation directory selected as your project folder.

For more on MkDocs / Read the Docs, see:

- [MkDocs - Getting Started](http://www.mkdocs.org/#getting-started)
- [ReadtheDocs Instructions](http://docs.readthedocs.io/en/latest/)

## Forking & Committing Changes

Follow the standard git commit process to request changes. For a full introduction see:

- [Fork a Repo](https://help.github.com/articles/fork-a-repo/)
- [Pull request tutorial](https://yangsu.github.io/pull-request-tutorial/)
- [Git Novie](http://swcarpentry.github.io/git-novice/)

In short, you'll want to create a fork of the repository from the terminal or the GitHub project's website. The repository includes markdown source files which previewing in html on your own machine. Once you've finished with your proposed changes, add & commit the changes to your fork & open a pull request to the master branch at TACC/EnviroTyping/docs.

## How it Works

MkDocs is a fast, simple and downright gorgeous static site generator that's geared towards building project documentation. Documentation source files are written in Markdown, and configured with a single YAML configuration file.


Contact azg5169@uncw.edu
