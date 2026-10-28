# WARNING

We are currently working on creating a new and improved CIFcurator, this version of the code is unstable. If you want to use CIFcurator, use our latest stable release [v.1.0.1](https://github.com/fatipardo/CIFcurator/releases/tag/v.1.0.1). 

# CIFcurator 

With CIFcurator, users can create their own self-consistent structure dataset, enabling more straightforward comparison among structures. The library creates separate files for each biological assembly present in a structure file and standardizes chain names and numbering. Our goal is to provide researchers with a consistent dataset that facilitates their analysis.

## Table of contents

* [CIFcurator workflow and tutorial](#cifcurator-workflow)
* [Other tools](#other-tools)
* [Installation](#installation)
* [CIFcurator team](#cifcurator-team)

## CIFcurator Workflow

We have created Jupyter Notebooks that provide a step-by-step guide for creating a curated ensemble of structures using CIFcurator. 

![Workflow flowchart](./images/FlowChart.png)

### [Step 1. Download structural ensemble form RCSB PDB.](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step1.DownloadStructuralEnsembleFromRCSBPDB.ipynb)

Download all structures that match the name and sequence of your molecule of interest.

> **Note:** This notebook sometimes does not display on the Github website, download and open in your browser.   

### [Step 2. Clean Structures and Create one CIF file per biological assembly.](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step2.CreateOneCIFFilePerBiologicalAssembly.ipynb)

A CIF file may contain multiple biological assemblies within one asymmetric unit. In this step we separate these biological assemblies, and create one CIF file for each one. We also reduce the number of data blocks included in the CIF file, keeping those that the coordinates refer to, so that the files stay consistent with the mmCIF dictionary. (This step used to be two steps, 'clean' and 'simplify'; they are now one.)

### [Step 3.1. Assign MOLID to the entities found in the CIF files, version 1](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step3.1.AssignMolIDToEntitiesFoundInCIFfiles1.ipynb)

The script goes over all the CIF files and collects all entities. The user can decide what Mol ID to assign them. In this example, we show the case in which we give a different ID to each entity found.
This step is also important because it lists all the entities that were found in your ensemble, so it allows you to identify if there is a structure that doesn't belong. We show an example of this in this notebook.

### [Step 3.2. Assign MOLID to the entities found in the CIF files, version 2](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step3.2.AssignMolIDToEntitiesFoundInCIFfiles2.ipynb)

Same as Step 3.1, but in our example, we give the same MOL ID to different entities. You may want to do this for example, if you want to give the same MOL ID to all ligands, or water molecules. Doing this will trigger a concatenation menu, which we show how to use.

### [Step 3.3. Assign MOLID to the entities found in the CIF files, version 3](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step3.3AssignMolIDToEntitiesFoundInCIFfiles3.ipynb)

In this notebook we show how to perform concatenations and conversions by using a conversion file (useful when a structure contains many entities). 
We also show an option that allows users to keep a record of the changes introduced in this step (track old chain names, new chain names, entity names and file names).

### [Step 4. Chain ID standardization](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step4.ChainIDStandardization.ipynb)

Step 4 allows us to name each entity with whatever name we want. Step 4 makes sure that the chains that are the same (we do sequence alignment to determine similarity) in different CIF files, have a consistent name. Sometimes entities/chains are mislabeled in deposited structures, this step is recommended to identify any such cases. This step can also be used to identify any possible outliers, by seeing how all chains score compared to our reference. 

We divide the tutorial for this step in two parts. The second part shows how to generate the reference sequences, as well as showing how to load them when running the script. Doing this could also help speed up this step, as it allows to run the script in parallel in batches. This is particularly important when working with large datasets, or with molecules with many chains. 

### [Step 4.2 Chain ID standardization: generate reference sequences and how to load them](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step4.2.ChainIDStandardization.ipynb)

In this tutorial, we show how the reference sequence is selected by our script, and show how the user can modify it. It also shows how to load the reference sequences, creating the opportunity for running this step in parallel, in batches, speeding up the whole process. 

### [Step 5. Residue ID Standardization](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Step5.ResidueIDStandardization.ipynb)

Following step 4, now that we have consistent chain (entity) naming among all structures in the ensemble, we want to make sure that the numbering is also consistent (that the same residue position has the same number in all structures).

This is also the last step! You have a curated dataset!


> **Note:** There are more advanced curation steps and analysis that we will cover in future releases.

## Other tools

[Check project mini tutorial](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/CheckProject_CheckCreateDelete.ipynb). This mini tutorial can be run after doing step 2. `Check_project` checks if a directory has been created, if not it creates the directory and an info.txt file with the creation date. 

[Dataset Summary](https://github.com/fatipardo/CIFcurator/blob/master/Notebooks/Analysis.SummaryPDBDataset.ipynb).
This notebook can be run after doing step 0. It creates plots that summarize important information from your dataset such as organism of origin, resolution, year, and method used to solve the structure. The notebook also creates a pandas dataframe so users can create their own personalized plots.

## Installation

We recommend installing CIFcurator inside a virtual environment. We provide an `environment.yml` with the libraries you will need. 
We have tested the installation on MacOS and Ubuntu 24.04.1 LTS.
Ensuring you have the pre-requisites will facilitate the installation process!

### Pre-requisites

- Xcode including Command Line tools (for MacOS)
- Anaconda
- Create SSH Key and link it to your GitHub account

### Installation steps

1. Download CIFcurator from GitHub and install environment from YML file

>git clone git@github.com:fatipardo/CIFcurator

>cd CIFcurator

>conda config --remove channels defaults

>conda env create -f environment.yml

2. Activate environment and install CIFcurator

>conda activate CIFcurator

>pip install .

3. Install Jupyter Notebook kernel

> python -m ipykernel install --user --name CIFcurator --display-name CIFcurator


4. Running notebook:

> cd Notebooks

> jupyter notebook

- Open any notebook you would like to run.
- If Jupyter does not recognize the kernel, select ‘CIFcurator’ from the drop down menu.


## CIFcurator team

The code in this repository is based on the code found [here](https://test.pypi.org/project/PDBClean/#files).
The code was originally written by Frédéric Poitevin and Nicholas Corsepius.
Fátima Pardo Avila and Liv Weiner created this repository. Paulina Cabral contributed to the code and documentation.
We all worked on this project while being part of the Levitt Lab at Stanford University. This project used to be known as PDBClean and PDBCleanV2.
