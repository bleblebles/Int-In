<img src="https://intein.biologie.uni-freiburg.de/images/logo.svg" data-canonical-src="https://intein.biologie.uni-freiburg.de/images/logo.svg" width="150" height="150" />

# Int&in website
The Int&in website can be accessed here: *https://intein.biologie.uni-freiburg.de/*

# Int&In standalone tool source code
This repo contains all files for replicating the results from the Int&in paper. 
The tool only works under Linux, but can easily be used under Windows through wsl.

This tool is not actively being maintained!

# Citation
If you use this software or any of it's data please cite: [Int&in: A machine learning-based web server for active split site identification in inteins](https://doi.org/10.1002/pro.4985)

# Required programs
For the functionning, the following tools need to be installed: hmmer (3.3), rate4site (3.0.0.), muscle (v3.8.1551), dssp (3.0.0), pdb2pqr (v2.1.1), and dotnet runtime (.NET 8.0)

For Ubuntu use:
```
sudo apt install hmmer
sudo apt install rate4site
sudo apt install muscle
sudo apt install dssp
```
For pdb2pqr version 2.1.1 can be downloaded from here: https://github.com/Electrostatics/pdb2pqr/releases/download/v2.1.1/pdb2pqr-linux-bin64-2.1.1.tar.gz

The dotnet runtime (.NET 8.0) can be installed from here: https://learn.microsoft.com/en-us/dotnet/core/install/linux-ubuntu#2110

For homology search Uniref90 2021_03 is needed (http://ftp.ebi.ac.uk/pub/databases/uniprot/previous_releases/release-2021_03/uniref/uniref2021_03.tar.gz), extract several times and convert the Uniref90.xml with the following perl script: https://proteininformationresource.org/download/uniref/xml2fasta/ to a fasta file.

Additionally python needs to be installed with the following libraries (Versions that were used for the paper given in brackets): biopython (1.81), matplotlib (3.7.0), pandas (1.5.3), numpy (1.23.5), scipy (1.10.0), scikit-learn (1.2.1), seaborn (0.12.2), xgboost (1.7.6), and imbalanced-learn (0.10.1)

# How to use
The standalone tool Int&in tool may be used under Linux from the terminal:
```
dotnet Backend/SplitProteinPrediction/bin/Release/net8.0/publish/Intin.dll InteinFile OutputFolder PrefixOutput
```
**InteinFile**: is the .pdb file to be investigated

**OutputFolder (optional)**: The folder where the output files will be generated

**PrefixOutput (optional)**: A prefix for the file names for the output files

You'll have to modify the config file too at `Backend/SplitProteinPrediction/bin/Release/net8.0/publish/Intin.dll.config`

It's contents look like this:
```<?xml version="1.0" encoding="utf-8" ?>
<configuration>
  <appSettings>
    <add key="Conservation_Indexes" value="-1;0;2" />
    <add key="RelASA_Indexes" value="-2;0" />
    <add key="CFragDocking_Residues" value="8" />
    <add key="Pdb2pqr" value="/mnt/d/linux_pdb2pqr/pdb2pqr" />
    <add key="HmmerUniref90" value="/mnt/d/hmmer/Uniref90/uniref90.fasta" />
    <add key="CurrentCulture" value="en-CA" />
  </appSettings>
</configuration>
```
Pdb2pqr and HmmerUniref90 need to be set to the correct paths.

# Replicating of the model & Data set

The model files and figures may be generated through the python notebook (`Int&In Model.ipynb`)

The file `final_training_data.csv` contains the output from the standalone tool, `ReadEfficiencies.xlsx` contains the efficiencies of the training data and the testing data is in `LiteratureDataset_final.csv`

# Westernblot data

All westernblot images are in the respective folder as well as the Supplementary Excel File compiling all images (`FileS2_Westernblots_SourceData.xlsx`)
