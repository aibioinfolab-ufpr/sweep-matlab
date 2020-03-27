readme.txt, v1.1, 2020-feb

SWeeP - Spaced words Projection

Roberto Tadeu Raittz (rraittz) - 2018-2019
Mariane Goncalves Kulik (mgkulik) - v1.0, 2019-jan
Mariane Goncalves Kulik (mgkulik) - v1.1, 2020-feb
------------------------------------------------------------
INSTALLING SWeeP

MATLAB Compiler

1. Prerequisites for Deployment 

. Verify the MATLAB Compiler Runtime (MCR) is installed and ensure you    
  have installed version 7.17 (R2012a).   

. If the MCR is not installed, do following:
  (1) enter
  
      >>mcrinstaller
      
      at MATLAB prompt. This MCR Installer command displays the 
      location of the MCR Installer.

  (2) run the MCR Installer.

Or download Windows 64bit version of MCR from the MathWorks website:

   http://www.mathworks.com/products/compiler/
   
   
For more information about the MCR and the MCR Installer, see 
“Working With the MCR” in the MATLAB Compiler User’s Guide.    


NOTE: You will need administrator rights to run MCRInstaller. 


2. Files to Deploy and Package

Files to package for Standalone 
================================
-SWeeP.exe
-MCRInstaller.exe 
   -include when building component by clicking "Add MCR" link 
    in deploytool
-This readme file 

3. Definitions

For information on deployment terminology, go to 
http://www.mathworks.com/help. Select your product and see 
the Glossary in the User’s Guide.


------------------------------------------------------------
RUNNING THE CONSOLE PROGRAM

Run the program with the following command:

SWeeP fasta_path additional_parameters

Windows Users: As Windows allows folders names with spaces, please use PowerShell and quotes to avoid errors.


Additional input parameters
===========================

The following parameters can be informed. fasta_path is the only required parameter:

fasta_path		The place where SWeeP will find the fasta file or files. Three 
				options can be provided:
				1 - The folder with all the multifasta files to pre-process
				in the SWeeP expected format (There must be at 
				least 3 files in the folder), OR;
				2 - A file previously formatted by SWeeP (extensions .fasta, 
				.fas, .faa, .fa or .fna), OR;
				3 - A multifasta composed by one organism genes (extensions 
				.fasta, .fas, .faa, .fa or .fna). In this case SWeeP will not
				transform your file, which will be used exactly as informed.
SeqType			The entry type to define the correct default mask. Espected
				values are AA (Default) or NT.				
SWeePMethod		binary (Default), prime or count. The default binary is
				faster and presented better results.
				Please reference our paper to learn more.
GenerateTree	Options 0 (Default) or 1. Change this option if you do not 
				want to generate the .tree file using SWeeP. As this option
				uses the NeighborJoin clustering algorithm by default, the 
				complete SWeeP execution will take several more minutes to 
				finish, so it is advisable to disable it if you do not intend
				to use SWeeP to generate the phylogenetic tree.
				* If you choose not to generate the .tree file, 
				you can use the .csv file in any other tool capable of
				calculating the distances from the trees and execute any
				phylogenetic algorithm.
				** You will be able to open the .tree file using tools such as 
				(dendroscope.org) or (itol.embl.de).
ClusteringType	Options 1 (Default - NeighborJoin) or 2 (Ward). Visualization
				tools (such as dendroscope.org and R) were unable to open
				large trees generated with NeighborJoin (with more than 10 
				thousand sequences). For these cases, we suggest using Ward's
				clustering method.
DistanceMethod	euclidean (Default), spearman or cosine. This parameter will be
				used only if the option generate_tree be true (1). Please reference 
				our paper to learn why the Euclidean distance was chosen.
SaveHVD			Options 0 (Default) or 1. Please change this option to 1 if
				you intend to use MATLAB to process your data later.
ProjMatrix		By default SWeeP uses our pre-loaded quasi Ortho-normal Projection
				Matrix with size 160k x 600 (file with extension .SWeePproj). 
				Before change it, please reference our paper to understand 
				better how those matrices are generated. To change it, it is 
				necessary to pass an integer K, starting at 10, which will generate 
				a new projection matrix with size 160k x K. For personal computers 
				with 16 GB of memory, we suggest a maximum projection size of 
				1,000 (one thousand). This value can be increased when executing 
				in server.
SWeepMask		Those are the available mask options to AA*: 11011 (Default), 10111, 
				11101, 1101, 1011, 1001, 101; to NT: 111100000001111 (Default),
				111100001111, 111000111, 11000011, 110011, 101.
				* Please refer our paper to know more about the k-mer
				mask selected as default.

Examples
===========================

# 1 - SWeeP with fasta folder
SWeeP '\User\fasta_folder'

# 2 - SWeeP with fasta file
SWeeP '\User\fasta_folder\SWeeP_fasta_format.faa'

# 3 - SWeeP for NT
SWeeP '\User\fasta_folder' SeqType NT

# 4 - Change the SWeeP method
SWeeP '\User\file.faa' SWeePMethod prime

# 5 - Disable the Neighbor Join/Ward Tree generation
SWeeP '\User\fasta_folder' GenerateTree 0

# 6 - Change the distance method
SWeeP '\User\fasta_folder' DistanceMethod spearman

# 7 - Save the HDV to disk (for MATLAB users only)
SWeeP '\User\fasta_folder' SaveLVD 1

# 8 - Choose another projection size
SWeeP '\User\file.faa' ProjMatrix 800

# 9 - Inform a previous generated projection file
SWeeP '\User\file.faa' ProjMatrix '\User\proj.SWeePproj'

Output files
===========================

A folder named output will be created in the SWeeP folder. Check that the SWeeP has 
adequate permission to do so. The following files will be generated:

1 - multifasta.faa
This will be the multifasta transformed by SWeeP to meet its format requirements.
This file will not be generated if you enter a fasta file in the fasta_path parameter.

2 - random_ortho-normal_matrix.SWeePproj
This is the quasi Ortho-normal matrix generated using random values. The file will
be generated only if you decide to generate a new one, otherwise the default from 
the source folder will be used.

3 - file_LDV.csv - MAIN OUTPUT
This is the SWeeP main transformation file. You can to use it as the entry for tools
capable of calculating distances between sequences and executing phylogeny algorithms
or as features for machine learning algorithms.

4 - SWeePlog.txt
This is the log file updated every time SWeeP is executed, with errors, warnings and
updates.

5 - file.tree - Requires the parameter GenerateTree = 1
This is the phylogenetic tree file, generated using Neighbor Join or Ward clustering
method. You can open the .tree file using tools such as (dendroscope.org),
(itol.embl.de) or any other phylogeny tree viewer you choose.

Optional outputs
================

5 - file_HDV.mat - Requires the parameter SaveHVD = 1
This is the high dimension vector sparse matrix (HDV) in MATLAB format. You can load it
in MATLAB to analyze the sequences information before the projection.