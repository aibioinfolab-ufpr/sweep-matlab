# SWeeP matlab
Source code for the tool SWeeP (Spaced Words Projection Tool) originally developed in Matlab.

The SWeeP method was developed to favor the analizes between amino acids sequences and to assist alignment free phylogenetic studies. This method is based on the concept of sparse words, which is applied in the scan of biological sequences and its the conversion in a  matrix of ocurrences. A quasi-orthonormal matrix is then projected onto it, generating a low dimensionality matrix capable of representing the original sequences.


## Citation
If you use this code please cite the corresponding paper where original method appeared:

De Pierri, C.R., Voyceik, R., Santos de Mattos, L.G.C. et al. SWeeP: representing large biological sequences datasets in compact vectors. Sci Rep 10, 91 (2020). <a href="https://doi.org/10.1038/s41598-019-55627-4">[paper]</a>


## Version
1.0.1


## Running the matlab source
To properly run the source code in the Matlab ID you must add all this structure to the path.
The "core code" is separate from the support code, generated to support all the user interfaces required by the executable version.

By default, **when running in the Matlab IDE, SWeeP does not generate any external files.**

If you expect to generate all the output files using the Matlab IDE, you must create a GLOBAL variable to set the developer mode to FALSE, as follows:

setGlobalDev(0);
global DEV_MODE

Please refer to the readme.txt file to understand the tool


## The SWeeP Standalone Executable
Please refer to our Source Forge repository <a href="https://sourceforge.net/projects/spacedwordsprojection/">SWeeP - Spaced Words Projection tool</a> to download the Linux and Windows executable versions and Default Projection Matrices.


## Dependencies
The Parallel Computing Toolbox may be used to increase SWeeP performance but it is not required.
