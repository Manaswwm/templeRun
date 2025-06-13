<p align="center">
<img src=https://github.com/user-attachments/assets/fd7990fb-7b85-4398-a808-3dc21be3ad8a width=600 height=140 align=center>
</p>

# TEMPLERUN
TEMPLERUN is a comparative and population genomics tool that performs in-depth analyses of the action of constraint and adaptation on the **non-coding TFBS** elements within the genomes. These elements are identified using [TEMPLE](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12535)

## Prerequisites

TEMPLRUN is exclusively written in R, with occasional internal calls to external tools (TEMPLE, vcftools, bcftools and vcfkit) that are to be preinstalled. The R libraries required for running TEMPLERUN are listed in **temple_run.R**.

## Structure

The tool's centre point is **temple_run.R**. This script is the import point for libraries and the species-specific input files that will be used later for the analyses. 
Additionally, this script is also the call point for scripts, stored in **includes** folder, that perform sequential processing of the input files.

## How to run

The input points for the species-specific files are documented in **temple_run.R**. These can be replaced with input files for other species. Example input files are listed in the folder __input_files__ . 
The resulting output files are listed in __example_output__ folder. Scripts used for performing analyses are listed in __analysis__ folder

## Questions?

If you have any questions, then you can reach me here - manaswwm@gmail.com :sunglasses:
