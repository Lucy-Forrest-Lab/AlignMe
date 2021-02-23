# AlignMe Installation Instructions

The main AlignMe code is quite straightforward to compile. However, to run in P, PS or PST modes require inputs from PSI-BLAST sequence searches, PSIPRED secondary structure predictions and/or OCTOPUS transmembrane span predictions. Installing those can be challenging. Please see section below on [Optional external code](#optional-external-code)

## Required external libraries
[Boost](https://www.boost.org) library headers, namely shared_ptr.hpp

## How to install
- change directory to the appropriate folder in your terminal 
- type `make`
- the `alignme` executable is created in this folder and ready to use 

## Optional external code
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - for generating position specific substitution matrices  
[PSIPRED](http://bioinf.cs.ucl.ac.uk/software_downloads/) - for generating secondary structure predictions  
[OCTOPUS](http://octopus.cbr.su.se/index.php?about=download) - for generating transmembrane predictions  
Refer to `doc/INSTALL_all.txt`
