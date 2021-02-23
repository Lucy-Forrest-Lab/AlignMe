# AlignMe
Pairwise alignment of membrane protein sequences   

## Contributors
René Staritzbichler  
Marcus Stamm  
Kamil Khafizov  
Edoardo Sarti  
Lucy R. Forrest  
Giacomo Fiorin

## Required external libraries
[Boost](https://www.boost.org) library headers, namely shared_ptr.hpp

## Optional external code
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - for generating position specific substitution matrices  
[PSIPRED](http://bioinf.cs.ucl.ac.uk/software_downloads/) - for generating secondary structure predictions  
[OCTOPUS](http://octopus.cbr.su.se/index.php?about=download) - for generating transmembrane predictions  
Refer to `doc/INSTALL_all.txt`

## How to install
- change directory to the appropriate folder in your terminal 
- type `make`
- the executable is created in this folder and ready to use 

## How to use
- have a look at the [online manual](https://lucy-forrest-lab.github.io/AlignMe/) and try out some examples
- to run AlignMePST jobs, use `scripts/use_best_parameters.pl`. This script assumes that PSSM, secondary structure predictions, and transmembrane predictions are available.

## Online server
Our online server is available at [http://www.bioinfo.mpg.de/AlignMe/](http://www.bioinfo.mpg.de/AlignMe/)

If you have questions about the usage of this program or discovered a 
bug, then please post a message on github or write to: <AlignMe@rzg.mpg.de>

## How can I cite AlignMe?
If you have used AlignMe in your work, please cite:  

For pairwise alignments:  
[Stamm M., Staritzbichler R., Khafizov K. and Forrest L.R. 2013 PLoS ONE](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0057731)  

For the alignment of two MSAs:  
[Khafizov K., Staritzbichler R., Stamm M. and Forrest L.R. 2010, Biochemistry](http://pubs.acs.org/doi/abs/10.1021/bi101256x)  

## Release Notes:

#### Version 1.2
Inclusion of anchors as constraints

#### Version 1.1
- Position Specific Subsititution Matrices (PSSMs) supported as an alignment input
- ClustalW and Fasta-Format output of the aligned sequences supported
- extraction of two sequences from the alignment of two averaged MSAs possible
- gaps are now "?0" in the profiles - no more confusion of gaps with profile values having a 0
- added a perl script with the best parameters for alpha-helical proteins to the package
- improved error feedback

#### Version 1.0 
- initial release 
