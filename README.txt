AlignMe 1.1 README 

Required external libraries:
- boost library headers (boost.org), namely shared_ptr.hpp 

How to install:
- enter this folder in your terminal
- type "make"
- the executable is created in this folder and ready to use 

How to use:
- have a look at the manual AlignMe_manual.pdf and try out some examples

Online server:
- our online server is available at http://www.forrestlab.org/AlignMe 

If you have questions about the usage of this program or discovered a 
bug, then please write to: AlignMe@rzg.mpg.de 


Release Notes:

Version 1.1
- Position Specific Subsititution Matrices (PSSMs) supported as an alignment input
- ClustalW and Fasta-Format output of the aligned sequences supported
- extraction of two sequences from the alignment of two averaged MSAs possible
- gaps are now "?0" in the profiles - no more confusion of gaps with profile values having a 0
- added a perl script with the best parameters for alpha-helical proteins to the package
- improved error feedback

Version 1.0 
- initial release 
