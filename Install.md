## AlignMe Installation Instructions

The main AlignMe code is quite straightforward to compile. However, to run in P, PS or PST modes require inputs from PSI-BLAST sequence searches, PSIPRED secondary structure predictions and/or OCTOPUS transmembrane span predictions. Installing those can be challenging. Please see section below on [Optional external code](#optional-external-code)

### Required external libraries
[Boost](https://www.boost.org) library headers, namely shared_ptr.hpp

### How to install
- change directory to the appropriate folder in your terminal 
- type `make`
- the `alignme` executable is created in this folder and ready to use 

### Optional external code
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - for generating position specific substitution matrices  
[PSIPRED](http://bioinf.cs.ucl.ac.uk/software_downloads/) - for generating secondary structure predictions  
[OCTOPUS](http://octopus.cbr.su.se/index.php?about=download) - for generating transmembrane predictions  

---

### Complete installation instructions
If you wish to install all the dependencies, the following instructions should be of some help. However, some of these links are external code, so might be broken or difficult to install. Please let us know if you have any problems. 

#### Step 1: Preparations
Open your terminal and enter bash mode  
`bash`

Check if you have these programs installed or not:  
```
which gengetopt
which cmake
which subversion
which xsltproc
which wget 
```

If one or more of these commands does not show a path including a file then you have install the corresponding program. 
Ubuntu users can use one of these commands:  
```
sudo apt-get install gengetopt
sudo apt-get install cmake
sudo apt-get install subversion
sudo apt-get install xsltproc
sudo apt-get install wget 
```

Get information about installation path for these programs  
```
cmake=`which cmake`; # Have cmake in PATH or specify path to binary here
gengetopt=`which gengetopt`; # Have gengetopt in path
gengetopt_install_dir=`dirname $gengetopt`; # or specify its location directly
```

Define an alias for a temporary folder into which files are to be downloaded, e.g.:
```
tmpdir=/home/me/tmp/
[[ -d $tmpdir ]] || mkdir -pv $tmpdir;
```

Define an alias for the folder into which all software will be installed, e.g. 
`prefix_path=/home/me/software` 

Define an alias for the folder that will hold the static sequence databases, e.g.:
`prefix_db_path=/home/me/db`

Define aliases for subfolders for the OCTOPUS software:
```
modhmm_install_dir=$prefix_path/modhmm/;
topology_predictors_install_dir=$prefix_path/topology_predictors;
```

#### Step 2: Install BLAST

Download BLAST from [NCBI ftp site](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/release/2.2.17/) to your $tmpdir. We tested v2.2.17:
```
cd $tmpdir 
tar xvzf blast-2.2.17-x64-linux.tar.gz 
mv blast-2.2.17 $prefix_path
rm  blast-2.2.17-x64-linux.tar.gz 
```

#### Step 3: Install PSIPRED 3.2 
Download psipred32.tar.gz from [UCL server](http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/old_versions/) to your $tmpdir:

```
psipred_install_dir=$prefix_path/psipred3.2/;
[[ -d $psipred_install_dir ]] || mkdir -pv $psipred_install_dir;
cd $tmpdir 
mv psipred32.tar.gz $psipred_install_dir
cd $psipred_install_dir
tar xvzf psipred32.tar.gz
cd src
make
make install 
```

Download databases for PSIPRED. Currently these instructions look for files on the AlignMe server. If this doesn't work, you may need to obtain your own versions of the databases from NCBI.
```
[[ -d $prefix_db_path/psipred   ]] || mkdir -pv $prefix_db_path/psipred;
wget  -O $prefix_db_path/psipred/uniref90.fasta "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta"
wget  -O $prefix_db_path/psipred/uniref90.fasta.pal "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.pal"
wget  -O $prefix_db_path/psipred/uniref90.fasta.00.phr "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.00.phr"
wget  -O $prefix_db_path/psipred/uniref90.fasta.00.pin "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.00.pin"
wget  -O $prefix_db_path/psipred/uniref90.fasta.00.psq "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.00.psq"
wget  -O $prefix_db_path/psipred/uniref90.fasta.01.phr "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.01.phr"
wget  -O $prefix_db_path/psipred/uniref90.fasta.01.pin "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.01.pin"
wget  -O $prefix_db_path/psipred/uniref90.fasta.01.psq "http://www.bioinfo.mpg.de/AlignMe/db/uniref90.fasta.01.psq"
```

Modify the file `runspipred` which is stored in `$psipred_install_dir` to define paths:
- change: set dbname = uniref90filt  to set dbname =  $prefix_db_path/psipred/uniref90.fasta  (example: dbname = /home/me/db/psipred/uniref90.fasta)
- change: set datadir = ./data to  $prefix_path/blast-2.2.17  
- change: set execdir = ./bin to $psipred_install_dir/bin
- change: set datadir = ./data to $psipred_install_dir/data
- change: $ncbidir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i $tmproot.fasta -C $tmproot.chk >& $tmproot.blast   to  $ncbidir/blastpgp -b 0 -j 3 -h 0.001 -d $dbname -i $tmproot.fasta -C $tmproot.chk -Q $basename.pssm >& $tmproot.blast


#### Step 4: Install OCTOPUS 
This install has two components: the HMM and the topology predictors.

##### Install MODHMM 
Download source code:
```
[[ -d $tmpdir/modhmm_src ]] || mkdir -pv $tmpdir/modhmm_src;
svn co https://svn.sbc.su.se/repos/cbr/modhmm-projects/trunk/modhmm $tmpdir/modhmm_src;
```

Build and install:
```
[[ -d $tmpdir/modhmm_build ]] || mkdir -pv $tmpdir/modhmm_build;
cd $tmpdir/modhmm_build;
[[ -d $modhmm_install_dir ]] || mkdir -pv $modhmm_install_dir;
cmake -D CMAKE_INSTALL_PREFIX=$modhmm_install_dir -D CMAKE_PREFIX_PATH=$gengetopt_install_dir $tmpdir/modhmm_src;
make
make install
```

##### Install Topology Predictors
Download source code:
```
[[ -d $tmpdir/topology_predictors_src ]] || mkdir -pv $tmpdir/topology_predictors_src;
svn co https://svn.sbc.su.se/repos/cbr/modhmm-projects/trunk/cmdline $tmpdir/topology_predictors_src;
```

Build and install:
```
[[ -d $tmpdir/topology_predictors_build ]] || mkdir -pv $tmpdir/topology_predictors_build;
cd $tmpdir/topology_predictors_build;
[[ -d $topology_predictors_install_dir ]] || mkdir -pv $topology_predictors_install_dir;
cmake -D CMAKE_INSTALL_PREFIX=$topology_predictors_install_dir -D TARGETS="topcons;spoctopus;prodiv_tmhmm;scampi;scampi-msa" -D CMAKE_PREFIX_PATH=$modhmm_install_dir $tmpdir/topology_predictors_src
make
make install
```

Clean up tmp folder:
```
cd $tmpdir 
rm * 
```

##### Get sequence databases for OCTOPUS
NB these are currently pointing to the AlignMe server; see instructions above if they are not found.
```
[[ -d $prefix_db_path/octopus  ]] || mkdir -pv $prefix_db_path/octopus;
wget  -O $prefix_db_path/octopus/uniref90.mem.fasta.phr "http://www.bioinfo.mpg.de/AlignMe/db/octopus/uniref90.mem.fasta.phr"
wget  -O $prefix_db_path/octopus/uniref90.mem.fasta.pin "http://www.bioinfo.mpg.de/AlignMe/db/octopus/uniref90.mem.fasta.pin"
wget  -O $prefix_db_path/octopus/uniref90.mem.fasta.psq "http://www.bioinfo.mpg.de/AlignMe/db/octopus/uniref90.mem.fasta.psq"

```

##### Download OCTOPUS fix
Download source code containing modifications that were need to obtain the same OCTOPUS predictions using the local version as those that can be obtained using their web server:
```
wget  -O $tmpdir/Octopus_for_AlignMe.tar.gz "http://www.bioinfo.mpg.de/AlignMe/download/Octopus_for_AlignMe.tar.gz"
```

Extract `Octopus_for_AlignMe.tar.gz` to `$topology_predictors_install_dir/spoctopus/`

Adjust folders in these files so that they fit to your local configuration:  
    - In `$topology_predictors_install_dir/spoctopus/BLOCTOPUS_modified.sh`  
        - Change: octopusdir=/home/me/software/topology_predictors/spoctopus  
        - Change: workingdir=`/bin/mktemp -d /home/me/software/topology_predictors/BLOCTOPUS_XXXXXXXXXX` || exit 1  
    - In `$topology_predictors_install_dir/spoctopus/modhmmblast_modified/run_psiblast.sh`  
        - Change: blastfolder=/home/me/software/blast-2.2.17/bin/  
        - Change: modhmmblast=/home/me/software/topology_predictors/spoctopus/modhmmblast_modified  


#### Step 5: Testing - Generate Inputs for AlignMePST 
Two fasta files are required as input for AlignMe. They have to be formatted like this:
```
>1KPL_A
TPLAILFMAAVVGTLTGLVGVAFEKAVSWVQNMRIGALVQVADHAFLLWPLAFILSALLAMVGYFLVRKFAPEAGGSGI
PEIEGALEELRPVRWWRVLPVKFIGGMGTLGAGMVLGREGPTVQIGGNLGRMVLDVFRMRSAEARHTLLATGAAAGLSA
```

Define shortcuts for directory, filenames and sequence identifiers:
```
fasta_dir=/home/me/fastas/
fasta1=1KPL_A.fa
fasta2=1OTS_B.fa 
fasta1_id=1KPL_A
fasta2_id=1OTS_B
```

Move to the correct directory:
```
AlignMe_input_dir=/home/me/AlignMe_inputs/
[[ -d $AlignMe_input_dir ]] || mkdir -pv $AlignMe_input_dir;
cd $AlignMe_input_dir
```

Make a PSIPRED prediction for each of these fasta files:
```
scp $fasta_dir/$fasta1 $AlignMe_input_dir/
scp $fasta_dir/$fasta2 $AlignMe_input_dir/ 
$psipred_install_dir/runpsipred $AlignMe_input_dir/$fasta1
$psipred_install_dir/runpsipred $AlignMe_input_dir/$fasta2
```

Make an OCTOPUS prediction:
You need to generate a file containing the names of your sequences: 
```
echo $fasta1_id > $AlignMe_input_dir/protnamefile.txt 
echo $fasta2_id >> $AlignMe_input_dir/protnamefile.txt 
```

Execute OCTOPUS:
```
echo $topology_predictors_install_dir/spoctopus/BLOCTOPUS_modified.sh  \
    $AlignMe_input_dir/protnamefile.txt  \
    $fasta_dir $AlignMe_input_dir \
    $prefix_path/blast-2.2.17/bin/blastall \
    $prefix_path/blast-2.2.17/bin/blastpgp \
    $prefix_db_path/octopus/uniref90.mem.fasta $prefix_path/blast-2.2.17/bin/makemat -N 
```

Now you should have the following inputs for AlignMe: 
- Fasta files: $AlignMe_input_dir/1KPL_A.fa and $AlignMe_input_dir/1OTS_B.fa 
- PSIPRED predictions: $AlignMe_input_dir/1KPL_A.ss2 and $AlignMe_input_dir/1OTS_B.ss2
- OCTOPUS predictions: $AlignMe_input_dir/NN_PRF_FILES/1KPL_A.prf and $AlignMe_input_dir/1OTS_B.prf
- PSSM files: $AlignMe_input_dir/1KPL_A.pssm and $AlignMe_input_dir/1OTS_B.pssm

The .prf files have to be modified manually so that they contain only the prediction. 
The prediction starts with " ALPHA:   M       L       G       I       -      <SPACE> <LABEL> <QUERY>"
and ends 1 line above "END1". Everything that is not(!) within that segment has to be deleted (incl. "END1")

Finally, define the output filenames and run AlignMe:
```
alignment_file=/home/me/aln.out
profile_file=/home/me/prof.out 

perl $AlignMe_install_dir/use_best_parameters.pl \
     -alignme_exe $AlignMe_install_dir/alignme \
     -fasta1 $fasta_dir/$fasta1 \
     -fasta2 $fasta_dir/$fasta2 \
     -sspred1 $AlignMe_input_dir/$fasta1_id.ss2 \
     -sspred2 $AlignMe_input_dir/$fasta2_id.ss2 \
     -tmpred1 $AlignMe_input_dir/NN_PRF_FILES/$fasta1_id.prf \
     -tmpred2 $AlignMe_input_dir/NN_PRF_FILES/$fasta2_id.prf \
     -pssm1 $AlignMe_input_dir/$fasta1_id.pssm \
     -pssm2 $AlignMe_input_dir/$fasta2_id.pssm \
     -output_alignment $alignment_file -output_profile $profile_file 
```


