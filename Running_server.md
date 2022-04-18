# Running AlignMe on the Web Server

Jump to section on:
- [Getting started](#getting-started)
- [Pairwise Sequence Alignment](#pairwise-sequence-alignment) 
   - [Alignment Inputs and Parameters](#alignment-inputs-and-parameters)
   - [Anchors](#anchors)
   - [Batch Mode](#batch-mode)
- [Visualizing alignments in structure space](#visualizing-alignments-in-structure-space)
- [Align Multiple Sequence Alignments (Hydropathy Profile Alignments)](#align-multiple-sequence-alignments)   
- Examples of running on the AlignMe Web Server
   - [Pairwise](#pairwise)
   - [Batch](#batch-mode)
   - [Aligning two MSAs to match two family-averaged hydropathy profiles](#aligning-two-msas-to-match-two-family-averaged-hydropathy-profiles)

   



### Getting started 

To access the AlignMe Web Server go to: [https://www.bioinfo.mpg.de/AlignMe/](https://www.bioinfo.mpg.de/AlignMe/)

There are two major versions of AlignMe available for convenience via the two tabs on the web server. The two versions, described in full in separate sections below, are:

1. Pair-wise alignments using [two primary sequences as input](#pairwise-sequence-alignment) 
2. Pair-wise alignment of [two families of sequences as input](#align-multiple-sequence-alignments), used to create hydropathy profile alignments.


---

### Pairwise Sequence Alignment

The pair-wise alignments [section of the AlignMe web server](https://www.bioinfo.mpg.de/AlignMe/AlignMePW.html) provides the most common type of alignment, in which two primary sequences are aligned using one of a range of different descriptors, depending on the problem at hand. 


#### Alignment Inputs and Parameters
There are three required inputs for every pairwise AlignMe calculation:
1. the first sequence (can be uploaded or pasted into the corresponding textbox) 
2. the second sequence (can be uploaded or pasted into the corresponding textbox) 
3. selection of an optimized parameter mode (or user-defined parameters)  

The input sequences should be provided in [fasta format](https://lucy-forrest-lab.github.io/AlignMe/Formats.html#Fasta).

The three above inputs can be combined to run AlignMe in a number of different modes. There are four modes with optimized parameters (including gap opening and extension parameters) that the user can choose:
1. Fast: based on an amino acid substitution matrix (VTML) and hydrophobicity scale (HWvH)
2. P Mode: based on a Position Specific Substitution Matrix (PSSM) which is generated during the AlignMe run 
3. PS Mode: based on a PSSM and secondary structure prediction which are both generated during the AlignMe run
4. PST Mode: based on a PSSM, secondary structure prediction, and a transmembrane prediction all of which are generated during the AlignMe run.  
 
As the input mode complexity increases, the speed of the calculation decreases, but the accuracy on more distantly-related sequences increases. We therefore recommend that the user first carry out a "Fast mode" run, to identify the approximate sequence identity between the two proteins ([Stamm et al 2014](https://doi.org/10.1093/nar/gku291)). Based on the sequence identity, the user may chose to rerun their alignment using a more complex, slower mode. Specifically, P mode is accurate for alignments >40% identity, PS mode gives more accurate aligmnents for proteins in the range 15-40% identity, and PST can be a good choice when the identity is <15% ([Stamm et al 2013](https://doi.org/10.1371/journal.pone.0057731)). 

In addition to these four pre-optimized modes, a user can also define their own parameters. This includes selecting gap penalties and choosing a PSSM, secondary structure prediction, and/or a transmembrane prediction to be generated during the AlignMe run. The user can also upload their own PSSM, secondary structure predictions, and/or transmembrane predictions. 

#### Anchors 

In some cases, initial alignments of  two proteins fail to match two positions that the expert user knows (from, e.g. biochemical or genetic data) should be related. AlignMe has the option to include position restraints, known as anchors, in the sequence alignment to remedy these errors. The user needs to specify which amino acids are to be matched in each sequence as well as an anchor constraint strength factor, known as an anchor weight.

There are two ways that an anchor can be input on the AlignMe Web Server:
1. Input into the form on the the server
2. Upload in a file

If the user chooses to upload a file containing the anchor positions, each anchor has to be on a separate line in the anchor file, following the scheme described in the [Formats section](#Formats.md):  

`position_in_first_sequence   position_in_second_sequence   strength`

e.g.: 
`25  36 1000`

The above example would align residue 25 in the first sequence with residue 36 of the second sequence with a relative strength of 1000. More information about anchors can be found in the corresponding article ([Staritzbichler et al 2020](PloS One 16(4): e0239881)).


#### Batch Mode

The AlignMe Web Server also allows the user to run many pairwise alignments at once, i.e., in batches. In order to run batches of alignments, instead of pasting individual sequences into each input textbox, the user can paste all of the sequences to be aligned into the textbox, or upload a file with multiple sequences. Sequences should be provied in [fasta format](https://lucy-forrest-lab.github.io/AlignMe/Formats.html#Fasta). When sets of sequences are provided, every sequence in the first input will be aligned with every sequence in the second input, up to a maximum of 1000 total alignments per submission. Every alignment will be pairwise and will use the same set of parameters. The results page will provide the alignment and plots (based on which mode is selected) for the first set of sequences that were aligned. The rest of the results can be downloaded via links provided on the results page. 

Batch mode also allows the input of anchors. To use anchors in Batch Mode, the input boxes cannot be used, and the user must upload a file containing the information about the anchors in the following format: 

`position of sequence in input file: position_in_first_sequence   position of sequence in input file: position_in_second_sequence   strength`

e.g.:
`1:18 7:48 1000`

This example would align residue 18 in the first sequence with residue 48 in the seventh sequence with a relative strength of 1000.

It is important to note that when running Batch Mode on the AlignMe Web Server the user must use the tools available on the Web Server when defining the parameters. Thus, the hydrophobicity scale, PSSM, secondary structure prediction, or the transmembrane prediction must be generated during the AlignMe run and cannot be uploaded by the user.


---

### Align Multiple Sequence Alignments

The second tab on the AlignMe Web Server allows the user to [align two hydropathy profiles based on multiple sequence alignments (MSAs)](https://www.bioinfo.mpg.de/AlignMe/AlignMe_MSA.html). This functionality then allows for the estimatmation of the structural similarity between two protein families rather than to align two multiple sequence alignments. 

There are three required inputs for every MSA-MSA AlignMe calculation:

   1. the first multiple sequence alignment (can be uploaded or pasted into the specified box)
   2. the second second multiple sequence alignment (can be uploaded or pasted in to the specified box)
   3. selection of an optimized parameter mode (or user-defined parameters)

In contrast to the pairwise sequence alignments, the MSA-MSA alignments only have the option to run in the parameter optimized Fast mode which is based on substitution matrix & hydrophobicity scale. In addition to Fast mode, the user can also choose/upload their own substitution matrix, hydrophobicity scale, and gap penalties. 

---

### Examples of running on the AlignMe Web Server
Below are some example inputs and outputs for runs on the AlignMe Web Server. 

#### Pairwise

For example, if we align two sequences pairwise in Fast Mode these could be our inputs:

Two sequences, in FASTA format:

```
>mELN
MGQMVPPRSIQNEDFWKNPWDVGGLTVIGLFTSTFLLFVLFAVVFGYVEKAVFEEE

>dsCLB
MNEAKSLFTTFLILAFLLFLLYAFYEAAF
```

The anchor input:
```
29	7	1000
32	10	1000
36	14	1000
40	18	1000
```

And our result alignment is:

```
CLUSTAL W formatted alignment obtained with AlignMe 1.2.2



mELN_________   MGQMVPPRSIQNEDFWKNPWDVGGLTVIGLFTSTFLLFVLFAVVFGYVEKAVFEEE
dsCLB________   MNEAKS----------------------LFTTFLILAFLLF-LLYAFYEAAF----

                *                           a  a   a * a*       * *     


Sequence identity: 8.93% 
Matched Positions: 51.79%
```

where the * (asterisk) indicates identical residues and the letter 'a' indicates an anchored position.


#### Batch Mode

If the user wants to run 9 pairwise alignments in the same submission the following example input can be used:

First file input:
```
>mALN
MEVSQAASGTDGVRERRGSFEAGRRNQDEAPQSGMNGLPKHSYWLDLWLFILFDLALFVFVYLLP
>mELN
MGQMVPPRSIQNEDFWKNPWDVGGLTVIGLFTSTFLLFVLFAVVFGYVEKAVFEEE
>hPLN
MEKVQYLTRSAIRRASTIEMPQQARQKLQNLFINFCLILICLLLICIIVMLL
```
Second file input:
```
>mELN
MGQMVPPRSIQNEDFWKNPWDVGGLTVIGLFTSTFLLFVLFAVVFGYVEKAVFEEE
>mALN
MEVSQAASGTDGVRERRGSFEAGRRNQDEAPQSGMNGLPKHSYWLDLWLFILFDLALFVFVYLLP
>hPLN
MEKVQYLTRSAIRRASTIEMPQQARQKLQNLFINFCLILICLLLICIIVMLL
```
Each sequence in the first input file will be aligned with every sequence in the second input file. Additionally, if we wanted to include anchors in the batch alignments, the following file can be uploaded (note: batch mode anchors cannot be input in the form on the server and must be input via file upload:

```
1:49 1:37 1000
1:49 2:49 1000
1:49 3:31 1000
```

where the 49th residue of the first sequence in the first input file would be aligned with the following residues in the second input file: the 37th residue of the first sequence, the 49th residue of the second sequence, and the 31st residue of the third sequence. 

So in this case, out of the 9 alignments produced, 3 of them would contain anchored residues. For example, the first alignment would look like:

```
CLUSTAL W formatted alignment obtained with AlignMe 1.2.2

mALN_________   MEVSQAASGTDGVRERRGSFEAGRRNQDEAPQSGMNGLPKHSYWLDLWLFILFDLALFVFVYLLP------
mELN_________   MGQMVPPRSIQNEDFWKNPWDVG------------GLTVIGLFTSTFLLFVLF---AVVFGYVEKAVFEEE

                *                     *                        a** **     ** *


Sequence identity: 12.68% 
Matched Positions: 70.42%
```

where we can see, indicated by the 'a', that the 49th residue in the mALN sequence is aligned to the 37th residue in the mELN sequence

---


### Visualizing alignments in structure space

Sequence alignments are very insightful for understanding how positions in two sequences are related, but sequence alignments only provide 2D information about these positions - this is where structural protein data becomes important as it allows for the 3D visualization residue positions. In order to understand the AlignMe pairwise sequence alignments in a structural context, the add-on MutationExplorer can be utilized. MutationExplorer allows the user to visualize the AlignMe calculation result mapped onto the PDB structure of one or more of the sequences being compared. 

The requirements for using MutationExplorer:

   1. At least one PDB structure from the sequences being aligned (use the coordinate file)
   2. The PDB sequence in the coordinate file must be identical to the aligned sequence

If one structure is provided, the pairwise sequence alignment will be mapped onto the provided structure with a color scheme corresponding to the alignment. For example:

- red: conserved residue
- yellow: similar residue (within one of the five amino acid groups)
- light blue: dissimilar residue
- dark blue: not aligned due to a gap or in a different chain 

If two structures are provided, the two structures will be superimposed according the the pairwise sequence alignment.

MutationExplorer can be accessed on the AlignMe results page which will then redirect the user into a web-based “NGL” molecular viewer on the MutationExplorer site hosted at http://proteinformatics.org/. 

---

#### Aligning two MSAs to match two family-averaged hydropathy profiles

If the user wishes to align two multiple sequence alignments they can use the following input:

First input file:
```
>ALP1_Saccharomyces_cerevisiae
MDETVNIQMSKEGQYEINSSSIIKEEEFVDEQYSGENVTKAITTERKVEDDAAKETESSPQERREVKRKLKQRHIGMIALGGTIGTGLIIGIGPPLAHAGPVGALISYLFMGTVIYSVTQSLGEMVTFIPVTSSFSVFAQRFLSP-ALGATNGYMYWLSWCFTFALELSVLGKVIQYWTEAVPLAA-------WIVIFWCLLTSMNMFPVKYYGEFEFCIASIKVIALLGFIIFS-FCVVCGA--GQSDGPIGFRYWRNPGAWGPGIISSDKNEGRFLGWVSSLINAAF-TYQGTELVGITAGEAANPRKALPRAIKKVVVRILVFYILSLFFIGLLVPYNDPKLDSDGIFVSSSPFMISIENSGTKVLPDIFNAVVLITILSAGNSNVYIGSRVLYSLSKNSLAPRFLSNVTRGGVPYF----SVLSTSVFG-FLAFLEVSAGSGKAFNWLLNITGVAGFFAWLLISFSHIRFMQAIRKRGISRDDLPYKAQMMPFLAYYASFFIALIVLIQGFTAFAPTFQPIDFVAAYISIFLFLAIWLSFQVWFKCRLLWKLQDIDIDSDRRQIEELVWIEPECKTRWQRVWDVLS
>PROY_Salmonella_typhimurium
-----------------------------------------------------------MESNNKLKRGLSTRHIRFMALGSAIGTGLFYGSADAIKMAGPS-VLLAYIIGGVAAYIIMRALGEMSVHNPAASSFSRYAQENLGP-LAGYITGWTYCFEILIVAIADVTAFGIYMGVWFPAVPHWI-------WVLSVVLIICAINLMSVKVFGELEFWFSFFKVATIIIMIVAGIGIIVWGI--GNGGQPTGIHNLWSNGGFFS---------NGWLGMIMSLQMVMF-AYGGIEIIGITAGEAKDPEKSIPRAINSVPMRILVFYVGTLFVIMSIYPWNQ-------VGTNGSPFVLTFQHMGITFAASILNFVVLTASLSAINSDVFGVGRMLHGMAEQGSAPKVFAKTSRRGIPWV----TVLVM-TIALLFAVYLNYIMPENVFLVIASLATFATVWVWIMILLSQIAFRRRLPPE--EVKALKFKV---P--GGVVTTIAGLIFLVFIIALIG--YHPDTRISLYVGFAWIV---LLLIGWIFKR----RRDRQLAQA--------------------------
>MMUP_Escherichia_coli
--------------------------------------------------------MQTTQQNAPLKRTMKTRHLIMLSLGGVIGTGLFFNTGYIISTTGAAGTLLAYLIGALVVWLVMQCLGELSVAMPETGAFHVYAARYLGP-ATGYTVAWLYWLTWTVALGSSFTAAGFCMQYWFPQVPVWV-------WCVVFCAIIFGLNVISTRFFAEGEFWFSLVKVVTIIAFIILG-GAAIFGFIPMQDGSPAPGLSNITAEGWFP---------HGGLPILMTMVAVNF-AFSGTELIGIAAGETENPRKVIPVAIRTTIARLIIFFIGTVFVLAALIPMQQ-------VGVEKSPFVLVFEKVGIPYAADIFNFVILTAILSAANSGLYASGRMLWSLSNERTLPACFARVTKNGVPLT----ALSVS-MLGGVLALFSSVVAPDTVFVALSAISGFAVVAVWLSICASHFVFRRRHLQQGKALSELHYRA---P--WYPLVPVLGFVLCLVACVGLA--FDPAQRIALWCGLPFVA---LCYGAYFLTQ----PRNAKQEPEHVAE----------------------
```
Second input file:
```
>LeuT_Aquifex_aeolicus
---------------------------------------------------------------------------MEVKREHWATRLGLILAMAGNAVGLGNFLRFPVQAAENGGGAFMIPYIIAFLLVGIPLMWIEWAMGRYGGAQGHGTTPAIFYLLWRNRFAKIL---GVFGLWIPLVVAIYYVYIESWTLGFAIKFLVGLVPEP-------PPNATDPD----------------SILRPFKEFLYSYIGVPKGDEPILK-PSLFAYIVFLITMFINVSILIRGISKGIERFAKIA-MPTLFILAVFLVIRVFLLETPNGTAADGLNFLWTPDFEKLKDPGVWIAAVGQIFFTLSLGFGAIITYASYVRKDQDIVLSGLTAATLNEKAEVILGGSISIPAAVAFF----GVANAVAIAKAGAFNLGFITLPAIFSQTAG-GTFLGFLWFFLLFFAGLTSSIAIMQPMIAFLEDEL-KLSRK--HAVLWTAAIVF----FSAHLVMFLNKS---LDEMDFWAG-TIGVVFF-GLTELIIFFWIFGADKAWEEINRGGIIKVPRIYYYVMRYITPAFLAVLLVVWAREYIPKIMEE----THWTVWITRFYIIGLFLFLTFLVFLA---------------------------------------------------ERRRNHESA----------------------
>AAT1_Aedes_aegypti
MPEIATISYPESKKNDEANSSHGNGNGVVQLNASQPENAAQNRPEWLELAESSNFLCHVFQCPLVNQLNASQPENAAQNRPEWSNKLEFLMSCISMSVGLGNVWRFPFTAYENGGGAFLIPYLIVLFVIGRPLYYLEMALGQFTSRSSVK--------IW--EVSPLFKGIGIGQLVGTTSVVSYYVSLIALTLHYIFASFASELPWATCKDNW-ADNCVDSSLVVNEMRDSNSTSGQ-KVSSSQIYFLD---IVLKEKDSIDDGIGAPDW-KLTLWLLLAWVVIFLVLVRGVKSSGKVAYFLAIFPYVVLIIILVRALTLEG--AVDGVIFFIKPQWGELLNPKVWYSATTQLFSSLSVGMGSIIMFSSYNNFHHNIYRDAMIVTTLDTFTSLLGG--MTIFSILGNLAHNLGIEDISKVVKS-GTGLAFISYPDAIAKFDIVPQLFSVLFFFMLFVLGVGSAVALHGAIITAFWDAFP-------KRKYWHLALILSTIGFCTGLVYITPGGQWILDLVDHYG--GTFLIYVLAIIEMVAIFWIYGLDNWCNDIEFMVQRRVGLYWRLCWGLITPLFMIAVFIYSLVEYKWPTYSG-QQYPGEALICG--IFVFIFGILQVLIWAVWTVTGHSTKESVWGKITKAAKPSSQWGPCTEHTRKAWLKYKEEAKSRRDDIIYAKNHSSIVRKLCVLFGMYDDFIRLDSTRL
>SerT_Homo_sapiens
-METTPLNSQKQLSACEDGEDCQENGVLQKVVPTPGDKVESGQISNGYSAVPSPGAGDDTRHSIPATTTTLVAELHQGERETWGKKVDFLLSVIGYAVDLGNVWRFPYICYQNGGGAFLLPYTIMAIFGGIPLFYMELALGQYHRNGCIS--------IWR-KICPIFKGIGYAICIIAFYIASYYNTIMAWALYYLISSFTDQLPWTSCKNSWNTGNCTNYF------SEDNITWTL-HSTSPAEEFYTRHVLQIHRSKGLQD-LGGISW-QLALCIMLIFTVIYFSIWKGVKTSGKVVWVTATFPYIILSVLLVRGATLPG--AWRGVLFYLKPNWQKLLETGVWIDAAAQIFFSLGPGFGVLLAFASYNKFNNNCYQDALVTSVVNCMTSFVSG--FVIFTVLGYMAEMRNE-DVSEVAKDAGPSLLFITYAEAIANMPA-STFFAIIFFLMLITLGLDSTFAGLEGVITAVLDEFPHVWAK--RRERFVLAVVITC--FFGSLVTLTFGGAYVVKLLEEYAT-GPAVLTV-ALIEAVAVSWFYGITQFCRDVKEMLGFSPGWFWRICWVAISPLFLLFIICSFLMSPPQLRLFQ-YNYPYWSIILG--YCIGTSSFICIPTYIAYRLI--ITPGTFKERIIKSITPE----------------------TPTE-IPCGDIRLNAV---------------------

```
The result will be a hydropathy plot showing the aligned hydropathy values averaged over all sequences in each of the two MSAs: 

<img width="795" alt="msa_msa_hydropathy_profile" src="https://user-images.githubusercontent.com/44010226/163870221-f14dd917-4639-4fbc-9e6e-f44f50858f55.png">

In addition to the hydropathy plot, a pairwise sequence alignment will be provided for the first sequence of each of the two MSAs, to allow the hydropathy profile alignment to be interpreted in more detail. 

---

[Return to home page](index.md)
