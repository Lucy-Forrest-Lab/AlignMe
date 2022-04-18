# Running AlignMe on the Web Server

Jump to section on:
- [Getting started](#getting-started)
- [Pairwise Sequence Alignment](#pairwise-sequence-alignment) 
   - [Alignment Inputs and Parameters](#alignment-inputs-and-parameters)
   - [Anchors](#anchors)
   - [Batch Mode](#batch-mode)
- [Align Multiple Sequence Alignments](#align-multiple-sequence-alignments)   
- Examples of running on the AlignMe Web Server
   - [Pairwise](#pairwise)
   - [Batch](#batch-mode)
   - [MSA](#msa)
   



### Getting started 

To access the AlignMe Web Server go to: 
[https://www.bioinfo.mpg.de/AlignMe/](https://www.bioinfo.mpg.de/AlignMe/]

There are two major versions of AlignMe available for convenience as two separate tabs on the web server:
1. Pair-wise alignments using [two primary sequences as input](#pairwise-sequence-alignment) 
2. Pair-wise alignment of [two families of sequences as input](#align-multiple-sequence-alignments). This option is often used for creating hydropathy profile alignments.

All of the output results can be visualized on the AlignMe server as well as via download to your computer.


---

### Pairwise Sequence Alignment

#### Alignment Inputs and Parameters
There are three required inputs for every pairwise AlignMe calculation:
1. the first sequence (can be uploaded or pasted into the corresponding window) 
2. the second sequence (can be uploaded or pasted into the corresponding window) 
3. selection of an optimized parameter mode (or user-defined parameters) 

These inputs can be combined to run AlignMe in a number of different modes. There are four modes with optimized parameters (including gap opening and extension parameters) that the user can choose:
1. Fast: based on a substituion matrix (VTML) and hydrophobicity scale (HWvH)
2. P Mode: based on a Position Specific Substitution Matrix (PSSM) which is generated during the AlignMe run 
3. PS Mode: based on a PSSM and secondary structure prediction which are both generated during the AlignMe run
4. PST Mode: based on a PSSM, secondary structure prediction, and a transmembrane prediction all of which are generated during the AlignMe run

In addition to these four modes, a user can also choose to define their own parameters. This includes choosing gap penalties and choosing a PSSM, secondary structure prediction, and/or a transmembrane prediction to be generated during the AlignMe run. The user can also upload their own PSSM, secondary structure predictions, and/or transmembrane prediction. 

#### Anchors 

AlignMe has the option to include position restraints, known as anchors, in the sequence alignment. The user needs to specify which amino acids to be matched in each sequence as well as an anchor constraint strength factor, known as an anchor weight.

There are two ways an anchor can be input on the AlignMe Web Server:
1. Input on the server page
2. Uploaded in a file

If the user chooses to upload a file containing the anchor positions, each anchor has to be given as a line in the anchor file, following the scheme described in the [Formats section](#Formats.md):  

`position_in_first_sequence   position_in_second_sequence   strength`

e.g.: 
`25  36 1000`

This would align residue 25 in the first sequence with residue 36 of the second sequence with a relative strength of 1000. 

More detail can be found in the section on [Anchors](#Anchors.md).

#### Batch Mode

The AlignMe Web Server also allows the user to run many pairwise alignments at once, this is called Batch Mode. In order to use Batch Mode, instead of pasting individual sequences into each input box, the user can paste all of the sequences to be aligned into the box, or upload a file with multiple FASTA sequences. For the provided sets of sequences, every sequence in the first input will be aligned with every sequence in the second input, up to a maximum of 1000 total alignments per submission. Every alignment will be pairwise and will use the same set of parameters. The result page will provide the alignment and plots (based on which mode is selected) for the first set of sequences that were aligned. The rest of the results can be downloaded via links on the results page. 

Batch mode also allows the input of anchors. To use anchors in Batch Mode, the user must upload a file (cannot use the input boxes on the server) in the following format: 

`position of sequence in input file: position_in_first_sequence   position of sequence in input file: position_in_second_sequence   strength`

e.g.:
`1:18 7:48 1000`

This would align residue 18 in the first sequence with residue 48 in the seventh sequence with a relative strength of 1000.


---

### Align Multiple Sequence Alignments

The AlignMe Web Server allows the user to align two hydropathy profiles based on multiple sequence alignments (MSAs). This functionality then allows for the estimatmation of the structural similarity between two protein families rather than to align two multiple sequence alignments. 

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

where the asterik indicates identical residues and the 'a' indicate an anchored position.


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
Each sequence in the first input file will be aligned with each sequence in the second input file. Additionally, if we wanted to include anchors in the batch alignments, the following file can be uploaded (note: batch mode anchors cannot be input in the boxes on the server; it must be a file upload):

```
1:49 1:37 1000
1:49 2:49 1000
1:49 3:31 1000
```

where the 49th residue of the first sequence in the first input file would be aligned with the 37th residue of the first sequence, 49th residue of the second sequence, and 31st residue of the third sequence in the second input file. 

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

#### MSA

If the use wishes to align two multiple sequence alignments they can use the following input:

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
The result will be a hydropathy plot showing the aligned hydropathy values of the two MSAs and a pairwise sequence alignment of the first sequence of each of the two MSAs.

