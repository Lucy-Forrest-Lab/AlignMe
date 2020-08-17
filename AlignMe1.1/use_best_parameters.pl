#!/usr/bin/perl -w
# This is a script to use the best gap penalty parameters of AlignMe in combination with the best input parameters. 
# Therefore, you need to have 
# - two protein sequences of alpha-helical proteins in fasta format
# - two corresponding secondary structure predictions of PSIPRED
# - two corresponding transmembrane prediction from OCTOPUS
# - two corresponding Position Specific Matrices 
# More information on how to obtain these inputs is available in the manual, which can be downloaded from the AlignMe website
# Alignments with optimal parameters for alpha-helical proteins can also be obtained online via http://www.forrestlab.org/AlignMe/ 


my %inputs = (
        "-alignme_exe" => '',
        "-fasta1" => '',
        "-fasta2" => '',
        "-sspred1" => '',
        "-sspred2" => '',
        "-tmpred1" => '',
        "-tmpred2" => '',
        "-pssm1" => '',
        "-pssm2" => '',
        "-output_alignment" => '',
        "-output_profile" => '',

);

my $best_gap_penalties = "-below_threshold_gap_opening_penalty 2.14 -above_threshold_gap_opening_penalty 2.96 -below_threshold_gap_extension_penalty 3.10 -above_threshold_gap_extension_penalty 3.06 -termini_gap_opening_penalty 0.07 -termini_gap_extension_penalty 1.18 -thresholds_for_penalties 0.5 ";


$j = 0;
$error = 0;
while ($j< @ARGV)
{
	if (exists $inputs{$ARGV[$j]})
	{
		$nr = $ARGV[$j];
		$j++;
		if (substr($ARGV[$j], 0, 1) ne "-")
		{
			$inputs{$nr} = $ARGV[$j];
			$j++;
		}
		else
		{
			$error = 1;
			print "ERROR: Value after the flag ".$inputs{$nr}." is missing \n";
		}
	}
	else
	{
		$error = 1;
		print "ERROR: Flag ".$ARGV[$j]." is not allowed! \n";
		$j++;
	}
}

foreach $key ( sort keys (%inputs))
{	
	if ($inputs{$key} eq '')
	{
		$error = 1;
		print "ERROR: The flag ".$key." was not defined! \n";
	}
	else
	{
		#print "key is ".$key."\t";
		#print  $inputs{$key}."\n";
	}
}

if (!$error)
{
	open(SIMSCORE, '>simscore.txt');
	print SIMSCORE "weight: 1.4 type: UniversalProfileSimilarity column: 5 headerlines: 1 profile1: ".$inputs{"-sspred1"}." profile2: ".$inputs{"-sspred2"}."\n";
	print SIMSCORE "weight: 0.2 type: PositionSpecificSimilarity PSSM1: ".$inputs{"-pssm1"}." PSSM2:  ".$inputs{"-pssm2"}."\n";
	print SIMSCORE "weight: 1.4 type: UniversalProfileSimilarity column: 4 headerlines: 1 profile1: ".$inputs{"-sspred1"}." profile2: ".$inputs{"-sspred2"}."\n";
	print SIMSCORE "weight: 1.4 type: UniversalProfileSimilarity column: 6 headerlines: 1 profile1: ".$inputs{"-sspred1"}." profile2: ".$inputs{"-sspred2"}."\n";
	print SIMSCORE "weight: 4.2 type: UniversalProfileSimilarity column: 3 headerlines: 1 profile1: ".$inputs{"-tmpred1"}." profile2: ".$inputs{"-tmpred2"}."\n";
	close(SIMSCORE);

	$commandline = $inputs{"-alignme_exe"}." -fasta_file1 ".$inputs{"-fasta1"}." -fasta_file2 ".$inputs{"-fasta2"}." -similarity_score_file ./simscore.txt ".$best_gap_penalties." -output_aligned_sequences ".$inputs{"-output_alignment"}." -output_aligned_profiles ".$inputs{"-output_profile"}."\n";

	system($commandline);
	#print $commandline;
}


