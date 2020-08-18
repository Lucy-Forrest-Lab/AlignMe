/////////////////////////////////////////////////////////////////////////
//!
//!  AlignMe - a program to identify and align similar proteins with or
//!  without sequence similarity. It allows to use a flexible number of
//!  scales or profiles for identifying similarity.
//!
//!  Copyright (C) 2010 by Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//!  AlignMe@rzg.mpg.de
//!
//!  AlignMe is free software; you can redistribute it and/or modify
//!  it under the terms of the GNU General Public License as published by
//!  the Free Software Foundation; either version 2 of the License, or
//!  (at your option) any later version.
//!
//!  AlignMe is distributed in the hope that it will be useful,
//!  but WITHOUT ANY WARRANTY; without even the implied warranty of
//!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//!  GNU General Public License for more details.
//!
//!  You should have received a copy of the GNU General Public License
//!  along with this program. If not, see <http://www.gnu.org/licenses/>.
//!
//!  @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//!  @date: 1.4.2010
/////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <map>
#include <cmath>
#include <set>


#include "../include/definitions.h"
#include "../include/macro_functions_read_write.h"
#include "../include/score_sequence_similarity.h"
#include "../include/score_profile_similarity.h"
#include "../include/amino_acid.h"
#include "../include/sum_function.t.h"
#include "../include/reader.h"
#include "../include/needleman_wunsch.h"
#include "../include/needleman_wunsch_affine_gaps.h"
#include "../include/smith_waterman.h"
#include "../include/command_line_manager.h"
#include "../include/string_functions.h"
#include "../include/alignment_write.h"
#include "../include/std_functions.h"
#include "../include/profile_step_function.h"
#include "../include/constant_function.h"
#include "../include/averager_msa.h"


// reading in the commandline. argc = number of parameters of the commandline, argv = string-array with the commandlineparameters
int main( const int argc, const char * argv[])
{
  ///////////////////////////////////
  ///  DECLARATIONS
  ///////////////////////////////////
  DebugWrite( argv[0]);
  DebugWrite( "Welcome to this sweetly aligned world");



  std::map<std::string,std::string>
	  allowed_flags;

  allowed_flags["fasta_file1"] = "allowed";
  allowed_flags["fasta_file2"] = "allowed";
  allowed_flags["similarity_score_file"] = "allowed";
  allowed_flags["gap_opening_penalty"] = "allowed";
  allowed_flags["gap_extension_penalty"] = "allowed";
  allowed_flags["msa_file1"] = "allowed";
  allowed_flags["msa_file2"] = "allowed";
  allowed_flags["fraction_allowed_gaps"] = "allowed";
  allowed_flags["algorithm"] = "allowed";
  allowed_flags["termini_gap_opening_penalty"] = "allowed";
  allowed_flags["termini_gap_extension_penalty"] = "allowed";
  allowed_flags["output_aligned_sequences"] = "allowed";
  allowed_flags["output_aligned_profiles"] = "allowed";
  allowed_flags["thresholds_for_penalties"] = "allowed";
  allowed_flags["below_threshold_gap_opening_penalty"] = "allowed";
  allowed_flags["below_threshold_gap_extension_penalty"] = "allowed";
  allowed_flags["above_threshold_gap_opening_penalty"] = "allowed";
  allowed_flags["above_threshold_gap_extension_penalty"] = "allowed";
  allowed_flags["profile_gap_value_for_plotting"] = "allowed";
  allowed_flags["help"] = "allowed";
  allowed_flags["-help"] = "allowed";

  // calls the class CommandLineManager in file command_line_manager.h up and reads in the commands that the user gives to the program
  CommandLineManager 
    cmd( argc, argv, allowed_flags);

  size_t
	  line_length = 60;  // for the output of the alignments

  double 
    gap_opening_penalty( 10.0),
    gap_extension_penalty( 1.0),
    termini_gap_opening_penalty,
    termini_gap_extension_penalty,
	below_threshold_gap_opening_penalty,
	above_threshold_gap_opening_penalty,
	below_threshold_gap_extension_penalty,
	above_threshold_gap_extension_penalty,
    last_element,
    fraction_allowed_gaps = 0.5,
    gap_value = 0.0;

  std::string
	first_fasta_id,
	second_fasta_id,
	algorithm = "global_affine",
    similarity_score_file,
    first_fasta_file(""),
    second_fasta_file(""),
    outputfile_aligned_sequences( ""),
    outputfile_aligned_profiles(""),
    first_msa_file(""),
    second_msa_file("");

  //calls definitions.h where AAsequence is per typedef std::vector< GeneralizedAminoAcid>
  // A GeneralizedAminoAcid has as initializations-list m_Type and m_Profiles
  AASequence 
    first_sequence,
    second_sequence,
    first_averaged_alignment,
    second_averaged_alignment;

  // to keep identity of profile for output
  std::vector<std::string>
	write_profiles_header1,
	write_profiles_header2;

  // MSA
  std::vector< AASequence>
    first_msa,
    second_msa;

  // for averaged MSA
  std::vector< int>
	  first_ignored_positions,
	  second_ignored_positions;

  // creates structure for the scores which will be a pair of AA combined with a certain value and the number
  // of the corresponding shared pointer...
  boost::shared_ptr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > 
    scores;

  // creates structure for alignment with a number combined with two other numbers... ???
  std::pair< double, std::list< std::pair< int, int> > >
    alignment;

  // collect amino acid types provided by user
  std::set<char>
	  defined_amino_acids;



//////////////////////////////////////////////////////////////////////
//////  LOGICAL FLAG ASSERTION   /////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#ifdef SECURE

  if( (cmd.IsFlagSet( "fasta_file1") && !cmd.IsFlagSet("fasta_file2")) || (!cmd.IsFlagSet( "fasta_file1") && cmd.IsFlagSet("fasta_file2")))
  {
	  std::cerr << "ERROR: If you provide the flags fasta_file1 or fasta_file2, you always have to provide both flags!" << "\n";
	  std::cerr << "ERROR: If you want to align profiles only, you don't need fasta files" << "\n";
	  return -1;
  }

  if( (cmd.IsFlagSet( "msa_file1") && !cmd.IsFlagSet("msa_file2")) || (!cmd.IsFlagSet( "msa_file1") && cmd.IsFlagSet("msa_file2")))
  {
	  std::cerr << "ERROR: It only makes sense to provide two multiple sequence alignment files, not only one!" << "\n";
	  return -1;
  }

  if
  (
		  ( cmd.IsFlagSet( "fasta_file1") && cmd.IsFlagSet( "msa_file1"))
		  || ( cmd.IsFlagSet( "fasta_file2") && cmd.IsFlagSet( "msa_file2"))
  )
  {
	  std::cerr << "ERROR: Please provide only MSA-files or fasta-files. Don't use both at the same time!" << "\n";
	  return -1;
  }

#endif




  ///////////////////////////////////////////
  /// READ VALUES FROM COMMANDLINE OPTIONS
  ///////////////////////////////////////////

 // ConvertStringToNumericalValue opens string_functions.h and does something... ???

  if( cmd.IsFlagSet( "help") || cmd.IsFlagSet( "-help"))
  {
	  std::cout << "\n";
	  std::cout << "This is the usage information for AlignMe Version 1.0\n";
	  std::cout << "\n";
	  std::cout << "Flags for input files:\n";
	  std::cout << "-similarity_score_file [file]         file containing information about the type \n";
	  std::cout << "                                      of alignment you want to do. See the\n";
	  std::cout << "                                      manual for detailed information!\n";
	  std::cout << "-fasta_file1 [file]                   file containing an amino acid sequence\n";
	  std::cout << "                                      in fasta format \n";
	  std::cout << "-fasta_file2 [file]                   file containing an amino acid sequence\n";
	  std::cout << "                                      in fasta format \n";
	  std::cout << "-msa_file1 [file]                     file containing a multiple sequence\n";
	  std::cout << "                                      in which all sequence are the\n";
      std::cout << "                                      same length\n";
	  std::cout << "-msa_file2  [file]                    file containing a multiple sequence\n";
	  std::cout << "                                      in which all sequence are the\n";
      std::cout << "                                      same length\n";
	  std::cout << "\n";
	  std::cout << "-similarity_score_file must always be provided!\n";
	  std::cout << "-fasta_file1 and -fasta_file2 must be provided together \n";
	  std::cout << " The same is true for -msa_file1 and msa_file2.\n";
	  std::cout << "\n";
	  std::cout << "Flags for output files:\n";
	  std::cout << "-output_aligned_sequences [file]      file to which the aligned amino acid\n";
	  std::cout << "                                      sequences are printed\n";
	  std::cout << "-output_aligned_profiles [file]       file to which the aligned profile values\n";
	  std::cout << "                                      are printed \n";
	  std::cout << "\n";
	  std::cout << "You can choose one of the following gap penalty sets\n";
	  std::cout << "A combination is not possible!\n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 2 gap penalties:\n";
	  std::cout << "-gap_opening_penalty [value]          penalty for opening gaps\n";
	  std::cout << "-gap_extension_penalty [value]        penalty for extending gaps \n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 4 gap penalties:\n";
	  std::cout << "-gap_opening_penalty [value]          penalty for opening gaps within\n";
	  std::cout << "                                      the sequence \n";
	  std::cout << "-gap_extension_penalty [value]        penalty for extending gaps within\n";
	  std::cout << "                                      the sequence \n";
	  std::cout << "-termini_gap_opening_penalty          additional penalty for opening gaps\n";
	  std::cout << " [value]                              at the end of the sequence\n";
	  std::cout << "-termini_gap_extension_penalty        additional penalty for extending gaps\n";
	  std::cout << " [value]                              at the end of the sequence\n";
	  std::cout << "\n";
	  std::cout << "Flags for a set of 6 gap penalties:	\n";
	  std::cout << "-thresholds_for_penalties [values]     defines a threshold according to which\n";
	  std::cout << "                                       gap penalties are assigned to a given\n";
	  std::cout << "                                       position of the sequence\n";
	  std::cout << "-below_threshold_gap_opening_penalty   penalty for opening gaps opposite to \n";
	  std::cout << " [value]                               residues with a value below the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-below_threshold_gap_extension_penalty penalty for extending gaps opposite to \n";
	  std::cout << " [value]                               residues with a value below the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-above_threshold_gap_opening_penalty   penalty for opening gaps opposite to \n";
	  std::cout << " [value]                               residues with a value above the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-above_threshold_gap_extension_penalty penalty for extending gaps opposite to \n";
	  std::cout << " [value]                               residues with a value above the\n";
	  std::cout << "                                       chosen threshold\n";
	  std::cout << "-termini_gap_opening_penalty           additional penalty for opening gaps\n";
	  std::cout << " [value]                               at the end of the sequence\n";
	  std::cout << "-termini_gap_extension_penalty         additional penalty for extending gaps\n";
	  std::cout << " [value]                               at the end of the sequence\n";
	  std::cout << "\n";
	  std::cout << "Other flags:";
	  std::cout << "\n";
	  std::cout << "-fraction_allowed_gaps [value]         columns of the MSA in which more than \n ";
	  std::cout << "                                       [value] positions are gaps are skipped  \n";
	  std::cout << "                                       and not considered in the alignment\n";
	  std::cout << "-profile_gap_value_for_plotting        define a value to be assigned to gaps\n";
	  std::cout << " [value]                               in the output aligned profiles \n";
	  std::cout << "\n";
	  exit (-1);
  }

  DebugWrite( "read values from options...");


  if( cmd.IsFlagSet( "thresholds_for_penalties"))
  {
	  if( cmd.IsFlagSet( "below_threshold_gap_opening_penalty"))
	  {
		  below_threshold_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "below_threshold_gap_opening_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
    	  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"below_threshold_gap_opening_penalty\" \n";
     	  return -1;
	  }
	  #endif
	  if( cmd.IsFlagSet( "above_threshold_gap_opening_penalty"))
	  {
		  above_threshold_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "above_threshold_gap_opening_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
		  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"above_threshold_gap_opening_penalty\" \n";
		  return -1;
	  }
	  #endif
	  if( cmd.IsFlagSet( "below_threshold_gap_extension_penalty"))
	  {
		  below_threshold_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "below_threshold_gap_extension_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
		  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"below_threshold_gap_extension_penalty\" \n";
		  return -1;
	  }
	  #endif
	  if( cmd.IsFlagSet( "above_threshold_gap_extension_penalty"))
	  {
		  above_threshold_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "above_threshold_gap_extension_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
		  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"above_threshold_gap_extension_penalty\" \n";
		  return -1;
	  }
	  #endif
	  if( cmd.IsFlagSet( "termini_gap_opening_penalty"))
	  {
		  termini_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "termini_gap_opening_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
		  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"termini_gap_opening_penalty\" \n";
		  return -1;
	  }
	  #endif
	  if( cmd.IsFlagSet( "termini_gap_extension_penalty"))
	  {
		  termini_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "termini_gap_extension_penalty"));
	  }
	  #ifdef SECURE
	  else
	  {
		  std::cerr << "ERROR: If you want to use thresholds, you have to provide a \"termini_gap_extension_penalty\" \n";
		  return -1;
	  }
	  #endif
  }
  else
  {
	  #ifdef SECURE
	  if( cmd.IsFlagSet( "above_threshold_gap_opening_penalty") || cmd.IsFlagSet( "above_threshold_gap_extension_penalty") || cmd.IsFlagSet( "below_threshold_gap_opening_penalty") || cmd.IsFlagSet( "below_threshold_gap_extension_penalty"))
	  {
		  std::cout << "WARNING: Gap penalties according to a threshold are not used because you did not provide a threshold-value \n";
	  }
      #endif
	  if( cmd.IsFlagSet( "gap_opening_penalty"))
	  {
		  gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "gap_opening_penalty"));
	  }
	#ifdef SECURE
	  else
	  {
		  if (!cmd.IsFlagSet( "below_threshold_gap_opening_penalty") && !cmd.IsFlagSet( "above_threshold_gap_opening_penalty"))
		  {
			  std::cout << "WARNING: No gap opening penalty provided. The default value will be used: " << gap_opening_penalty  <<"\n";
		  }
	  }
	#endif

	  if( cmd.IsFlagSet( "gap_extension_penalty"))
	  {
		  gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "gap_extension_penalty"));
	  }
	#ifdef SECURE
	  else
	  {
		  if (!cmd.IsFlagSet( "below_threshold_gap_extension_penalty") && !cmd.IsFlagSet( "above_threshold_gap_extension_penalty"))
		  {
			  std::cout << "WARNING: No gap extension penalty provided. The default value will be used: " << gap_extension_penalty << "\n";
		  }
	  }
	#endif

	  if( cmd.IsFlagSet( "termini_gap_opening_penalty"))
	  {
		  termini_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "termini_gap_opening_penalty"));
	  }
	  else
	  {
		  termini_gap_opening_penalty = gap_extension_penalty;
	#ifdef SECURE
		  std::cout << "WARNING: No termini gap opening penalty provided. The default value will be used: " << termini_gap_opening_penalty << "\n";
	#endif
	  }

	  if( cmd.IsFlagSet( "termini_gap_extension_penalty"))
	  {
		  termini_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "termini_gap_extension_penalty"));
	  }
	  else
	  {
		  termini_gap_extension_penalty = gap_extension_penalty;
	#ifdef SECURE
		  std::cout << "WARNING: No termini gap extension penalty provided. The default value will be used: " << termini_gap_extension_penalty  << "\n";
	#endif
	 }
  }



  if( cmd.IsFlagSet( "similarity_score_file"))
  {
      similarity_score_file = cmd.GetFirstArgument( "similarity_score_file");
#ifdef SECURE
      if( similarity_score_file.length() == 0)
      {
    	  std::cerr << "ERROR: No file name passed to '-similarity_score_file' flag!" << "\n";
     	  return -1;
      }
#endif
  }
  else
  {
      std::cerr << "\n" << "ERROR: The flag '-similarity_score_file' has to be set!!" << "\n";
      return -1;
  }


  if( cmd.IsFlagSet( "fasta_file1"))
  {
      first_fasta_file = cmd.GetFirstArgument( "fasta_file1");
	  outputfile_aligned_sequences =  "aligned_sequences.aln";
  }
  else if(cmd.IsFlagSet( "msa_file1"))
  {
#ifdef SECURE
	  // consistency check:
	  if( cmd.IsFlagSet( "fasta_file2"))
	  {
		  std::cout << "ERROR: fasta_file2 can be used as flag only in combination with fasta_file1" << "\n";
		  exit(-1);
	  }
	  if( cmd.IsFlagSet( "fasta_file1"))
	  {
		  std::cout << "ERROR: fasta_file1 and msa_file1 can not be used at the same time" << "\n";
		  exit(-1);
	  }
#endif

	  first_msa_file  = cmd.GetFirstArgument( "msa_file1");
	  outputfile_aligned_profiles =  "aligned_profiles.aln";

	  // when averaging an multiple sequence alignment, this fraction determines whether a position in the alignment is ignored
	  if(cmd.IsFlagSet( "fraction_allowed_gaps"))
	  {
		  fraction_allowed_gaps = mystr::ConvertStringToNumericalValue< double> (cmd.GetFirstArgument( "fraction_allowed_gaps"));
	  }
#ifdef SECURE
	  else
	  {
		std::cout << "WARNING: No value for \"fraction of allow gaps\" provided. The default value will be used: " << fraction_allowed_gaps  << "\n";
	  }
#endif
  }

  if( cmd.IsFlagSet( "fasta_file2"))
  {
      second_fasta_file = cmd.GetFirstArgument( "fasta_file2");
  }
  else if (cmd.IsFlagSet( "msa_file2"))
  {
#ifdef SECURE
	  // consistency check:
	  if( cmd.IsFlagSet( "fasta_file1"))
	  {
		  std::cout << "ERROR: fasta_file1 can be used as flag only in combination with fasta_file2" << "\n";
		  exit(-1);
	  }
	  if( cmd.IsFlagSet( "fasta_file2"))
	  {
		  std::cout << "ERROR: fasta_file2 and msa_file2 can not be used at the same time" << "\n";
		  exit(-1);
	  }
#endif

   	  second_msa_file  = cmd.GetFirstArgument( "msa_file2");
  }

  if( cmd.IsFlagSet( "algorithm"))
  {
	  algorithm = cmd.GetFirstArgument( "algorithm");
#ifdef SECURE
	  if( algorithm == "")
	  {
		  std::cerr << "ERROR: If '-algorithm' is set, an algorithm type has to be provided! " << "\n";
		  return -1;
	  }
#endif
  }
  else
  {
	  std::cout << "WARNING: You did not provide an algorithm type. The default algorithm " << algorithm <<" will be used for the alignment" << "\n";
  }

  DebugWrite( "values read");

  if( cmd.IsFlagSet( "output_aligned_sequences"))
  {
	  outputfile_aligned_sequences = cmd.GetFirstArgument( "output_aligned_sequences");
#ifdef SECURE
	  if( outputfile_aligned_sequences == "")
	  {
		  std::cerr << "ERROR: if '-output_aligned_sequences' is set, a filename has to be provided! " << "\n";
 		  return -1;
	  }
#endif
  }
  else
  {
	  if( outputfile_aligned_sequences != "")
	  {
		  std::cout << "WARNING: You did not provide a filename for the output of the sequence alignment. It will be written to " << outputfile_aligned_sequences << "\n";
	  }
  }


  if( cmd.IsFlagSet( "output_aligned_profiles"))
  {
	  outputfile_aligned_profiles = cmd.GetFirstArgument( "output_aligned_profiles");
#ifdef SECURE
	  if( outputfile_aligned_profiles == "")
	  {
		  std::cerr << "ERROR: if '-output_aligned_profiles' is set, a filename has to be provided! " << "\n";
 		  return -1;
	  }
#endif
  }
  else
  {
	  if( outputfile_aligned_profiles != "")
	  {
		  std::cout << "WARNING: You did not provide a filename for the output of the aligned profiles. They will be written to " << outputfile_aligned_profiles << "\n";
	  }
  }




  /////////////////////////////////////////////
  ///  READ PROFILES AND SEQUENCES
  /////////////////////////////////////////////


  if( first_fasta_file != "")
  {
	  DebugWrite( "read first sequence:...");
      first_fasta_id = ReadSequence( first_sequence, first_fasta_file, defined_amino_acids);
#ifdef SECURE
      if (first_sequence.size() == 0)
      {
    	  std::cerr << "ERROR: fasta_file1 "<< first_fasta_file <<" only contains a header and misses a sequence \n";
    	  std::cerr << "ERROR: A fasta-file has to contain a header starting with '>' followed by a sequence in a new line \n";
   	  exit(-1);
      }
#endif
      DebugWrite( "sequence is read: " << first_sequence);
  }
  else if( first_msa_file != "")
  {
	  DebugWrite( "read first MSA:...");
	  ReadMultipleSequences( first_msa, first_msa_file, defined_amino_acids);
#ifdef SECURE
      if (first_msa.size() == 0)
      {
    	  std::cerr << "ERROR: msa_file1 "<< first_msa_file <<" only contains a header and misses a sequence \n";
       	  std::cerr << "ERROR: A msa-file has at least to contain one sequence with a header starting with '>' followed by a sequence in a new line \n";
    	  exit(-1);
      }
#endif
	  DebugWrite( "First MSA is read: " << first_msa);
  }


  if( second_fasta_file != "")
  {
	  DebugWrite( "read second sequence:...");
      second_fasta_id = ReadSequence( second_sequence, second_fasta_file, defined_amino_acids);
#ifdef SECURE
      if (second_sequence.size() == 0)
      {
    	  std::cerr << "ERROR: fasta_file2 "<< second_fasta_file <<" only contains a header and misses a sequence \n";
    	  std::cerr << "ERROR: A fasta-file has to contain a header starting with '>' followed by a sequence in a new line \n";

      }
#endif
      DebugWrite( "sequence is read" << second_sequence);
  }
  else if( second_msa_file != "")
  {
	  DebugWrite( "read second MSA:...");
      ReadMultipleSequences( second_msa, second_msa_file, defined_amino_acids);
#ifdef SECURE
      if (second_msa.size() == 0)
      {
    	  std::cerr << "ERROR: msa_file1 "<< second_msa_file <<" only contains a header and misses a sequence \n";
       	  std::cerr << "ERROR: A msa-file has at least to contain one sequence with a header starting with '>' followed by a sequence in a new line \n";
    	  exit(-1);
      }
#endif
      DebugWrite( "Second MSA is read" << second_msa);
  }

  // opens the file "profile and scales" and gets the values

  if( first_msa_file == "" && second_msa_file == "")
  {
	  DebugWrite( "read profiles and scales...");
	  ReadProfilesAndScales( first_sequence, second_sequence, similarity_score_file, defined_amino_acids, write_profiles_header1, write_profiles_header2, first_fasta_id, second_fasta_id);
	  DebugWrite( "first sequence: " << first_sequence << "\nsecond sequence: " << second_sequence);
  }

  if (first_msa_file != "")
  {
	  std::cout << "WARNING: At the moment, if you provide multiple sequences alignments, it will only work if only scales (values per AA-type) are involved and no profiles (values per AA)" << "\n";
	  write_profiles_header1.push_back( " is based on the msa_file: " + first_msa_file + "\n");

	  for (size_t i = 0; i < first_msa.size( ); i++)
	  {
		  ReadProfilesAndScalesMsa( first_msa[i], similarity_score_file, defined_amino_acids);  // TODO: This will only work if only scales are involved and no profiles !!!
	  }
	  first_sequence = AverageMsa( first_msa, fraction_allowed_gaps, first_ignored_positions);
  }

  if (second_msa_file != "")
  {
	  write_profiles_header2.push_back( " is based on the msa_file: " + second_msa_file+ "\n");
	  for (size_t i = 0; i < second_msa.size( ); i++)
	  {
		  ReadProfilesAndScalesMsa( second_msa[i], similarity_score_file, defined_amino_acids);
	  }

	  second_sequence = AverageMsa( second_msa, fraction_allowed_gaps, second_ignored_positions);
  }

  DebugWrite( "read scores, including substitution matrices");
  ReadScores( scores, similarity_score_file, defined_amino_acids);
  DebugWrite( "scores read: " << scores);


  /////////////////////////////////////////////
  ////  PERFORM ALIGNMENT
  ///////////////////////////////////////////

  Matrix< DynamicProgrammingMatrixElement>
	  matrix( first_sequence.size() + 1, second_sequence.size() + 1);


  std::ofstream write;

  if( algorithm == "global_linear")
  {
	  if( cmd.IsFlagSet( "thresholds_for_penalties"))
	  {
		  std::cout << "Warning: the 'thresholds_for_penalties' flag has no meaning when the algorithm is set to 'global_linear', it only makes sense with 'global_affine'\n" << "\n";
	  }

	  NeedlemanWunsch
		needleman_wunsch
		(
				gap_opening_penalty,
				gap_extension_penalty,
				termini_gap_opening_penalty,
				termini_gap_extension_penalty,
				first_sequence,
				second_sequence,
				scores,
				matrix
		);

	  DebugWrite( "needleman wunsch: calculate matrix:...");
	  needleman_wunsch.CalculateMatrix();

	// needleman_wunsch.Write( write);

	  DebugWrite( "needleman wunsch: traceback...");
	  alignment = needleman_wunsch.TraceBack();
  }
  else if( algorithm == "global_affine")
  {
	  boost::shared_ptr< Function< std::vector< double>, double> >
		  gap_opening_penalty_function,
		  gap_extension_penalty_function;

	  std::vector< double>
		  thresholds;

	  if( cmd.IsFlagSet( "thresholds_for_penalties"))
	  {
		  thresholds =  cmd.GetConvertedArguments< double>( "thresholds_for_penalties");

		  DebugWrite( "thresholds_for_penalties: " << thresholds);

		  gap_opening_penalty_function
			  = boost::shared_ptr< Function< std::vector< double>, double> >
			  (
					  new ProfileStepFunction
					  (
							  below_threshold_gap_opening_penalty,
							  above_threshold_gap_opening_penalty,
							  thresholds
					  )
			  );
		  gap_extension_penalty_function
			  = boost::shared_ptr< Function< std::vector< double>, double> >
			  (
					  new ProfileStepFunction
					  (
							  below_threshold_gap_extension_penalty,
							  above_threshold_gap_extension_penalty,
							  thresholds
					  )
			  );
	  }
	  else  // simple gap penalties
	  {
		  gap_opening_penalty_function
			  = boost::shared_ptr< Function< std::vector< double>, double> >( new ConstantFunction< std::vector< double>, double>( gap_opening_penalty));
		  gap_extension_penalty_function = boost::shared_ptr< Function< std::vector< double>, double> >( new ConstantFunction< std::vector< double>, double>( gap_extension_penalty));
	  }

	  NeedlemanWunschAffineGaps
		needleman_wunsch_affine_gaps
		(
				gap_opening_penalty_function,
				gap_extension_penalty_function,
				termini_gap_opening_penalty,
				termini_gap_extension_penalty,
				first_sequence,
				second_sequence,
				scores,
				matrix
		);

	  DebugWrite( "needleman wunsch affine gaps: calculate matrix:...");
	  needleman_wunsch_affine_gaps.CalculateMatrix();

	  DebugWrite( "needleman wunsch affine gaps: trace back:...");
	  alignment = needleman_wunsch_affine_gaps.TraceBack();

	  const Matrix< DynamicProgrammingMatrixElement> * matrix_ptr = &needleman_wunsch_affine_gaps.GetMatrix();
	  std::vector< double> values = ( *matrix_ptr)( matrix_ptr->GetNumberOfRows() - 1, matrix_ptr->GetNumberOfColumns() - 1).GetAffinePathWays();
	  last_element = *std::max_element( values.begin(), values.end());
	  DebugWrite ("alignment score is " <<  alignment.first);
  }

  //write.close();
  //write.clear();

  /////////////////////////////////////
  /// WRITE ALIGNMENT
  ////////////////////////////////////

  if (outputfile_aligned_sequences != "")
  {
	  DebugWrite( "write aligned sequences");
	  WriteAlignedSequences( alignment, first_sequence, second_sequence, first_fasta_id, second_fasta_id, outputfile_aligned_sequences, line_length, gap_extension_penalty);
	  DebugWrite( "done");
  }

  if( outputfile_aligned_profiles != "")
  {
	  DebugWrite( "write aligned profiles");
	  if( cmd.IsFlagSet( "profile_gap_value_for_plotting"))
	  {
		  gap_value = mystr::ConvertStringToNumericalValue< double>( cmd.GetFirstArgument( "profile_gap_value_for_plotting"));
	  }
	  WriteAlignedProfiles( alignment, first_sequence, second_sequence, outputfile_aligned_profiles, gap_extension_penalty, gap_value, last_element, write_profiles_header1, write_profiles_header2 );
	  DebugWrite( "done");
  }
  else
  {
	  std::cout << "Warning: The flag -output_aligned_profiles was not set. A file containing aligned profiles is not created " << "\n";
  }

  return 0;
}






