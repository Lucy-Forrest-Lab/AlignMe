/////////////////////////////////////////////////////////////////////////
//!
//  This file is part of the AlignMe program
//
//  Copyright (C) 2010 by Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//  AlignMe@rzg.mpg.de
//
//  AlignMe is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  AlignMe is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//!
//!
//!
//!
//!
//!
//! @author: Rene Staritzbichler
//! @date: 9.6.2011
/////////////////////////////////////////////////////////////////////////

#ifndef COMMANDS_DEFINED_H_
#define COMMANDS_DEFINED_H_



void ReadGeneralCommands( CommandLineManager &CMD, AlignmentVariables &VARS, bool &inputs_correct)
{
	  DebugWrite( "read values from options...");


	  if( CMD.IsFlagSet( "similarity_score_file"))
	  {
		  VARS.similarity_score_file = CMD.GetFirstArgument( "similarity_score_file");
#ifdef SECURE
	      if( VARS.similarity_score_file.length() == 0)
	      {
	    	  std::cerr << "ERROR: A filename is missing after the flag '-similarity_score_file'!" << "\n";
	     	  inputs_correct = false;
	      }
#endif
	  }
	  else
	  {
	      std::cerr << "\n" << "ERROR: The flag '-similarity_score_file' has to be set!!" << "\n";
     	  inputs_correct = false;
	  }


	  if( CMD.IsFlagSet( "algorithm"))
	  {
		  VARS.algorithm = CMD.GetFirstArgument( "algorithm");
#ifdef SECURE
		  if( VARS.algorithm == "")
		  {
			  std::cerr << "ERROR: if the flag '-algorithm' is set, an algorithm type has to be provided! " << "\n";
			  inputs_correct = false;
		  }
#endif
	  }
//	  else
//	  {
//		  if (inputs_correct != false)
//		  {
//			  std::cout << "NOTE: You did not specify an algorithm-type - default algorithm " << VARS.algorithm <<" is used for the alignment" << "\n";
//		  }
//	  }

	  DebugWrite( "values read");


//	  if( CMD.IsFlagSet( "norm"))
//	  {
//		  VARS.normalization_factor = CMD.GetConvertedArguments< double>( "norm")[0];
//	  }


	  if( CMD.IsFlagSet( "thresholds_for_penalties"))
	  {
		  if( CMD.IsFlagSet( "below_threshold_gap_opening_penalty"))
		  {
			  VARS.below_threshold_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "below_threshold_gap_opening_penalty"));
		  }
#ifdef SECURE
		  else
		  {
	    	  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"below_threshold_gap_opening_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "above_threshold_gap_opening_penalty"))
		  {
			  VARS.above_threshold_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "above_threshold_gap_opening_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"above_threshold_gap_opening_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "below_threshold_gap_extension_penalty"))
		  {
			  VARS.below_threshold_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "below_threshold_gap_extension_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"below_threshold_gap_extension_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "above_threshold_gap_extension_penalty"))
		  {
			  VARS.above_threshold_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "above_threshold_gap_extension_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"above_threshold_gap_extension_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "termini_gap_opening_penalty"))
		  {
			  VARS.termini_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "termini_gap_opening_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"termini_gap_opening_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "termini_gap_extension_penalty"))
		  {
			  VARS.termini_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "termini_gap_extension_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  std::cerr << "ERROR: If you want to use one or more thresholds, you have to provide a \"termini_gap_extension_penalty\" \n";
			  inputs_correct = false;
		  }
#endif
		  VARS.thresholds =  CMD.GetConvertedArguments< double>( "thresholds_for_penalties");

		  VARS.gap_opening_penalty_function
			  = ShPtr< Function< std::vector< double>, double> >
			  (
					  new ProfileStepFunction
					  (
							  VARS.below_threshold_gap_opening_penalty,
							  VARS.above_threshold_gap_opening_penalty,
							  VARS.thresholds
					  )
			  );
		  VARS.gap_extension_penalty_function
			  = ShPtr< Function< std::vector< double>, double> >
			  (
					  new ProfileStepFunction
					  (
							  VARS.below_threshold_gap_extension_penalty,
							  VARS.above_threshold_gap_extension_penalty,
							  VARS.thresholds
					  )
			  );
	  }
	  else
	  {
#ifdef SECURE
		  if( CMD.IsFlagSet( "above_threshold_gap_opening_penalty") || CMD.IsFlagSet( "above_threshold_gap_extension_penalty") || CMD.IsFlagSet( "below_threshold_gap_opening_penalty") || CMD.IsFlagSet( "below_threshold_gap_extension_penalty"))
		  {
			  std::cout << "ERROR: Gap penalties according to a threshold are not used because you did not provide a threshold-value. Please provide a threshold value or use only threshold-independent gap penalties! \n";
			  inputs_correct = false;
		  }
#endif
		  if( CMD.IsFlagSet( "gap_opening_penalty"))
		  {
			  VARS.gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "gap_opening_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  if (!CMD.IsFlagSet( "below_threshold_gap_opening_penalty") && !CMD.IsFlagSet( "above_threshold_gap_opening_penalty") && inputs_correct != false)
			  {
				  std::cout << "WARNING: No gap opening penalty provided. The default value will be used: " << VARS.gap_opening_penalty  <<"\n";
			  }
		  }
#endif

		  if( CMD.IsFlagSet( "gap_extension_penalty"))
		  {
			  VARS.gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "gap_extension_penalty"));
		  }
#ifdef SECURE
		  else
		  {
			  if (!CMD.IsFlagSet( "below_threshold_gap_extension_penalty") && !CMD.IsFlagSet( "above_threshold_gap_extension_penalty") && inputs_correct != false)
			  {
				  std::cout << "WARNING: No gap extension penalty provided. The default value will be used: " << VARS.gap_extension_penalty << "\n";
			  }
		  }
#endif

		  if( CMD.IsFlagSet( "termini_gap_opening_penalty"))
		  {
			  VARS.termini_gap_opening_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "termini_gap_opening_penalty"));
		  }
		  else
		  {
			  VARS.termini_gap_opening_penalty = VARS.gap_extension_penalty;
#ifdef SECURE
			  if (inputs_correct != false)
			  {
				  std::cout << "WARNING: No termini gap opening penalty provided. The default value will be used: " << VARS.termini_gap_opening_penalty << "\n";
			  }
				  #endif
		  }

		  if( CMD.IsFlagSet( "termini_gap_extension_penalty"))
		  {
			  VARS.termini_gap_extension_penalty = mystr::ConvertStringToNumericalValue< double>( CMD.GetFirstArgument( "termini_gap_extension_penalty"));
		  }
		  else
		  {
			  VARS.termini_gap_extension_penalty = VARS.gap_extension_penalty;
#ifdef SECURE
			  if (inputs_correct != false)
			  {
				  std::cout << "WARNING: No termini gap extension penalty provided. The default value will be used: " << VARS.termini_gap_extension_penalty  << "\n";
			  }
#endif
		  }

		  VARS.gap_opening_penalty_function
			  = ShPtr< Function< std::vector< double>, double> >( new ConstantFunction< std::vector< double>, double>( VARS.gap_opening_penalty));
		  VARS.gap_extension_penalty_function
			  = ShPtr< Function< std::vector< double>, double> >( new ConstantFunction< std::vector< double>, double>( VARS.gap_extension_penalty));

	  }
	///////////  WRITE VARs    ///////////




	  if( CMD.IsFlagSet( "output_aligned_sequences"))
	  {
		  VARS.outputfile_aligned_sequences = CMD.GetFirstArgument( "output_aligned_sequences");
#ifdef SECURE
		  if( VARS.outputfile_aligned_sequences == "")
		  {
			  std::cerr << "ERROR: if '-output_aligned_sequences' is set, a filename has to be provided! " << "\n";
	 		  inputs_correct = false;
		  }
#endif
	  }
	  else
	  {
		  if( VARS.outputfile_aligned_sequences != "" && inputs_correct != false)
		  {
			  std::cout << "WARNING: You did not provide a filename for the output of the sequence alignment. It will be written to " << VARS.outputfile_aligned_sequences << "\n";
		  }
	  }


	  if( CMD.IsFlagSet( "output_aligned_profiles"))
	  {
		  VARS.outputfile_aligned_profiles = CMD.GetFirstArgument( "output_aligned_profiles");
#ifdef SECURE
		  if( VARS.outputfile_aligned_profiles == "")
		  {
			  std::cerr << "ERROR: if '-output_aligned_profiles' is set, a filename has to be provided! " << "\n";
	 		  inputs_correct = false;
		  }
#endif
	  }
	  else
	  {
		  if( VARS.outputfile_aligned_profiles != "" && inputs_correct != false)
		  {
			  std::cout << "WARNING: You did not provide a filename for the output of the aligned profiles. They will be written to " << VARS.outputfile_aligned_profiles << "\n";
		  }
	  }

	  if( CMD.IsFlagSet( "alignment_output_format"))
	  {
		  VARS.alignment_output_format = CMD.GetFirstArgument( "alignment_output_format");
	  }


	  if( !inputs_correct)
	  {
		  exit( -1);
	  }

} // END ReadVariables( CMD, VARS)




void ReadCommandsForAlignPairs( CommandLineManager &CMD, AlignmentVariables &VARS)
{
	bool
		inputs_correct = true;

	#ifdef SECURE

		if( (CMD.IsFlagSet( "fasta_file1") && !CMD.IsFlagSet("fasta_file2")) || (!CMD.IsFlagSet( "fasta_file1") && CMD.IsFlagSet("fasta_file2")))
		{
			std::cerr << "ERROR: If you provide the flags fasta_file1 or fasta_file2, you always have to provide both flags! If you want to align profiles only, you don't need fasta files \n ";
			exit(-1);
		}

		if( (CMD.IsFlagSet( "msa_file1") && !CMD.IsFlagSet("msa_file2")) || (!CMD.IsFlagSet( "msa_file1") && CMD.IsFlagSet("msa_file2")))
		{
			std::cerr << "ERROR: It only makes sense to provide two multiple sequence alignment files, not only one!" << "\n";
			exit(-1);
		}

		if
		(
				( CMD.IsFlagSet( "fasta_file1") && CMD.IsFlagSet( "msa_file1"))
				|| ( CMD.IsFlagSet( "fasta_file2") && CMD.IsFlagSet( "msa_file2"))
		)
		{
			std::cerr << "ERROR: Please provide only MSA-files or fasta-files. Don't use both at the same time!" << "\n";
			exit(-1);
		}

	#endif


	  if( CMD.IsFlagSet( "fasta_file1"))
	  {
		  VARS.first_fasta_file = CMD.GetFirstArgument( "fasta_file1");
	  }
	  else if(CMD.IsFlagSet( "msa_file1"))
	  {
		  VARS.outputfile_aligned_sequences = "";
		  VARS.first_msa_file  = CMD.GetFirstArgument( "msa_file1");

		  // when averaging an multiple sequence alignment, this fraction determines whether a position in the alignment is ignored
		  if(CMD.IsFlagSet( "fraction_allowed_gaps"))
		  {
			  VARS.fraction_allowed_gaps = mystr::ConvertStringToNumericalValue< double> (CMD.GetFirstArgument( "fraction_allowed_gaps"));
		  }
		  #ifdef SECURE
			  else
			  {
				  std::cout << "WARNING: No value for the flag \"fraction of allow gaps\" provided. The default value will be used: " << VARS.fraction_allowed_gaps  << "\n";
			  }
		  #endif
	  }

	  if( CMD.IsFlagSet( "fasta_file2"))
	  {
		  VARS.second_fasta_file = CMD.GetFirstArgument( "fasta_file2");
	  }
	  else if (CMD.IsFlagSet( "msa_file2"))
	  {
		  VARS.second_msa_file  = CMD.GetFirstArgument( "msa_file2");
	  }


	ReadGeneralCommands( CMD, VARS, inputs_correct);
}

void ReadSequencesProfilesAndScalesForAlignPairs( AlignmentVariables &VARS)
{



	  // read in the fasta sequence(s) that are provided to AlignMe
	  if( VARS.first_fasta_file != "")
	  {
		  DebugWrite( "read first sequence:...");
		  VARS.first_fasta_id = ReadSequence( VARS.first_sequence, VARS.first_fasta_file, VARS.defined_amino_acids);
	      DebugWrite( "sequence is read: " << VARS.first_sequence);
	  }
	  else if( VARS.first_msa_file != "")
	  {
		  DebugWrite( "read first MSA:...");
		  ReadMultipleSequences( VARS.first_msa, VARS.first_msa_file, VARS.defined_amino_acids);
		  DebugWrite( "First MSA is read: " << VARS.first_msa);
	  }

	  if( VARS.second_fasta_file != "")
	  {
		  DebugWrite( "read second sequence:...");
		  VARS.second_fasta_id = ReadSequence( VARS.second_sequence, VARS.second_fasta_file, VARS.defined_amino_acids);
	      DebugWrite( "sequence is read" << VARS.second_sequence);
	  }
	  else if( VARS.second_msa_file != "")
	  {
		  DebugWrite( "read second MSA:...");
	      ReadMultipleSequences( VARS.second_msa, VARS.second_msa_file, VARS.defined_amino_acids);
	      DebugWrite( "Second MSA is read" << VARS.second_msa);
	  }

	  //reading in all information stored in the similarity_score_file
	  //generation of the scoring object based on all these inputs

	  if( VARS.first_msa_file == "" && VARS.second_msa_file == "")
	  {
		  DebugWrite( "read profiles and scales...");
		  ReadSimilarityScoreFile( VARS.scores, VARS.first_sequence, VARS.second_sequence, VARS.similarity_score_file, VARS.defined_amino_acids, VARS.write_profiles_header1, VARS.write_profiles_header2, VARS.first_fasta_id, VARS.second_fasta_id);
		  DebugWrite( "first sequence: " << VARS.first_sequence << "\nsecond sequence: " << VARS.second_sequence);
	  }

	  if (VARS.first_msa_file != "" && VARS.second_msa_file != "")
	  {
		  std::cout << "NOTE: At the moment, if you provide multiple sequences alignments, it will only work if only scales (values per AA-type) are involved and no profiles (values per AA) \n";
		  VARS.write_profiles_header1.push_back( " is based on the msa_file: " + VARS.first_msa_file + "\n");
		  VARS.write_profiles_header2.push_back( " is based on the msa_file: " + VARS.second_msa_file + "\n");
		  ReadSimilarityScoreFileMsa( VARS.scores, VARS.first_msa, VARS.second_msa, VARS.similarity_score_file, VARS.defined_amino_acids);
		  VARS.first_sequence = AverageMsa( VARS.first_msa, VARS.fraction_allowed_gaps, VARS.first_ignored_positions, VARS.first_msa_file);
		  VARS.second_sequence = AverageMsa( VARS.second_msa, VARS.fraction_allowed_gaps, VARS.second_ignored_positions, VARS.second_msa_file);
	  }

	 // DebugWrite( "read scores, including substitution matrices");
	 // ReadScores( VARS.scores, VARS.similarity_score_file, VARS.defined_amino_acids);
	 // DebugWrite( "scores read: " << VARS.scores);
}

#endif /* COMMANDS_DEFINED_H_ */
