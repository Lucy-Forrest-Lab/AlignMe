/*
 * allowed_flags.h
 *
 *  Created on: Dec 21, 2011
 *      Author: marcus
 */

#ifndef ALLOWED_FLAGS_H_
#define ALLOWED_FLAGS_H_

void DefineAllowedFlagsGeneral( std::map<std::string,std::string> &ALLOWED_FLAGS)
{
	  ALLOWED_FLAGS["similarity_score_file"] = "allowed";
	  ALLOWED_FLAGS["gap_opening_penalty"] = "allowed";
	  ALLOWED_FLAGS["gap_extension_penalty"] = "allowed";
	  ALLOWED_FLAGS["algorithm"] = "allowed";
	  ALLOWED_FLAGS["termini_gap_opening_penalty"] = "allowed";
	  ALLOWED_FLAGS["termini_gap_extension_penalty"] = "allowed";
	  ALLOWED_FLAGS["output_aligned_sequences"] = "allowed";
	  ALLOWED_FLAGS["output_aligned_profiles"] = "allowed";
	  ALLOWED_FLAGS["thresholds_for_penalties"] = "allowed";
	  ALLOWED_FLAGS["below_threshold_gap_opening_penalty"] = "allowed";
	  ALLOWED_FLAGS["below_threshold_gap_extension_penalty"] = "allowed";
	  ALLOWED_FLAGS["above_threshold_gap_opening_penalty"] = "allowed";
	  ALLOWED_FLAGS["above_threshold_gap_extension_penalty"] = "allowed";
	  ALLOWED_FLAGS["profile_gap_value_for_plotting"] = "allowed";
	  ALLOWED_FLAGS["alignment_output_format"] = "allowed";
	  ALLOWED_FLAGS["anchors"] = "allowed";
	  ALLOWED_FLAGS["help"] = "allowed";
	  ALLOWED_FLAGS["-help"] = "allowed";
	#ifdef SECURE
	  ALLOWED_FLAGS["completed"] = "allowed";
	#endif
}



void DefineAllowedFlagsForAlignPairs( std::map<std::string,std::string> &ALLOWED_FLAGS)
{
	  DefineAllowedFlagsGeneral( ALLOWED_FLAGS);
	  ALLOWED_FLAGS["fasta_file1"] = "allowed";
	  ALLOWED_FLAGS["fasta_file2"] = "allowed";
	  ALLOWED_FLAGS["msa_file1"] = "allowed";
	  ALLOWED_FLAGS["msa_file2"] = "allowed";
	  ALLOWED_FLAGS["fraction_allowed_gaps"] = "allowed";
	  ALLOWED_FLAGS["extract_from_MSA_sequences_with_ids"] = "allowed";
	  ALLOWED_FLAGS["output_extracted_sequences"] = "allowed";
}



#endif /* ALLOWED_FLAGS_H_ */
