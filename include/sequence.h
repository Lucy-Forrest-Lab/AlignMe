/*
 * sequence.h
 *
 *  Created on: Jun 21, 2010
 *      Author: marcus, rene
 */

#ifndef Sequence_H_
#define Sequence_H_

#include <vector>
#include <string>

#include "amino_acid.h"


class Sequence
: public std::vector< GeneralizedAminoAcid>
{
private:
	std::string        m_Header;    //!< fasta header
	std::string		   m_FileName;  //!< file name

public:
	//! default constructor
	Sequence( const std::string &FILE = "", const std::string &HEADER = "")
	: std::vector< GeneralizedAminoAcid>(),
	  m_Header( HEADER),
	  m_FileName( FILE)
	  {}

	//! copy constructor
	Sequence( const Sequence &ORIG)
	: std::vector< GeneralizedAminoAcid>( ORIG),
	  m_Header( ORIG.m_Header),
	  m_FileName( ORIG.m_FileName)
	  {/*std::cout << __FUNCTION__ << std::endl;*/}

	Sequence( const int &SIZE, const GeneralizedAminoAcid &AA)
	: std::vector< GeneralizedAminoAcid>( SIZE, AA),
	  m_Header(),
	  m_FileName()
	{}

	//! virtual destructor
	virtual ~Sequence(){/*std::cout << __FUNCTION__ << std::endl;*/}

	//! change fasta header
	void SetFastaHeader( const std::string &HEADER)
	{
		m_Header = HEADER;
	}

	//! get fasta header
	const std::string &GetFastaHeader() const
	{
		return m_Header;
	}

	void SetFileName( const std::string &FILE)
	{
		m_FileName = FILE;
	}

	const std::string &GetFileName() const
	{
		return m_FileName;
	}


	void ReadPositionSpecificSimilarityMatrix( const std::string &FILE);

	std::ostream &WriteAsFasta( std::ostream &STREAM, const int &LINE_LENGTH = 60) const;

}; // end class Sequence


void
WriteAsFasta( std::ostream &STREAM, const std::vector<Sequence> &MS, const int &LINE_LENGTH = 60);

#endif /* Sequence_H_ */
