/*
 * sequence.h
 *
 *  Created on: Jun 21, 2010
 *      Author: marcus, rene
 */


#include "../include/sequence.h"

#include "../include/string_functions.h"


	std::ostream &Sequence::WriteAsFasta( std::ostream &STREAM, const int &LINE_LENGTH) const
	{
		STREAM << ">" << m_Header << " " << size() << std::endl;
		int i = 1;
		for( std::vector< GeneralizedAminoAcid>::const_iterator itr = this->begin(); itr != this->end(); ++itr, ++i)
		{
			STREAM << itr->GetType();
			if( i % LINE_LENGTH == 0)
			{
				STREAM << std::endl;
			}
		}
		if( i % LINE_LENGTH != 0)
		{
			STREAM << std::endl;
		}
		return STREAM;
	}


void
WriteAsFasta( std::ostream &STREAM, const std::vector< Sequence> &MULTSEQALG, const int &LINE_LENGTH)
{
	  for( std::vector< Sequence>::const_iterator itr = MULTSEQALG.begin(); itr != MULTSEQALG.end(); ++itr)
	  {
		  itr->WriteAsFasta( STREAM, LINE_LENGTH);
	  }
}

