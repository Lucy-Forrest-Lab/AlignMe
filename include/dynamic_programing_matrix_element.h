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
//!  NEEDS CLEANUP
//! Beside the similarity value in a matrix element the DynamicProgrammingMatrixElement
//! contains links to previous elements for determining the gap penalty.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////

#ifndef DYNAMIC_PROGRAMMING_MATRIX_ELEMENT_H
#define DYNAMIC_PROGRAMMING_MATRIX_ELEMENT_H

#include <cstddef>

#include <boost/shared_ptr.hpp>

#include "pearson_correlation.h"

class DynamicProgrammingMatrixElement
{
private:
	double                                  m_Value;                                 //!<  the total similarity score at this element of the dynamic programming matrix
	std::pair< size_t, size_t>              m_PreviousElementIndices;                //!<  indices of the previous element, the one with the best score
	ShPtr< PearsonCorrelation>  m_Correlation;
	std::vector<double>                     m_AffinePathWays;                          //!<  needed for Needleman Wunsch Affine Gap (2 shells)

public:

	// default constructor
	DynamicProgrammingMatrixElement()
	: m_Value( 0.0),
	m_PreviousElementIndices(),
	m_Correlation(),
    m_AffinePathWays( 9, 0.0)
	{}

    // construct from data
    DynamicProgrammingMatrixElement( const double &VALUE, const std::pair< size_t, size_t> &INDICES)
    : m_Value( VALUE),
    m_PreviousElementIndices( INDICES),
	m_Correlation(),
	m_AffinePathWays( 9, 0.0)
    {}
    
    // construct from data (for convenience)
    DynamicProgrammingMatrixElement( const double &VALUE, const size_t &INDEX_I, const size_t &INDEX_J)
	: m_Value( VALUE),
	m_PreviousElementIndices( std::make_pair( INDEX_I, INDEX_J)),
	m_Correlation(),
	m_AffinePathWays( 9, 0.0)
	{}

    // construct from data (for convenience)
    DynamicProgrammingMatrixElement( const double &VALUE, const size_t &INDEX_I, const size_t &INDEX_J, const ShPtr< PearsonCorrelation> &CORRELATION)
	: m_Value( VALUE),
	m_PreviousElementIndices( std::make_pair( INDEX_I, INDEX_J)),
	m_Correlation( CORRELATION),
	m_AffinePathWays( 9, 0.0)
	{}


    // copy constructor
    DynamicProgrammingMatrixElement( const DynamicProgrammingMatrixElement &ORIGINAL)
    : m_Value( ORIGINAL.m_Value),
    m_PreviousElementIndices( ORIGINAL.m_PreviousElementIndices),
	m_Correlation( ORIGINAL.m_Correlation),
	m_AffinePathWays( ORIGINAL.m_AffinePathWays)
    {}

    
	// destructor
	~DynamicProgrammingMatrixElement(){/*std::cout << __FUNCTION__ << std::endl;*/}

	//! add the total similarity VALUE to the path way given by INDEX
	void AddValue(const int &INDEX, const double &VALUE)
	{
		m_AffinePathWays[INDEX] += VALUE;
	}

	//! add the total similarity VALUE to all path ways
	void AddValue(const double &VALUE)
	{
		for( std::vector<double>::iterator itr = m_AffinePathWays.begin(); itr != m_AffinePathWays.end(); ++itr)
		{
			*itr += VALUE;
		}
		m_Value += VALUE; // todo: own function???
	}

	const std::vector<double> &GetAffinePathWays() const
	{
		return m_AffinePathWays;
	}

	const double &GetAffinePathWay( const int &INDEX) const
	{
		return m_AffinePathWays[INDEX];
	}

	//! from the current element search one of the three neighbor elements determined by MAIN_PATH_WAY_INDEX (0..2) for the best total path way (0..8)
	const double BestSubPathWayScore( const size_t &MAIN_PATH_WAY_INDEX) const
	{
		double best = m_AffinePathWays[ 3 * MAIN_PATH_WAY_INDEX];
		for( size_t i = 3 * MAIN_PATH_WAY_INDEX + 1; i < 3 * MAIN_PATH_WAY_INDEX + 3; i++)
		{
			if( m_AffinePathWays[i] > best)
			{
				best = m_AffinePathWays[i];
			}
		}
		return best;
	}

	//! return the total similarity value of this element
	const double &GetValue() const
	{
		return m_Value;
	}

	//! set the total similarity value of this element
	void SetValue( const double &VALUE)
	{
		m_Value = VALUE;
	}

	// returns the indices of the previous matrix element (where it comes from)
	const std::pair< size_t, size_t> &GetIndicesOfPreviousElement() const
	  {
	    return m_PreviousElementIndices;
	  }

	std::ostream &Write( std::ostream &STREAM) const
	  {
	    STREAM  << m_Value;
//	    STREAM << "(" << m_Value << "   " << m_PreviousElementIndices.first << " " << m_PreviousElementIndices.second << " " << m_AffinePathWays << ")" << std::endl;
	    return STREAM;
	  }

	const ShPtr< PearsonCorrelation> &GetCorrelation() const
	{
		return m_Correlation;
	}

	void SetCorrelation( const ShPtr< PearsonCorrelation> &CORRELATION)
	{
		m_Correlation = CORRELATION;
	}

	double operator *= ( const double &VALUE)
	{
		return (m_Value *= VALUE);
	}
};


inline
std::ostream &operator << ( std::ostream &STREAM, const DynamicProgrammingMatrixElement &ELEMENT)
{
  return ELEMENT.Write( STREAM);
}

inline
double operator *( const double &VALUE, const DynamicProgrammingMatrixElement &ELEMENT)
{
	return VALUE * ELEMENT.GetValue();
}

inline
double operator *( const DynamicProgrammingMatrixElement &ELEMENT, const double &VALUE)
{
	return VALUE * ELEMENT.GetValue();
}

#endif
