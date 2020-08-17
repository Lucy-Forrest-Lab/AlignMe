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
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 10.9.2009
/////////////////////////////////////////////////////////////////////////


#ifndef PEARSON_CORRELATION_H_
#define PEARSON_CORRELATION_H_

#include <cmath>

class PearsonCorrelation
{
private:
	double  m_SumX;                                                          //!<  for the calculation of the correlation: the sum over the values of the first sequence
	double  m_SumY;                                                          //!<  for the calculation of the correlation: the sum over the profile values of the second sequence
	double  m_SumXSquare;                                                    //!<  for the calculation of the correlation:
	double  m_SumYSquare;                                                    //!<  for the calculation of the correlation:
	double  m_SumXTimesY;                                                    //!<  for the calculation of the correlation:
	double  m_Correlation;                                                   //!<  the resulting correlation of the sequence until this element
	size_t  m_Counter;


public:
	PearsonCorrelation()
	: m_SumX( 0.0),
	m_SumY( 0.0),
	m_SumXSquare( 0.0),
	m_SumYSquare( 0.0),
	m_SumXTimesY( 0.0),
	m_Correlation( 0.0),
	m_Counter( 0)
	{}


	PearsonCorrelation( const PearsonCorrelation &ORIGINAL)
	: m_SumX(      ORIGINAL.m_SumX),
	m_SumY(        ORIGINAL.m_SumY),
	m_SumXSquare(  ORIGINAL.m_SumXSquare),
	m_SumYSquare(  ORIGINAL.m_SumYSquare),
	m_SumXTimesY(  ORIGINAL.m_SumXTimesY),
	m_Correlation( ORIGINAL.m_Correlation),
	m_Counter(     ORIGINAL.m_Counter)
	{}


	//! return the total similarity value of this element
	const double &GetSumX() const
	{
		return m_SumX;
	}

	//! set the total similarity SumX of this element
	void SetSumX( const double &VALUE)
	{
		m_SumX = VALUE;
	}

	//! return the total similarity SumY of this element
	const double &GetSumY() const
	{
		return m_SumY;
	}

	//! set the total similarity SumY of this element
	void SetSumY( const double &VALUE)
	{
		m_SumY = VALUE;
	}

	//! return the total similarity SumXSquare of this element
	const double &GetSumXSquare() const
	{
		return m_SumXSquare;
	}

	//! set the total similarity SumXSquare of this element
	void SetSumXSquare( const double &VALUE)
	{
		m_SumXSquare = VALUE;
	}

	//! return the total similarity SumXSquare of this element
	const double &GetSumYSquare() const
	{
		return m_SumYSquare;
	}

	//! set the total similarity value of this element
	void SetSumYSquare( const double &VALUE)
	{
		m_SumYSquare = VALUE;
	}

	//! return the total similarity value of this element
	const double &GetSumXTimesY() const
	{
		return m_SumXTimesY;
	}

	//! set the total similarity value of this element
	void SetSumXTimesY( const double &VALUE)
	{
		m_SumXTimesY = VALUE;
	}

	//! return the total similarity value of this element
	const double &GetCorrelation() const
	{
		return m_Correlation;
	}

	//! set the total similarity value of this element
	void SetCorrelation( const double &VALUE)
	{
		m_Correlation = VALUE;
	}

	void Add( const double &X, const double &Y)
	{
		m_SumX += X;
		m_SumY += Y;
		m_SumXSquare += X * X;
		m_SumYSquare += Y * Y;
		m_SumXTimesY += X * Y;
		++m_Counter;
	}


	void Clear()
	{
		m_SumX        = 0.0;
		m_SumY        = 0.0;
		m_SumXSquare  = 0.0;
		m_SumYSquare  = 0.0;
		m_SumXTimesY  = 0.0;
		m_Counter     = 0;
		m_Correlation = 0.0;
	}

	double operator()( const std::vector< double> &V1, const std::vector< double> &V2)
	{
		assert( V1.size() == V2.size());
		Clear();
		for( std::vector< double>::const_iterator itr_1( V1.begin()), itr_2( V2.begin()); itr_1 != V1.end(); ++itr_1, ++itr_2)
		{
			m_SumX += *itr_1;
			m_SumY += *itr_2;
			m_SumXSquare += *itr_1 * *itr_1;
			m_SumYSquare += *itr_2 * *itr_2;
			m_SumXTimesY += *itr_1 * *itr_2;
			++m_Counter;
		}
		return UpdateCorrelation();
	}

	double UpdateCorrelation()
	{
		return m_Correlation = ( m_Counter * m_SumXTimesY - m_SumX * m_SumY)
				/ ( sqrt( m_Counter * m_SumXSquare - m_SumX * m_SumX) * sqrt( m_Counter * m_SumYSquare - m_SumY * m_SumY));
	}

}; // end class PearsonCorrelation


#endif /* PEARSON_CORRELATION_H_ */
