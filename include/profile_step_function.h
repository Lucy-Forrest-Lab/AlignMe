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
//!  This function handles the min/max gap penalties
//! which allow to conserve e.g. secondary structure elements.
//!
//!
//! @author: Kamil Khafizov, Rene Staritzbichler, Marcus Stamm
//! @date: 21.10.2009
/////////////////////////////////////////////////////////////////////////



#ifndef PROFILE_STEP_FUNCTION_H_
#define PROFILE_STEP_FUNCTION_H_

class ProfileStepFunction
: public Function< std::vector< double>, double>
{
private:
	std::vector< double>  m_Thresholds;
	double 				  m_UpperValue;
	double                m_LowerValue;


public:

	ProfileStepFunction
	(
			  const double &LOWER_VALUE,
			  const double &UPPER_VALUE,
			  const std::vector< double> &THRESHOLDS
//			  const int &NR_OF_PROFILES
	)
	: m_Thresholds( THRESHOLDS),
	  m_UpperValue( UPPER_VALUE),
	  m_LowerValue( LOWER_VALUE)
	{}

	virtual ~ProfileStepFunction(){}

	virtual double operator()( const std::vector< double> &PROFILES)
	{
		size_t size = m_Thresholds.size();
#ifdef SECURE
		if( size > PROFILES.size())
		{
			std::cerr << "ERROR: Number of thresholds (" << size << ") is larger than the number of submitted profiles (" << PROFILES.size() << ")\n";
			exit( 1);
		}
#endif
		for( size_t i = 0; i < size; ++i)
		{
			if( PROFILES[i] < m_Thresholds[i])
			{
				return m_LowerValue;
			}
		}
		return m_UpperValue;
	}

};



#endif /* PROFILE_STEP_FUNCTION_H_ */
