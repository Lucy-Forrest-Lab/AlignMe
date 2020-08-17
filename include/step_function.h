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
//!  Simple step function
//! f(x) = a if x < threshold
//! f(x) = b if x >= threshold
//!
//!
//! @author: Kamil Khafizov, Rene Staritzbichler, Marcus Stamm
//! @date: 22.10.2009
/////////////////////////////////////////////////////////////////////////


#ifndef STEP_FUNCTION_H_
#define STEP_FUNCTION_H_

template< typename t_INPUT, typename t_RETURN>
class StepFunction
: public Function< t_INPUT, t_RETURN>
{
protected:
	t_RETURN   m_LowerValue;
	t_RETURN   m_UpperValue;
	t_INPUT    m_Threshold;

public:
	StepFunction
	(
			  const double &LOWER_VALUE,
			  const double &UPPER_VALUE,
			  const double &THRESHOLD
	)
	: m_LowerValue( LOWER_VALUE),
	m_UpperValue( UPPER_VALUE),
	m_Threshold( THRESHOLD)
	{}

	virtual ~StepFunction(){}

	virtual t_RETURN operator()( const t_INPUT &VALUE)
	{
		if( VALUE < m_Threshold)
		{
			return m_LowerValue;
		}
		return m_UpperValue;
	}
};


#endif /* STEP_FUNCTION_H_ */
