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
//!  A function that returns always its stored value,
//! regardless of input value
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 21.10.2009
/////////////////////////////////////////////////////////////////////////

#ifndef CONSTANT_FUNCTION_H_
#define CONSTANT_FUNCTION_H_


template< typename t_INPUT, typename t_RETURN>
class ConstantFunction
: public Function< t_INPUT, t_RETURN>
{
protected:
	t_RETURN   m_ReturnValue;
public:
	ConstantFunction( const t_RETURN &VALUE)
	: m_ReturnValue( VALUE)
	{}

	virtual ~ConstantFunction(){}

	virtual t_RETURN operator()( const t_INPUT &VALUE)
	{
		return m_ReturnValue;
	}
};

#endif /* CONSTANT_FUNCTION_H_ */
