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
//!  The abstract base class for all functors.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////

#ifndef FUNCTION_H
#define FUNCTION_H

#include <limits>
#include <iostream>

#include "definitions.h"

/////////////////////////////////////////////////
// a class to derive all sorts of functions from
/////////////////////////////////////////////////

template< typename t_INPUT, typename t_RETURN>
class Function
{
 public:
  virtual ~Function(){}

	virtual t_RETURN operator()(const t_INPUT &DATA) = 0; // = 0 means this class is pure virtual, means one cannot create an object of this !! it forces to implement this function in any derived class!!!

	virtual std::ostream &Write( std::ostream &STREAM) const {
		STREAM << "Function::Write() ... why?" << "\n";
		return STREAM;
	}

	virtual int GetClassID() const
	{
		return std::numeric_limits< int>::max();
	}

};

template<typename t_INPUT, typename t_RETURN>
inline std::ostream & operator << ( std::ostream &STREAM, ShPtr< Function<t_INPUT, t_RETURN> > &FUNCT)
{
	STREAM << __FUNCTION__ << "\n";
	FUNCT->Write(STREAM);
	return STREAM;
}

#endif
