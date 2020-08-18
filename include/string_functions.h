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
//!  Collection of string functions.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef STRING_FUNCTIONS
#define STRING_FUNCTIONS

#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>

#include "definitions.h"


namespace mystr
{

	//! trims spaces from beginning and end of a copied string
	std::string TrimString( const std::string &STRING);


	//! test whether string is numerical (double, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
	bool IsNumerical( const std::string &STRING);


//    std::string &RemoveSubstr( std::string &STR, const std::string SUBSTR)
//    { exit( -1); return STR;}

	//! converts std::string to object of template class T1
	template< class T1>
	T1 ConvertStringToNumericalValue( const std::string &STRING)
	{
		// remove tabs, spaces, new lines, etc.
		const std::string string( TrimString( STRING));

		// if string is empty an undefined is returned
		if( string.empty())
		{ return std::numeric_limits< T1>::max();}

//		if( !IsNumerical( string))
//		{
//			std::cerr << " ERROR: " << __FUNCTION__ << ": you provided (" << string << ") which is not numerical. Only integers or fractions are allowed (no commas)!" << "\n";
//			exit(-1);
//		}
		T1 new_t;
		std::stringstream iss( string);
		iss >> new_t;
		return new_t;
	}

	template< class T>
	std::vector< T>
	ConvertStringVector( const std::vector< std::string> &VEC)
	{
		std::vector< T>
			transformed( VEC.size());
		std::vector< std::string>::const_iterator
			in_itr = VEC.begin();
		typename std::vector< T>::iterator
			out_itr = transformed.begin();
		for( ; in_itr != VEC.end(); ++in_itr, ++out_itr)
		{
			*out_itr = ConvertStringToNumericalValue< T>( *in_itr);
		}
		return transformed;
	}

	std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER = " ");

	template< class T>
	std::vector< T> SplitAndConvertString( const std::string &STRING, const std::string &SPLITTER = " ")
	{
		return ConvertStringVector< T>( SplitString( STRING, SPLITTER));
	}

	template<class T>
	std::string NumericalValueToString(const T& VALUE)
	{
	    std::stringstream strm;
	    strm << VALUE;
	    return strm.str();
	}

//	template<class T>
//	std::string NumericalValueToString(const T& VALUE, const size_t &PRECISION)
//	{
//	    std::stringstream strm;
//	    strm << std::setprecision( PRECISION);
////	    T value( T( size_t( VALUE * pow( 10, PRECISION)) / pow( 10, PRECISION))
//	    strm << VALUE;
//	    return strm.str();
//	}

	//! transforms given string into lower case
	std::string &AllToLowerCase( std::string &STRING);

	//! returns a lower case copy of string
	std::string ToLower( const std::string &STRING);

	//! transforms given string to upper case
	std::string &AllToUpperCase( std::string &STRING);

	//! returns an upper case copy of string
	std::string ToUpper( const std::string &STRING);

	//! translates strings up to size 9 to size_t; needed for example for sending messages among processes -msgsnd, msgrcv cannot handle strings correctly on all OS
	size_t StringToSizeT( const std::string &STRING);

	//! translates size_t to string up to size 9
	std::string SizeTToString( const size_t &VALUE);

	//! searches for spaces and removes them from the string
	inline
	std::string RemoveSpacesFromString( const std::string &STRING)
	{
		std::string cleaned_string;
		for( size_t i = 0 ;i < STRING.size(); i++ )
		{
			if( STRING[i] != ' ') cleaned_string.push_back( STRING[i]);
		}
		return cleaned_string;
	  }

	inline
	bool
	IsCapitolLetter( const char &CHAR)
	{
		size_t nr( CHAR);
		return nr >= 65 && nr <= 90;
	}

	inline
	bool
	IsSmallLetter( const char &CHAR)
	{
		size_t nr( CHAR);
		return nr >= 97 && nr <= 122;
	}

	inline
	bool
	IsLetter( const char &CHAR)
	{
		size_t nr( CHAR);
		return ( nr >= 65 && nr <= 90) || ( nr >= 97 && nr <= 122);
	}

	inline
	bool
	IsNumber( const char &CHAR)
	{
		size_t nr( CHAR);
		return nr >= 48 && nr <= 57;
	}

	inline
	std::string
	RemoveEnding( const std::string &ORIG)
	{
		return ORIG.substr( 0, ORIG.find( "."));
	}


} // end namespace mystr

#endif
