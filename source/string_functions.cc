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
//!  Functions on strings.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#include "../include/string_functions.h"
#include <cmath>

namespace mystr
{

    //! trims spaces from beginning and end of a copied string
    std::string TrimString( const std::string &STRING)
    {
      //searches for last leading space
      const std::string::size_type pos1( STRING.find_first_not_of( " \n\t\r\0"));
      //searches for first tailing space
      const std::string::size_type pos2( STRING.find_last_not_of( " \n\t\r\0"));

      //returns substring from pos1 of length pos2 - pos1 + 1
      return STRING.substr( ( pos1 == std::string::npos ? 0 : pos1), ( pos2 == std::string::npos ? 0 : pos2 - pos1 + 1));
    }


   //! test whether string is numerical (double, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    bool IsNumerical( const std::string &STRING)
    {
      //remove leading and tailing spaces
      const std::string trimmed_string( TrimString( STRING));

      //check that string is not empty
      if( trimmed_string.empty())
      {
        return false;
      }
      bool contains_point( false);
      bool contains_eE( false);

      if( trimmed_string.size() == 1)
      {
    	  if( trimmed_string[0] < '0' || trimmed_string[0] > '0' + 9)
    	  {
    		  return false;
    	  }
      }


      // test all characters
      for( size_t i( 0); i < trimmed_string.length(); ++i)
      {
        // character is between 0 and 9 -> ok
        if( trimmed_string[i] >= '0' && trimmed_string[i] <= '0' + 9)
        {
          continue;
        }
        //character is '.' and the only point-> ok
        else if( trimmed_string[i] == '.')
        {
        	if( !contains_point)
        	{
			  contains_point = true;
			  continue;
        	}
        	else
        	{
            	return false;
        	}
        }
        else if( trimmed_string[i] == 'e' || trimmed_string[i] == 'E')
        {
        	if( !contains_eE)
        	{
        		contains_eE = true;
        		continue;
        	}
        	else
        	{
        		return false;
        	}
        }
        //is charater '+' or '-' -> ok
        else if( trimmed_string[i] == '-' || trimmed_string[i] == '+')
        {
#ifdef SECURE
        	if( i == 0 || trimmed_string[ i-1] == 'e' || trimmed_string[ i-1] == 'E')
        	{
        		continue;
        	}
        	else
        	{
        		return false;
        	}
#else
        	continue;
#endif
        }
        //if neither of above cases -> not a numerical value
        else
        {
          return false;
        }

      }

      // return true
      return true;
    }


    std::vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER)
    {
      std::vector< std::string> result;
      size_t a(0),b(0);
      do
        {
          a = STRING.find_first_not_of( SPLITTER,b);
          b = STRING.find_first_of( SPLITTER, a);
          if( b != std::string::npos)
          { result.push_back( STRING.substr( a, b-a));}
          else
          { result.push_back( STRING.substr( a));}
        } while( a != std::string::npos && b != std::string::npos && STRING.find_first_not_of( SPLITTER,b) != std::string::npos);

      return result;
    }


	std::string &AllToLowerCase( std::string &STRING)  // BE AWARE: this always mutates string, non const
	{
		std::transform( STRING.begin(), STRING.end(), STRING.begin(), ::tolower);
		return STRING;
	}

	std::string ToLower( const std::string &STRING)
	{
		std::string str( STRING);
		std::transform( str.begin(), str.end(), str.begin(), ::tolower);
		return str;
	}

	std::string &AllToUpperCase( std::string &STRING)
	{
		std::transform( STRING.begin(), STRING.end(), STRING.begin(), toupper);
		return STRING;
	}

	std::string ToUpper( const std::string &STRING)
	{
		std::string str( STRING);
		std::transform( str.begin(), str.end(), str.begin(), toupper);
		return str;
	}

	size_t StringToSizeT( const std::string &STRING) // no encryption
	{
		assert( STRING.size() < 10);
		int i(0);
		size_t result(0);
		for( std::string::const_reverse_iterator itr = STRING.rbegin(); itr != STRING.rend(); ++itr, ++i)
		{
			result += ( size_t( *itr) - 32) * size_t( pow( 100.0, i));
		}
		return result;
	}

	std::string SizeTToString( const size_t &VALUE)
	{
		size_t value( VALUE);
		std::string str( "");
		while( value > 0)
		{
			size_t rest = value % 100;
			char c ( rest + 32);
//			std::cout << c;
			str.insert( str.begin(), c);
			value = size_t( value / 100);
		}
//		std::cout << "\n";
		return str;
	}



} // end namespace mystr
