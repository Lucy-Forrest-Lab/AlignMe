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
//!  A collections of functions and operators using std containers.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef EXTERNAL_STD_FUNCTIONS_H
#define EXTERNAL_STD_FUNCTIONS_H

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <map>
#include <list>
#include <boost/shared_ptr.hpp>


template< typename t_DATA>
inline
std::vector< t_DATA> &operator *= ( std::vector< t_DATA> &VEC, const double &FACTOR)
{
	DebugWrite( __PRETTY_FUNCTION__);
	for( typename std::vector< t_DATA>::iterator itr = VEC.begin(); itr != VEC.end(); ++itr)
	{
		*itr *= FACTOR;
	}
	return VEC;
}

template< typename t_DATA>
inline
std::vector< t_DATA> operator * ( const std::vector< t_DATA> &VEC, const double &FACTOR)
{
	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_DATA> tmp( VEC);
	for( typename std::vector< t_DATA>::iterator itr = tmp.begin(); itr != tmp.end(); ++itr)
	{
		*itr *= FACTOR;
	}
	return tmp;
}

template< typename t_DATA>
inline
std::vector< t_DATA> operator * ( const double &FACTOR, const std::vector< t_DATA> &VEC)
{
//	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_DATA> tmp( VEC);
	for( typename std::vector< t_DATA>::iterator itr = tmp.begin(); itr != tmp.end(); ++itr)
	{
		*itr *= FACTOR;
	}
	return tmp;
}

template< typename t_DATA>
inline
std::string GetStdVectorClassName( const std::vector< t_DATA> &VEC)
{ return "std::vector<t_DATA>";}

template< typename t_DATA>
inline
std::istream& operator >> ( std::istream &STREAM, std::vector< t_DATA> &VEC)
{
  std::string str;
  STREAM >> str;
  assert( str == GetStdVectorClassName< t_DATA>( VEC));
  size_t size;
  STREAM >> size;
  std::vector< t_DATA> tmp( size);
  for( typename std::vector< t_DATA>::iterator itr( tmp.begin()); itr != tmp.end(); ++itr)
    { STREAM >> *itr;}
  VEC = tmp;
  return STREAM;
}

template< typename t_DATA>
inline
std::vector< t_DATA> MinPerElement(const std::vector< t_DATA> &VEC1, const std::vector< t_DATA> &VEC2)
{
	std::vector< t_DATA> vec (VEC1.size());

	typename std::vector< t_DATA>::iterator vec_itr = vec.begin();
	typename std::vector< t_DATA>::const_iterator vec1_itr = VEC1.begin(), vec2_itr = VEC2.begin();

	for (; vec_itr != vec.end(); ++vec_itr, ++vec1_itr, ++vec2_itr)
	{
		*vec_itr = std::min( *vec1_itr, *vec2_itr);
	}
	return vec;
}

//template< typename t_DATA>
//inline
//std::ostream& operator << ( std::ostream &STREAM, const std::vector< t_DATA> &VEC)
//{
//	STREAM << GetStdVectorClassName<t_DATA>( VEC) << "\n";
//	STREAM << VEC.size() << "\n";
//	if( VEC.size() > 0)
//    {
//      for( typename std::vector< t_DATA>::const_iterator itr( VEC.begin()); itr != VEC.end(); ++itr)
//	{ STREAM << "<" << *itr << ">  ";}
//      STREAM << "\n";
//    }
//	return STREAM;
//}

template< typename t_FIRST, typename t_SECOND>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::pair< t_FIRST, t_SECOND> &PAIR)
{
	STREAM << "std::pair" << "\n";
	STREAM << PAIR.first << "   " << PAIR.second << "\n";
	return STREAM;
}

template< typename t_DATA>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< t_DATA> &LIST)
{
  STREAM << "std::vector" << "\n";
  STREAM << LIST.size() << "\n";
  if( LIST.size() == 0)
  {
	  return STREAM;
  }
  for( typename std::vector< t_DATA>::const_iterator itr( LIST.begin()); itr != LIST.end(); ++itr)
  {
	  STREAM  << *itr << "  ";
  }
  STREAM << "\n";
  return STREAM;
}


template< typename t_KEY, typename t_VALUE>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::map< t_KEY, t_VALUE> &MAP)
{
	STREAM << "std::map" << "\n";
	for( typename std::map< t_KEY, t_VALUE>::const_iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
	{
		STREAM << itr->first << "   " << itr->second << "\n";
	}
	return STREAM;
}

template< typename t_KEY, typename t_VALUE>
inline
std::ostream& operator << ( std::ostream &STREAM, const std::multimap< t_KEY, t_VALUE> &MAP)
{
	STREAM << "std::multimap" << "\n";
	for( typename std::multimap< t_KEY, t_VALUE>::const_iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
	{
		STREAM << itr->first << "   " << itr->second << "\n";
	}
	return STREAM;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::ostream& operator >> ( std::ostream &STREAM, std::pair< t_FIRST, t_SECOND> &PAIR)
{
	t_FIRST first;
	t_SECOND second;
	STREAM >> first >> second;
	PAIR.first = first;
	PAIR.second = second;
	return STREAM;
}

inline
std::istream& operator >> ( std::istream &STREAM, std::vector< double> &VEC)
{
  std::string str;
  STREAM >> str;
  assert( str == "std::vector< double>");
  size_t size;
  STREAM >> size;
  std::vector< double> tmp( size);
  for( size_t i = 0; i < size; ++i)
    { STREAM >> tmp[i];}
  VEC = tmp;
  return STREAM;
}

inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< double> &VEC)
{
  STREAM << "std::vector<double>" << "\n";
  STREAM << VEC.size() << "\n";
  if( VEC.size() > 0)
    {
      for( size_t i = 0; i < VEC.size(); ++i)
	{ STREAM << VEC[i] << "  ";}
      STREAM << "\n";
    }
  return STREAM;
}


inline
std::istream& operator >> ( std::istream &STREAM, std::vector< size_t> &VEC)
{
  std::string str;
  STREAM >> str;
  assert( str == "std::vector< size_t>");
  size_t size;
  STREAM >> size;
  std::vector< size_t> tmp( size);
  for( size_t i = 0; i < size; ++i)
    { STREAM >> tmp[i];}
  VEC = tmp;
  return STREAM;
}

inline
std::ostream& operator << ( std::ostream &STREAM, const std::vector< size_t> &VEC)
{
  STREAM << "std::vector< size_t>" << "\n";
  STREAM << VEC.size() << "\n";
  for( size_t i = 0; i < VEC.size(); ++i)
    { STREAM << VEC[i] << "   ";}
  STREAM << "\n";
  return STREAM;
}


template< typename t_TYPE>
inline
std::vector< t_TYPE> operator + ( const std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
//	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_TYPE> vec( V1);
	return vec += V2;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> operator += ( std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
	if( V1.size() != V2.size())
	{
		std::cout << "ERROR: adding two vectors makes only sense if the two vectors are of equal size! " << std::endl;
		exit(-1);
	}
	typename std::vector< t_TYPE>::const_iterator vitr = V2.begin();
	typename std::vector< t_TYPE>::iterator itr = V1.begin();
	for(  ;itr != V1.end(); ++itr, ++vitr)
	{
		*itr += *vitr;
	}
	return V1;
}



template< typename t_TYPE>
inline
std::vector< t_TYPE> operator - ( const std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
//	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_TYPE> vec( V1);
	return vec -= V2;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> operator -= ( std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
	if( V1.size() != V2.size())
	{
		std::cout << "ERROR: adding two vectors makes only sense if the two vectors are of equal size! " << std::endl;
		exit(-1);
	}
	typename std::vector< t_TYPE>::const_iterator vitr = V2.begin();
	typename std::vector< t_TYPE>::iterator itr = V1.begin();
	for(  ;itr != V1.end(); ++itr, ++vitr)
	{
		*itr -= *vitr;
	}
	return V1;
}


template< typename t_FIRST, typename t_SECOND>
inline
std::pair< t_FIRST, t_SECOND>  operator += ( const std::pair< t_FIRST, t_SECOND> &PAIR_A, const std::pair< t_FIRST, t_SECOND> &PAIR_B)
{
	DebugWrite( __PRETTY_FUNCTION__);
	std::pair< t_FIRST, t_SECOND> pair;
	pair.first = PAIR_A.first + PAIR_B.first;
	pair.second = PAIR_B.second + PAIR_B.second;
	return pair;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::vector< t_SECOND> CastVector( const std::vector< t_FIRST> &SECOND)
{
	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_SECOND> vector( SECOND.size());
	typename std::vector< t_SECOND>::iterator first_itr( vector.begin());
	typename std::vector< t_FIRST>::const_iterator second_itr( SECOND.begin());
	for( ; second_itr != SECOND.end(); ++first_itr, ++second_itr)
	{
		*first_itr = t_SECOND( **second_itr);
	}
	return vector;
}

template< typename t_FIRST, typename t_SECOND>
inline
std::map< t_FIRST, t_SECOND> &operator /= ( std::map< t_FIRST, t_SECOND> &MAP, const t_SECOND &OBJECT)
{
	DebugWrite( __PRETTY_FUNCTION__);
	for( typename std::map< t_FIRST, t_SECOND>::iterator itr = MAP.begin(); itr != MAP.end(); ++itr)
	{
		itr->second /= OBJECT;
	}
	return MAP;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> Merge( const std::vector< t_TYPE> &V1, const std::vector< t_TYPE> &V2)
{
	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_TYPE> v( V1);
	v.insert( v.end(), V2.begin(), V2.end());
	return v;
}

template< typename t_TYPE>
inline
std::vector< t_TYPE> RemoveDuplicates( const std::vector< t_TYPE> &VECTOR)
{
	DebugWrite( __PRETTY_FUNCTION__);
	std::vector< t_TYPE> v( VECTOR), result;
	typename std::vector< t_TYPE>::iterator new_end = unique(v.begin(), v.end());
	result.insert( result.end(), v.begin(), new_end);
	return result;
}

template< typename t_TYPE>
inline
bool IsElement( const std::vector< t_TYPE> &VECTOR, const t_TYPE &VALUE)
{
	return(  std::find( VECTOR.begin(), VECTOR.end(), VALUE) != VECTOR.end());
}

#endif
