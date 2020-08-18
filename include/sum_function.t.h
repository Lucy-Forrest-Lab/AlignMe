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
//!  A linear combination of functors, itself a functor.
//!
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef SUM_FUNCTION_T_H
#define SUM_FUNCTION_T_H

#include "function.t.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <iostream>



template< typename t_INPUT, typename t_RETURN>
class SumFunction
  : public Function< t_INPUT, t_RETURN>
{
 protected:
  std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > > m_FunctionsAndWeights;   //!< vector where each element holds a weight and a function


 public:

  //! construct a sum function from two pointers to functions
  SumFunction( const ShPtr< Function< t_INPUT, t_RETURN> > &F1, const ShPtr< Function< t_INPUT, t_RETURN> > &F2)
  : m_FunctionsAndWeights( 2) // construct vector of size 2
  {
 	   m_FunctionsAndWeights[0] = std::make_pair( 1.0, F1);
 	   m_FunctionsAndWeights[1] = std::make_pair( 1.0, F2);
  }


  //! construct a sum function from a weight and a pointer to a function
  SumFunction( const double &WEIGHT, const ShPtr< Function< t_INPUT, t_RETURN> > &FUNCTION)
  : m_FunctionsAndWeights( 1) // construct vector of size 2
  {
	  m_FunctionsAndWeights[0] = std::make_pair( WEIGHT, FUNCTION);
  }


    
  virtual t_RETURN operator()( const t_INPUT &DATA)
  {
	  t_RETURN result( 0);

	  // loop/iterate through the container
	  for( typename std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > >::const_iterator itr = m_FunctionsAndWeights.begin(); itr != m_FunctionsAndWeights.end(); ++itr)
	  {
		  // multiply weight with result of operator() call of each function and sum it to total result
		  result += itr->first * itr->second->operator()( DATA); // calls the operator () for each of the functions stored in m_FunctionsAndWeights
	  }

	  return result;
  }


  virtual void AddNewElement( const std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > &PAIR)
  {
    m_FunctionsAndWeights.push_back( PAIR);
  }


  virtual void AddSumFunction( const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION)
  {
    m_FunctionsAndWeights.insert( m_FunctionsAndWeights.end(), ( *SUM_FUNCTION).m_FunctionsAndWeights.begin(), ( *SUM_FUNCTION).m_FunctionsAndWeights.end());
  }


  virtual double SumOfWeights() const
  {
	  double
		  sum = 0.0;
	  for( typename std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > >::const_iterator itr = m_FunctionsAndWeights.begin(); itr != m_FunctionsAndWeights.end(); ++itr)
	  {
		  sum += itr->first;
	  }
	  return sum;
  }

  virtual void NormalizeWeights()
  {
	  double
		  inverse = 1.0 / SumOfWeights();

	  for( typename std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > >::iterator itr = m_FunctionsAndWeights.begin(); itr != m_FunctionsAndWeights.end(); ++itr)
	  {
		  itr->first *= inverse;
	  }
  }

  virtual const std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > > &
  GetData() const
  {
	  return m_FunctionsAndWeights;
  }

  virtual std::ostream &Write( std::ostream &STREAM) const
  {
	  STREAM << "SumFunction::Write()" << "\n";
	  for( typename std::vector< std::pair< double, ShPtr< Function< t_INPUT, t_RETURN> > > >::const_iterator itr = m_FunctionsAndWeights.begin(); itr != m_FunctionsAndWeights.end(); ++itr)
	  {
		  itr->second->Write( STREAM);
	  }
	  return STREAM;
  }

}; // end class SumFunction


//! operator + that takes two pointer to functions and returns a pointer to a SumFunction object
template< typename t_INPUT, typename t_RETURN>
ShPtr< SumFunction< t_INPUT, t_RETURN> > 
operator + ( const ShPtr< Function< t_INPUT, t_RETURN> > &F1, const ShPtr< Function< t_INPUT, t_RETURN> > &F2)
{
//	DebugWrite( __FUNCTION__);

	return ShPtr< SumFunction< t_INPUT, t_RETURN> >( new SumFunction< t_INPUT, t_RETURN>( F1, F2));
}


// is called by weight * a pointer to a pair of AA and a double and a pointer_counter <(X,X),double>(0)
template< typename t_INPUT, typename t_RETURN>
ShPtr< SumFunction< t_INPUT, t_RETURN> > 
operator * ( const double &WEIGHT, const ShPtr< Function< t_INPUT, t_RETURN> > &FUNCTION)
{
//	DebugWrite( __FUNCTION__);

  return ShPtr< SumFunction< t_INPUT, t_RETURN> >( new SumFunction< t_INPUT, t_RETURN>( WEIGHT, FUNCTION));
  // ptr< classX> name( new classX( value));
}



//! operator + that takes two pointer to functions and returns a pointer to a SumFunction object
template< typename t_INPUT, typename t_RETURN>
ShPtr< SumFunction< t_INPUT, t_RETURN> > 
operator + ( const ShPtr< Function< t_INPUT, t_RETURN> > &FUNCTION, const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION)
{
//  DebugWrite( __FUNCTION__);

  SUM_FUNCTION->AddNewElement( std::make_pair( 1.0, FUNCTION)); // add new element to existing SumFunction
  
  return SUM_FUNCTION;
}

//! operator + that takes two pointer to functions and returns a pointer to a SumFunction object
template< typename t_INPUT, typename t_RETURN>
ShPtr< SumFunction< t_INPUT, t_RETURN> > 
  operator + (  const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION, const ShPtr< Function< t_INPUT, t_RETURN> > &FUNCTION)
{
//  DebugWrite( __FUNCTION__);

  SUM_FUNCTION->AddNewElement( std::make_pair( 1.0, FUNCTION)); // add new element to existing SumFunction
  
  return SUM_FUNCTION;
}

//! operator + that takes two pointer to functions and returns a pointer to a SumFunction object
template< typename t_INPUT, typename t_RETURN>
ShPtr< SumFunction< t_INPUT, t_RETURN> > 
  operator + (  const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION_1, const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION_2)
{
 // DebugWrite( __FUNCTION__);

  SUM_FUNCTION_1->AddSumFunction( SUM_FUNCTION_2); // add new element to existing SumFunction
  
  return SUM_FUNCTION_1;
}

//template< typename t_INPUT, typename t_RETURN>
//ShPtr< Function< t_INPUT, t_RETURN> > 
//  operator += (  ShPtr< Function< t_INPUT, t_RETURN> > &FUNCTION, const ShPtr< SumFunction< t_INPUT, t_RETURN> > &SUM_FUNCTION)
//{
//  DebugWrite( __FUNCTION__);
//
//  ShPtr< SumFunction< t_INPUT, t_RETURN> > copy( SUM_FUNCTION);
//
//  copy->AddNewElement( std::make_pair( 1.0, FUNCTION));
//
//
//  SUM_FUNCTION_1->AddSumFunction( SUM_FUNCTION_2); // add new element to existing SumFunction
//  
//  return SUM_FUNCTION_1;
//}
//



#endif
