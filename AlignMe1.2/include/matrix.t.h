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
//!  A very simple matrix class based on a dynamically allocated
//! pointer to pointer implementation.
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>

template< typename T>
class Matrix
{
 protected:
  T ** m_Matrix;
  size_t m_NrRows;
  size_t m_NrCols;

 public:
  Matrix( const size_t &NR_ROWS = 0, const size_t &NR_COLS = 0)
    : m_Matrix( new T*[NR_ROWS]),
    m_NrRows( NR_ROWS),
    m_NrCols( NR_COLS)
  {
	  for( size_t i = 0; i < NR_ROWS; ++i)
	  {
		  m_Matrix[i] = new T[NR_COLS];
	  }
  }

    
  Matrix( const Matrix &ORIGINAL)
  : m_Matrix( new T*[ORIGINAL.m_NrRows]), // no soft copy because of destructor!!!
    m_NrRows( ORIGINAL.m_NrRows),
    m_NrCols( ORIGINAL.m_NrCols)
  {
	  for( size_t i = 0; i < m_NrRows; ++i)
	  {
		  m_Matrix[i] = new T[ m_NrCols];
		  for( size_t j = 0; j < m_NrCols; ++j)
			{
			  m_Matrix[i][j] = ORIGINAL( i, j);
			}
	  }
  }


    ~Matrix()
    {
    	for( size_t i = 0; i < m_NrRows; ++i)
    	{
			delete [] m_Matrix[i];
		}
    	delete [] m_Matrix;
    }


    const T& operator()( const size_t &I, const size_t &J) const
    {
    	return m_Matrix[I][J];
    }


    T& operator()( const size_t &I, const size_t &J)
    {
    	return m_Matrix[I][J];
    }


    const size_t &GetNumberOfColumns() const
    {
    	return m_NrCols;
    }


    const size_t &GetNumberOfRows() const
    {
    	return m_NrRows;
    }


    Matrix<T> & operator = ( const Matrix<T> &MATRIX)
    {
    	std::cout << __FUNCTION__ << " better avoid this, seems slow!" << std::endl;
    	for( size_t i = 0; i < m_NrRows; ++i)
    	{
			delete [] m_Matrix[i];
		}
    	delete [] m_Matrix;

    	m_Matrix = new T*[MATRIX.m_NrRows]; // no soft copy because of destructor!!!
    	m_NrRows = MATRIX.m_NrRows;
    	m_NrCols = MATRIX.m_NrCols;
    	for( size_t i = 0; i < m_NrRows; ++i)
    	{
    		m_Matrix[i] = new T[ m_NrCols];
    		for( size_t j = 0; j < m_NrCols; ++j)
    		{
    			m_Matrix[i][j] = MATRIX( i, j);
    		}
    	}
    	return *this;
    }


    void Set( const size_t &I, const size_t &J)
    {
    	if( m_NrRows > 0 && m_NrCols > 0)
    	{
        	for( size_t i = 0; i < m_NrRows; ++i)
        	{
    			delete [] m_Matrix[i];
    		}
        	delete [] m_Matrix;
    	}

    	m_NrRows = I;
    	m_NrCols = J;
    	m_Matrix = new T*[m_NrRows];
       	for( size_t i = 0; i < m_NrRows; ++i)
       	{
       		m_Matrix[i] = new T[m_NrCols];
       	}
    }

    std::ostream &Write( std::ostream &STREAM )const
    {
  	  for( size_t i = 0; i < m_NrRows; ++i)
  	  {
   		  for( size_t j = 0; j < m_NrCols; ++j)
  			{
  			  STREAM << m_Matrix[i][j] << " ";
  			}
   		  STREAM << std::endl;
  	  }
  	  return STREAM;
    }

}; // end class Matrix //



template< typename T>
inline
Matrix< T>  operator *( const T &VALUE, const Matrix<T> &MATRIX)
{
	Matrix<T> matrix( MATRIX.GetNumberOfRows(), MATRIX.GetNumberOfColumns());
	for( int i = 0; i < MATRIX.GetNumberOfRows(); ++i)
		for( int j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
		{
			matrix(i,j) = VALUE * MATRIX( i, j);
		}
	return matrix;
}



#endif
