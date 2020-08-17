/////////////////////////////////////////////////////////////////////////
//!
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
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 1.4.2010
/////////////////////////////////////////////////////////////////////////


#ifndef TRIPLET_H_
#define TRIPLET_H_





  template< class T1, class T2, class T3>
  class Triplet
  {

  private:

  //////////
  // data //
  //////////

    T1 m_First;    //!< First
    T2 m_Second;   //!< Second
    T3 m_Third;    //!< Third

  public:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Triplet() :
      m_First(), m_Second(), m_Third()
    {};

    //! construct Triplet from three values
    Triplet( const T1 &FIRST, const T2 &SECOND, const T3 &THIRD) :
      m_First( FIRST),
      m_Second( SECOND),
      m_Third( THIRD)
    {};

    //! copy constructor
    Triplet( const Triplet &TRIPLET) :
      m_First( TRIPLET.m_First),
      m_Second( TRIPLET.m_Second),
      m_Third( TRIPLET.m_Third)
    {};

    //! virtual destructor
    virtual ~Triplet()
    {};

    //! virtual copy constructor, needed to make hard copies from pointer to it
    virtual Triplet* Clone() const
    { return new Triplet( *this);}


  /////////////////
  // data access //
  /////////////////


    //! data manipultation, first value
    virtual T1       &First()
    { return m_First;}

    //! data access, first value
    virtual T1 const &First() const
    { return m_First;}

    //! data manipultation, second value
    virtual T2       &Second()
    { return m_Second;}

    //! data access, second value
    virtual T2 const &Second() const
    { return m_Second;}

    //! data manipultation, third value
    virtual T3       &Third()
    { return m_Third;}

    //! data access, third value
    virtual T3 const &Third() const
    { return m_Third;}

//    virtual std::string GetClassName() const
//    {
////            std::cout << __FUNCTION__ << std::endl;
//        return mystr::GetClassName( __PRETTY_FUNCTION__);
//    }
  ///////////////
  // operators //
  ///////////////

    //! equal operator Triplet = TRIPLET
    Triplet< T1, T2, T3> &operator = ( const Triplet< T1, T2, T3> &TRIPLET)
    {
      m_First  = TRIPLET.First();
      m_Second = TRIPLET.Second();
      m_Third  = TRIPLET.Third();
      return *this;
    }

    //! equal operator Triplet = TRIPLET
    Triplet< T1, T2, T3> &operator += ( const Triplet< T1, T2, T3> &TRIPLET)
    {
      m_First  += TRIPLET.First();
      m_Second += TRIPLET.Second();
      m_Third  += TRIPLET.Third();
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Triple to std::ostream
    virtual std::ostream& Write( std::ostream &STREAM) const
    {
//        STREAM << GetClassName() << std::endl;
        STREAM << "Triplet( "
            <<  m_First
            << ", "
            << m_Second
            << ", "
            << m_Third
            << " )"
            << std::endl;
      return STREAM;
    }

    //! read Triple from std::ifstream
    virtual std::istream& Read( std::istream &STREAM)
    {
      //read member
      std::string identify;
      STREAM >> identify;
      STREAM >> m_First
             >> identify;
      STREAM >> m_Second
             >> identify;
      STREAM >> m_Third
             >> identify;
      return STREAM;
    }

  };  // end class Triplet

  template< class T1, class T2, class T3>
  std::ostream & operator << ( std::ostream & STREAM, const Triplet< T1, T2, T3> &TRIP)
  {
	  return TRIP.Write( STREAM);
  }

  //! this is a helper struct for ordering sets or other standard containers of triplets
  template< class T1, class T2, class T3>
  struct TripletCompare
  {
    bool operator()( const Triplet< T1, T2, T3> &TRIPLET1, const Triplet< T1, T2, T3> &TRIPLET2) const
    {
      if( TRIPLET1.First() < TRIPLET2.First())
      {
        return true;
      }
      else if( TRIPLET1.First() > TRIPLET2.First())
      {
        return false;
      }
      else if( TRIPLET1.Second() < TRIPLET2.Second())
      {
        return true;
      }
      else if( TRIPLET1.Second() > TRIPLET2.Second())
      {
        return false;
      }
      else if( TRIPLET1.Third() < TRIPLET2.Third())
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  };

  //! returns bool, true if two std::pairs are same
  template< class T1, class T2>
    inline bool operator == ( std::pair<T1, T2> &PAIR_A, std::pair<T1, T2> &PAIR_B)
  { return ( PAIR_A.first == PAIR_B.first && PAIR_A.second == PAIR_B.second);}

  //! tests wether two Triplets are the same
  template< class T1, class T2, class T3>
    inline bool operator == ( const Triplet< T1, T2, T3> &TRIPLET_A, const Triplet< T1, T2, T3> &TRIPLET_B)
  { return ( TRIPLET_A.First() == TRIPLET_B.First() && TRIPLET_A.Second() == TRIPLET_B.Second() && TRIPLET_A.Third() == TRIPLET_B.Third());}


#endif // TRIPLET_H_
