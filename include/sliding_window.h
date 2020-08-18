///////////////////////////////////////////////////////////////////
//!  A collection of functions for calculating different sliding
//!  window algorithms for data smoothing.
//!
//!  @date Aug 26, 2009
//!  @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
///////////////////////////////////////////////////////////////////

#ifndef SLIDING_WINDOW_H_
#define SLIDING_WINDOW_H_

#include <cmath>

inline
void RectangularSlidingWindow( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum;
  size_t
    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
    i( 0),
    number_inwindow_counter( 0);
  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif


  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;

	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
	  }

	  if( i < SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  number_inwindow_counter = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter)
	  {
		  sum += SCALE_MAP.find( win_itr->GetType())->second;  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
	  }
	  sum /= double (number_inwindow_counter);
	  // assigns to a position in the sequence a value for each command-line
	  // so if you have two scales there would be a list of two values
      seq_itr->AddNewProfile( sum);
    }
}




inline
void TriangularSlidingWindow( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum,
    multiplier,
    integral;

  int
	  // division between integers
	  counter_in_window( 0),
	  half_window( int( double( WINDOW_SIZE) / 2.0)),
	  i( 0);
#ifdef SECURE
	  if( WINDOW_SIZE < 2)
	  {
		  std::cerr << "ERROR: window size for triangular sliding window needs to be at least of size 2" << "\n";
		  exit( -1);
	  }
#endif

  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif

// this for-loop goes through all positions of a sequence and computes a value for every position
// according to the window (here triangular) and window_size
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;

	  // check if the window starts within sequence or if a part of the left half of the window
	  // is not in the sequence
	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
		  counter_in_window = i - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
		  counter_in_window = 0;
	  }
	  // check if the window ends within sequence or if a part of the right half of the window
	  // is not in the sequence

	  if( i < (signed) SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  DebugWrite( " i: " << i);

	  // the calculated sum will be divided through the integral
	  // the integral is the sum of all multipliers that are involved in scoring
	  // and is necessary for normalization of the score, so that one can compare
	  // scores of different windowsizes or shapes
	  integral = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++counter_in_window)
	  {
		  DebugWriteNoFlush(  win_itr->GetType());

		  multiplier =  1.0 - fabs( double (i - counter_in_window) / double (half_window));
		  DebugWrite( "multiplier:  " << multiplier);
		  integral += multiplier;
		  // because SCALE_MAP is passed as const to function,
		  // map<>::find( X) has to be used that returns a pointer to a key/value pair;
		  // '->second' returns the value
		  sum += multiplier * ( SCALE_MAP.find( win_itr->GetType())->second);
	  }
	  ///here put that thing later
	  sum /= integral;
      DebugWrite( "sum:  " << sum);
      DebugWrite( "\n" << "\n");
      seq_itr->AddNewProfile( sum);
    }
}


//!  for sliding window flag 'triangular_msa'
inline
void TriangularSlidingWindowIgnoringGaps( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{

  double
	  sum,		//sum of hydrophobicities over all positions, normalized by integral
	  multiplier, //multiplier of the triangular window
	  integral( 0.0); //to divide the sum and to get a hydrophobicity value after window-averaging

  size_t
	  i( 0), //iterator
	  j, //iterator within a window
	  window_counter,
	  window_start, //windows will start from these indices
	  window_end, //window will end at these indices
	  first_counter,
	  second_counter,
	  half_window( int( double( WINDOW_SIZE) / 2.0)); //this is obvious

  int
	  tmp2;

  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    start,
    end;

#ifdef SECURE
	  if( WINDOW_SIZE < 2)
	  {
		  std::cerr << "ERROR: window size for triangular sliding window needs to be at least of size 2" << "\n";
		  exit( -1);
	  }
#endif

  for (std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;

      if( seq_itr->GetType() != '-')
      {
    	  // determine window-start
          first_counter = 1;
          second_counter = 1;
          window_start = i;

          while( first_counter <= half_window && (int (i) - int( second_counter)) >= 0)
          {
        	  if( ( seq_itr - second_counter)->GetType() != '-')
              {
                  ++first_counter;
              }
              ++second_counter;
              --window_start;
          }

          // determine window end
          first_counter = 1;
          second_counter = 1;
          window_end = i;
          while( ( first_counter <= half_window) and (( i + second_counter) < SEQUENCE.size()))
          {
              if( ( seq_itr + second_counter)->GetType() != '-')
              {
                  ++first_counter;
              }
              ++second_counter;
              ++window_end;
          }

          // sliding window calculation (sum of values and integral of sliding function)
          window_counter = 0;
          integral = 0;
          for( j = window_start; j <= window_end; j++)  //goes along "non-gaps" inside of the window (which may contain gaps)
          {
        	  if( (SEQUENCE.begin() + j) -> GetType()  != '-')
        	  {   // if not a gap
        		  ++window_counter;
                  tmp2 = ( half_window + 1) - window_counter; //walker (or position) with respect to the window size
                  if( tmp2 < 0)
                  {
                	  tmp2  = -1 * tmp2; // for negative values
                  }
                  multiplier = ( double( half_window) - double( tmp2)) / double( half_window); //triangular value
                  integral += multiplier; // a sum of multipliers to normalize the sum later
                  sum = sum + multiplier*( SCALE_MAP.find( (SEQUENCE.begin() + j) ->GetType())->second); // win_itr j
              }
          }

          sum = sum/ integral;
          seq_itr->AddNewProfile( sum);
      }
  }
}


inline
void ZigZagSlidingWindow( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
	 DebugWrite( __FUNCTION__);
	  double
	    sum;

	  size_t
	    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
	    i( 0),
	    number_inwindow_counter( 0);
	  std::vector< GeneralizedAminoAcid>::iterator
	    win_itr,
	    start,
	    end;

	  DebugWrite( "half_window: " << half_window);

	#ifdef DEBUG
	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
	  {
		  std::cout << "<" << seq_itr->GetType() << "> ";
	  }
	  std::cout << "\n" << "\n";
	#endif


	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
	  {
		  sum = 0.0;

		  if( i >= half_window)
		  {
			  DebugWrite( "window starts within sequence");
			  start = seq_itr - half_window;
		  }
		  else
		  {
			  DebugWrite( "sequence window from beginning");
			  start = SEQUENCE.begin();
//			  std::cerr << half_window - i  << "\n";
//			  std::cerr << (half_window - i) % 2  << "\n";
			  if (((half_window - i) % 2) == 1)
				  start++;
		  }

		  if( i < SEQUENCE.size() - half_window - 1)
		  {
			  DebugWrite( "window ends within sequence");
			  end = seq_itr + half_window + 2;
			  DebugWrite( *end << " the end");
		  }
		  else
		  {
			  DebugWrite( "sequence window till end");
			  end = SEQUENCE.end();
		  }

		  number_inwindow_counter = 0;

		  for( win_itr = start; win_itr != end && win_itr <= SEQUENCE.end(); ++number_inwindow_counter)
		  {
			  sum += SCALE_MAP.find( win_itr->GetType())->second;  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
			  win_itr++;
			  win_itr++;
		  }
		  sum /= double (number_inwindow_counter);
		  // assigns to a position in the sequence a value for each command-line
		  // so if you have two scales there would be a list of two values
	      seq_itr->AddNewProfile( sum);
	    }
}

inline
void SinoidSlidingWindow(Sequence &SEQUENCE,const std::map<char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{

	DebugWrite(__FUNCTION__);
	double
		sum,
		multiplier,
		integral;
	size_t
		half_window(size_t(double(WINDOW_SIZE) / 2.0)), // devision between integers
		i(0),
		shift (0),
		number_inwindow_counter(0);
	std::vector<GeneralizedAminoAcid>::iterator
		win_itr,
		start,
		end;
  const double
	  pi( acos( -1.0));

	DebugWrite("half_window: " << half_window);

//	std::ofstream out( "gnuplot.txt");

#ifdef DEBUG
	for( std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
	{
		std::cout << "<" << seq_itr->GetType() << "> ";
	}
	std::cout << "\n" << "\n";
#endif

	for (std::vector<GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr!= SEQUENCE.end(); ++seq_itr, ++i)
	{
		sum = 0.0;
		shift = 0;

		if (i >= half_window)
		{
			DebugWrite("window starts within sequence");
			start = seq_itr - half_window;
		}
		else
		{
			DebugWrite("sequence window from beginning");
			start = SEQUENCE.begin();
			shift = half_window - i;
		}

		if (i < SEQUENCE.size() - half_window - 1)
		{
			DebugWrite("window ends within sequence");
			end = seq_itr + half_window + 1;
			DebugWrite(*end << " the end");
		}
		else {
			DebugWrite("sequence window till end");
			end = SEQUENCE.end();
			shift = (half_window - (SEQUENCE.size()-1 - i));
		}

		number_inwindow_counter = 0;
		integral = 0;
		for (win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter, ++shift)
		{
			multiplier = 0.5 * ( cos (2.0 / 3.6 * pi * fabs( double( shift) - double( half_window))) + 1);
			integral += multiplier;
			sum += multiplier * ( SCALE_MAP.find( win_itr->GetType())->second);
		}
		sum /= integral;

		// assigns to a position in the sequence a value for each command-line
		// so if you have two scales there would be a list of two values
		seq_itr->AddNewProfile(sum);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////// windows to detect amphiphaticity patterns within the sequence //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

inline
void AmphiphilicityPure( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum,
    sum2;
  size_t
    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
    i( 0),
    number_inwindow_counter( 0);
  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    win_itr2,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif


  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;
	  sum2 = 0.0;

	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
	  }

	  if( i < SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  number_inwindow_counter = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter)
	  {
		  win_itr2 = win_itr;
		  win_itr2++;
		  if (win_itr2 != end)
		  {
			  sum += fabs(SCALE_MAP.find( win_itr->GetType())->second - SCALE_MAP.find( win_itr2->GetType())->second);  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
		  }


	  }
	  sum /= double (number_inwindow_counter);
	 // assigns to a position in the sequence a value for each command-line
	  // so if you have two scales there would be a list of two values
      seq_itr->AddNewProfile( sum);
    }
}


inline
void AmphiphilicityPureRectangular( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum,
    sum2;
  size_t
    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
    i( 0),
    number_inwindow_counter( 0);
  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    win_itr2,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif


  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;
	  sum2 = 0.0;

	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
	  }

	  if( i < SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  number_inwindow_counter = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter)
	  {
		  win_itr2 = win_itr;
		  win_itr2++;
		  if (win_itr2 != end)
		  {
			  sum += fabs(SCALE_MAP.find( win_itr->GetType())->second - SCALE_MAP.find( win_itr2->GetType())->second);  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
		  }

		  sum2 += SCALE_MAP.find( win_itr->GetType())->second;
		  if ( win_itr ==seq_itr)
		  {
			  sum2 -= SCALE_MAP.find( win_itr->GetType())->second;
		  }

	  }
	  sum /= double (number_inwindow_counter);
	  sum2 /= double (number_inwindow_counter);

	  sum =  sum - fabs(sum2);
	 // assigns to a position in the sequence a value for each command-line
	  // so if you have two scales there would be a list of two values
      seq_itr->AddNewProfile( sum);
    }
}



inline
void AmphiphilicitySquared( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum,
    sum2;
  size_t
    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
    i( 0),
    number_inwindow_counter( 0);
  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    win_itr2,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif


  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;
	  sum2 = 0.0;

	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
	  }

	  if( i < SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  number_inwindow_counter = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter)
	  {
		  win_itr2 = win_itr;
		  win_itr2++;
		  if (win_itr2 != end)
		  {
			  sum += pow(fabs(SCALE_MAP.find( win_itr->GetType())->second - SCALE_MAP.find( win_itr2->GetType())->second), 2);  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
		  }
	  }
	  sum /= double (number_inwindow_counter);
	  sum  =  sqrt(sum);
	 // assigns to a position in the sequence a value for each command-line
	  // so if you have two scales there would be a list of two values
      seq_itr->AddNewProfile( sum);
    }
}



inline
void AmphiphilicitySquaredRectangular( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE)
{
  DebugWrite( __FUNCTION__);
  double
    sum,
    sum2;
  size_t
    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
    i( 0),
    number_inwindow_counter( 0);
  std::vector< GeneralizedAminoAcid>::iterator
    win_itr,
    win_itr2,
    start,
    end;

  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
  {
	  std::cout << "<" << seq_itr->GetType() << "> ";
  }
  std::cout << "\n" << "\n";
#endif


  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
  {
	  sum = 0.0;
	  sum2 = 0.0;

	  if( i >= half_window)
	  {
		  DebugWrite( "window starts within sequence");
		  start = seq_itr - half_window;
	  }
	  else
	  {
		  DebugWrite( "sequence window from beginning");
		  start = SEQUENCE.begin();
	  }

	  if( i < SEQUENCE.size() - half_window - 1)
	  {
		  DebugWrite( "window ends within sequence");
		  end = seq_itr + half_window + 1;
		  DebugWrite( *end << " the end");
	  }
	  else
	  {
		  DebugWrite( "sequence window till end");
		  end = SEQUENCE.end();
	  }

	  number_inwindow_counter = 0;
	  for( win_itr = start; win_itr != end; ++win_itr, ++number_inwindow_counter)
	  {
		  win_itr2 = win_itr;
		  win_itr2++;
		  if (win_itr2 != end)
		  {
			  sum += pow(fabs(SCALE_MAP.find( win_itr->GetType())->second - SCALE_MAP.find( win_itr2->GetType())->second), 2);  // because SCALE_MAP is passed as const to function, map<>::find( X) has to be used that returns a pointer to a key/value pair; '->second' returns the value
		  }

		  sum2 += SCALE_MAP.find( win_itr->GetType())->second;
		  if ( win_itr == seq_itr)
		  {
			  sum2 -= SCALE_MAP.find( win_itr->GetType())->second;
		  }

	  }
	  sum /= double (number_inwindow_counter);
	  sum  =  sqrt(sum);

	  sum2 /= double (number_inwindow_counter);

	  sum =  sum - fabs(sum2);
	 // assigns to a position in the sequence a value for each command-line
	  // so if you have two scales there would be a list of two values
      seq_itr->AddNewProfile( sum);
    }
}


inline
void AmphiphaticSlidingWindow( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE, const double &DELTA, const double &THRESHOLD)
{
	 DebugWrite( __FUNCTION__);
	  double
		  first = 0,
		  second = 0,
		  sum;

	  std::vector< std::vector< double> > values;

	  size_t
	    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
	    i( 0),
	    number_inwindow_counter( 0);
	  std::vector< GeneralizedAminoAcid>::iterator
	    win_itr,
	    start,
	    end;

	  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
	  {
		  std::cout << "<" << seq_itr->GetType() << "> ";
	  }
	  std::cout << "\n" << "\n";
#endif


	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
	  {
		  sum = 0.0;
		  values = std::vector< std::vector< double> >( 2, std::vector< double>());

		  if( i >= half_window)
		  {
			  DebugWrite( "window starts within sequence");
			  start = seq_itr - half_window;
		  }
		  else
		  {
			  DebugWrite( "sequence window from beginning");
			  start = SEQUENCE.begin();
//			  std::cerr << half_window - i  << "\n";
//			  std::cerr << (half_window - i) % 2  << "\n";
			  if (((half_window - i) % 2) == 1)
				  start++;
		  }

		  if( i < SEQUENCE.size() - half_window - 1)
		  {
			  DebugWrite( "window ends within sequence");
			  end = seq_itr + half_window + 2;
			  DebugWrite( *end << " the end");
		  }
		  else
		  {
			  DebugWrite( "sequence window till end");
			  end = SEQUENCE.end();
		  }


		  for( win_itr = start; win_itr != end && win_itr <= SEQUENCE.end(); ++number_inwindow_counter, ++win_itr)
		  {
			  std::cout << win_itr - SEQUENCE.begin() << " " << (win_itr - SEQUENCE.begin()) % 2 << " " << win_itr->GetType() << " " << SCALE_MAP.find( win_itr->GetType())->second << std::endl;
			  values[ (win_itr - SEQUENCE.begin()) % 2].push_back( SCALE_MAP.find( win_itr->GetType())->second);
		  }

		  for( int j = 0; j < (int) values[0].size(); ++j)
		  {
			  first += values[0][j];
		  }
		  first /= values[0].size();

		  for( int j = 0; j < (signed) values[1].size(); ++j)
		  {
			  second += values[1][j];
		  }
		  second /= values[1].size();

		  std::cout << first << " " << second << " " << std::abs( first - second) / DELTA;
		  if( std::abs( first - second) / DELTA > THRESHOLD)
		  {
			  std::cout << " => 1.0 " << std::endl;
			  sum = 1.0;
		  }
		  else
		  {
			  std::cout << " => 0.0 " << std::endl;
			  sum = 0.0;
		  }
		  std::cout << std::endl;
		  // assigns to a position in the sequence a value for each command-line
		  // so if you have two scales there would be a list of two values
	      seq_itr->AddNewProfile( std::abs( first - second) / DELTA);
	      seq_itr->AddNewProfile( sum);
	    }
}




inline
void AmphiphaticNextNeighborSlidingWindow( Sequence &SEQUENCE, const std::map< char, double> &SCALE_MAP, const size_t &WINDOW_SIZE, const double &DELTA, const double &THRESHOLD)
{
	 DebugWrite( __FUNCTION__);
	  double
		  last,
		  value;

//	  double
//		  sum;

	  std::vector< std::vector< double> > values;

	  size_t
	    i( 0);

//	  size_t
//	    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
//	    number_inwindow_counter( 0);
	  std::vector< GeneralizedAminoAcid>::iterator
	    win_itr,
	    start,
	    end;

	  DebugWrite( "half_window: " << half_window);

#ifdef DEBUG
	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
	  {
		  std::cout << "<" << seq_itr->GetType() << "> ";
	  }
	  std::cout << "\n" << "\n";
#endif


//	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
//	  {
//		  sum = 0.0;
//		  number_inwindow_counter = 0;
//
//		  if( i >= half_window)
//		  {
//			  DebugWrite( "window starts within sequence");
//			  start = seq_itr - half_window;
//		  }
//		  else
//		  {
//			  DebugWrite( "sequence window from beginning");
//			  start = SEQUENCE.begin();
////			  std::cerr << half_window - i  << "\n";
////			  std::cerr << (half_window - i) % 2  << "\n";
//			  if (((half_window - i) % 2) == 1)
//				  start++;
//		  }
//
//		  if( i < SEQUENCE.size() - half_window - 1)
//		  {
//			  DebugWrite( "window ends within sequence");
//			  end = seq_itr + half_window + 2;
//			  DebugWrite( *end << " the end");
//		  }
//		  else
//		  {
//			  DebugWrite( "sequence window till end");
//			  end = SEQUENCE.end();
//		  }
//
//
//		  last = SCALE_MAP.find( start->GetType())->second;
//		  float vel = 0.0, acc = 0.0, prev = 0, diff;
//		  for( win_itr = start + 1; win_itr != end && win_itr <= SEQUENCE.end(); ++number_inwindow_counter, ++win_itr)
//		  {
//			  std::cout << win_itr - SEQUENCE.begin() << " " << (win_itr - SEQUENCE.begin()) % 2 << " " << win_itr->GetType() << " " << SCALE_MAP.find( win_itr->GetType())->second << std::endl;
//
//			  value = SCALE_MAP.find( win_itr->GetType())->second;
//			  diff = value - last;
//			  last = value;
//
//			  if( std::abs( diff) / DELTA > THRESHOLD)
//			  {
//				  sum += 1.0;
//			  }
//
//			  vel += std::abs( diff);
//
//			  acc += std::abs( prev - diff);
//			  prev = diff;
//
//		  }
//		  sum /= (double) number_inwindow_counter;
//		  acc /= (double) number_inwindow_counter;
//		  vel /= (double) number_inwindow_counter;
//
//		  std::cout << std::endl;
//		  // assigns to a position in the sequence a value for each command-line
//		  // so if you have two scales there would be a list of two values
//	      seq_itr->AddNewProfile( vel);
//	      seq_itr->AddNewProfile( acc);
//	      seq_itr->AddNewProfile( sum);
//	    }

//	  float vel = 0.0;

	  float acc = 0.0, prev = 0, diff;

	  last = 0.0;

	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
	  {
		  value = SCALE_MAP.find( seq_itr->GetType())->second;
		  diff = value - last;
		  last = value;

		  acc = diff - prev;
		  prev = diff;

	      seq_itr->AddNewProfile( value);
	      seq_itr->AddNewProfile( diff);
	      seq_itr->AddNewProfile( acc);
	  }

}

inline
void ReadBetaStrandIntoSequence
(
		Sequence &SEQUENCE,
		const std::map< char, double> &PORE,
		const std::map< char, double> &MEMBRANE,
		const size_t &WINDOW_SIZE
)
{
	 DebugWrite( __FUNCTION__);
	  double
	    poremem_score,
	    mempore_score;
	  size_t
	    half_window( size_t( double( WINDOW_SIZE) / 2.0)),  // devision between integers
	    i( 0),
	    number_inwindow_counter( 0);
	  std::vector< GeneralizedAminoAcid>::iterator
	    win_itr,
	    win_itr2,
	    start,
	    end;

	  DebugWrite( "half_window: " << half_window);

	#ifdef DEBUG
	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr)
	  {
		  std::cout << "<" << seq_itr->GetType() << "> ";
	  }
	  std::cout << "\n" << "\n";
	#endif


	  for(  std::vector< GeneralizedAminoAcid>::iterator seq_itr = SEQUENCE.begin(); seq_itr != SEQUENCE.end(); ++seq_itr, ++i)
	  {
		  poremem_score = 0.0;
		  mempore_score = 0.0;

		  if( i >= half_window)
		  {
			  DebugWrite( "window starts within sequence");
			  start = seq_itr - half_window;
		  }
		  else
		  {
			  DebugWrite( "sequence window from beginning");
			  start = SEQUENCE.begin();
		  }

		  if( i < SEQUENCE.size() - half_window - 1)
		  {
			  DebugWrite( "window ends within sequence");
			  end = seq_itr + half_window + 1;
			  DebugWrite( *end << " the end");
		  }
		  else
		  {
			  DebugWrite( "sequence window till end");
			  end = SEQUENCE.end();
		  }

		  number_inwindow_counter = 0;

		  win_itr = start;
		  while (win_itr != end)
		  {
			poremem_score += PORE.find( win_itr->GetType())->second;
			mempore_score += MEMBRANE.find( win_itr->GetType())->second;
			++number_inwindow_counter;
			win_itr++;
			if (win_itr != end)
			{
				win_itr++;
				win_itr2 = win_itr;
				if (win_itr2 != end)
				{
					++number_inwindow_counter;
					poremem_score += MEMBRANE.find( win_itr2->GetType())->second;
					mempore_score += PORE.find( win_itr2->GetType())->second;
				}
			}
		  }

		  poremem_score /= double (number_inwindow_counter);
		  mempore_score /= double (number_inwindow_counter);
		  // assigns to a position in the sequence a value for each command-line
		  // so if you have two scales there would be a list of two values
	      seq_itr->AddNewProfile( std::max(poremem_score,mempore_score));
	    }

}



#endif /* SLIDING_WINDOW_H_ */
