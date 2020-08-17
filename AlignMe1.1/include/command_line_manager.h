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
//!  The simple CommandLineManager class reads all arguments passed to the
//! program and identifies the flags and their arguments.
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////

#ifndef COMMANDLINE_MANAGER_H
#define COMMANDLINE_MANAGER_H

#include <string>
#include <map>
#include <vector>
#include <cassert>
#include "string_functions.h"

class CommandLineManager
{
 private:
	std::map< std::string, std::vector< std::string> >  m_FlagsAndValues;

 public:
	CommandLineManager( const size_t ARGC, const char *ARGV[], std::map<std::string,std::string> ALLOWED_FLAGS)
	: m_FlagsAndValues()
	  {
		ReadCommandLine( ARGC, ARGV, ALLOWED_FLAGS);
	  }

	CommandLineManager( const size_t ARGC, const char *ARGV[])
	: m_FlagsAndValues()
	  {
		std::map<std::string,std::string>  ALLOWED_FLAGS;
		ReadCommandLine( ARGC, ARGV, ALLOWED_FLAGS);
	  }


	~CommandLineManager(){}

	void ReadCommandLine( const size_t ARGC, const char *ARGV[],std::map<std::string,std::string>  ALLOWED_FLAGS)
	{
		std::vector< std::string>
			values;

		std::string
			option,
			flag;

		if( ARGC == 1)
		{
			flag = "-help";
		}
		else
		{
			flag = ARGV[1];
		}
	
		// program aborts if no - is set... so the first parameter is always a flag
#ifdef SECURE
		if(  flag.substr( 0, 1) != "-")
		{
			std::cerr << "ERROR: You have to provide a flag starting with - after the executable file" << "\n";
			exit(-1);
		}
#endif
		bool
			is_numerical,
			test_for_allowed_flags;

		if (ALLOWED_FLAGS.size() > 0 )
		{
			test_for_allowed_flags = true;
		}

		//deletes element to which flag.begin is pointing to. The returned iterator then points to the following element if one exists otherwise to the end-iterator
		// "-" will be removed from the flagname
		flag.erase( flag.begin());

		if( ARGC > 2)
		{
			for( size_t i = 2; i < ARGC; ++i)
			{
		    	option = std::string( ARGV[i]);
		
		    	is_numerical = mystr::IsNumerical( option);

		    	// if no flag is set, then the commandline parameter will be added at the end of a vector called values
				if( option.substr( 0, 1) != "-" || is_numerical)
				{
					values.push_back( option);
				}

				// this loop will be entered if the loop reaches the next "-" in the command line
				// searches through all already found flags and checks if the flag is now or already inserted
				// then it makes a pair of the last found flag and all parameters that followed this flag (stored in the vector called values) until the next "-"
				// to flag the new name which appears after "-" will be assigned and deletes the "-" from the flagname
				if( (option.substr( 0, 1) == "-" && !is_numerical) || i == ARGC - 1)
				{
#ifdef SECURE
					if (test_for_allowed_flags && ALLOWED_FLAGS.find( flag) == ALLOWED_FLAGS.end())
					{
						std::cerr << "ERROR: Flag <" << flag << "> is not provided by this program!" << "\n";
						exit( -1);
					}
#endif
					std::map< std::string, std::vector< std::string> >::const_iterator itr = m_FlagsAndValues.find( flag);

#ifdef SECURE
					if( itr != m_FlagsAndValues.end())
					{
						std::cerr << "ERROR: Flag <" << flag << "> was already given - please use each flag only once!" << "\n";
						exit( -1);
					}
#endif
					m_FlagsAndValues.insert( make_pair( flag, values));
					flag = option;
					flag.erase( flag.begin());
					values.clear();
				}

				if( option.substr( 0, 1) == "-" && i == ARGC - 1 && !is_numerical)
				{
#ifdef SECURE
					if (test_for_allowed_flags && ALLOWED_FLAGS.find( flag) == ALLOWED_FLAGS.end())
					{
						std::cerr << "ERROR: Flag <" << flag << "> is not provided by this program!" << "\n";
						exit( -1);
					}
#endif
					m_FlagsAndValues.insert( make_pair( flag, values));
				}
			}
		}

		// this loop will only be entered if only one flag without any values is set
		// in such a case the program will crash ----> sense of this loop????
		else
		{
#ifdef SECURE
			if (test_for_allowed_flags && ALLOWED_FLAGS.find( flag) == ALLOWED_FLAGS.end())
			{
				std::cerr << "ERROR: Flag <" << flag << "> is not provided by this program2!" << "\n";
				exit( -1);
			}
#endif
			m_FlagsAndValues.insert( make_pair( flag, values));
		}
	}

	// checks the list of flags until its end if a certain flag has been inserted
	bool IsFlagSet( const std::string &FLAG) const
	{
		return  m_FlagsAndValues.find( FLAG) != m_FlagsAndValues.end();
	}

	// finds a flag and returns all the values that belong to this flag
	std::vector< std::string> GetArguments( const std::string &FLAG) const
	{
		return m_FlagsAndValues.find( FLAG)->second;
	}



	// finds a flag and returns all the values that belong to this flag
	template< typename TYPE>
	std::vector< TYPE> GetConvertedArguments( const std::string &FLAG) const
	{
		std::vector< std::string> strings = m_FlagsAndValues.find( FLAG)->second;
		std::vector< TYPE> values( strings.size());
		std::vector< std::string>::const_iterator str_itr = strings.begin();
		typename std::vector< TYPE>::iterator val_itr = values.begin();

		size_t count = 1;
		for( ; str_itr != strings.end(); ++str_itr, ++val_itr, ++count)
		{
#ifdef SECURE
			if( !mystr::IsNumerical( *str_itr))
			{
				std::cout << "ERROR: The " << count << "-th value after the flag: "  << FLAG << " has to be a numerical value but is " << *str_itr << " \n";
				exit(-1);
			}
#endif
			*val_itr = mystr::ConvertStringToNumericalValue< TYPE>( *str_itr);
		}
#ifdef SECURE
		if (values.size() == 0)
		{
			std::cout << "ERROR: The flag " << FLAG << " is provided but a corresponding value or filename after it is missing. Flags can not be used without corresponding inputs! \n";
			exit(-1);
		}
#endif
		return values;
	}



	// finds a flag and returns the first added value belonging to this flag (which stands directly after the flag)
	std::string GetFirstArgument( const std::string &FLAG) const
	{
		std::vector< std::string>
			strs = m_FlagsAndValues.find( FLAG)->second;
#ifdef SECURE

		if( strs.size() == 0)
		{
			std::cout << "ERROR: The flag " << FLAG << " is provided but a corresponding value or filename after it is missing. Flags can not be used without corresponding inputs! \n";
			exit(-1);
		}
#endif
		return strs.front();
	}
	

	std::ostream &Write( std::ostream &STREAM) const
	{
		for( std::map< std::string, std::vector< std::string> >::const_iterator itr = m_FlagsAndValues.begin(); itr != m_FlagsAndValues.end(); ++itr)
		{
			STREAM << itr->first << ":  <";
			for( std::vector< std::string>::const_iterator str_itr = itr->second.begin(); str_itr != itr->second.end(); ++str_itr)
			{
				STREAM << *str_itr;
				if( str_itr + 1 != itr->second.end())
				{
					STREAM << ",";
				}
			}
			STREAM << ">" << "\n";
		}
		return STREAM;
	}

};
#endif
