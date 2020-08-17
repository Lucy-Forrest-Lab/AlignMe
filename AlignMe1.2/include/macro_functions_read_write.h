////////////////////////////////////////////////
//!  macro_functions_read_write.h
//!  read/write functions that depend on a
//!  message level, based on macro definitions.
//!
//!  @date May 13, 2009
//!  @author: Rene Staritzbichler
//!  @example:
////////////////////////////////////////////////

#ifndef MACRO_FUNCTIONS_READ_WRITE_H_
#define MACRO_FUNCTIONS_READ_WRITE_H_

// these macro definitions allow to control output according to DEBUG or STANDARD message level
// this implementation allows to
// switching to another message level requires recompilation
// ===>  DO NOT write stuff like:  DebugWrite( arg++)  it can cause a double increment!  <===

#ifdef DEBUG
#define DebugWriteNoFlush( ARGUMENT) std::cout << ARGUMENT << " "
#define DebugWrite( ARGUMENT) std::cout << ARGUMENT << "\n"
#else
#define DebugWriteNoFlush( ARGUMENT)
#define DebugWrite( ARGUMENT)
#endif

#endif /* MACRO_FUNCTIONS_READ_WRITE_H_ */
