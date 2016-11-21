/*
 *  Number type selection header 	
 *
 *  Created on: Jan 30, 2015
 *  Author: Jingjin Yu
 */
#ifndef _O_CGAL_MS_RATIONAL_NT_H
#define _O_CGAL_MS_RATIONAL_NT_H

#ifndef PI
#define PI           3.14159
#endif

#ifndef BOOST_ALL_NO_LIB 
#if _MSC_VER == 1600
#define BOOST_ALL_NO_LIB
#endif 
#endif

#include <CGAL/basic.h>

#ifdef CGAL_USE_GMP

  // GMP is installed. Use the GMP rational number-type.
  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Number_type;

#else

  // GMP is not installed. Use CGAL's exact rational number-type.
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif

#endif
