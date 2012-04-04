/**
 * @file    data_access.h
 * \author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
 * \author  Math Dept, PUC-Rio
 * \author  Matmidia Labs
 * \date    14/02/2006
 *
 * @brief   data accesss
 */
//________________________________________________


#pragma once

#ifndef WIN32
#pragma interface
#endif // WIN32

#include <string>
#include <stdio.h>
#include "cube.h"


//_____________________________________________________________________________
/// data access abstract class
class data_access
//_____________________________________________________________________________
{
  public :
  Level  _max_level;  ///< maximal level of a cell
  real   _iso_val;  ///< iso value
  float  _sx ;  ///< grid x proportion
  float  _sy ;  ///< grid y proportion
  float  _sz ;  ///< grid z proportion
  
  
  public :
  /// constructor
  data_access( Level max_level_ = 4, real iso_val_ = 0.0 ) : _max_level(max_level_), _iso_val(iso_val_), _sx(1.0f), _sy(1.0f), _sz(1.0f) {}
  /// destructor
  virtual ~data_access() {}
  /// test function
  virtual bool need_refine( const Cube &c ) = 0 ;
  /// data function
  virtual real value_at   ( const Cube &c ) = 0 ;
  /// interpolation function
  virtual void interpolate( real ci, real cj, const Point &pi, const Point &pj, Point &p )
  {
    real u = ( ci ) / ( ci - cj ) ;
    p = ((1-u)*pi) + (u*pj) ;
  }
} ;
//-----------------------------------------------------------------------------



// may be overridden by the preprocessor
#ifndef RAND_THRES
# define RAND_THRES 0.2
#endif // RAND_THRES

//_____________________________________________________________________________
/// basic data access
class data_rand : public data_access
//-----------------------------------------------------------------------------
{
public:
  /// test function
  bool need_refine( const Cube &c )
  {
    int l = c.lv() ;
    if( l < 2 ) return true ;
    if( l >= _max_level ) return false ;
    return value_at(c) > RAND_THRES ;
  }
  
  /// data function
  real value_at( const Cube &c ) { return ( (real)rand() ) / ( RAND_MAX >> 1 ) - 1.0 ; }
  
  /// interpolation function
  void interpolate( real ci, real cj, const Point &pi, const Point &pj, Point &p )
  {
    real u = ( (real)rand() ) / ( RAND_MAX ) ;
    p = ((1-u)*pi) + (u*pj) ;
  }

  //-----------------------------------------------------------------------------
  // Basics
public:
  /// default constructor
  data_rand( Level max_level_ ) : data_access(max_level_) {}
  /// destructor
  virtual ~data_rand() {}
} ;
//_____________________________________________________________________________



#include "fptypes.h"
#include "fparser.h"

//_____________________________________________________________________________
/// function crossing data access
class data_func : public data_access
//-----------------------------------------------------------------------------
{
public:
  /// implicit function
  std::string _formula ; 
  
protected:
  double          _v[3]   ;
  FunctionParser  _parser ;
  
protected:
  void parse() ;
  
public:
  /// test function
  bool need_refine( const Cube &c ) ;
  
  /// data function
  real value_at( const Cube &c )
  {
    if( _parser.GetParseErrorType() != FunctionParser::FP_NO_ERROR ) return false ;
    _v[X] = c.cx() ;
    _v[Y] = c.cy() ;
    _v[Z] = c.cz() ;
    real v = _parser.Eval(_v) - _iso_val ;
    if( _parser.EvalError() ) return -1.0 ;
    return v ;
  }
  
  /// interpolation function
  void interpolate( real ci, real cj, const Point &pi, const Point &pj, Point &p ) ;  
  
  //-----------------------------------------------------------------------------
  // Basics
public:
  /// default constructor
  data_func( Level max_level_, real iso_val_, const std::string &formula_ ) : data_access(max_level_,iso_val_), _formula(formula_) { parse() ; }
  /// destructor
  virtual ~data_func() {}  
} ;
//_____________________________________________________________________________




//_____________________________________________________________________________
/// iso file data access
class data_file : public data_access
//-----------------------------------------------------------------------------
{
public:
  uint   _size_x     ;  ///< width  of the grid
  uint   _size_y     ;  ///< depth  of the grid
  uint   _size_z     ;  ///< height of the grid
  
protected:
  float *_data       ;  /**< implicit function values sampled on the grid */
  uint   _size_yz    ;  /**< k skip size */
  uint   _max_size   ;  /**< bigger side to fit into the cube */
  
  
public:
  bool  load_data  ( const char *fn ) ;
  bool  has_data   ( uint i, uint j, uint k ) { return (i < _size_x) && (j < _size_y) && (k < _size_z) ; }
  float get_data   ( uint i, uint j, uint k ) { return _data[ k + j*_size_z + i*_size_yz] - _iso_val ; }
  
  /// test function
  bool  need_refine( const Cube &c ) ;
  /// data function
  real value_at( const Cube &c )
  {
    uint i = (uint)(c.cx()*_max_size) ;
    uint j = (uint)(c.cy()*_max_size) ;
    uint k = (uint)(c.cz()*_max_size) ;
    if( has_data(i,j,k) ) return (real)get_data(i,j,k) ;
    return -1.0 ;
  }
  // linear interpolation: no specification
  
  //-----------------------------------------------------------------------------
  // Basics
public:
  /// default constructor
  data_file( Level max_level_, real iso_val_, const char *fn ) : data_access(max_level_,iso_val_), _data((float*)NULL) { load_data( fn ) ; }
  /// destructor
  virtual ~data_file() {}  
} ;
//_____________________________________________________________________________
