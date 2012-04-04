/**
 * \file    point.h
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 * \brief  (Point/Vector Class).
 */
/**------------------------------------------------------------------------------------*/


#pragma once

#include <stdlib.h>
#include <float.h>   // FLT_EPSILON
#include <math.h>   // nan
#ifdef UNIX
#include <GL/gl.h>
#else
#include <OpenGL/GL.h>
#endif



//------------------------------------------------
/** real number type */
typedef float           real ;
typedef const real     creal ;
#define R_PI      ((creal)3.1415926535897932384626433832795)
#define R_EPSILON ((real)FLT_EPSILON)
inline bool is_inv( real x ) { return isnan(x) ; }
extern real R_INV ;

/** integer number aliases */
typedef   signed char  schar ;
typedef unsigned char  uchar ;
typedef unsigned  int   uint ;
typedef const     int   cint ;
typedef const    uint  cuint ;

/** \brief Axis Definition */
enum Axis { X = 0, Y = 1, Z = 2 } ;
//------------------------------------------------


/** \brief point with normal */
class Point
  {
    protected :
    real  _x, _y, _z ;
    
    public :
    Point() : _x(R_INV), _y(R_INV), _z(R_INV) {}
    Point( creal x_, creal y_, creal z_ ) : _x(x_), _y(y_), _z(z_) {}
    Point( const Point &p_ ) : _x(p_._x), _y(p_._y), _z(p_._z) {}
    Point  &operator = ( const Point &p_ )
    { _x=p_._x ;  _y=p_._y ;  _z=p_._z ;  return *this; }
    
    
    public :
    inline real &x()       { return _x ; }
    inline real &y()       { return _y ; }
    inline real &z()       { return _z ; }
    inline creal x() const { return _x ; }
    inline creal y() const { return _y ; }
    inline creal z() const { return _z ; }
    
    inline creal &operator[] ( Axis a ) const
    { switch( a ) {
      case X : return _x ;
      case Y : return _y ;
      case Z : return _z ;
      default: return R_INV ;
    } }
    inline real  &operator[] ( Axis a )
    { switch( a ) {
      case X : return _x ;
      case Y : return _y ;
      case Z : return _z ;
      default: return R_INV ;
    } }    
    
    /**------------------------------------------------------------------------------------*/
    inline bool valid()
    {
      return
      ( !is_inv( x() ) ) &&
      ( !is_inv( y() ) ) &&
      ( !is_inv( z() ) ) ;
    }
    inline bool invalid()
    {
      return
      ( is_inv( x() ) ) ||
      ( is_inv( y() ) ) ||
      ( is_inv( z() ) ) ;
    }
    
    
    /**------------------------------------------------------------------------------------*/
    // Basic operations
    
    /// +=
    inline Point &operator+= ( const Point &p_ )
    {
      x() += p_.x() ;  y() += p_.y() ;  z() += p_.z() ;
      return *this ;
    }
    /// unary -
    inline const Point operator- ()
    {
      x() = -x() ;  y() = -y() ;  z() = -z() ;
      return *this ;
    }
    /// -=
    inline Point &operator-= ( const Point &p_ )
    {
      x() -= p_.x() ;  y() -= p_.y() ;  z() -= p_.z() ;
      return *this ;
    }
    /// *= scalar
    inline Point &operator*= ( const real l )
    {
      x() *= l ;  y() *= l ;  z() *= l ;
      return *this ;
    }
    /// /= scalar
    inline Point &operator/= ( const real l )
    {
      if( fabs(l) < R_EPSILON ) return *this ;
      real s = 1.0/l ;
      return *this *= s ;
    }

    
    /**------------------------------------------------------------------------------------*/
    // Norms
    inline real length () const
    {
      return hypot( hypot( x(), y() ), z() ) ;
    }
    inline real norm   () const { return length() ; }
    inline real dist   ( const Point &p_ ) const { Point r = *this ; r -= p_ ;  return r.length() ; }

    inline bool operator== ( const Point &p_ ) const { return dist( p_ ) < R_EPSILON ; }
    inline bool operator!= ( const Point &p_ ) const { return !( *this == p_ ) ; }
    
    inline bool normalize()
    {
      real l = length() ;
      if( l < R_EPSILON ) return false ;
      *this /= l ;
      return true ;
    }
    
    /** replaces p by the by-coordinate minimum of p and p_ */
    inline void pmin( const Point &p_ )
    { 
      if( x() > p_.x() ) x() = p_.x() ; 
      if( y() > p_.y() ) y() = p_.y() ; 
      if( z() > p_.z() ) z() = p_.z() ; 
    }
    
    /** replaces p by the by-coordinate maximum of p and p_ */
    inline void pmax( const Point &p_ )
    { 
      if( x() < p_.x() ) x() = p_.x() ; 
      if( y() < p_.y() ) y() = p_.y() ; 
      if( z() < p_.z() ) z() = p_.z() ; 
    }
    
    //-----------------------------------------------------------------------------
    // Drawing
  public:
    /// Draws the point with opengl
    void draw () const { ::glVertex3f( (float)x(), (float)y(), (float)z() ) ; }
  } ;

// Valid point
extern Point P_INV ;


/**------------------------------------------------------------------------------------*/
// Geometric Operations

/// +
inline const Point operator+  ( const Point &p, const Point &p_ )
{
  Point r = p ;
  return (r += p_) ;
}
/// -
inline const Point operator-  ( const Point &p, const Point &p_ )
{
  Point r = p ;
  return r -= p_ ;
}
/// * scalar
inline const Point  operator*  ( const Point &p, const real l )
{
  Point r = p ;
  return r *= l ;
}
/// * scalar
inline const Point  operator*  ( const real l, const Point &p )
{
  return p * l ;
}
/// / scalar
inline const Point  operator/  ( const Point &p, const real l )
{
  Point r = p ;
  return r /= l ;
}

/// scalar (dot) product
inline real operator* ( const Point &p, const Point &p_ )
{
  return p.x()*p_.x() + p.y()*p_.y() + p.z()*p_.z() ;
}
/// vector (cross) product 
static inline const Point operator^  ( const Point &p, const Point &p_ )
{
  Point r ;
  r.x() = p.y() * p_.z()  -  p.z() * p_.y() ;
  r.y() = p.z() * p_.x()  -  p.x() * p_.z() ;
  r.z() = p.x() * p_.y()  -  p.y() * p_.x() ;
  return r ;
}
static inline real det( const Point &p, const Point &q, const Point &r )
{
  return
  p.x()*q.y()*r.z() + p.y()*q.z()*r.x() + p.z()*q.x()*r.y() -
  p.z()*q.y()*r.x() - p.x()*q.z()*r.y() - p.y()*q.x()*r.z() ;
}


/**------------------------------------------------------------------------------------*/
// Geometric Operations

inline real length ( const Point &p ) { return p.length() ; }
inline real norm   ( const Point &p ) { return p.norm() ; }
inline real dist   ( const Point &p, const Point &p_ ) { return p.dist(p_) ; }
inline real sqdist ( const Point &p, const Point &p_ ) { Point r = p-p_ ; return r*r ; }

/** \brief Computes the normal of a triangle
 * \param const Point &v0
 * \param const Point &v1
 * \param const Point &v2 */
inline const Point normal (const Point & v0, const Point & v1, const Point & v2)
{
  Point n = ( v1-v0 ) ^ ( v2 - v0 ) ;
  n.normalize() ;
  return n ;
}

/** \brief Computes the area of a triangle
 * \param const Point &v0
 * \param const Point &v1
 * \param const Point &v2 */
inline real area (const Point & v0, const Point & v1, const Point & v2)
{
  return 0.5 * norm( ( v1-v0 ) ^ ( v2 - v0 ) );
}

/** \brief Computes the cotangent of ange v v1,v v2
 * \param const Point &v
 * \param const Point &v1
 * \param const Point &v2 */
inline const real cotan (const Point & v, const Point & v1, const Point & v2)
{
  Point u = v1-v ;
  Point w = v2-v ;
  return (u*w) / norm( u ^ w ) ;
}

/** \brief Computes the middle of an edge
 * \param const Point &v0
 * \param const Point &v1 */
inline const Point middle (const Point & v0, const Point & v1)
{
  return 0.5 * (v0+v1) ;
}


inline void pmin( Point &p, const Point &p_ ) { p.pmin( p_ ) ; }
inline void pmax( Point &p, const Point &p_ ) { p.pmax( p_ ) ; }


