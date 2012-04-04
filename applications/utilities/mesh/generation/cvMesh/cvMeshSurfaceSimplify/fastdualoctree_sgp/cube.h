/**
 * \file    cube.h
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 * \brief  Simple cube class
 *
 * Simple cube class
 */
//_____________________________________________________________________________



#pragma once


#include "point.h"

typedef uchar Level ;
extern Level L_INV ;

//_____________________________________________________________________________
/// Simple cube class
class Cube : public Point
{
  // Elements
protected :
  Level    _lv ;  ///< cube level (from 0 on): cube corners have coordinates c(xyz) (+/-) (1/2^(lv+1))

  //-----------------------------------------------------------------------------
  // Constructors
public:
  /// Default constructor
  Cube( const real &cx_ = R_INV, const real &cy_ = R_INV, const real &cz_ = R_INV, Level lv_ = L_INV )
    : Point(cx_,cy_,cz_), _lv(lv_) {}

  /// Default constructor
  Cube( const Point& c_, Level lv_ = L_INV )
    : Point(c_), _lv(lv_) {}
  
  /// Default constructor
  Cube( const Cube& c_ )
    : Point((const Point&)c_), _lv(c_._lv) {}
  
  /// Destructor
  ~Cube() {}
    
  /// Assignment operator
  Cube &operator= ( const Cube &c )
  {  Point::operator=(c) ;  _lv = c._lv ;  return *this; }
  
  //-----------------------------------------------------------------------------
  // Accessors
public:
  /// const accessor for point center: x coordinate
  real  cx() const { return Point::x() ; }
  /// accessor for point center: x coordinate
  real &cx() { return Point::x() ; }

  /// const accessor for point center: y coordinate
  real  cy() const { return Point::y() ; }
  /// accessor for point center: y coordinate
  real &cy() { return Point::y() ; }

  /// const accessor for point center: z coordinate
  real  cz() const { return Point::z() ; }
  /// accessor for point center: z coordinate
  real &cz() { return Point::z() ; }

  /// const accessor for a coordinate
  real  coord( Axis a ) const { return Point::operator[](a) ; }
  /// accessor for a coordinate
  real &coord( Axis a ) { return Point::operator[](a) ; }


  /// const accessor for the cube half size
  real   sz() const { return (real)1.0 / (2<<_lv) ; }
  /// const accessor for cube level
  Level  lv() const { return _lv ; }
  /// accessor for cube level
  Level &lv() { return _lv ; }

  /// const accessor for cube center: min x coordinate
  real  xmin() const { return cx() - sz() ; }
  /// const accessor for cube center: max x coordinate
  real  xmax() const { return cx() + sz() ; }

  /// const accessor for cube center: min y coordinate
  real  ymin() const { return cy() - sz() ; }
  /// const accessor for cube center: max y coordinate
  real  ymax() const { return cy() + sz() ; }

  /// const accessor for cube center: min z coordinate
  real  zmin() const { return cz() - sz() ; }
  /// const accessor for cube center: max z coordinate
  real  zmax() const { return cz() + sz() ; }

  /// const accessor for cube center: min coordinate
  real  coordmin( Axis a ) const { return coord(a) - this->sz() ; }
  /// const accessor for cube center: max coordinate
  real  coordmax( Axis a ) const { return coord(a) + this->sz() ; }

  //------------------------------------------------
  // Test
public :
  /// conained
  bool contains( const real &x_, const real &y_, const real &z_ ) const
  {
    real sz_ = sz() ;
    return
      ( fabs(x_-cx()) <= sz_ ) &&
      ( fabs(y_-cy()) <= sz_ ) &&
      ( fabs(z_-cz()) <= sz_ ) ;
  }

  /// conained
  bool contains( const Point &p ) const { return contains( p.x(), p.y(), p.z() ) ; }

  /// epsilon-conained
  bool contains( const real &x_, const real &y_, const real &z_, const real &epsilon ) const
  {
    real sz_ = sz() + epsilon ;
    return
      ( fabs(x_-cx()) <= sz_ ) &&
      ( fabs(y_-cy()) <= sz_ ) &&
      ( fabs(z_-cz()) <= sz_ ) ;
  }

  /// epsilon-conained
  bool contains( const Point &p, const real &epsilon ) const { return contains( p.x(), p.y(), p.z(), epsilon ) ; }

  //-----------------------------------------------------------------------------
  // Drawing
public:
  /// Draws the cube center with opengl
  void draw_center () const { ::glVertex3f( (float)x(), (float)y(), (float)z() ) ; }

  /// Draws the cube corners with opengl
  void draw_corners( real fact = 0.9 ) const
  {
    const float ix = (const float) x() ;
    const float iy = (const float) y() ;
    const float iz = (const float) z() ;
    const float is = (const float) (fact*sz()) ;

    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix+is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;
    ::glVertex3f( ix-is, iy-is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;
  }

  /// Draws the cube wire with opengl
  void draw_wire  ( real fact = 0.9 ) const
  {
    const float ix = (const float) x() ;
    const float iy = (const float) y() ;
    const float iz = (const float) z() ;
    const float is = (const float) (fact*sz()) ;
    
    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz-is ) ;

    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix+is, iy-is, iz-is ) ;

    ::glVertex3f( ix+is, iy-is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;

    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;


    ::glVertex3f( ix-is, iy-is, iz+is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;

    ::glVertex3f( ix-is, iy-is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;

    ::glVertex3f( ix+is, iy-is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;

    ::glVertex3f( ix-is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;


    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy-is, iz+is ) ;

    ::glVertex3f( ix+is, iy-is, iz-is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;

    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;

    ::glVertex3f( ix+is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;
  }



  /// Draws the cube wire with opengl
  void draw_fill  ( real fact = 0.9 ) const
  {
    const float ix = (const float) x() ;
    const float iy = (const float) y() ;
    const float iz = (const float) z() ;
    const float is = (const float) (fact*sz()) ;
    
    ::glNormal3f(0,0,-1);
    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy-is, iz-is ) ;

    ::glNormal3f(0,0,1);
    ::glVertex3f( ix-is, iy-is, iz+is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;

    ::glNormal3f(-1,0,0);
    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;
    ::glVertex3f( ix-is, iy-is, iz+is ) ;

    ::glNormal3f(1,0,0);
    ::glVertex3f( ix+is, iy-is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;

    ::glNormal3f(0,-1,0);
    ::glVertex3f( ix-is, iy-is, iz-is ) ;
    ::glVertex3f( ix-is, iy-is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz+is ) ;
    ::glVertex3f( ix+is, iy-is, iz-is ) ;

    ::glNormal3f(0,1,0);
    ::glVertex3f( ix-is, iy+is, iz-is ) ;
    ::glVertex3f( ix-is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz+is ) ;
    ::glVertex3f( ix+is, iy+is, iz-is ) ;
  }
};
//_____________________________________________________________________________

// Invalid cube
extern Cube C_INV ;


