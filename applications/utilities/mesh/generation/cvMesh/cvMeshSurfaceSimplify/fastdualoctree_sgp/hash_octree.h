/**
 * \file    hash_octree.h
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree structure with hashtable and hierarchical operations
 */
//_____________________________________________________________________________


#pragma once

#ifndef WIN32
#pragma interface
#endif // WIN32


#include "mlist.h"
#include "cube.h"
#include "hash.h"
#include "MarchingCubes.h"
#include "mc_draw.h"
#include "data_access.h"



//_____________________________________________________________________________
// HashOctree
/// \class HashOctree HashOctree.h
class HashOctree
//-----------------------------------------------------------------------------
{
// forward declaration
public:
  class geom_cell  ;
  class cell_iterator  ;
  class leaf_iterator  ;

// Elements
protected:
  /// Octree cell data structure
  typedef Hash<real> HashField ;

  /// Hashtable of the octree nodes
  HashField _hash ;

  /// Maximal level of the octree (for find_adjacent)
  Level _max_level ;
  
  /// Maximal field of the octree
  real _max_field ;
  
  /// Isosurface
  MarchingCubes _mc ;

  /// Isosurface Draw
  MC_Draw _mc_draw ;
  
  /// Dual run memory consumption
  uint _dual_temp_memory ;
  
//-----------------------------------------------------------------------------
// Constructors
public:
  /// Default constructor: Constructs a octree with and empty root
  HashOctree() { init() ; }

  /// Destructor: Free memory
  ~HashOctree() { clear() ; }


  /// Create the root
  void init() ;

  /// Memory cleaning
  void clear() { clear_octree() ;  _mc.clean_all() ; }

  /// Delete the content of the octree
  void clear_octree() ;

  /// Check that only the leaves have data
  bool check () ;

  /// Prints statistics about the octree
  void stats() ;
  
  
  /// Return the maximal level of the octree
  Level  max_level() const { return _max_level ; }
  
  /// Return the maximal field of the octree
  real max_field() const { return _max_field ; }
   
  /// Return the isosurface
  MarchingCubes &mc() { return _mc ; }

  /// Return the isosurface draw
  MC_Draw &mc_draw() { return _mc_draw ; }
  
  /// set the values of each leaf from the implicit function
  bool set_impl( data_access *ref = NULL ) ;
  
  /// Refine the octree according to the data access
  bool refine( data_access *ref = NULL ) ;

  /// Adapt the octree according to the data access
  bool adapt( data_access *ref = NULL ) ;
  
  /// Draw the octree with wireframe
  bool draw_wire() ;
  
  /// Draw the octree with dots
  bool draw_centers() ;
  
  /// Dual function type
  typedef bool hash_dual_walker( HashOctree &fo, Key *keys ) ;

  /// Walk on the dual cubes
  bool dual_cubes_walk( hash_dual_walker &walker ) ;

  /// Build the isosurface using dual marching cubes
  bool build_isosurface( data_access *ref = NULL ) ;
  
  /// Draw the isosurface on-the-fly using dual marching cubes
  bool direct_draw_isosurface( data_access *ref = NULL ) ;
  
  /// Draw the dual octree with wireframe
  bool draw_dual() ;

  /// Do nothing, just for dual generation timing
  bool dual_timing() ;
    
//-----------------------------------------------------------------------------
// iterators
public:
  /// Create an iterator traversing the tree from the root
  inline cell_iterator cells_begin() { return cell_iterator( *this ) ; }

  /// Create an iterator traversing the leaves of the tree from the root
  inline leaf_iterator leaves_begin() { return leaf_iterator( *this ) ; }

//-----------------------------------------------------------------------------
// search operations
public:
  /// Find cells of the octree at a given position
  bool find_leaf( real x, real y, real z, geom_cell &cell ) const ;

  /// Find cells of the octree inside a given box of center x,y,z and half side r
  bool find_radius( real x, real y, real z, real r, List<geom_cell> &cells ) const ;

  /// Find adjacent cells of the octree to a given cell
  bool adjacent( const geom_cell &cell, List<geom_cell> &cells ) const ;

  /// Leaf test based on the field value
  bool is_leaf( Key k ) const { return !is_inv( _hash[k].data ) ; }

//-----------------------------------------------------------------------------
// I/O
public:

  /// Draws the intersecting cells of a plane with the octree with color
  void draw_plane ( real nx, real ny, real nz, real d	  ) ;

  /// Draws the intersection of a plane with the octree with color
  void draw_slice ( real nx, real ny, real nz, real d, float alpha ) ;

  /// Draws the isosurface of level l inside the dual graph
  void draw_iso   () { _mc.draw_surf() ; }

//_____________________________________________________________________________
// Iterator Cell
public :
  /// Auxiliary structure to traverse the octree
  class geom_cell : public Cube
  //---------------------------------------------------------------------------
  {
    friend class HashOctree ;

  protected:
    Key   _key   ;  ///< octree cell
    real  _field ;  ///< field associated to the cell

  //---------------------------------------------------------------------------
  // Constructors
  public:
    /// Default constructor: Constructs an iterator from a cell
    geom_cell( Key key_ = KEY_INV, real field_ = R_INV ) : Cube(key2cube(key_)), _key(key_), _field(field_) {}

    /// Destructor
    ~geom_cell() {}

    /// Copy constructor
    geom_cell( const geom_cell &i ) : Cube(i), _key(i._key), _field(i._field) {}

    /// Assignment operator
    geom_cell &operator = ( const geom_cell &i )
    {  Cube::operator=(i) ;  _key = i._key ;  _field = i._field ;  return *this; }
    

  //---------------------------------------------------------------------------
  // Public constant accessors
  public  :
    /// key const accessor
    inline Key  key() const { return _key ; }
    /// key accessor
    inline Key &key() { return _key ; }

    /// id const accessor
    inline real  operator*() const { return _field ; }
    
    //---------------------------------------------------------------------------
  // Tests
  public :
    /// equality operator
    inline bool operator ==( const geom_cell &i ) const { return key() == i.key() ; }

    /// inequality operator
    inline bool operator !=( const geom_cell &i ) const { return key() != i.key() ; }

    /// leaf test
    inline bool is_leaf() const { return !is_inv(*(*this)) ; }

    /// validation operator
    inline bool      operator ()() const { return key() != KEY_INV ; }

  //---------------------------------------------------------------------------
  // Operations
  public  :
    /// sons
    inline bool sons( geom_cell *s /*[8]*/, const HashField &hash )
    {
      if( is_leaf() ) return false ;

      Key k = _key << 3 ;
      for( int i = 0 ; i < 8 ; ++i )
      {
        s[i]._key   = k | i ;
        s[i]._field = hash[k|i].data ;
        (Cube&)s[i] = key2cube(s[i]._key) ;
      }
      return true ;
    }

    /// get son from side i
    inline bool son( int i , geom_cell &s, const HashField &hash )
    {
      if( is_leaf() ) return false ;
      s._key   = (key() << 3) | i ;
      s._field = hash[s._key].data ;
      (Cube&)s = key2cube(s._key) ;
      return true ;
    }
  };

  const geom_cell geom_root() const { const Key root_key = 1 ; return geom_cell( root_key, _hash[root_key].data ) ; }
  const geom_cell geom_key ( Key k ) const { return geom_cell( k, _hash[k].data ) ; }



//_____________________________________________________________________________
// Cell Iterator
public :
  /// Octree cell iterator : Traverse the octree returning basic information on the cells
  class cell_iterator
  //---------------------------------------------------------------------------
  {
    friend class HashOctree ;

  protected:
    /// Octree traversal iterator
    HashField::iterator _it ;

  //---------------------------------------------------------------------------
  // Constructors
  public:
    /// Default constructor : Constructs an iterator from a cell
    cell_iterator( HashOctree &o ) : _it( o._hash.begin() ) {}

    /// Destructor
    ~cell_iterator() {}

    /// Copy constructor
    cell_iterator( const cell_iterator &i ) : _it(i._it) {}

    /// Assignment operator
    cell_iterator &operator = ( const cell_iterator &i )
    { _it = i._it; return *this; }

  //---------------------------------------------------------------------------
  // Operations
  public  :
    /// equality operator
    inline bool operator ==( const cell_iterator &i ) const { return _it == i._it ; }

    /// inequality operator
    inline bool operator !=( const cell_iterator &i ) const { return _it != i._it ; }

    /// validation operator
    inline bool      operator ()() const { return _it() ; }

    /// next position
    inline cell_iterator &operator ++() { ++_it ; return *this ; }

  //---------------------------------------------------------------------------
  // Accessors
  public  :
    // cell accessor
    inline geom_cell top() const { return geom_cell( key(), *_it ) ; }

    /// id accessor
    inline real &operator*() { return *_it ; }

    /// level accessor
    inline Level lv() { return key_level(_it.key()) ; }

    /// size accessor
    inline real  sz() { return Cube(0,0,0,lv()).sz() ; }
    
    /// key accessor
    inline Key   key() const { return _it.key() ; }
    
    /// points accessor
    inline bool is_leaf() const { return !is_inv(*_it) ; }

    /// Draws the cell wire with opengl
    void draw_wire  () const { top().draw_wire  () ; }

  };


  /// Octree leaf iterator : Traverse the octree returning basic information on the leaves
  class leaf_iterator : public cell_iterator
  //---------------------------------------------------------------------------
  {
    public :
    leaf_iterator( HashOctree &o ) : cell_iterator( o )
    { if( (*this)() && !this->is_leaf() ) ++(*this) ; }

    
    /// next position
    inline leaf_iterator &operator ++()
    {
      cell_iterator &it = *this ;
      do ++it ; while ( it() && !it.is_leaf() ) ;
      return *this ;
    }
  } ;
} ;
//_____________________________________________________________________________


