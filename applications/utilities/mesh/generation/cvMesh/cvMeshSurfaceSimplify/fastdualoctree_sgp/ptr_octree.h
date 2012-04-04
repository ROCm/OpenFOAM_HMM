/**
 * \file    ptr_octree.h
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree structure with son-brother pointers
 */
//_____________________________________________________________________________


#pragma once


#ifndef WIN32
#pragma interface
#endif // WIN32


#include <stack>
#include "mlist.h"
#include "cube.h"
#include "MarchingCubes.h"
#include "mc_draw.h"
#include "data_access.h"


//_____________________________________________________________________________
// PtrOctree
/// \class PtrOctree PtrOctree.h
class PtrOctree
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
  typedef struct cell
  {
    struct cell *brother   ;  ///< next brother of the cell
    struct cell *first_son ;  ///< first son of the cell
    real         field     ;  ///< field supported by the cell
  } cell ;

  /// Root of the octree
  cell _root ;

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
  PtrOctree() { init() ; }

  /// Destructor: Free memory
  ~PtrOctree() { clear() ; }


  /// Create the root
  void init() ;

  /// Memory cleaning
  void clear() { clear_octree() ;  _mc.clean_all() ; }

  /// Delete the branch below a given node
  void clear_branch( cell *b ) ;

  /// Delete the content of the octree
  void clear_octree() { clear_branch( &_root ) ; }

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
  typedef bool ptr_dual_walker( PtrOctree &fo, geom_cell *cells ) ;

  /// Walk on the dual cubes
  bool dual_cubes_walk( ptr_dual_walker &walker ) ;

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
  inline cell_iterator cells_begin() { return cell_iterator( &_root ) ; }

  /// Create an iterator traversing the leaves of the tree from the root
  inline leaf_iterator leaves_begin() { return leaf_iterator( &_root ) ; }

//-----------------------------------------------------------------------------
// search operations
public:
  /// Find cells of the octree at a given position
  bool find_leaf( real x, real y, real z, geom_cell &cell ) ;

  /// Find cells of the octree inside a given box of center x,y,z and half side r
  bool find_radius( real x, real y, real z, real r, List<geom_cell> &cells ) ;

  /// Find adjacent cells of the octree to a given cell
  bool adjacent( const geom_cell &cell, List<geom_cell> &cells ) ;


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
    friend class PtrOctree ;

  protected:
    PtrOctree::cell *_cell   ;  ///< octree cell

  //---------------------------------------------------------------------------
  // Constructors
  public:
    /// Default constructor: Constructs an iterator from a cell
    geom_cell( PtrOctree::cell *cell_ = NULL, real cx_ = 0.5, real cy_ = 0.5, real cz_ = 0.5, Level lv_ = 0 )
      : Cube(cx_,cy_,cz_,lv_), _cell(cell_) {}

    /// Destructor
    ~geom_cell() {}


    /// Copy constructor
    geom_cell( const geom_cell &i )
      : Cube(i), _cell(i._cell) {}

    /// Assignment operator
    geom_cell &operator = ( const geom_cell &i )
    { Cube::operator=(i); _cell=i._cell; return *this; }

  //---------------------------------------------------------------------------
  // Public constant accessors
  public  :
    /// cell const accessor
    inline const PtrOctree::cell * cell() const { return _cell ; }
    /// cell accessor
    inline PtrOctree::cell *&cell() { return _cell ; }

    /// first son const accessor
    inline const PtrOctree::cell * first_son() const { return _cell->first_son ; }
    /// first son accessor
    inline PtrOctree::cell *&first_son() { return _cell->first_son ; }
    
    /// brother const accessor
    inline const PtrOctree::cell * brother() const { return _cell->brother ; }
    /// brother accessor
    inline PtrOctree::cell *&brother() { return _cell->brother ; }
    
    /// id const accessor
    inline real  operator*() const { return cell()->field ; }
    /// id accessor
    inline real &operator*() { return cell()->field ; }
    
    //---------------------------------------------------------------------------
  // Tests
  public :
    /// equality operator
    inline bool operator ==( const geom_cell &i ) const { return cell() == i.cell() ; }

    /// inequality operator
    inline bool operator !=( const geom_cell &i ) const { return cell() != i.cell() ; }

    /// leaf test
    inline bool is_leaf() const { return first_son() == NULL ; }

    /// validation operator
    inline bool      operator ()() const { return cell() != NULL ; }

  //---------------------------------------------------------------------------
  // Operations
  public  :
    /// sons
    inline bool sons( geom_cell *s /*[8]*/ )
    {
      if( !first_son() ) return false ;

      int    level_ = lv() + 1 ;
      real sz_    = sz() / 2.0 ;
      for( int i = 0 ; i < 8 ; ++i )
      {
        s[i].cell () = (i==0) ? first_son() : s[i-1].brother() ;
        s[i].lv   () = level_ ;
        s[i].cx   () = (i&1) ? cx() + sz_ : cx() - sz_ ;
        s[i].cy   () = (i&2) ? cy() + sz_ : cy() - sz_ ;
        s[i].cz   () = (i&4) ? cz() + sz_ : cz() - sz_ ;
      }
      return true ;
    }

    /// get son from side i
    inline bool son( int i , geom_cell &s )
    {
      if( is_leaf() ) return false ;

      real sz_ = s.sz() / 2.0 ;
      ++s.lv() ;
      s.cx() = (i&1) ? cx() + sz_ : cx() - sz_ ;
      s.cy() = (i&2) ? cy() + sz_ : cy() - sz_ ;
      s.cz() = (i&4) ? cz() + sz_ : cz() - sz_ ;

      s.cell() = first_son() ;
      while( --i > -1 ) s.cell() = s.brother() ;

      return true ;
    }
  };


//_____________________________________________________________________________
// Cell Iterator
public :
  /// Octree cell iterator : Traverse the octree returning basic information on the cells
  class cell_iterator
  //---------------------------------------------------------------------------
  {
    friend class PtrOctree ;

  protected:
    /// Octree traversal stack
    std::stack<geom_cell> _s ;

  //---------------------------------------------------------------------------
  // Constructors
  public:
    /// Default constructor : Constructs an iterator from a cell
    cell_iterator( PtrOctree::cell *root = NULL, real cx_ = 0.5, real cy_ = 0.5, real cz_ = 0.5, Level lv_ = 0 )
    { if( root ) _s.push( geom_cell( root, cx_, cy_, cz_, lv_ ) ) ; }

    /// Destructor
    ~cell_iterator() {}

    /// Copy constructor
    cell_iterator( const cell_iterator &i ) : _s(i._s) {}

    /// Assignment operator
    cell_iterator &operator = ( const cell_iterator &i )
    { _s = i._s; return *this; }

  //---------------------------------------------------------------------------
  // Operations
  public  :
    /// validation operator
    inline bool      operator ()() const { return !_s.empty() ; }

    /// next position
    inline cell_iterator &operator ++()
    {
      if( _s.empty() ) return *this ;
      geom_cell n = _s.top() ;  _s.pop() ;
      if( n.is_leaf() ) return *this ;

      // depth first search
      geom_cell sons[8] ;
      n.sons( sons ) ;
      for( int i = 0 ; i < 8 ; ++i )
        _s.push( sons[i] ) ;

      return *this ;
    }

  //---------------------------------------------------------------------------
  // Accessors
  public  :
    /// geom_cell accessor
    inline geom_cell &top() { return _s.top() ; }

    /// id accessor
    inline real &operator*() { return * _s.top() ; }

    /// level accessor
    inline Level  &lv() { return _s.top().lv() ; }

    /// center x position accessor
    inline real &cx() { return _s.top().cx() ; }

    /// center y position accessor
    inline real &cy() { return _s.top().cy() ; }

    /// center y position accessor
    inline real &cz() { return _s.top().cz() ; }

    /// size accessor
    inline real  sz() { return _s.top().sz() ; }

    /// points accessor
    inline bool is_leaf() const { return _s.top().is_leaf() ; }

    /// inside test
    inline bool contains( real x, real y, real z ) const { return _s.top().contains( x,y,z ) ; }

    /// get son from its side:
    inline bool son( int i , geom_cell &s_ )       { return _s.top().son(i,s_) ; }

    /// Draws the cell wire with opengl
    void draw_wire  () const { _s.top().draw_wire  () ; }

  protected :
    /// stack accessor for ad hoc iterations
    inline std::stack< geom_cell > &s() { return _s ; }
  };


  /// Octree leaf iterator : Traverse the octree returning basic information on the leaves
  class leaf_iterator : public cell_iterator
  //---------------------------------------------------------------------------
  {
    public :
    leaf_iterator( PtrOctree::cell *root = NULL, real cx_ = 0.5, real cy_ = 0.5, real cz_ = 0.5, Level lv_ = 0 ) : cell_iterator( root, cx_,cy_,cz_,lv_ )
    { if( root && !this->is_leaf() ) ++(*this) ; }

    
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


