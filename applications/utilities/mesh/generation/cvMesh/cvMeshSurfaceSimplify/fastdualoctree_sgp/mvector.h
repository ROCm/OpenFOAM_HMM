/**
 * \file    mvector.h
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 * \brief  light array container
 *
 * Light array container
 */
//_____________________________________________________________________________


#pragma once

#include <stdlib.h>
#include <string.h>  // memcpy, memset
#include <mlist.h>

/// unsigned int alias
typedef unsigned int uint ;
/// const unsigned int alias
typedef const uint cuint ;

/// default data size
#define INITIAL_MVECTOR_SIZE ((uint)1024)


//_____________________________________________________________________________
// mvector
/// \class mvector mvector.h
template< typename Data > class mvector
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /// Constructor
  mvector() : _data((Data*)NULL), _allocated_size(0), _used_size(0) {} ;

  /// Destructor
  ~mvector() { clear() ; }

  // no copy on purpose

//------------------------------------------------
// Operations
public :
  /// delete the data
  void clear() { free(_data) ; _data = (Data*)NULL ;  _allocated_size = 0 ;  _used_size = 0 ;  _deleted.clear() ; }

  /// estimation of the memory that will be used
  bool reserve( uint n )
  {
    if( n < INITIAL_MVECTOR_SIZE ) n = INITIAL_MVECTOR_SIZE ;
    if( n <= _allocated_size ) return true ;

    Data *temp = _data ;
    _data = (Data*)malloc( n*sizeof(Data) ) ;
    if( temp ) memcpy( _data, temp, _used_size*sizeof(Data) ) ;
    memset( _data + _used_size, 0, (n-_used_size)*sizeof(Data) ) ;
    free( temp ) ;

    _allocated_size = n ;
    return _data != (Data*)NULL ;
  }

  /// reserve space for n new elements
  bool reserve_more( uint n )
  {
    if( _used_size + n > _allocated_size + _deleted.size() )
    {
      uint new_size = 2*_allocated_size ;
      if( new_size < _used_size + n ) new_size = _used_size + n ;
      //return reserve( 2*_allocated_size ) ;
      return reserve( new_size ) ;
    }
    return true ;
  }

  /// insert n new elements
  bool allocate_more( uint n )
  {
    if( !reserve_more( n ) ) return false ;
    _used_size += n ;
    return true ;
  }

  /// add one more element and
  Data &add( uint &id )
  {
    if( _deleted.empty() )
    {
      id = _used_size ;
      allocate_more( 1 ) ;
    }
    else
    {
      id = _deleted.first() ;
      _deleted.remove_first() ;
    }
    return at(id) ;
  }

  /// lazy remotion
  bool remove( uint id )
  {
    if( id+1 == _used_size )
      --_used_size ;
    else
      _deleted.insert(id) ;
    return true ;
  }


//------------------------------------------------
// Accessor
public :
  /// size accessor
  uint size() {  return _used_size - _deleted.size() ; }

  /// size accessor
  uint used_size() {  return _used_size ; }

  /// non checked const accessor
  const Data &at( uint id ) const { return _data[id] ; }
  /// non checked accessor
  Data &at( uint id ) { return _data[id] ; }

  /// non checked const accessor
  const Data &operator[]( uint id ) const { return at(id) ; }
  /// non checked accessor
  Data &operator[]( uint id ) { return at(id) ; }

  /// non checked pointer accessor
  Data *operator +( uint id ) { return _data + id ; }

  /// id range validity check
  bool is_valid( uint id ) {  return (id < size() && id>=0) ; }

  /// id validity check
  bool is_deleted( uint id ) {  return is_valid(id) && _deleted.find(id) () ; }


//------------------------------------------------
// Elements
private :
  /// main data array
  Data *_data ;

  /// real data size
  uint  _allocated_size ;

  /// filled data size
  uint  _used_size ;

  /// deleted data
  List<uint> _deleted ;


//------------------------------------------------
// iterator
public :
  /// iterator
  class iterator
  {
  public :
  // constructors
    /// default constructors
    iterator( mvector<Data> &vec_, uint id_ = 0 ) : _vec(vec_), _id(id_) {}

    /// destructor
    ~iterator() {}

    /// copy constructor
    iterator( const iterator &it ) : _vec(it._vec), _id(it._id) {}

    /// assignment operator
    iterator &operator = ( const iterator &it )
    { _vec = it._vec; _id = it._id;  return *this; }

  //-----------------------------------------------------------------------------
  // Operations
  public  :
    /// equality operator
    inline bool        operator ==( const iterator &it ) const { return &_vec == &it._vec && _id == it._id ; }
    /// inequality operator
    inline bool        operator !=( const iterator &it ) const { return &_vec != &it._vec || _id != it._id ; }

    /// validation operator
    inline bool        operator ()() const { return _vec.is_valid( _id ) ; }

    /// id accessor
    inline uint        id         () const { return _id                  ; }

    /// value accessor
    inline const Data &operator * () const { return _vec[_id]            ; }

    /// value accessor
    inline Data       &operator * ()       { return _vec[_id]            ; }

    /// accesses the next position
    inline iterator   &operator ++()       { do ++_id ; while( _vec.is_valid(_id) && _vec.is_deleted(_id) ) ;  return *this ; }

  //-----------------------------------------------------------------------------
  // Elements
  private :
    mvector<Data> &_vec ; ///< data vector
    uint           _id  ; ///< node position
  };

  //------------------------------------------------
  // const_iterator
  public :
  /// const_iterator
  class const_iterator
  {
    public :
    // constructors
    /// default constructors
    const_iterator( const mvector<Data> &vec_, uint id_ = 0 ) : _vec(vec_), _id(id_) {}
    
    /// destructor
    ~const_iterator() {}
    
    //-----------------------------------------------------------------------------
    // Operations
    public  :
    /// equality operator
    inline bool        operator ==( const const_iterator &it ) const { return &_vec == &it._vec && _id == it._id ; }
    /// inequality operator
    inline bool        operator !=( const const_iterator &it ) const { return &_vec != &it._vec || _id != it._id ; }
    
    /// validation operator
    inline bool        operator ()() const { return _vec.is_valid( _id ) ; }
    
    /// id accessor
    inline uint        id         () const { return _id                  ; }
    
    /// value accessor
    inline const Data &operator * () const { return _vec[_id]            ; }
    
    /// accesses the next position
    inline const_iterator   &operator ++()       { do ++_id ; while( _vec.is_valid(_id) && _vec.is_deleted(_id) ) ;  return *this ; }
    
    //-----------------------------------------------------------------------------
    // Elements
    private :
    const mvector<Data> &_vec ; ///< data vector
    uint                 _id  ; ///< node position
  };
  
  public :
  /// Node iterator creation
  iterator begin( uint id = 0 ) { return iterator( *this, id ) ; }
  /// Node const iterator creation
  const_iterator cbegin( uint id = 0 ) { return const_iterator( *this, id ) ; }

};
//_____________________________________________________________________________


