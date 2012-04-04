/**
 * @file    hash_noptr.h
 * @author  Thomas Lewiner 
 * @version 
 * @date    
 *
 * @brief Hash table without pointers
 *
 */
//________________________________________________


#ifndef _HASH_H_
# error This header should not be included directly! include hash.h instead
#endif // _HASH_H_

#pragma once


#include <stdlib.h> // memset

#include "morton.h"

//_____________________________________________________________________________
// Class hash
template < typename Data > class Hash
//-----------------------------------------------------------------------------
{
public:
  typedef uint cntr ;
//  static const cntr HASH_BITS = 16 ;
//  static const cntr NON_HASH_BITS = 64-HASH_BITS ;
  static const cntr HASH_SIZE = ((cntr)1<<HASH_BITS) ;
  static const cntr HASH_MASK = HASH_SIZE-1 ;
  static const cntr SKIP_MINI = 11 ;
  static const cntr SKIP_BITS = 3 ;
  static const cntr SKIP_MASK = ((cntr)1<<SKIP_BITS)-1 ;

  typedef struct { Key key ; Data data ; } KeyData ;
  static KeyData KD_INV ;

protected:
//  KeyData _hash[HASH_SIZE] ;
  KeyData *_hash ; // allocate dynamically since stack memory is limited
  
public:
	// Constructor & Desctructor  
  Hash   () { _hash = new KeyData[HASH_SIZE] ; reset() ; } 
  ~Hash  () { delete [] _hash ; }

  void reset() { memset( _hash, KEY_INV, HASH_SIZE*sizeof(KeyData) ) ; KD_INV.key = KEY_INV ; }
  
  inline cntr hash( Key k ) const ; 
  
  cntr size() const
  {
    cntr s = 0 ;
    cntr n = HASH_SIZE ;
    KeyData *ptr = _hash ;
    for( ; n > 0 ; --n, ++ptr )
      s += ( ptr->key != KEY_INV ) ;
    return s ;
  }
  
  bool insert( KeyData d )
  {
    cntr h = hash( d.key ) ;
    KeyData *p = _hash + h ;
    const cntr skip = ((d.key>>HASH_BITS) & SKIP_MASK) + SKIP_MINI ;
    cntr overflow   = skip ;
    while( p->key != KEY_INV && p->key != d.key )
    { // collision
      h += skip ;
      if( h < HASH_SIZE )
        p += skip ;
      else
      {
        h &= HASH_MASK ;
        p = _hash + h ;
        if( !(overflow--) )
          return false ;
      }
    }
    *p = d ;
    return true ;
  }

  const KeyData operator[] ( Key k ) const
  {
    cntr h = hash( k ) ;
    const KeyData *p = _hash + h ;
    const cntr skip = ((k>>HASH_BITS) & SKIP_MASK) + SKIP_MINI ;
    cntr overflow   = skip ;
    while( p->key != k && p->key != KEY_INV )
    { // collision
      h += skip ;
      if( h < HASH_SIZE )
        p += skip ;
      else
      {
        h &= HASH_MASK ;
        p = _hash + h ;
        if( !(overflow--) )
          return KD_INV ;
      }
    }
    return *p ;
  }
  
  KeyData &operator[] ( Key k )
  {
    cntr h = hash( k ) ;
    KeyData *p = _hash + h ;
    const cntr skip = ((k>>HASH_BITS) & SKIP_MASK) + SKIP_MINI ;
    cntr overflow   = skip ;
    while( p->key != k && p->key != KEY_INV )
    { // collision
      h += skip ;
      if( h < HASH_SIZE )
        p += skip ;
      else
      {
        h &= HASH_MASK ;
        p = _hash + h ;
        if( !(overflow--) )
          return KD_INV ;
      }
    }
    return *p ;
  }
  

  bool contains( Key k ) const { return (*this)[k].key == k ; }
  
  
  
  //------------------------------------------------
  // iterator
  public :
  /// iterator
  class iterator
  {
    public :
    // constructors
    /// default constructors
    iterator( Hash<Data> &hash_, cntr id_ = 0 ) : _ptr(hash_._hash+id_), _last(hash_._hash+HASH_SIZE) { while( _ptr != _last && _ptr->key == KEY_INV ) ++_ptr ; }
    
    /// destructor
    ~iterator() {}
    
    /// copy constructor
    iterator( const iterator &it ) : _ptr(it._ptr), _last(it._last) {}
    
    /// assignment operator
    iterator &operator = ( const iterator &it )
    { _ptr = it._ptr; _last = it._last;  return *this; }
    
    //-----------------------------------------------------------------------------
    // Operations
    public  :
    /// equality operator
    inline bool        operator ==( const iterator &it ) const { return _ptr == it._ptr ; }
    /// inequality operator
    inline bool        operator !=( const iterator &it ) const { return _ptr != it._ptr ; }
    
    /// validation operator
    inline bool        operator ()() const { return _ptr != _last ; }
    
    /// key accessor
    inline const Key  &key() const { return _ptr->key ; }
    
    /// value accessor
    inline const Data &operator * () const { return _ptr->data ; }
    
    /// value accessor
    inline Data       &operator * ()       { return _ptr->data ; }
    
    /// accesses the next position
    inline iterator   &operator ++()       { while( (++_ptr != _last) && _ptr->key == KEY_INV ) ;  return *this ; }
    
    //-----------------------------------------------------------------------------
    // Elements
    private :
    KeyData       *_ptr  ; ///< position pointer
    const KeyData *_last ; ///< last position pointer
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
    const_iterator( const Hash<Data> &hash_, cntr id_ = 0 ) : _ptr(hash_._hash+id_), _last(hash_._hash+HASH_SIZE) { while( _ptr != _last && _ptr->key == KEY_INV ) ++_ptr ; }
    
    /// destructor
    ~const_iterator() {}
    
    /// copy constructor
    const_iterator( const const_iterator &it ) : _ptr(it._ptr), _last(it._last) {}
    
    /// assignment operator
    const_iterator &operator = ( const const_iterator &it )
    { _ptr = it._ptr; _last = it._last;  return *this; }
    
    //-----------------------------------------------------------------------------
    // Operations
    public  :
    /// equality operator
    inline bool        operator ==( const const_iterator &it ) const { return _ptr == it._ptr ; }
    /// inequality operator
    inline bool        operator !=( const const_iterator &it ) const { return _ptr != it._ptr ; }
    
    /// validation operator
    inline bool        operator ()() const { return _ptr != _last ; }
    
    /// key accessor
    inline const Key  &key() const { return _ptr->key ; }
    
    /// value accessor
    inline const Data &operator * () const { return _ptr->data ; }
    
    /// accesses the next position
    inline const_iterator   &operator ++() { while( (++_ptr != _last) && _ptr->key == KEY_INV ) ;  return *this ; }
    
    //-----------------------------------------------------------------------------
    // Elements
    private :
    const KeyData * _ptr  ; ///< position pointer
    const KeyData * const _last ; ///< last position pointer
  };
  
  public :
  /// Node iterator creation
  iterator begin( cntr id = 0 ) { return iterator( *this, id ) ; }
  /// Node const iterator creation
  const_iterator cbegin( cntr id = 0 ) { return const_iterator( *this, id ) ; }
  
};
//_____________________________________________________________________________



/// hash function specification: no inversion
template <> inline Hash<real >::cntr Hash<real >::hash( Key k ) const { return k & HASH_MASK ; }

/// hash function specification: inversion
template <> inline Hash<Level>::cntr Hash<Level>::hash( Key k ) const { return (k>>NON_HASH_BITS) & HASH_MASK ; }
