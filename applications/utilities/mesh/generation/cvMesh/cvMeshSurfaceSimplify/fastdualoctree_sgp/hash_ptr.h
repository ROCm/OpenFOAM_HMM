/**
 * @file    hash_ptr.h
 * @author  Thomas Lewiner 
 * @version 
 * @date    
 *
 * @brief
 *
 */
//________________________________________________


#ifndef _HASH_H_
# error This header should not be included directly! include hash.h instead
#endif // _HASH_H_

#pragma once

#define HASH_HAS_ERASE 1

#include <stdlib.h> // memset
#include <string.h>
#include <stdio.h> //printf
#include "mlist.h" // memset

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

  // Bucket associative data
  typedef struct KeyData
  {
    Key key ;
    Data data ;

    /// comparison based on the key only
    bool operator==( const struct KeyData& rhs) const { return key == rhs.key ; }

    /// comparison based on the key only
    bool operator!=( const struct KeyData& rhs) const { return key != rhs.key ; }
    
  } KeyData ;
  static KeyData KD_INV ;
  
  class iterator ;

public:
  typedef List<KeyData>                      KeyBucket ;
  typedef typename KeyBucket::iterator       keybucket_iterator ;
  typedef typename KeyBucket::const_iterator keybucket_const_iterator ;

protected:
//  KeyBucket _hash[HASH_SIZE] ;
  KeyBucket *_hash ; // allocate dynamically since stack memory is limited
  
public:
	// Constructor & Desctructor  
  Hash   () { _hash = new KeyBucket[HASH_SIZE] ; } 
  ~Hash  () { delete [] _hash ; }

  void reset()
  {
    cntr n = HASH_SIZE ;
    KeyBucket *ptr = _hash ;
    for( ; n > 0 ; --n, ++ptr )
    {
      (*ptr).clear() ;
    }
    KD_INV.key = KEY_INV ;
  }
  
  cntr size() const
  {
    cntr s = 0 ;
    cntr n = HASH_SIZE ;
    KeyBucket *ptr = _hash ;
    for( ; n > 0 ; --n, ++ptr )
      s += ptr->size() ;
    return s ;
  }
  
  void stats() const
  {
    static const cntr MAX_COLLISION = 32 ;
    cntr colls[MAX_COLLISION] ;
    memset( colls, 0, MAX_COLLISION*sizeof(cntr) ) ;
    KeyBucket *ptr = _hash ;
    for( cntr n = HASH_SIZE ; n > 0 ; --n, ++ptr )
    {
      cntr s = ptr->size() ;
      if( s >= MAX_COLLISION ) s = MAX_COLLISION-1 ;
      ++colls[s] ;
    }

    printf( "Hashtable collisions: [\t") ;
    for( cntr n = 0 ; n < MAX_COLLISION ; ++n )
      printf( "%d\t", colls[n] ) ;
    printf( "\t]\n") ;
  }
  
  inline cntr hash( Key k ) const ; 

  bool insert( KeyData d )
  {
    cntr h = hash( d.key ) ;
    KeyBucket &l = _hash[h] ;
    return l.insert_unique( d ) ;
  }

  const KeyData operator[] ( Key k ) const
  {
    cntr h = hash( k ) ;
    KeyData kd ; kd.key = k ;
    keybucket_const_iterator it = _hash[h].cfind( kd ) ;
    if( !it() ) return KD_INV ;
    return *it ;
  }
  
  KeyData &operator[] ( Key k )
  {
    cntr h = hash( k ) ;
    KeyData kd ; kd.key = k ;
    keybucket_iterator it = _hash[h].find( kd ) ;
    if( !it() ) return KD_INV ;
    return *it ;
  }
  
  KeyData erase( Key k )
  {
    cntr h = hash( k ) ;
    KeyData kd ; kd.key = k ;
    if( _hash[h].remove( kd ) )
      return kd ;
    return KD_INV ;
  }
  
  bool contains( Key k ) const { return (*this)[k].key == k ; }
  
  
  //------------------------------------------------
  // iterators
  public :
  /// iterator
  class iterator
  {
    public :
    // constructors
    /// default constructors
    iterator( Hash<Data> &hash_, Key id_ = 0 ) : _ptr(hash_._hash+id_), _last(hash_._hash+HASH_SIZE)
    {
      while( _ptr != _last && _ptr->empty() )
        ++_ptr ;
      if( _ptr != _last )
        _l_it = _ptr->begin() ;
    }
    
    /// destructor
    ~iterator() {}
    
    /// copy constructor
    iterator( const iterator &it ) : _ptr(it._ptr), _last(it._last), _l_it(it._l_it) {}
    
    /// assignment operator
    iterator &operator = ( const iterator &it )
    { _ptr = it._ptr; _last = it._last;  _l_it = it._l_it ;  return *this; }
    
    //-----------------------------------------------------------------------------
    // Operations
    public  :
    /// equality operator
    inline bool        operator ==( const iterator &it ) const { return _ptr == it._ptr && _l_it == it._l_it ; }
    /// inequality operator
    inline bool        operator !=( const iterator &it ) const { return _ptr != it._ptr || _l_it != it._l_it ; }
    
    /// validation operator
    inline bool        operator ()() const { return _ptr != _last ; }
    
    /// key accessor
    inline const Key  &key() const { return (*_l_it).key ; }
    
    /// value accessor
    inline const Data &operator * () const { return (*_l_it).data ; }
    
    /// value accessor
    inline Data       &operator * ()       { return (*_l_it).data ; }
    
    /// accesses the next position
    inline iterator   &operator ++()       {
      if( !(++_l_it)() )
      {
        while( _ptr != _last && (*(++_ptr)).empty() ) ;
        if( _ptr != _last )
          _l_it = _ptr->begin() ;
      }
      return *this ;
    }
    
    //-----------------------------------------------------------------------------
    // Elements
    private :
    KeyBucket           *_ptr  ; ///< position pointer
    const KeyBucket     *_last ; ///< last position pointer
    keybucket_iterator   _l_it ; ///< list iterator inside a bucket
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
    const_iterator( Hash<Data> &hash_, Key id_ = 0 ) : _ptr(hash_._hash+id_), _last(hash_._hash+HASH_SIZE)
    {
      while( _ptr != _last && _ptr->empty() )
        ++_ptr ;
      if( _ptr != _last )
        _l_it = _ptr->cbegin() ;
    }
    
    /// destructor
    ~const_iterator() {}
    
    /// copy constructor
    const_iterator( const const_iterator &it ) : _ptr(it._ptr), _last(it._last), _l_it(it._l_it) {}
    
    /// assignment operator
    const_iterator &operator = ( const const_iterator &it )
    { _ptr = it._ptr; _last = it._last;  _l_it = it._l_it ;  return *this; }
    
    //-----------------------------------------------------------------------------
    // Operations
    public  :
    /// equality operator
    inline bool        operator ==( const const_iterator &it ) const { return _ptr == it._ptr && _l_it == it._l_it ; ; }
    /// inequality operator
    inline bool        operator !=( const const_iterator &it ) const { return _ptr != it._ptr || _l_it != it._l_it ; ; }
    
    /// validation operator
    inline bool        operator ()() const { return _ptr != _last ; }
    
    /// key accessor
    inline const Key  &key() const { return (*_l_it).key ; }
    
    /// value accessor
    inline const Data &operator * () const { return (*_l_it).data ; }
    
    /// accesses the next position
    inline const_iterator   &operator ++()       {
      if( !(++_l_it)() )
      {
        while( _ptr != _last && (*(++_ptr)).empty() ) ;
        if( _ptr != _last )
          _l_it = _ptr->cbegin() ;
      }
      return *this ;
    }
    
    //-----------------------------------------------------------------------------
    // Elements
    private :
    const KeyBucket  *         _ptr  ; ///< position pointer
    const KeyBucket  * const   _last ; ///< last position pointer
    keybucket_const_iterator   _l_it ; ///< list iterator inside a bucket
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
