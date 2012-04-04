/**
 * @file    mlist.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @version 0.3
 * @date    28/06/2004
 *
 * @brief   template list class
 */
//_____________________________________________________________________________



#pragma once

#include <stdlib.h> // qsort

typedef unsigned int  lsize ;
typedef const lsize   clsize ;

//_____________________________________________________________________________
/// template list class
template <typename T> class List
//-----------------------------------------------------------------------------
{
// types
public :
  class node ;
  class iterator ;
  class const_iterator ;

//-----------------------------------------------------------------------------
// default constructors
public  :
  /// default constructors
  List() : _first((node*)NULL) {}

  /// one-element constructor
  List( const T &val ) : _first((node*)NULL) { insert(val) ; }

  /// copy constructor (reverse order)
  List( const List<T> &l ) : _first((node*)NULL)
    {
      for( const_iterator it = l.cbegin() ; it() ; ++it )
        insert( *it ) ;
    }

  /// assignment operator
  List<T> &operator = ( const List<T> &l ) { if( this != &l) { clear() ;  for( const_iterator it = l.cbegin() ; it() ; ++it ) insert( *it ) ; }  return *this; }

  /// destructor
  ~List() { clear() ; }

  /// number of elements of the list
  clsize size() const { lsize n = 0 ; for( const_iterator it = cbegin() ; it() ; ++it ) ++n ; return n ; }

  /// number of elements of the list
  bool size_greater( clsize m ) const { lsize n = 0 ; for( const_iterator it = cbegin() ; it() ; ++it ) { if(++n > m) return true ; } return false ; }

  /// number of elements of the list
  bool size_greater_or_equal( clsize m ) const { if( m == 0 ) return true ; lsize n = 0 ; for( const_iterator it = cbegin() ; it() ; ++it ) { if(++n >= m) return true ; } return false ; }

  /// tests if the list has one element
  bool empty () const { return _first == (node*)NULL ; }

  /// tests if the list has two elements
  bool single() const { return !empty() && _first->next() == (node*)NULL ; }

  /// tests if the list has three elements
  bool pair  () const { return !empty() && _first->next() != (node*)NULL && _first->next()->next() == (node*)NULL ; }

  /// tests if the list has more than three elements
  bool triple() const { return !empty() && _first->next() != (node*)NULL && _first->next()->next() != (node*)NULL ; }

  const T  first () const { return _first->val() ; }  ///< first element accessor
  T       &first ()       { return _first->val() ; }  ///< first element accessor
  const T  second() const { return _first->next()->val() ; }  ///< second element accessor
  T       &second()       { return _first->next()->val() ; }  ///< second element accessor
  const T  third () const { return _first->next()->next()->val() ; }  ///< third element accessor
  T       &third ()       { return _first->next()->next()->val() ; }  ///< third element accessor

  /// get the first two elements
  bool firsts( T &t1, T &t2 ) const { if( !empty() ) { t1 = first() ; if( !single() ) { t2 = second() ; return true ; } } return false ; }
  /// get the first three elements
  bool firsts( T &t1, T &t2, T &t3 ) const { if( firsts( t1,t2 ) ) { if( !pair() ) t3 = third() ; return true ; } return false ; }

  const T  top () const { return _first->val() ; }  ///< first element accessor
  T       &top ()       { return _first->val() ; }  ///< first element accessor

//-----------------------------------------------------------------------------
// clear
public  :
  /// destructor
  void clear()
    {
      node *p = _first ;
      while( p != (node*)NULL )
      {
        node *n = p->next() ;
        delete p ;
        p = n ;
      }
      _first = (node*)NULL ;
    }

  /// emptyness test
  inline bool empty() { return _first == (node*)NULL ; }

//-----------------------------------------------------------------------------
// find and replace
public  :
  /// membership test
  const_iterator cfind( const T &val ) const { for( const_iterator it = cbegin() ; it() ; ++it )  if( (*it) == val ) return it ; return cend() ; }

  /// membership test
  iterator find( const T &val ) { for( iterator it = begin() ; it() ; ++it )  if( (*it) == val ) return it ; return end() ; }

  /// membership test (return also the last pointer for eventual remotion)
  iterator find( const T &val, iterator &lst ) { lst = end() ; for( iterator it = begin() ; it() ; ++it ) { if( *it == val ) return it ;  lst = it ; } return end() ; }

  /// membership test
  clsize pos( const T &val ) const { lsize i = 0 ; for( const_iterator it = cbegin() ; it() ; ++it, ++i )  if( *it == val ) return i ; return (clsize)-1 ; }

  /// get the minimal element
  const_iterator min_val() const { if( empty() ) return cend() ; const_iterator it = cbegin(), m = it ; for( ++it ; it() ; ++it ) { if( *it < *m  ) m = it ; } return m ; }

  /// get the minimal element
  iterator min_val() { if( empty() ) return end() ; iterator it = begin(), m = it ; for( ++it ; it() ; ++it ) { if( *it < *m  ) m = it ; } return m ; }

  /// get the maximal element
  const_iterator max_val() const { if( empty() ) return cend() ; const_iterator it = cbegin(), m = it ; for( ++it ; it() ; ++it ) { if( *it > *m  ) m = it ; } return m ; }

  /// get the maximal element
  iterator max_val() { if( empty() ) return end() ; iterator it = begin(), m = it ; for( ++it ; it() ; ++it ) { if( *it > *m  ) m = it ; } return m ; }

  /// replacing the elements of a list
  clsize replace( const T &old_val, const T &new_val )
    {
      lsize n = 0 ;
      for( iterator it = begin() ; it() ; ++it )
        if( *it == old_val ) { *it = new_val ;  ++n ; }
      return n ;
    }

//-----------------------------------------------------------------------------
// element insertion
public  :
  /// insertion at the beginning
  inline void push_front( const T &val ) { _first = new node( val, _first ) ; }

  /// insertion in the middle
  iterator push_after( const T &val, iterator i = begin() )
    {
      if( i == end() ) { push_front( val ) ; return begin() ; }
      i._node->next() = new node( val, i._node->next() ) ;
      return ++i ;
    }

  /// insertion at the end
  iterator push_back( const T &val )
    {
      return push_after( val, last() ) ;
    }

  /// insertion
  inline void insert( const T &val ) { push_front( val ) ; }
  
  /// insertion
  inline void push( const T &val ) { push_front( val ) ; }
  
  /// insertion without duplicates
  inline bool insert_unique( const T &val ) { if( !find(val)() ) {push_front( val ); return true;} return false; }


  /// list insertion (reverse, at the begining)
  void insert( const List<T> &l )
    { for( const_iterator lit = l.cbegin() ; lit() ; ++lit ) insert( *lit ) ; }

  /// list insertion without duplicates (reverse, at the begining)
  void insert_unique( const List<T> &l )
    { for( const_iterator lit = l.cbegin() ; lit() ; ++lit ) insert_unique( *lit ) ; }

  /// list insertion (ordered)
  void insert_ordered( const List<T> &l )
    {
      const_iterator lit = l.cbegin() ;
      iterator it = push_back( *lit ) ;
      for( ; lit() ; ++lit )
        it = push_after( *lit, it ) ;
    }


  /// copy a list and clear its content
  void insert_and_clear( List<T> &l )
    {
      if( this == &l ) return ;
      iterator lst = last() ;
      if( lst() )
        lst._node->next() = l._first ;
      else
        _first = l._first ;
      l._first = NULL ;
    }

//-----------------------------------------------------------------------------
// element remotion
public  :
  /// removing the first element of the list
  void remove_first()
    {
      node *p = _first->next() ;
      delete _first ;
      _first = p ;
    }

  /// removing the first element of the list
  void pop() { remove_first() ; }

  /// removing elements from the list
  void remove_next( iterator &prev )
    {
      node *p ;
      if( prev == end() )
      { p = _first->next() ; delete _first ; _first = p ; }
      else if( prev._node )
      { p = prev._node->next() ; prev._node->next() = p->next() ; delete p ; }
    }

  /// removing elements from the list
  bool remove( T &val )
  {
    iterator prev = end() ;
    for( iterator it = begin() ; it() ; )
    {
      if( *it != val ) { prev = it ; ++it ; continue ; }
      val = *it ;
      node *p ;
      if( prev == end() )
      { p = _first->next() ; delete _first ; _first = p ; it = begin() ; }
      else
      { p = it._node ; prev._node->next() = it._node->next() ; delete p ; it = prev ; ++it ; }
      return true ;
    }
    return false ;
  }

  /// removing elements from the list
  clsize remove_all( const T &val )
    {
      lsize n = 0 ;
      iterator prev = end() ;
      for( iterator it = begin() ; it() ; )
      {
        if( *it != val ) { prev = it ; ++it ; continue ; }
        node *p ;
        if( prev == end() )
        { p = _first->next() ; delete _first ; _first = p ; it = begin() ; }
        else
        { p = it._node ; prev._node->next() = it._node->next() ; delete p ; it = prev ; ++it ; }
        ++n ;
        if( !it() ) return n ;
      }
      return n ;
    }


//-----------------------------------------------------------------------------
// sort
public  :
  /// array transform
  void get_array( T *a, clsize sz /*= size()*/ )
    {
      lsize i = 0 ;
      for( const_iterator it = cbegin() ; it() && i < sz ; ++it, ++i )
        a[i] = *it ;
    }

  /// array transform
  T *get_array( clsize sz /*= size()*/ )
    {
      T  *a = new T[sz] ;
      get_array( a,sz ) ;
      return a ;
    }

  /// array read
  void set_array( clsize sz = 0, const T *a = (const T*)NULL )
    {
      clear() ;
      if( sz == 0 ) return ;
      push_front( a[0] ) ;
      iterator it = begin() ;
      for( lsize i = 1 ; i < sz ; ++i )
        it = push_after( a[i], it ) ;
    }

  /// sort
  void sort( int compare(const void*,const void*) )
    {
      clsize sz = size() ;
      T  *a  = get_array( sz ) ;
      qsort( (void*)a, sz, sizeof(T), compare ) ;
      set_array( sz, a ) ;
      delete [] a ;
    }

//-----------------------------------------------------------------------------
// Elements
private :
  node  *_first; ///< link to the first element



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// node representation
public :
  /// node of a list
  class node
  {
  public :
  // constructors
    /// default constructors
    node( T val_, node *next_ = (node*)NULL ) : _val (val_) , _next (next_) {}

    /// destructor
    ~node() {}

    /// copy constructor
    node( const node &n ) : _val(n.val()), _next(n.next()) {}

    /// assignment operator
    node &operator = ( const node &n )
    { _val = n.val(); _next = n.next(); return this; }

  //-----------------------------------------------------------------------------
  // Accessors
  public  :
    inline const T    &val () const { return _val  ; }  ///< accesses the value of the node
    inline const node *next() const { return _next ; }  ///< accesses the next element

    inline T     &val () { return _val  ; }  ///< accesses the value of the node
    inline node *&next() { return _next ; }  ///< accesses the next element

  //-----------------------------------------------------------------------------
  // Elements
  private :
    T      _val ; ///< value
    node  *_next; ///< link to the next element
  };

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// const_iterator
public :
  /// list const iterator
  class const_iterator
  {
  public :
  // constructors
    /// default constructors
    const_iterator( node *n = (node*)NULL ) : _node(n) {}

    /// destructor
    ~const_iterator() {}

    /// copy constructor
    const_iterator( const const_iterator &i ) : _node(i._node) {}

    /// assignment operator
    const_iterator &operator = ( const const_iterator &i )
    { _node = i._node; return *this; }

    friend class List<T> ;

  //-----------------------------------------------------------------------------
  // Operations
  public  :
    inline bool            operator ==( const const_iterator &i ) const { return _node == i._node ; }  ///< equality operator
    inline bool            operator !=( const const_iterator &i ) const { return _node != i._node ; }  ///< inequality operator
    inline bool            operator ()() const { return _node != (node*)NULL ; }                       ///< validation operator
    inline const T        &operator * () const { return _node->val() ; }                               ///< value accessor
    inline const_iterator &operator ++() { _node = _node->next() ; return *this ; }              ///< accesses the next position

  //-----------------------------------------------------------------------------
  // Elements
  private :
    node  *_node; ///< node position
  };

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// iterator
public :
  /// list iterator
  class iterator
  {
  public :
  // constructors
    /// default constructors
    iterator( node *n = (node*)NULL ) : _node(n) {}

    /// destructor
    ~iterator() {}

    /// copy constructor
    iterator( const iterator &i ) : _node(i._node) {}

    /// assignment operator
    iterator &operator = ( const iterator &i )
    { _node = i._node; return *this; }

    friend class List <T> ;

  //-----------------------------------------------------------------------------
  // Operations
  public  :
    inline bool      operator ==( const iterator &i ) const { return _node == i._node ; }  ///< equality operator
    inline bool      operator !=( const iterator &i ) const { return _node != i._node ; }  ///< inequality operator
    inline bool      operator ()() const { return _node != (node*)NULL ; }                 ///< validation operator
    inline T        &operator * () const { return _node->val() ; }                         ///< value accessor
    inline iterator &operator ++() { _node = _node->next() ; return *this ; }        ///< accesses the next position

  //-----------------------------------------------------------------------------
  // Elements
  private :
    node  *_node; ///< node position
  };

//-----------------------------------------------------------------------------
// iterator creation
public  :
  inline const_iterator cbegin() const { return const_iterator( _first ) ; }      ///< iterator creation
  inline const_iterator cend  () const { return const_iterator( (node*)NULL ) ; }  ///< iterator end
  inline iterator begin() { return iterator( _first ) ; }      ///< iterator creation
  inline iterator end  () { return iterator( (node*)NULL ) ; }  ///< iterator end

  /// last position
  iterator last () { iterator lst = end(), it = begin() ; while( it() ) { lst = it ; ++it ; } return lst ; }
};
//_____________________________________________________________________________




//_____________________________________________________________________________
/// template list class
template <typename T> class Queue : public List<T>
//-----------------------------------------------------------------------------
{
//-----------------------------------------------------------------------------
// default constructors
public  :
  /// default constructors
  Queue() : List<T>() { _last = this->end() ; }

  /// one-element constructor
  Queue( const T &val ) : List<T>( val ) { _last = this->begin() ; }

  /// copy constructor (reverse order)
  Queue( const Queue<T> &q ) : List<T>()
    { for( typename List<T>::const_iterator it = q.cbegin() ; it() ; ++it ) Queue<T>::insert( *it ) ; }

  /// assignment operator
  Queue<T> &operator = ( const Queue<T> &q )
   { if( this != &q) { this->clear() ;  for( typename List<T>::const_iterator it = q.cbegin() ; it() ; ++it ) insert( *it ) ; }  return *this; }

  /// destructor
  ~Queue() { this->clear() ; }


//-----------------------------------------------------------------------------
// element insertion
public  :
  /// insertion at the end
  inline void insert( const T &val )
    { _last = push_after( val, _last ) ; }

  /// queue insertion
  void insert( const Queue<T> &l )
    { for( typename List<T>::const_iterator lit = l.cbegin() ; lit() ; ++lit ) insert( *lit ) ; }

//-----------------------------------------------------------------------------
// stqndqrd function
public  :
  /// first element
  inline const T front() const { return this->first() ; }

  /// insertion at the end
  inline void push( const T &val ) { this->insert( val ) ; }

  /// deletion at the beginning
  inline void pop () { this->remove_first() ;  if( this->empty() ) _last = this->end() ; }


//-----------------------------------------------------------------------------
// Elements
private :
  typename List<T>::iterator _last; ///< link to the last element
};
//_____________________________________________________________________________


