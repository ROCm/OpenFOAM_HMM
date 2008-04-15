/*
 *  MICO --- an Open Source CORBA implementation
 *  Copyright (c) 1997-2006 by The Mico Team
 *
 *  This file was automatically generated. DO NOT EDIT!
 */

#include <FoamX.H>


using namespace std;

//--------------------------------------------------------
//  Implementation of stubs
//--------------------------------------------------------
namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXType;
}

void operator<<=( CORBA::Any &_a, const FoamXServer::FoamXType &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXType, &_e);
  _a.from_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::FoamXType &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXType, &_e);
  return _a.to_static_any (_sa);
}

class _Marshaller_FoamXServer_FoamXType : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::FoamXType _MICO_T;
  public:
    ~_Marshaller_FoamXServer_FoamXType();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_FoamXType::~_Marshaller_FoamXServer_FoamXType()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_FoamXType::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_FoamXType::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_FoamXType::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_FoamXType::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong ul;
  if( !dc.enumeration( ul ) )
    return FALSE;
  *(_MICO_T*) v = (_MICO_T) ul;
  return TRUE;
}

void _Marshaller_FoamXServer_FoamXType::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.enumeration( (::CORBA::ULong) *(_MICO_T *) v );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_FoamXType::typecode()
{
  return FoamXServer::_tc_FoamXType;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXType;

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXAny;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXAny::FoamXAny()
{
}

FoamXServer::FoamXAny::FoamXAny( const FoamXAny& _s )
{
  type = ((FoamXAny&)_s).type;
  value = ((FoamXAny&)_s).value;
}

FoamXServer::FoamXAny::~FoamXAny()
{
}

FoamXServer::FoamXAny&
FoamXServer::FoamXAny::operator=( const FoamXAny& _s )
{
  type = ((FoamXAny&)_s).type;
  value = ((FoamXAny&)_s).value;
  return *this;
}
#endif

class _Marshaller_FoamXServer_FoamXAny : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::FoamXAny _MICO_T;
  public:
    ~_Marshaller_FoamXServer_FoamXAny();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_FoamXAny::~_Marshaller_FoamXServer_FoamXAny()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_FoamXAny::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_FoamXAny::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_FoamXAny::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_FoamXAny::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    _marshaller_FoamXServer_FoamXType->demarshal( dc, &((_MICO_T*)v)->type ) &&
    CORBA::_stc_any->demarshal( dc, &((_MICO_T*)v)->value ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_FoamXAny::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  _marshaller_FoamXServer_FoamXType->marshal( ec, &((_MICO_T*)v)->type );
  CORBA::_stc_any->marshal( ec, &((_MICO_T*)v)->value );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_FoamXAny::typecode()
{
  return FoamXServer::_tc_FoamXAny;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXAny;

void operator<<=( CORBA::Any &_a, const FoamXServer::FoamXAny &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXAny, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::FoamXAny *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::FoamXAny &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXAny, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::FoamXAny *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_FoamXAny, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXAnyList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_StringList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_TypeDescriptorList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_ErrorCode;
}

void operator<<=( CORBA::Any &_a, const FoamXServer::ErrorCode &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ErrorCode, &_e);
  _a.from_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::ErrorCode &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ErrorCode, &_e);
  return _a.to_static_any (_sa);
}

class _Marshaller_FoamXServer_ErrorCode : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::ErrorCode _MICO_T;
  public:
    ~_Marshaller_FoamXServer_ErrorCode();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_ErrorCode::~_Marshaller_FoamXServer_ErrorCode()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_ErrorCode::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_ErrorCode::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_ErrorCode::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_ErrorCode::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong ul;
  if( !dc.enumeration( ul ) )
    return FALSE;
  *(_MICO_T*) v = (_MICO_T) ul;
  return TRUE;
}

void _Marshaller_FoamXServer_ErrorCode::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.enumeration( (::CORBA::ULong) *(_MICO_T *) v );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_ErrorCode::typecode()
{
  return FoamXServer::_tc_ErrorCode;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_ErrorCode;

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXError;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXError::FoamXError()
{
}

FoamXServer::FoamXError::FoamXError( const FoamXError& _s )
{
  errorCode = ((FoamXError&)_s).errorCode;
  errorMessage = ((FoamXError&)_s).errorMessage;
  methodName = ((FoamXError&)_s).methodName;
  fileName = ((FoamXError&)_s).fileName;
  lineNo = ((FoamXError&)_s).lineNo;
}

FoamXServer::FoamXError::~FoamXError()
{
}

FoamXServer::FoamXError&
FoamXServer::FoamXError::operator=( const FoamXError& _s )
{
  errorCode = ((FoamXError&)_s).errorCode;
  errorMessage = ((FoamXError&)_s).errorMessage;
  methodName = ((FoamXError&)_s).methodName;
  fileName = ((FoamXError&)_s).fileName;
  lineNo = ((FoamXError&)_s).lineNo;
  return *this;
}
#endif

#ifndef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXError::FoamXError()
{
}

#endif

FoamXServer::FoamXError::FoamXError( FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2, const char* _m3, CORBA::Long _m4 )
{
  errorCode = _m0;
  errorMessage = _m1;
  methodName = _m2;
  fileName = _m3;
  lineNo = _m4;
}

class _Marshaller_FoamXServer_FoamXError : public ::CORBA::StaticTypeInfo {
    typedef ::FoamXServer::FoamXError _MICO_T;
  public:
    ~_Marshaller_FoamXServer_FoamXError();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_FoamXError::~_Marshaller_FoamXServer_FoamXError()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_FoamXError::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_FoamXError::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_FoamXError::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_FoamXError::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  string repoid;
  return
    dc.except_begin( repoid ) &&
    _marshaller_FoamXServer_ErrorCode->demarshal( dc, &((_MICO_T*)v)->errorCode ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->errorMessage._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->methodName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->fileName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->lineNo ) &&
    dc.except_end();
}

void _Marshaller_FoamXServer_FoamXError::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.except_begin( "IDL:FoamXServer/FoamXError:1.0" );
  _marshaller_FoamXServer_ErrorCode->marshal( ec, &((_MICO_T*)v)->errorCode );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->errorMessage.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->methodName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->fileName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->lineNo );
  ec.except_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_FoamXError::typecode()
{
  return FoamXServer::_tc_FoamXError;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXError;

void operator<<=( CORBA::Any &_a, const FoamXServer::FoamXError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXError, &_e);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::FoamXError *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::FoamXError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXError, &_e);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::FoamXError *&_e )
{
  return _a.to_static_any (_marshaller_FoamXServer_FoamXError, (void *&)_e);
}

void FoamXServer::FoamXError::_throwit() const
{
  #ifdef HAVE_EXCEPTIONS
  #ifdef HAVE_STD_EH
  throw *this;
  #else
  throw FoamXError_var( (FoamXServer::FoamXError*)_clone() );
  #endif
  #else
  CORBA::Exception::_throw_failed( _clone() );
  #endif
}

const char *FoamXServer::FoamXError::_repoid() const
{
  return "IDL:FoamXServer/FoamXError:1.0";
}

void FoamXServer::FoamXError::_encode( CORBA::DataEncoder &_en ) const
{
  _marshaller_FoamXServer_FoamXError->marshal( _en, (void*) this );
}

void FoamXServer::FoamXError::_encode_any( CORBA::Any &_a ) const
{
  _a <<= *this;
}

CORBA::Exception *FoamXServer::FoamXError::_clone() const
{
  return new FoamXError( *this );
}

FoamXServer::FoamXError *FoamXServer::FoamXError::_downcast( CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXError:1.0" ) )
    return (FoamXError *) _ex;
  return NULL;
}

const FoamXServer::FoamXError *FoamXServer::FoamXError::_downcast( const CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXError:1.0" ) )
    return (FoamXError *) _ex;
  return NULL;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXIOError;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXIOError::FoamXIOError()
{
}

FoamXServer::FoamXIOError::FoamXIOError( const FoamXIOError& _s )
{
  errorMessage = ((FoamXIOError&)_s).errorMessage;
  ioFileName = ((FoamXIOError&)_s).ioFileName;
  ioStartLineNo = ((FoamXIOError&)_s).ioStartLineNo;
  ioEndLineNo = ((FoamXIOError&)_s).ioEndLineNo;
  methodName = ((FoamXIOError&)_s).methodName;
  fileName = ((FoamXIOError&)_s).fileName;
  lineNo = ((FoamXIOError&)_s).lineNo;
}

FoamXServer::FoamXIOError::~FoamXIOError()
{
}

FoamXServer::FoamXIOError&
FoamXServer::FoamXIOError::operator=( const FoamXIOError& _s )
{
  errorMessage = ((FoamXIOError&)_s).errorMessage;
  ioFileName = ((FoamXIOError&)_s).ioFileName;
  ioStartLineNo = ((FoamXIOError&)_s).ioStartLineNo;
  ioEndLineNo = ((FoamXIOError&)_s).ioEndLineNo;
  methodName = ((FoamXIOError&)_s).methodName;
  fileName = ((FoamXIOError&)_s).fileName;
  lineNo = ((FoamXIOError&)_s).lineNo;
  return *this;
}
#endif

#ifndef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXIOError::FoamXIOError()
{
}

#endif

FoamXServer::FoamXIOError::FoamXIOError( const char* _m0, const char* _m1, CORBA::Long _m2, CORBA::Long _m3, const char* _m4, const char* _m5, CORBA::Long _m6 )
{
  errorMessage = _m0;
  ioFileName = _m1;
  ioStartLineNo = _m2;
  ioEndLineNo = _m3;
  methodName = _m4;
  fileName = _m5;
  lineNo = _m6;
}

class _Marshaller_FoamXServer_FoamXIOError : public ::CORBA::StaticTypeInfo {
    typedef ::FoamXServer::FoamXIOError _MICO_T;
  public:
    ~_Marshaller_FoamXServer_FoamXIOError();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_FoamXIOError::~_Marshaller_FoamXServer_FoamXIOError()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_FoamXIOError::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_FoamXIOError::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_FoamXIOError::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_FoamXIOError::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  string repoid;
  return
    dc.except_begin( repoid ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->errorMessage._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->ioFileName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->ioStartLineNo ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->ioEndLineNo ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->methodName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->fileName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->lineNo ) &&
    dc.except_end();
}

void _Marshaller_FoamXServer_FoamXIOError::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.except_begin( "IDL:FoamXServer/FoamXIOError:1.0" );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->errorMessage.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->ioFileName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->ioStartLineNo );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->ioEndLineNo );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->methodName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->fileName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->lineNo );
  ec.except_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_FoamXIOError::typecode()
{
  return FoamXServer::_tc_FoamXIOError;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXIOError;

void operator<<=( CORBA::Any &_a, const FoamXServer::FoamXIOError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXIOError, &_e);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::FoamXIOError *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::FoamXIOError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXIOError, &_e);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::FoamXIOError *&_e )
{
  return _a.to_static_any (_marshaller_FoamXServer_FoamXIOError, (void *&)_e);
}

void FoamXServer::FoamXIOError::_throwit() const
{
  #ifdef HAVE_EXCEPTIONS
  #ifdef HAVE_STD_EH
  throw *this;
  #else
  throw FoamXIOError_var( (FoamXServer::FoamXIOError*)_clone() );
  #endif
  #else
  CORBA::Exception::_throw_failed( _clone() );
  #endif
}

const char *FoamXServer::FoamXIOError::_repoid() const
{
  return "IDL:FoamXServer/FoamXIOError:1.0";
}

void FoamXServer::FoamXIOError::_encode( CORBA::DataEncoder &_en ) const
{
  _marshaller_FoamXServer_FoamXIOError->marshal( _en, (void*) this );
}

void FoamXServer::FoamXIOError::_encode_any( CORBA::Any &_a ) const
{
  _a <<= *this;
}

CORBA::Exception *FoamXServer::FoamXIOError::_clone() const
{
  return new FoamXIOError( *this );
}

FoamXServer::FoamXIOError *FoamXServer::FoamXIOError::_downcast( CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXIOError:1.0" ) )
    return (FoamXIOError *) _ex;
  return NULL;
}

const FoamXServer::FoamXIOError *FoamXServer::FoamXIOError::_downcast( const CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXIOError:1.0" ) )
    return (FoamXIOError *) _ex;
  return NULL;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_ValidationError;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::ValidationError::ValidationError()
{
}

FoamXServer::ValidationError::ValidationError( const ValidationError& _s )
{
  errorCode = ((ValidationError&)_s).errorCode;
  errorMessage = ((ValidationError&)_s).errorMessage;
  itemPath = ((ValidationError&)_s).itemPath;
}

FoamXServer::ValidationError::~ValidationError()
{
}

FoamXServer::ValidationError&
FoamXServer::ValidationError::operator=( const ValidationError& _s )
{
  errorCode = ((ValidationError&)_s).errorCode;
  errorMessage = ((ValidationError&)_s).errorMessage;
  itemPath = ((ValidationError&)_s).itemPath;
  return *this;
}
#endif

#ifndef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::ValidationError::ValidationError()
{
}

#endif

FoamXServer::ValidationError::ValidationError( FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2 )
{
  errorCode = _m0;
  errorMessage = _m1;
  itemPath = _m2;
}

class _Marshaller_FoamXServer_ValidationError : public ::CORBA::StaticTypeInfo {
    typedef ::FoamXServer::ValidationError _MICO_T;
  public:
    ~_Marshaller_FoamXServer_ValidationError();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_ValidationError::~_Marshaller_FoamXServer_ValidationError()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_ValidationError::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_ValidationError::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_ValidationError::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_ValidationError::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  string repoid;
  return
    dc.except_begin( repoid ) &&
    _marshaller_FoamXServer_ErrorCode->demarshal( dc, &((_MICO_T*)v)->errorCode ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->errorMessage._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->itemPath._for_demarshal() ) &&
    dc.except_end();
}

void _Marshaller_FoamXServer_ValidationError::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.except_begin( "IDL:FoamXServer/ValidationError:1.0" );
  _marshaller_FoamXServer_ErrorCode->marshal( ec, &((_MICO_T*)v)->errorCode );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->errorMessage.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->itemPath.inout() );
  ec.except_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_ValidationError::typecode()
{
  return FoamXServer::_tc_ValidationError;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_ValidationError;

void operator<<=( CORBA::Any &_a, const FoamXServer::ValidationError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ValidationError, &_e);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::ValidationError *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::ValidationError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ValidationError, &_e);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::ValidationError *&_e )
{
  return _a.to_static_any (_marshaller_FoamXServer_ValidationError, (void *&)_e);
}

void FoamXServer::ValidationError::_throwit() const
{
  #ifdef HAVE_EXCEPTIONS
  #ifdef HAVE_STD_EH
  throw *this;
  #else
  throw ValidationError_var( (FoamXServer::ValidationError*)_clone() );
  #endif
  #else
  CORBA::Exception::_throw_failed( _clone() );
  #endif
}

const char *FoamXServer::ValidationError::_repoid() const
{
  return "IDL:FoamXServer/ValidationError:1.0";
}

void FoamXServer::ValidationError::_encode( CORBA::DataEncoder &_en ) const
{
  _marshaller_FoamXServer_ValidationError->marshal( _en, (void*) this );
}

void FoamXServer::ValidationError::_encode_any( CORBA::Any &_a ) const
{
  _a <<= *this;
}

CORBA::Exception *FoamXServer::ValidationError::_clone() const
{
  return new ValidationError( *this );
}

FoamXServer::ValidationError *FoamXServer::ValidationError::_downcast( CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/ValidationError:1.0" ) )
    return (ValidationError *) _ex;
  return NULL;
}

const FoamXServer::ValidationError *FoamXServer::ValidationError::_downcast( const CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/ValidationError:1.0" ) )
    return (ValidationError *) _ex;
  return NULL;
}


/*
 * Base interface for class ITypeDescriptor
 */

FoamXServer::ITypeDescriptor::~ITypeDescriptor()
{
}

void *
FoamXServer::ITypeDescriptor::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/ITypeDescriptor:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::ITypeDescriptor_ptr
FoamXServer::ITypeDescriptor::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::ITypeDescriptor_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/ITypeDescriptor:1.0" )))
      return _duplicate( (FoamXServer::ITypeDescriptor_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/ITypeDescriptor:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/ITypeDescriptor:1.0")) {
      _o = new FoamXServer::ITypeDescriptor_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::ITypeDescriptor_ptr
FoamXServer::ITypeDescriptor::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_ITypeDescriptor;
}
class _Marshaller_FoamXServer_ITypeDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::ITypeDescriptor_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_ITypeDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_ITypeDescriptor::~_Marshaller_FoamXServer_ITypeDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_ITypeDescriptor::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_ITypeDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::ITypeDescriptor::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_ITypeDescriptor::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_ITypeDescriptor::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_ITypeDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::ITypeDescriptor::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_ITypeDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_ITypeDescriptor::typecode()
{
  return FoamXServer::_tc_ITypeDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_ITypeDescriptor;

void
operator<<=( CORBA::Any &_a, const FoamXServer::ITypeDescriptor_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ITypeDescriptor, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::ITypeDescriptor_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ITypeDescriptor, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::ITypeDescriptor_ptr &_obj )
{
  FoamXServer::ITypeDescriptor_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_ITypeDescriptor, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class ITypeDescriptor
 */

FoamXServer::ITypeDescriptor_stub::~ITypeDescriptor_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::ITypeDescriptor::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/ITypeDescriptor:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::ITypeDescriptor *
POA_FoamXServer::ITypeDescriptor::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/ITypeDescriptor:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::ITypeDescriptor *) p;
  }
  return NULL;
}

FoamXServer::ITypeDescriptor_stub_clp::ITypeDescriptor_stub_clp ()
{
}

FoamXServer::ITypeDescriptor_stub_clp::ITypeDescriptor_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::ITypeDescriptor_stub_clp::~ITypeDescriptor_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::FoamXType FoamXServer::ITypeDescriptor_stub::type()
{
  FoamXServer::FoamXType _res;
  CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXType, &_res );

  CORBA::StaticRequest __req( this, "_get_type" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::FoamXType
FoamXServer::ITypeDescriptor_stub_clp::type()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::FoamXType __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->type();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::type();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::type( FoamXServer::FoamXType _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXType, &_par__value );
  CORBA::StaticRequest __req( this, "_set_type" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::type( FoamXServer::FoamXType _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->type(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::type(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::isPrimitiveType()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_isPrimitiveType" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::isPrimitiveType()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->isPrimitiveType();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::isPrimitiveType();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::isCompoundType()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_isCompoundType" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::isCompoundType()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->isCompoundType();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::isCompoundType();
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::path()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_path" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::path()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->path();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::path();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::path( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_path" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::path( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->path(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::path(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::displayName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_displayName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::displayName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->displayName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::displayName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::displayName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_displayName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::displayName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->displayName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::displayName(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::comment()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_comment" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::comment()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->comment();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::comment();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::comment( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_comment" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::comment( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->comment(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::comment(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::category()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_category" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::category()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->category();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::category();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::category( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_category" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::category( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->category(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::category(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::helpURL()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_helpURL" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::helpURL()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->helpURL();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::helpURL();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::helpURL( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_helpURL" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::helpURL( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->helpURL(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::helpURL(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::iconURL()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_iconURL" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::iconURL()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->iconURL();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::iconURL();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::iconURL( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_iconURL" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::iconURL( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->iconURL(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::iconURL(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::optional()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_optional" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::optional()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->optional();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::optional();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::optional( CORBA::Boolean _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );
  CORBA::StaticRequest __req( this, "_set_optional" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::optional( CORBA::Boolean _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->optional(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::optional(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::visible()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_visible" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::visible()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->visible();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::visible();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::visible( CORBA::Boolean _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );
  CORBA::StaticRequest __req( this, "_set_visible" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::visible( CORBA::Boolean _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->visible(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::visible(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::editable()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_editable" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::editable()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->editable();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::editable();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::editable( CORBA::Boolean _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );
  CORBA::StaticRequest __req( this, "_set_editable" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::editable( CORBA::Boolean _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->editable(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::editable(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::FoamXAny* FoamXServer::ITypeDescriptor_stub::minValue()
{
  CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );

  CORBA::StaticRequest __req( this, "_get_minValue" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::FoamXAny*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::FoamXAny*
FoamXServer::ITypeDescriptor_stub_clp::minValue()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::FoamXAny* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->minValue();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::minValue();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::minValue( const FoamXServer::FoamXAny& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );
  CORBA::StaticRequest __req( this, "_set_minValue" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::minValue( const FoamXServer::FoamXAny& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->minValue(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::minValue(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::FoamXAny* FoamXServer::ITypeDescriptor_stub::maxValue()
{
  CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );

  CORBA::StaticRequest __req( this, "_get_maxValue" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::FoamXAny*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::FoamXAny*
FoamXServer::ITypeDescriptor_stub_clp::maxValue()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::FoamXAny* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->maxValue();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::maxValue();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::maxValue( const FoamXServer::FoamXAny& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );
  CORBA::StaticRequest __req( this, "_set_maxValue" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::maxValue( const FoamXServer::FoamXAny& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->maxValue(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::maxValue(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::lookupDict()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_lookupDict" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::lookupDict()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->lookupDict();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::lookupDict();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::lookupDict( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_lookupDict" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::lookupDict( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->lookupDict(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::lookupDict(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::FoamXAnyList* FoamXServer::ITypeDescriptor_stub::valueList()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_FoamXAny );

  CORBA::StaticRequest __req( this, "_get_valueList" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::FoamXAnyList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::FoamXAnyList*
FoamXServer::ITypeDescriptor_stub_clp::valueList()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::FoamXAnyList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->valueList();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::valueList();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::valueList( const FoamXServer::FoamXAnyList& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller__seq_FoamXServer_FoamXAny, &_par__value );
  CORBA::StaticRequest __req( this, "_set_valueList" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::valueList( const FoamXServer::FoamXAnyList& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->valueList(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::valueList(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::ITypeDescriptor_stub::dictionaryPath()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_dictionaryPath" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::ITypeDescriptor_stub_clp::dictionaryPath()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->dictionaryPath();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::dictionaryPath();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::dictionaryPath( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_dictionaryPath" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::dictionaryPath( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->dictionaryPath(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::dictionaryPath(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::ITypeDescriptor_stub::numElements()
{
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "_get_numElements" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::ITypeDescriptor_stub_clp::numElements()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->numElements();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::numElements();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::numElements( CORBA::Long _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_long, &_par__value );
  CORBA::StaticRequest __req( this, "_set_numElements" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::numElements( CORBA::Long _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->numElements(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::numElements(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::ITypeDescriptor_stub::elementLabels()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_elementLabels" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::ITypeDescriptor_stub_clp::elementLabels()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->elementLabels();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::elementLabels();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::elementLabels( const FoamXServer::StringList& _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stcseq_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_elementLabels" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::elementLabels( const FoamXServer::StringList& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->elementLabels(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::elementLabels(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::TypeDescriptorList* FoamXServer::ITypeDescriptor_stub::subTypes()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ITypeDescriptor );

  CORBA::StaticRequest __req( this, "_get_subTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::TypeDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::TypeDescriptorList*
FoamXServer::ITypeDescriptor_stub_clp::subTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::TypeDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->subTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::subTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr FoamXServer::ITypeDescriptor_stub::elementType()
{
  FoamXServer::ITypeDescriptor_ptr _res = FoamXServer::ITypeDescriptor::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );

  CORBA::StaticRequest __req( this, "_get_elementType" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr
FoamXServer::ITypeDescriptor_stub_clp::elementType()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ITypeDescriptor_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->elementType();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::elementType();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::ITypeDescriptor_stub::hasDefaultValue()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "hasDefaultValue" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::ITypeDescriptor_stub_clp::hasDefaultValue()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->hasDefaultValue();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::ITypeDescriptor_stub::hasDefaultValue();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::getDefaultValue( FoamXServer::IDictionaryEntry_out _par_defaultValue )
{
  CORBA::StaticAny _sa_defaultValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_defaultValue.ptr() );
  CORBA::StaticRequest __req( this, "getDefaultValue" );
  __req.add_out_arg( &_sa_defaultValue );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::getDefaultValue( FoamXServer::IDictionaryEntry_out _par_defaultValue )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getDefaultValue(_par_defaultValue);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::getDefaultValue(_par_defaultValue);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::addSubType( FoamXServer::FoamXType _par_type, FoamXServer::ITypeDescriptor_out _par_subEntry )
{
  CORBA::StaticAny _sa_type( _marshaller_FoamXServer_FoamXType, &_par_type );
  CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_ITypeDescriptor, &_par_subEntry.ptr() );
  CORBA::StaticRequest __req( this, "addSubType" );
  __req.add_in_arg( &_sa_type );
  __req.add_out_arg( &_sa_subEntry );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::addSubType( FoamXServer::FoamXType _par_type, FoamXServer::ITypeDescriptor_out _par_subEntry )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addSubType(_par_type, _par_subEntry);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::addSubType(_par_type, _par_subEntry);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::removeSubType( FoamXServer::ITypeDescriptor_ptr _par_subEntry )
{
  CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_ITypeDescriptor, &_par_subEntry );
  CORBA::StaticRequest __req( this, "removeSubType" );
  __req.add_in_arg( &_sa_subEntry );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::removeSubType( FoamXServer::ITypeDescriptor_ptr _par_subEntry )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->removeSubType(_par_subEntry);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::removeSubType(_par_subEntry);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::ITypeDescriptor_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::ITypeDescriptor_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::ITypeDescriptor * _myserv = POA_FoamXServer::ITypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::ITypeDescriptor_stub::validate();
}

#endif // MICO_CONF_NO_POA

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_DictionaryEntryList;
}


/*
 * Base interface for class IDictionaryEntry
 */

FoamXServer::IDictionaryEntry::~IDictionaryEntry()
{
}

void *
FoamXServer::IDictionaryEntry::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/IDictionaryEntry:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::IDictionaryEntry_ptr
FoamXServer::IDictionaryEntry::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::IDictionaryEntry_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/IDictionaryEntry:1.0" )))
      return _duplicate( (FoamXServer::IDictionaryEntry_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/IDictionaryEntry:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/IDictionaryEntry:1.0")) {
      _o = new FoamXServer::IDictionaryEntry_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::IDictionaryEntry_ptr
FoamXServer::IDictionaryEntry::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_IDictionaryEntry;
}
class _Marshaller_FoamXServer_IDictionaryEntry : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::IDictionaryEntry_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_IDictionaryEntry();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_IDictionaryEntry::~_Marshaller_FoamXServer_IDictionaryEntry()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_IDictionaryEntry::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_IDictionaryEntry::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::IDictionaryEntry::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_IDictionaryEntry::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_IDictionaryEntry::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_IDictionaryEntry::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::IDictionaryEntry::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_IDictionaryEntry::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_IDictionaryEntry::typecode()
{
  return FoamXServer::_tc_IDictionaryEntry;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_IDictionaryEntry;

void
operator<<=( CORBA::Any &_a, const FoamXServer::IDictionaryEntry_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_IDictionaryEntry, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::IDictionaryEntry_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_IDictionaryEntry, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::IDictionaryEntry_ptr &_obj )
{
  FoamXServer::IDictionaryEntry_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_IDictionaryEntry, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IDictionaryEntry
 */

FoamXServer::IDictionaryEntry_stub::~IDictionaryEntry_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::IDictionaryEntry::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/IDictionaryEntry:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::IDictionaryEntry *
POA_FoamXServer::IDictionaryEntry::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/IDictionaryEntry:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::IDictionaryEntry *) p;
  }
  return NULL;
}

FoamXServer::IDictionaryEntry_stub_clp::IDictionaryEntry_stub_clp ()
{
}

FoamXServer::IDictionaryEntry_stub_clp::IDictionaryEntry_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::IDictionaryEntry_stub_clp::~IDictionaryEntry_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr FoamXServer::IDictionaryEntry_stub::typeDescriptor()
{
  FoamXServer::ITypeDescriptor_ptr _res = FoamXServer::ITypeDescriptor::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );

  CORBA::StaticRequest __req( this, "_get_typeDescriptor" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr
FoamXServer::IDictionaryEntry_stub_clp::typeDescriptor()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ITypeDescriptor_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->typeDescriptor();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::typeDescriptor();
}

#endif // MICO_CONF_NO_POA

FoamXServer::FoamXAny* FoamXServer::IDictionaryEntry_stub::value()
{
  CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );

  CORBA::StaticRequest __req( this, "_get_value" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::FoamXAny*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::FoamXAny*
FoamXServer::IDictionaryEntry_stub_clp::value()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      FoamXServer::FoamXAny* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->value();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::value();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::value( const FoamXServer::FoamXAny& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );
  CORBA::StaticRequest __req( this, "_set_value" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::value( const FoamXServer::FoamXAny& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->value(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::value(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::DictionaryEntryList* FoamXServer::IDictionaryEntry_stub::subElements()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_IDictionaryEntry );

  CORBA::StaticRequest __req( this, "_get_subElements" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::DictionaryEntryList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::DictionaryEntryList*
FoamXServer::IDictionaryEntry_stub_clp::subElements()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      FoamXServer::DictionaryEntryList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->subElements();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::subElements();
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::IDictionaryEntry_stub::selection()
{
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "_get_selection" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::IDictionaryEntry_stub_clp::selection()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->selection();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::selection();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::selection( CORBA::Long _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_long, &_par__value );
  CORBA::StaticRequest __req( this, "_set_selection" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::selection( CORBA::Long _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->selection(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::selection(_par__value);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::setValue( const FoamXServer::FoamXAny& _par_value )
{
  CORBA::StaticAny _sa_value( _marshaller_FoamXServer_FoamXAny, &_par_value );
  CORBA::StaticRequest __req( this, "setValue" );
  __req.add_in_arg( &_sa_value );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::setValue( const FoamXServer::FoamXAny& _par_value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->setValue(_par_value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::setValue(_par_value);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::IDictionaryEntry_stub::nSubElements()
{
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "nSubElements" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::IDictionaryEntry_stub_clp::nSubElements()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->nSubElements();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::nSubElements();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::IDictionaryEntry_stub::packedList()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "packedList" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::IDictionaryEntry_stub_clp::packedList()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->packedList();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::packedList();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::addElement( FoamXServer::IDictionaryEntry_out _par_subEntry )
{
  CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_IDictionaryEntry, &_par_subEntry.ptr() );
  CORBA::StaticRequest __req( this, "addElement" );
  __req.add_out_arg( &_sa_subEntry );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::addElement( FoamXServer::IDictionaryEntry_out _par_subEntry )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addElement(_par_subEntry);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::addElement(_par_subEntry);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::removeElement( FoamXServer::IDictionaryEntry_ptr _par_subEntry )
{
  CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_IDictionaryEntry, &_par_subEntry );
  CORBA::StaticRequest __req( this, "removeElement" );
  __req.add_in_arg( &_sa_subEntry );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::removeElement( FoamXServer::IDictionaryEntry_ptr _par_subEntry )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->removeElement(_par_subEntry);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::removeElement(_par_subEntry);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::validate();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::IDictionaryEntry_stub::modified()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "modified" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::IDictionaryEntry_stub_clp::modified()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->modified();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::IDictionaryEntry_stub::modified();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::IDictionaryEntry_stub::save()
{
  CORBA::StaticRequest __req( this, "save" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::IDictionaryEntry_stub_clp::save()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::IDictionaryEntry * _myserv = POA_FoamXServer::IDictionaryEntry::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->save();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::IDictionaryEntry_stub::save();
}

#endif // MICO_CONF_NO_POA

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_JobStatus;
}

void operator<<=( CORBA::Any &_a, const FoamXServer::JobStatus &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobStatus, &_e);
  _a.from_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::JobStatus &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobStatus, &_e);
  return _a.to_static_any (_sa);
}

class _Marshaller_FoamXServer_JobStatus : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::JobStatus _MICO_T;
  public:
    ~_Marshaller_FoamXServer_JobStatus();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_JobStatus::~_Marshaller_FoamXServer_JobStatus()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_JobStatus::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_JobStatus::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_JobStatus::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_JobStatus::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong ul;
  if( !dc.enumeration( ul ) )
    return FALSE;
  *(_MICO_T*) v = (_MICO_T) ul;
  return TRUE;
}

void _Marshaller_FoamXServer_JobStatus::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.enumeration( (::CORBA::ULong) *(_MICO_T *) v );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_JobStatus::typecode()
{
  return FoamXServer::_tc_JobStatus;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobStatus;

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_MessageType;
}

void operator<<=( CORBA::Any &_a, const FoamXServer::MessageType &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_MessageType, &_e);
  _a.from_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::MessageType &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_MessageType, &_e);
  return _a.to_static_any (_sa);
}

class _Marshaller_FoamXServer_MessageType : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::MessageType _MICO_T;
  public:
    ~_Marshaller_FoamXServer_MessageType();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_MessageType::~_Marshaller_FoamXServer_MessageType()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_MessageType::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_MessageType::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_MessageType::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_MessageType::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong ul;
  if( !dc.enumeration( ul ) )
    return FALSE;
  *(_MICO_T*) v = (_MICO_T) ul;
  return TRUE;
}

void _Marshaller_FoamXServer_MessageType::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.enumeration( (::CORBA::ULong) *(_MICO_T *) v );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_MessageType::typecode()
{
  return FoamXServer::_tc_MessageType;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_MessageType;

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_DoubleList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FloatList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_LongList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_Point3;
}

void operator<<=( CORBA::Any &_a, const FoamXServer::Point3_forany &_const_s )
{
  _const_s.alloc();
  CORBA::StaticAny _sa (_marshaller__a3_float, _const_s.in());
  _a.from_static_any (_sa);
  _const_s.clear();
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::Point3_forany &_s )
{
  FoamXServer::Point3_slice *p;
  if (_a.to_static_any (_marshaller__a3_float, (void *&)p)) {
    _s = p;
    return TRUE;
  }
  return FALSE;
}


namespace FoamXServer
{
CORBA::TypeCodeConst _tc_StringPair;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::StringPair::StringPair()
{
}

FoamXServer::StringPair::StringPair( const StringPair& _s )
{
  name = ((StringPair&)_s).name;
  value = ((StringPair&)_s).value;
}

FoamXServer::StringPair::~StringPair()
{
}

FoamXServer::StringPair&
FoamXServer::StringPair::operator=( const StringPair& _s )
{
  name = ((StringPair&)_s).name;
  value = ((StringPair&)_s).value;
  return *this;
}
#endif

class _Marshaller_FoamXServer_StringPair : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::StringPair _MICO_T;
  public:
    ~_Marshaller_FoamXServer_StringPair();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_StringPair::~_Marshaller_FoamXServer_StringPair()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_StringPair::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_StringPair::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_StringPair::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_StringPair::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->name._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->value._for_demarshal() ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_StringPair::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->name.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->value.inout() );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_StringPair::typecode()
{
  return FoamXServer::_tc_StringPair;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_StringPair;

void operator<<=( CORBA::Any &_a, const FoamXServer::StringPair &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_StringPair, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::StringPair *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::StringPair &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_StringPair, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::StringPair *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_StringPair, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_StringPairList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_HostDescriptor;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::HostDescriptor::HostDescriptor()
{
}

FoamXServer::HostDescriptor::HostDescriptor( const HostDescriptor& _s )
{
  name = ((HostDescriptor&)_s).name;
  alive = ((HostDescriptor&)_s).alive;
}

FoamXServer::HostDescriptor::~HostDescriptor()
{
}

FoamXServer::HostDescriptor&
FoamXServer::HostDescriptor::operator=( const HostDescriptor& _s )
{
  name = ((HostDescriptor&)_s).name;
  alive = ((HostDescriptor&)_s).alive;
  return *this;
}
#endif

class _Marshaller_FoamXServer_HostDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::HostDescriptor _MICO_T;
  public:
    ~_Marshaller_FoamXServer_HostDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_HostDescriptor::~_Marshaller_FoamXServer_HostDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_HostDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_HostDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_HostDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_HostDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->name._for_demarshal() ) &&
    CORBA::_stc_boolean->demarshal( dc, &((_MICO_T*)v)->alive ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_HostDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->name.inout() );
  CORBA::_stc_boolean->marshal( ec, &((_MICO_T*)v)->alive );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_HostDescriptor::typecode()
{
  return FoamXServer::_tc_HostDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_HostDescriptor;

void operator<<=( CORBA::Any &_a, const FoamXServer::HostDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_HostDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::HostDescriptor *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::HostDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_HostDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::HostDescriptor *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_HostDescriptor, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_HostDescriptorList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_ApplicationDescriptor;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::ApplicationDescriptor::ApplicationDescriptor()
{
}

FoamXServer::ApplicationDescriptor::ApplicationDescriptor( const ApplicationDescriptor& _s )
{
  name = ((ApplicationDescriptor&)_s).name;
  category = ((ApplicationDescriptor&)_s).category;
  path = ((ApplicationDescriptor&)_s).path;
  systemClass = ((ApplicationDescriptor&)_s).systemClass;
}

FoamXServer::ApplicationDescriptor::~ApplicationDescriptor()
{
}

FoamXServer::ApplicationDescriptor&
FoamXServer::ApplicationDescriptor::operator=( const ApplicationDescriptor& _s )
{
  name = ((ApplicationDescriptor&)_s).name;
  category = ((ApplicationDescriptor&)_s).category;
  path = ((ApplicationDescriptor&)_s).path;
  systemClass = ((ApplicationDescriptor&)_s).systemClass;
  return *this;
}
#endif

class _Marshaller_FoamXServer_ApplicationDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::ApplicationDescriptor _MICO_T;
  public:
    ~_Marshaller_FoamXServer_ApplicationDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_ApplicationDescriptor::~_Marshaller_FoamXServer_ApplicationDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_ApplicationDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_ApplicationDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_ApplicationDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_ApplicationDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->name._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->category._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->path._for_demarshal() ) &&
    CORBA::_stc_boolean->demarshal( dc, &((_MICO_T*)v)->systemClass ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_ApplicationDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->name.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->category.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->path.inout() );
  CORBA::_stc_boolean->marshal( ec, &((_MICO_T*)v)->systemClass );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_ApplicationDescriptor::typecode()
{
  return FoamXServer::_tc_ApplicationDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_ApplicationDescriptor;

void operator<<=( CORBA::Any &_a, const FoamXServer::ApplicationDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ApplicationDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::ApplicationDescriptor *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::ApplicationDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_ApplicationDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::ApplicationDescriptor *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_ApplicationDescriptor, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_ApplicationDescriptorList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_CaseDescriptor;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::CaseDescriptor::CaseDescriptor()
{
}

FoamXServer::CaseDescriptor::CaseDescriptor( const CaseDescriptor& _s )
{
  rootDir = ((CaseDescriptor&)_s).rootDir;
  rawRootDir = ((CaseDescriptor&)_s).rawRootDir;
  caseName = ((CaseDescriptor&)_s).caseName;
  app = ((CaseDescriptor&)_s).app;
  nProcs = ((CaseDescriptor&)_s).nProcs;
  managed = ((CaseDescriptor&)_s).managed;
  locked = ((CaseDescriptor&)_s).locked;
  error = ((CaseDescriptor&)_s).error;
}

FoamXServer::CaseDescriptor::~CaseDescriptor()
{
}

FoamXServer::CaseDescriptor&
FoamXServer::CaseDescriptor::operator=( const CaseDescriptor& _s )
{
  rootDir = ((CaseDescriptor&)_s).rootDir;
  rawRootDir = ((CaseDescriptor&)_s).rawRootDir;
  caseName = ((CaseDescriptor&)_s).caseName;
  app = ((CaseDescriptor&)_s).app;
  nProcs = ((CaseDescriptor&)_s).nProcs;
  managed = ((CaseDescriptor&)_s).managed;
  locked = ((CaseDescriptor&)_s).locked;
  error = ((CaseDescriptor&)_s).error;
  return *this;
}
#endif

class _Marshaller_FoamXServer_CaseDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseDescriptor _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseDescriptor::~_Marshaller_FoamXServer_CaseDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_CaseDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_CaseDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->rootDir._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->rawRootDir._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->caseName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->app._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->nProcs ) &&
    CORBA::_stc_boolean->demarshal( dc, &((_MICO_T*)v)->managed ) &&
    CORBA::_stc_boolean->demarshal( dc, &((_MICO_T*)v)->locked ) &&
    CORBA::_stc_boolean->demarshal( dc, &((_MICO_T*)v)->error ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_CaseDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->rootDir.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->rawRootDir.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->caseName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->app.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->nProcs );
  CORBA::_stc_boolean->marshal( ec, &((_MICO_T*)v)->managed );
  CORBA::_stc_boolean->marshal( ec, &((_MICO_T*)v)->locked );
  CORBA::_stc_boolean->marshal( ec, &((_MICO_T*)v)->error );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseDescriptor::typecode()
{
  return FoamXServer::_tc_CaseDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseDescriptor;

void operator<<=( CORBA::Any &_a, const FoamXServer::CaseDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::CaseDescriptor *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::CaseDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::CaseDescriptor *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_CaseDescriptor, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_CaseDescriptorList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_JobID;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::JobID::JobID()
{
}

FoamXServer::JobID::JobID( const JobID& _s )
{
  hostName = ((JobID&)_s).hostName;
  processID = ((JobID&)_s).processID;
}

FoamXServer::JobID::~JobID()
{
}

FoamXServer::JobID&
FoamXServer::JobID::operator=( const JobID& _s )
{
  hostName = ((JobID&)_s).hostName;
  processID = ((JobID&)_s).processID;
  return *this;
}
#endif

class _Marshaller_FoamXServer_JobID : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::JobID _MICO_T;
  public:
    ~_Marshaller_FoamXServer_JobID();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_JobID::~_Marshaller_FoamXServer_JobID()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_JobID::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_JobID::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_JobID::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_JobID::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->hostName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->processID ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_JobID::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->hostName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->processID );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_JobID::typecode()
{
  return FoamXServer::_tc_JobID;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobID;

void operator<<=( CORBA::Any &_a, const FoamXServer::JobID &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobID, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::JobID *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::JobID &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobID, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::JobID *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_JobID, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_JobIDList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_JobDescriptor;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::JobDescriptor::JobDescriptor()
{
}

FoamXServer::JobDescriptor::JobDescriptor( const JobDescriptor& _s )
{
  jobID = ((JobDescriptor&)_s).jobID;
  ppid = ((JobDescriptor&)_s).ppid;
  pgid = ((JobDescriptor&)_s).pgid;
  startDate = ((JobDescriptor&)_s).startDate;
  startTime = ((JobDescriptor&)_s).startTime;
  userName = ((JobDescriptor&)_s).userName;
  foamVersion = ((JobDescriptor&)_s).foamVersion;
  code = ((JobDescriptor&)_s).code;
  argList = ((JobDescriptor&)_s).argList;
  currentDir = ((JobDescriptor&)_s).currentDir;
  rootDir = ((JobDescriptor&)_s).rootDir;
  caseName = ((JobDescriptor&)_s).caseName;
  nProcs = ((JobDescriptor&)_s).nProcs;
  slaves = ((JobDescriptor&)_s).slaves;
  nCountedProcs = ((JobDescriptor&)_s).nCountedProcs;
  cpuTime = ((JobDescriptor&)_s).cpuTime;
  endDate = ((JobDescriptor&)_s).endDate;
  endTime = ((JobDescriptor&)_s).endTime;
  status = ((JobDescriptor&)_s).status;
}

FoamXServer::JobDescriptor::~JobDescriptor()
{
}

FoamXServer::JobDescriptor&
FoamXServer::JobDescriptor::operator=( const JobDescriptor& _s )
{
  jobID = ((JobDescriptor&)_s).jobID;
  ppid = ((JobDescriptor&)_s).ppid;
  pgid = ((JobDescriptor&)_s).pgid;
  startDate = ((JobDescriptor&)_s).startDate;
  startTime = ((JobDescriptor&)_s).startTime;
  userName = ((JobDescriptor&)_s).userName;
  foamVersion = ((JobDescriptor&)_s).foamVersion;
  code = ((JobDescriptor&)_s).code;
  argList = ((JobDescriptor&)_s).argList;
  currentDir = ((JobDescriptor&)_s).currentDir;
  rootDir = ((JobDescriptor&)_s).rootDir;
  caseName = ((JobDescriptor&)_s).caseName;
  nProcs = ((JobDescriptor&)_s).nProcs;
  slaves = ((JobDescriptor&)_s).slaves;
  nCountedProcs = ((JobDescriptor&)_s).nCountedProcs;
  cpuTime = ((JobDescriptor&)_s).cpuTime;
  endDate = ((JobDescriptor&)_s).endDate;
  endTime = ((JobDescriptor&)_s).endTime;
  status = ((JobDescriptor&)_s).status;
  return *this;
}
#endif

class _Marshaller_FoamXServer_JobDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::JobDescriptor _MICO_T;
  public:
    ~_Marshaller_FoamXServer_JobDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_JobDescriptor::~_Marshaller_FoamXServer_JobDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_JobDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_JobDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_JobDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_JobDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    _marshaller_FoamXServer_JobID->demarshal( dc, &((_MICO_T*)v)->jobID ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->ppid ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->pgid ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->startDate._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->startTime._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->userName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->foamVersion._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->code._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->argList._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->currentDir._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->rootDir._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->caseName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->nProcs ) &&
    _marshaller__seq_FoamXServer_JobID->demarshal( dc, &((_MICO_T*)v)->slaves ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->nCountedProcs ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->cpuTime ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->endDate._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->endTime._for_demarshal() ) &&
    _marshaller_FoamXServer_JobStatus->demarshal( dc, &((_MICO_T*)v)->status ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_JobDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  _marshaller_FoamXServer_JobID->marshal( ec, &((_MICO_T*)v)->jobID );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->ppid );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->pgid );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->startDate.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->startTime.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->userName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->foamVersion.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->code.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->argList.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->currentDir.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->rootDir.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->caseName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->nProcs );
  _marshaller__seq_FoamXServer_JobID->marshal( ec, &((_MICO_T*)v)->slaves );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->nCountedProcs );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->cpuTime );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->endDate.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->endTime.inout() );
  _marshaller_FoamXServer_JobStatus->marshal( ec, &((_MICO_T*)v)->status );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_JobDescriptor::typecode()
{
  return FoamXServer::_tc_JobDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobDescriptor;

void operator<<=( CORBA::Any &_a, const FoamXServer::JobDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::JobDescriptor *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::JobDescriptor &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_JobDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::JobDescriptor *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_JobDescriptor, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_JobDescriptorList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_DimensionSet;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::DimensionSet::DimensionSet()
{
}

FoamXServer::DimensionSet::DimensionSet( const DimensionSet& _s )
{
  mass = ((DimensionSet&)_s).mass;
  length = ((DimensionSet&)_s).length;
  time = ((DimensionSet&)_s).time;
  temperature = ((DimensionSet&)_s).temperature;
  moles = ((DimensionSet&)_s).moles;
  current = ((DimensionSet&)_s).current;
  luminousIntensity = ((DimensionSet&)_s).luminousIntensity;
}

FoamXServer::DimensionSet::~DimensionSet()
{
}

FoamXServer::DimensionSet&
FoamXServer::DimensionSet::operator=( const DimensionSet& _s )
{
  mass = ((DimensionSet&)_s).mass;
  length = ((DimensionSet&)_s).length;
  time = ((DimensionSet&)_s).time;
  temperature = ((DimensionSet&)_s).temperature;
  moles = ((DimensionSet&)_s).moles;
  current = ((DimensionSet&)_s).current;
  luminousIntensity = ((DimensionSet&)_s).luminousIntensity;
  return *this;
}
#endif

class _Marshaller_FoamXServer_DimensionSet : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::DimensionSet _MICO_T;
  public:
    ~_Marshaller_FoamXServer_DimensionSet();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_DimensionSet::~_Marshaller_FoamXServer_DimensionSet()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_DimensionSet::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_DimensionSet::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_DimensionSet::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_DimensionSet::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  return
    dc.struct_begin() &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->mass ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->length ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->time ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->temperature ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->moles ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->current ) &&
    CORBA::_stc_double->demarshal( dc, &((_MICO_T*)v)->luminousIntensity ) &&
    dc.struct_end();
}

void _Marshaller_FoamXServer_DimensionSet::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.struct_begin();
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->mass );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->length );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->time );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->temperature );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->moles );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->current );
  CORBA::_stc_double->marshal( ec, &((_MICO_T*)v)->luminousIntensity );
  ec.struct_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_DimensionSet::typecode()
{
  return FoamXServer::_tc_DimensionSet;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_DimensionSet;

void operator<<=( CORBA::Any &_a, const FoamXServer::DimensionSet &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_DimensionSet, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::DimensionSet *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::DimensionSet &_s )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_DimensionSet, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::DimensionSet *&_s )
{
  return _a.to_static_any (_marshaller_FoamXServer_DimensionSet, (void *&)_s);
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_DimensionSetList;
}

namespace FoamXServer
{
CORBA::TypeCodeConst _tc_FoamXSYSError;
}

#ifdef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXSYSError::FoamXSYSError()
{
}

FoamXServer::FoamXSYSError::FoamXSYSError( const FoamXSYSError& _s )
{
  errorCode = ((FoamXSYSError&)_s).errorCode;
  errorMessage = ((FoamXSYSError&)_s).errorMessage;
  hostName = ((FoamXSYSError&)_s).hostName;
  methodName = ((FoamXSYSError&)_s).methodName;
  fileName = ((FoamXSYSError&)_s).fileName;
  lineNo = ((FoamXSYSError&)_s).lineNo;
}

FoamXServer::FoamXSYSError::~FoamXSYSError()
{
}

FoamXServer::FoamXSYSError&
FoamXServer::FoamXSYSError::operator=( const FoamXSYSError& _s )
{
  errorCode = ((FoamXSYSError&)_s).errorCode;
  errorMessage = ((FoamXSYSError&)_s).errorMessage;
  hostName = ((FoamXSYSError&)_s).hostName;
  methodName = ((FoamXSYSError&)_s).methodName;
  fileName = ((FoamXSYSError&)_s).fileName;
  lineNo = ((FoamXSYSError&)_s).lineNo;
  return *this;
}
#endif

#ifndef HAVE_EXPLICIT_STRUCT_OPS
FoamXServer::FoamXSYSError::FoamXSYSError()
{
}

#endif

FoamXServer::FoamXSYSError::FoamXSYSError( FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2, const char* _m3, const char* _m4, CORBA::Long _m5 )
{
  errorCode = _m0;
  errorMessage = _m1;
  hostName = _m2;
  methodName = _m3;
  fileName = _m4;
  lineNo = _m5;
}

class _Marshaller_FoamXServer_FoamXSYSError : public ::CORBA::StaticTypeInfo {
    typedef ::FoamXServer::FoamXSYSError _MICO_T;
  public:
    ~_Marshaller_FoamXServer_FoamXSYSError();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_FoamXSYSError::~_Marshaller_FoamXServer_FoamXSYSError()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_FoamXSYSError::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller_FoamXServer_FoamXSYSError::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller_FoamXServer_FoamXSYSError::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller_FoamXServer_FoamXSYSError::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  string repoid;
  return
    dc.except_begin( repoid ) &&
    _marshaller_FoamXServer_ErrorCode->demarshal( dc, &((_MICO_T*)v)->errorCode ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->errorMessage._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->hostName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->methodName._for_demarshal() ) &&
    CORBA::_stc_string->demarshal( dc, &((_MICO_T*)v)->fileName._for_demarshal() ) &&
    CORBA::_stc_long->demarshal( dc, &((_MICO_T*)v)->lineNo ) &&
    dc.except_end();
}

void _Marshaller_FoamXServer_FoamXSYSError::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.except_begin( "IDL:FoamXServer/FoamXSYSError:1.0" );
  _marshaller_FoamXServer_ErrorCode->marshal( ec, &((_MICO_T*)v)->errorCode );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->errorMessage.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->hostName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->methodName.inout() );
  CORBA::_stc_string->marshal( ec, &((_MICO_T*)v)->fileName.inout() );
  CORBA::_stc_long->marshal( ec, &((_MICO_T*)v)->lineNo );
  ec.except_end();
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_FoamXSYSError::typecode()
{
  return FoamXServer::_tc_FoamXSYSError;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXSYSError;

void operator<<=( CORBA::Any &_a, const FoamXServer::FoamXSYSError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXSYSError, &_e);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, FoamXServer::FoamXSYSError *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, FoamXServer::FoamXSYSError &_e )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_FoamXSYSError, &_e);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const FoamXServer::FoamXSYSError *&_e )
{
  return _a.to_static_any (_marshaller_FoamXServer_FoamXSYSError, (void *&)_e);
}

void FoamXServer::FoamXSYSError::_throwit() const
{
  #ifdef HAVE_EXCEPTIONS
  #ifdef HAVE_STD_EH
  throw *this;
  #else
  throw FoamXSYSError_var( (FoamXServer::FoamXSYSError*)_clone() );
  #endif
  #else
  CORBA::Exception::_throw_failed( _clone() );
  #endif
}

const char *FoamXServer::FoamXSYSError::_repoid() const
{
  return "IDL:FoamXServer/FoamXSYSError:1.0";
}

void FoamXServer::FoamXSYSError::_encode( CORBA::DataEncoder &_en ) const
{
  _marshaller_FoamXServer_FoamXSYSError->marshal( _en, (void*) this );
}

void FoamXServer::FoamXSYSError::_encode_any( CORBA::Any &_a ) const
{
  _a <<= *this;
}

CORBA::Exception *FoamXServer::FoamXSYSError::_clone() const
{
  return new FoamXSYSError( *this );
}

FoamXServer::FoamXSYSError *FoamXServer::FoamXSYSError::_downcast( CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXSYSError:1.0" ) )
    return (FoamXSYSError *) _ex;
  return NULL;
}

const FoamXServer::FoamXSYSError *FoamXServer::FoamXSYSError::_downcast( const CORBA::Exception *_ex )
{
  if( _ex && !strcmp( _ex->_repoid(), "IDL:FoamXServer/FoamXSYSError:1.0" ) )
    return (FoamXSYSError *) _ex;
  return NULL;
}


/*
 * Base interface for class ICaseServer
 */

FoamXServer::CaseServer::ICaseServer::~ICaseServer()
{
}

void *
FoamXServer::CaseServer::ICaseServer::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/ICaseServer:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::ICaseServer_ptr
FoamXServer::CaseServer::ICaseServer::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::ICaseServer_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/ICaseServer:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::ICaseServer_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/ICaseServer:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/ICaseServer:1.0")) {
      _o = new FoamXServer::CaseServer::ICaseServer_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::ICaseServer_ptr
FoamXServer::CaseServer::ICaseServer::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_ICaseServer;
}
}
class _Marshaller_FoamXServer_CaseServer_ICaseServer : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::ICaseServer_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_ICaseServer();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_ICaseServer::~_Marshaller_FoamXServer_CaseServer_ICaseServer()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_ICaseServer::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_ICaseServer::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::ICaseServer::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_ICaseServer::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_ICaseServer::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_ICaseServer::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::ICaseServer::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_ICaseServer::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_ICaseServer::typecode()
{
  return FoamXServer::CaseServer::_tc_ICaseServer;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_ICaseServer;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::ICaseServer_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_ICaseServer, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::ICaseServer_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_ICaseServer, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::ICaseServer_ptr &_obj )
{
  FoamXServer::CaseServer::ICaseServer_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_ICaseServer, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class ICaseServer
 */

FoamXServer::CaseServer::ICaseServer_stub::~ICaseServer_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::ICaseServer::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/ICaseServer:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::ICaseServer *
POA_FoamXServer::CaseServer::ICaseServer::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/ICaseServer:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::ICaseServer *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::ICaseServer_stub_clp::ICaseServer_stub_clp ()
{
}

FoamXServer::CaseServer::ICaseServer_stub_clp::ICaseServer_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::ICaseServer_stub_clp::~ICaseServer_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseServer::ICaseServer_stub::managed()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_managed" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseServer::ICaseServer_stub_clp::managed()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->managed();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::managed();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::managed( CORBA::Boolean _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );
  CORBA::StaticRequest __req( this, "_set_managed" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::managed( CORBA::Boolean _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->managed(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::managed(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::ICaseServer_stub::caseRoot()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_caseRoot" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::ICaseServer_stub_clp::caseRoot()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->caseRoot();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::caseRoot();
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::ICaseServer_stub::caseName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_caseName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::ICaseServer_stub_clp::caseName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->caseName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::caseName();
}

#endif // MICO_CONF_NO_POA

FoamXServer::CaseServer::IApplication_ptr FoamXServer::CaseServer::ICaseServer_stub::application()
{
  FoamXServer::CaseServer::IApplication_ptr _res = FoamXServer::CaseServer::IApplication::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IApplication, &_res );

  CORBA::StaticRequest __req( this, "_get_application" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::CaseServer::IApplication_ptr
FoamXServer::CaseServer::ICaseServer_stub_clp::application()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::CaseServer::IApplication_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->application();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::application();
}

#endif // MICO_CONF_NO_POA

FoamXServer::CaseServer::IFoamProperties_ptr FoamXServer::CaseServer::ICaseServer_stub::foamProperties()
{
  FoamXServer::CaseServer::IFoamProperties_ptr _res = FoamXServer::CaseServer::IFoamProperties::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IFoamProperties, &_res );

  CORBA::StaticRequest __req( this, "_get_foamProperties" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::CaseServer::IFoamProperties_ptr
FoamXServer::CaseServer::ICaseServer_stub_clp::foamProperties()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::CaseServer::IFoamProperties_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->foamProperties();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::foamProperties();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::ICaseServer_stub::availableTimeSteps()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_availableTimeSteps" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::ICaseServer_stub_clp::availableTimeSteps()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->availableTimeSteps();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::availableTimeSteps();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseServer::ICaseServer_stub::meshDefined()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_meshDefined" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseServer::ICaseServer_stub_clp::meshDefined()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->meshDefined();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::meshDefined();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::ICaseServer_stub::patchNames()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_patchNames" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::ICaseServer_stub_clp::patchNames()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchNames();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::patchNames();
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::ICaseServer_stub::getTime()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "getTime" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::ICaseServer_stub_clp::getTime()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getTime();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::getTime();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::setTime( const char* _par_timeName, CORBA::Long _par_timeIndex )
{
  CORBA::StaticAny _sa_timeName( CORBA::_stc_string, &_par_timeName );
  CORBA::StaticAny _sa_timeIndex( CORBA::_stc_long, &_par_timeIndex );
  CORBA::StaticRequest __req( this, "setTime" );
  __req.add_in_arg( &_sa_timeName );
  __req.add_in_arg( &_sa_timeIndex );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::setTime( const char* _par_timeName, CORBA::Long _par_timeIndex )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->setTime(_par_timeName, _par_timeIndex);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::setTime(_par_timeName, _par_timeIndex);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::readMesh()
{
  CORBA::StaticRequest __req( this, "readMesh" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::readMesh()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->readMesh();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::readMesh();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::importMesh( const char* _par_hostName, const char* _par_rootDir, const char* _par_caseName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticRequest __req( this, "importMesh" );
  __req.add_in_arg( &_sa_hostName );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::importMesh( const char* _par_hostName, const char* _par_rootDir, const char* _par_caseName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->importMesh(_par_hostName, _par_rootDir, _par_caseName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::importMesh(_par_hostName, _par_rootDir, _par_caseName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::getFieldValues( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricField_out _par_fieldValues )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_fieldValues( _marshaller_FoamXServer_CaseServer_IGeometricField, &_par_fieldValues.ptr() );
  CORBA::StaticRequest __req( this, "getFieldValues" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_out_arg( &_sa_fieldValues );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::getFieldValues( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricField_out _par_fieldValues )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getFieldValues(_par_fieldName, _par_fieldValues);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::getFieldValues(_par_fieldName, _par_fieldValues);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::addPatch( const char* _par_patchName, const char* _par_patchPhysicalType )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType );
  CORBA::StaticRequest __req( this, "addPatch" );
  __req.add_in_arg( &_sa_patchName );
  __req.add_in_arg( &_sa_patchPhysicalType );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::addPatch( const char* _par_patchName, const char* _par_patchPhysicalType )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addPatch(_par_patchName, _par_patchPhysicalType);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::addPatch(_par_patchName, _par_patchPhysicalType);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::deletePatch( const char* _par_patchName )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticRequest __req( this, "deletePatch" );
  __req.add_in_arg( &_sa_patchName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::deletePatch( const char* _par_patchName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deletePatch(_par_patchName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::deletePatch(_par_patchName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::deleteAllPatches()
{
  CORBA::StaticRequest __req( this, "deleteAllPatches" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::deleteAllPatches()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteAllPatches();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::deleteAllPatches();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::setPatchPhysicalType( const char* _par_patchName, const char* _par_patchPhysicalType )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType );
  CORBA::StaticRequest __req( this, "setPatchPhysicalType" );
  __req.add_in_arg( &_sa_patchName );
  __req.add_in_arg( &_sa_patchPhysicalType );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::setPatchPhysicalType( const char* _par_patchName, const char* _par_patchPhysicalType )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->setPatchPhysicalType(_par_patchName, _par_patchPhysicalType);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::setPatchPhysicalType(_par_patchName, _par_patchPhysicalType);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::getPatchPhysicalType( const char* _par_patchName, CORBA::String_out _par_patchPhysicalType )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType.ptr() );
  CORBA::StaticRequest __req( this, "getPatchPhysicalType" );
  __req.add_in_arg( &_sa_patchName );
  __req.add_out_arg( &_sa_patchPhysicalType );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::getPatchPhysicalType( const char* _par_patchName, CORBA::String_out _par_patchPhysicalType )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchPhysicalType(_par_patchName, _par_patchPhysicalType);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::getPatchPhysicalType(_par_patchName, _par_patchPhysicalType);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::getDictionary( const char* _par_dictionaryName, CORBA::Boolean _par_forceRead, FoamXServer::IDictionaryEntry_out _par_dictRoot )
{
  CORBA::StaticAny _sa_dictionaryName( CORBA::_stc_string, &_par_dictionaryName );
  CORBA::StaticAny _sa_forceRead( CORBA::_stc_boolean, &_par_forceRead );
  CORBA::StaticAny _sa_dictRoot( _marshaller_FoamXServer_IDictionaryEntry, &_par_dictRoot.ptr() );
  CORBA::StaticRequest __req( this, "getDictionary" );
  __req.add_in_arg( &_sa_dictionaryName );
  __req.add_in_arg( &_sa_forceRead );
  __req.add_out_arg( &_sa_dictRoot );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::getDictionary( const char* _par_dictionaryName, CORBA::Boolean _par_forceRead, FoamXServer::IDictionaryEntry_out _par_dictRoot )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getDictionary(_par_dictionaryName, _par_forceRead, _par_dictRoot);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::getDictionary(_par_dictionaryName, _par_forceRead, _par_dictRoot);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::readFile( const char* _par_name, CORBA::String_out _par_contents )
{
  CORBA::StaticAny _sa_name( CORBA::_stc_string, &_par_name );
  CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents.ptr() );
  CORBA::StaticRequest __req( this, "readFile" );
  __req.add_in_arg( &_sa_name );
  __req.add_out_arg( &_sa_contents );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::readFile( const char* _par_name, CORBA::String_out _par_contents )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->readFile(_par_name, _par_contents);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::readFile(_par_name, _par_contents);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::writeFile( const char* _par_name, const char* _par_contents )
{
  CORBA::StaticAny _sa_name( CORBA::_stc_string, &_par_name );
  CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents );
  CORBA::StaticRequest __req( this, "writeFile" );
  __req.add_in_arg( &_sa_name );
  __req.add_in_arg( &_sa_contents );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::writeFile( const char* _par_name, const char* _par_contents )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->writeFile(_par_name, _par_contents);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::writeFile(_par_name, _par_contents);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::CaseServer::ICaseServer_stub::fileModificationDate( const char* _par_fileName )
{
  CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName );
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "fileModificationDate" );
  __req.add_in_arg( &_sa_fileName );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::CaseServer::ICaseServer_stub_clp::fileModificationDate( const char* _par_fileName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->fileModificationDate(_par_fileName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::fileModificationDate(_par_fileName);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::CaseServer::ICaseServer_stub::runCase( const char* _par_arguments )
{
  CORBA::StaticAny _sa_arguments( CORBA::_stc_string, &_par_arguments );
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "runCase" );
  __req.add_in_arg( &_sa_arguments );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::CaseServer::ICaseServer_stub_clp::runCase( const char* _par_arguments )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->runCase(_par_arguments);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::runCase(_par_arguments);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::killCase()
{
  CORBA::StaticRequest __req( this, "killCase" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::killCase()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->killCase();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::killCase();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::validate();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseServer::ICaseServer_stub::modified()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "modified" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseServer::ICaseServer_stub_clp::modified()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->modified();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::ICaseServer_stub::modified();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::save()
{
  CORBA::StaticRequest __req( this, "save" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::save()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->save();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::save();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::ICaseServer_stub::close()
{
  CORBA::StaticRequest __req( this, "close" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::ICaseServer_stub_clp::close()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::ICaseServer * _myserv = POA_FoamXServer::CaseServer::ICaseServer::_narrow (_serv);
    if (_myserv) {
      _myserv->close();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::ICaseServer_stub::close();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IFoamProperties
 */

FoamXServer::CaseServer::IFoamProperties::~IFoamProperties()
{
}

void *
FoamXServer::CaseServer::IFoamProperties::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IFoamProperties:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IFoamProperties_ptr
FoamXServer::CaseServer::IFoamProperties::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IFoamProperties_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IFoamProperties:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IFoamProperties_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IFoamProperties:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IFoamProperties:1.0")) {
      _o = new FoamXServer::CaseServer::IFoamProperties_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IFoamProperties_ptr
FoamXServer::CaseServer::IFoamProperties::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IFoamProperties;
}
}
class _Marshaller_FoamXServer_CaseServer_IFoamProperties : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IFoamProperties_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IFoamProperties();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IFoamProperties::~_Marshaller_FoamXServer_CaseServer_IFoamProperties()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IFoamProperties::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IFoamProperties::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IFoamProperties::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IFoamProperties::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IFoamProperties::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IFoamProperties::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IFoamProperties::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IFoamProperties::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IFoamProperties::typecode()
{
  return FoamXServer::CaseServer::_tc_IFoamProperties;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IFoamProperties;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IFoamProperties_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IFoamProperties, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IFoamProperties_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IFoamProperties, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IFoamProperties_ptr &_obj )
{
  FoamXServer::CaseServer::IFoamProperties_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IFoamProperties, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IFoamProperties
 */

FoamXServer::CaseServer::IFoamProperties_stub::~IFoamProperties_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IFoamProperties::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IFoamProperties:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IFoamProperties *
POA_FoamXServer::CaseServer::IFoamProperties::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IFoamProperties:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IFoamProperties *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IFoamProperties_stub_clp::IFoamProperties_stub_clp ()
{
}

FoamXServer::CaseServer::IFoamProperties_stub_clp::IFoamProperties_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IFoamProperties_stub_clp::~IFoamProperties_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::availableModules()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_availableModules" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::availableModules()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->availableModules();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::availableModules();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::rootDirectories()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_rootDirectories" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::rootDirectories()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->rootDirectories();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::rootDirectories();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::rawRootDirectories()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_rawRootDirectories" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::rawRootDirectories()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->rawRootDirectories();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::rawRootDirectories();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::foamTypes()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_foamTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::foamTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->foamTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::foamTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::geometryTypes()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_geometryTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::geometryTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->geometryTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::geometryTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::patchTypes()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_patchTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::patchTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::patchTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IFoamProperties_stub::patchFieldTypes()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_patchFieldTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::patchFieldTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchFieldTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::patchFieldTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::ApplicationDescriptorList* FoamXServer::CaseServer::IFoamProperties_stub::applicationes()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ApplicationDescriptor );

  CORBA::StaticRequest __req( this, "_get_applicationes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::ApplicationDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ApplicationDescriptorList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::applicationes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ApplicationDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->applicationes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::applicationes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::ApplicationDescriptorList* FoamXServer::CaseServer::IFoamProperties_stub::utilities()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ApplicationDescriptor );

  CORBA::StaticRequest __req( this, "_get_utilities" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::ApplicationDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ApplicationDescriptorList*
FoamXServer::CaseServer::IFoamProperties_stub_clp::utilities()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ApplicationDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->utilities();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IFoamProperties_stub::utilities();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::addRootDirectory( const char* _par_rawRootDir )
{
  CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir );
  CORBA::StaticRequest __req( this, "addRootDirectory" );
  __req.add_in_arg( &_sa_rawRootDir );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::addRootDirectory( const char* _par_rawRootDir )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addRootDirectory(_par_rawRootDir);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::addRootDirectory(_par_rawRootDir);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::deleteRootDirectory( const char* _par_rawRootDir )
{
  CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir );
  CORBA::StaticRequest __req( this, "deleteRootDirectory" );
  __req.add_in_arg( &_sa_rawRootDir );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::deleteRootDirectory( const char* _par_rawRootDir )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteRootDirectory(_par_rawRootDir);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::deleteRootDirectory(_par_rawRootDir);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getFoamType( const char* _par_foamTypeName, FoamXServer::ITypeDescriptor_out _par_typeDesc )
{
  CORBA::StaticAny _sa_foamTypeName( CORBA::_stc_string, &_par_foamTypeName );
  CORBA::StaticAny _sa_typeDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_typeDesc.ptr() );
  CORBA::StaticRequest __req( this, "getFoamType" );
  __req.add_in_arg( &_sa_foamTypeName );
  __req.add_out_arg( &_sa_typeDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getFoamType( const char* _par_foamTypeName, FoamXServer::ITypeDescriptor_out _par_typeDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getFoamType(_par_foamTypeName, _par_typeDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getFoamType(_par_foamTypeName, _par_typeDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getGeometryType( const char* _par_geometryTypeName, FoamXServer::CaseServer::IGeometryDescriptor_out _par_geometryDesc )
{
  CORBA::StaticAny _sa_geometryTypeName( CORBA::_stc_string, &_par_geometryTypeName );
  CORBA::StaticAny _sa_geometryDesc( _marshaller_FoamXServer_CaseServer_IGeometryDescriptor, &_par_geometryDesc.ptr() );
  CORBA::StaticRequest __req( this, "getGeometryType" );
  __req.add_in_arg( &_sa_geometryTypeName );
  __req.add_out_arg( &_sa_geometryDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getGeometryType( const char* _par_geometryTypeName, FoamXServer::CaseServer::IGeometryDescriptor_out _par_geometryDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getGeometryType(_par_geometryTypeName, _par_geometryDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getGeometryType(_par_geometryTypeName, _par_geometryDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getPatchType( const char* _par_patchTypeName, FoamXServer::CaseServer::IPatchDescriptor_out _par_patchDesc )
{
  CORBA::StaticAny _sa_patchTypeName( CORBA::_stc_string, &_par_patchTypeName );
  CORBA::StaticAny _sa_patchDesc( _marshaller_FoamXServer_CaseServer_IPatchDescriptor, &_par_patchDesc.ptr() );
  CORBA::StaticRequest __req( this, "getPatchType" );
  __req.add_in_arg( &_sa_patchTypeName );
  __req.add_out_arg( &_sa_patchDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getPatchType( const char* _par_patchTypeName, FoamXServer::CaseServer::IPatchDescriptor_out _par_patchDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchType(_par_patchTypeName, _par_patchDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getPatchType(_par_patchTypeName, _par_patchDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::findPatchType( const char* _par_patchTypeName, FoamXServer::CaseServer::IPatchDescriptor_out _par_patchDesc )
{
  CORBA::StaticAny _sa_patchTypeName( CORBA::_stc_string, &_par_patchTypeName );
  CORBA::StaticAny _sa_patchDesc( _marshaller_FoamXServer_CaseServer_IPatchDescriptor, &_par_patchDesc.ptr() );
  CORBA::StaticRequest __req( this, "findPatchType" );
  __req.add_in_arg( &_sa_patchTypeName );
  __req.add_out_arg( &_sa_patchDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::findPatchType( const char* _par_patchTypeName, FoamXServer::CaseServer::IPatchDescriptor_out _par_patchDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->findPatchType(_par_patchTypeName, _par_patchDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::findPatchType(_par_patchTypeName, _par_patchDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getPatchFieldType( const char* _par_patchFieldTypeName, FoamXServer::ITypeDescriptor_out _par_patchFieldDesc )
{
  CORBA::StaticAny _sa_patchFieldTypeName( CORBA::_stc_string, &_par_patchFieldTypeName );
  CORBA::StaticAny _sa_patchFieldDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_patchFieldDesc.ptr() );
  CORBA::StaticRequest __req( this, "getPatchFieldType" );
  __req.add_in_arg( &_sa_patchFieldTypeName );
  __req.add_out_arg( &_sa_patchFieldDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getPatchFieldType( const char* _par_patchFieldTypeName, FoamXServer::ITypeDescriptor_out _par_patchFieldDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchFieldType(_par_patchFieldTypeName, _par_patchFieldDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getPatchFieldType(_par_patchFieldTypeName, _par_patchFieldDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::findPatchFieldType( const char* _par_patchFieldTypeName, FoamXServer::ITypeDescriptor_out _par_patchFieldDesc )
{
  CORBA::StaticAny _sa_patchFieldTypeName( CORBA::_stc_string, &_par_patchFieldTypeName );
  CORBA::StaticAny _sa_patchFieldDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_patchFieldDesc.ptr() );
  CORBA::StaticRequest __req( this, "findPatchFieldType" );
  __req.add_in_arg( &_sa_patchFieldTypeName );
  __req.add_out_arg( &_sa_patchFieldDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::findPatchFieldType( const char* _par_patchFieldTypeName, FoamXServer::ITypeDescriptor_out _par_patchFieldDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->findPatchFieldType(_par_patchFieldTypeName, _par_patchFieldDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::findPatchFieldType(_par_patchFieldTypeName, _par_patchFieldDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getFoamControlDict( FoamXServer::IDictionaryEntry_out _par_controlDict )
{
  CORBA::StaticAny _sa_controlDict( _marshaller_FoamXServer_IDictionaryEntry, &_par_controlDict.ptr() );
  CORBA::StaticRequest __req( this, "getFoamControlDict" );
  __req.add_out_arg( &_sa_controlDict );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getFoamControlDict( FoamXServer::IDictionaryEntry_out _par_controlDict )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getFoamControlDict(_par_controlDict);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getFoamControlDict(_par_controlDict);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getApplication( const char* _par_appName, FoamXServer::CaseServer::IApplication_out _par_app )
{
  CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName );
  CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app.ptr() );
  CORBA::StaticRequest __req( this, "getApplication" );
  __req.add_in_arg( &_sa_appName );
  __req.add_out_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getApplication( const char* _par_appName, FoamXServer::CaseServer::IApplication_out _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getApplication(_par_appName, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getApplication(_par_appName, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::addApplication( const char* _par_appName, FoamXServer::CaseServer::IApplication_out _par_app )
{
  CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName );
  CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app.ptr() );
  CORBA::StaticRequest __req( this, "addApplication" );
  __req.add_in_arg( &_sa_appName );
  __req.add_out_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::addApplication( const char* _par_appName, FoamXServer::CaseServer::IApplication_out _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addApplication(_par_appName, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::addApplication(_par_appName, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::deleteApplication( const char* _par_appName )
{
  CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName );
  CORBA::StaticRequest __req( this, "deleteApplication" );
  __req.add_in_arg( &_sa_appName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::deleteApplication( const char* _par_appName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteApplication(_par_appName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::deleteApplication(_par_appName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::cloneApplication( const char* _par_appNameSrc, const char* _par_appNameDest, const char* _par_appDestPath, FoamXServer::CaseServer::IApplication_out _par_app )
{
  CORBA::StaticAny _sa_appNameSrc( CORBA::_stc_string, &_par_appNameSrc );
  CORBA::StaticAny _sa_appNameDest( CORBA::_stc_string, &_par_appNameDest );
  CORBA::StaticAny _sa_appDestPath( CORBA::_stc_string, &_par_appDestPath );
  CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app.ptr() );
  CORBA::StaticRequest __req( this, "cloneApplication" );
  __req.add_in_arg( &_sa_appNameSrc );
  __req.add_in_arg( &_sa_appNameDest );
  __req.add_in_arg( &_sa_appDestPath );
  __req.add_out_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::cloneApplication( const char* _par_appNameSrc, const char* _par_appNameDest, const char* _par_appDestPath, FoamXServer::CaseServer::IApplication_out _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->cloneApplication(_par_appNameSrc, _par_appNameDest, _par_appDestPath, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::cloneApplication(_par_appNameSrc, _par_appNameDest, _par_appDestPath, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::getUtilityControlDict( const char* _par_utilityName, const char* _par_rootDir, const char* _par_caseName, FoamXServer::IDictionaryEntry_out _par_controlDict )
{
  CORBA::StaticAny _sa_utilityName( CORBA::_stc_string, &_par_utilityName );
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_controlDict( _marshaller_FoamXServer_IDictionaryEntry, &_par_controlDict.ptr() );
  CORBA::StaticRequest __req( this, "getUtilityControlDict" );
  __req.add_in_arg( &_sa_utilityName );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_out_arg( &_sa_controlDict );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::getUtilityControlDict( const char* _par_utilityName, const char* _par_rootDir, const char* _par_caseName, FoamXServer::IDictionaryEntry_out _par_controlDict )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getUtilityControlDict(_par_utilityName, _par_rootDir, _par_caseName, _par_controlDict);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::getUtilityControlDict(_par_utilityName, _par_rootDir, _par_caseName, _par_controlDict);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::validate();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::saveSystemProperties()
{
  CORBA::StaticRequest __req( this, "saveSystemProperties" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::saveSystemProperties()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->saveSystemProperties();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::saveSystemProperties();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IFoamProperties_stub::saveUserProperties()
{
  CORBA::StaticRequest __req( this, "saveUserProperties" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IFoamProperties_stub_clp::saveUserProperties()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IFoamProperties * _myserv = POA_FoamXServer::CaseServer::IFoamProperties::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->saveUserProperties();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IFoamProperties_stub::saveUserProperties();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IApplication
 */

FoamXServer::CaseServer::IApplication::~IApplication()
{
}

void *
FoamXServer::CaseServer::IApplication::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IApplication:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IApplication_ptr
FoamXServer::CaseServer::IApplication::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IApplication_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IApplication:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IApplication_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IApplication:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IApplication:1.0")) {
      _o = new FoamXServer::CaseServer::IApplication_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IApplication_ptr
FoamXServer::CaseServer::IApplication::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IApplication;
}
}
class _Marshaller_FoamXServer_CaseServer_IApplication : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IApplication_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IApplication();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IApplication::~_Marshaller_FoamXServer_CaseServer_IApplication()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IApplication::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IApplication::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IApplication::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IApplication::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IApplication::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IApplication::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IApplication::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IApplication::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IApplication::typecode()
{
  return FoamXServer::CaseServer::_tc_IApplication;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IApplication;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IApplication_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IApplication, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IApplication_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IApplication, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IApplication_ptr &_obj )
{
  FoamXServer::CaseServer::IApplication_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IApplication, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IApplication
 */

FoamXServer::CaseServer::IApplication_stub::~IApplication_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IApplication::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IApplication:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IApplication *
POA_FoamXServer::CaseServer::IApplication::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IApplication:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IApplication *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IApplication_stub_clp::IApplication_stub_clp ()
{
}

FoamXServer::CaseServer::IApplication_stub_clp::IApplication_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IApplication_stub_clp::~IApplication_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IApplication_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IApplication_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IApplication_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IApplication_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IApplication_stub::category()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_category" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IApplication_stub_clp::category()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->category();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::category();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::category( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_category" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::category( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->category(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::category(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IApplication_stub::modules()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_modules" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IApplication_stub_clp::modules()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->modules();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::modules();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::modules( const FoamXServer::StringList& _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stcseq_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_modules" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::modules( const FoamXServer::StringList& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->modules(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::modules(_par__value);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseServer::IApplication_stub::systemClass()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "_get_systemClass" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseServer::IApplication_stub_clp::systemClass()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->systemClass();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::systemClass();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IApplication_stub::fields()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_fields" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IApplication_stub_clp::fields()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->fields();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::fields();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IApplication_stub::patchPhysicalTypes()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_patchPhysicalTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IApplication_stub_clp::patchPhysicalTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchPhysicalTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::patchPhysicalTypes();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CaseServer::IApplication_stub::dictionaries()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_dictionaries" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CaseServer::IApplication_stub_clp::dictionaries()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->dictionaries();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IApplication_stub::dictionaries();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::getField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "getField" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_out_arg( &_sa_fieldDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::getField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getField(_par_fieldName, _par_fieldDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::getField(_par_fieldName, _par_fieldDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::findField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "findField" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_out_arg( &_sa_fieldDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::findField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->findField(_par_fieldName, _par_fieldDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::findField(_par_fieldName, _par_fieldDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::addField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "addField" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_out_arg( &_sa_fieldDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::addField( const char* _par_fieldName, FoamXServer::CaseServer::IGeometricFieldDescriptor_out _par_fieldDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addField(_par_fieldName, _par_fieldDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::addField(_par_fieldName, _par_fieldDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::deleteField( const char* _par_fieldName )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticRequest __req( this, "deleteField" );
  __req.add_in_arg( &_sa_fieldName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::deleteField( const char* _par_fieldName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteField(_par_fieldName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::deleteField(_par_fieldName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::getPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName );
  CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "getPatchPhysicalType" );
  __req.add_in_arg( &_sa_patchPhysicalTypeName );
  __req.add_out_arg( &_sa_patchPhysicalTypeDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::getPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::getPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::findPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName );
  CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "findPatchPhysicalType" );
  __req.add_in_arg( &_sa_patchPhysicalTypeName );
  __req.add_out_arg( &_sa_patchPhysicalTypeDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::findPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->findPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::findPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::addPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName );
  CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "addPatchPhysicalType" );
  __req.add_in_arg( &_sa_patchPhysicalTypeName );
  __req.add_out_arg( &_sa_patchPhysicalTypeDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::addPatchPhysicalType( const char* _par_patchPhysicalTypeName, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out _par_patchPhysicalTypeDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::addPatchPhysicalType(_par_patchPhysicalTypeName, _par_patchPhysicalTypeDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::deletePatchPhysicalType( const char* _par_patchPhysicalTypeName )
{
  CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName );
  CORBA::StaticRequest __req( this, "deletePatchPhysicalType" );
  __req.add_in_arg( &_sa_patchPhysicalTypeName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::deletePatchPhysicalType( const char* _par_patchPhysicalTypeName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deletePatchPhysicalType(_par_patchPhysicalTypeName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::deletePatchPhysicalType(_par_patchPhysicalTypeName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::getDictionary( const char* _par_dictName, FoamXServer::ITypeDescriptor_out _par_dictTypeDescriptor )
{
  CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName );
  CORBA::StaticAny _sa_dictTypeDescriptor( _marshaller_FoamXServer_ITypeDescriptor, &_par_dictTypeDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "getDictionary" );
  __req.add_in_arg( &_sa_dictName );
  __req.add_out_arg( &_sa_dictTypeDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::getDictionary( const char* _par_dictName, FoamXServer::ITypeDescriptor_out _par_dictTypeDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getDictionary(_par_dictName, _par_dictTypeDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::getDictionary(_par_dictName, _par_dictTypeDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::addDictionary( const char* _par_dictName, FoamXServer::ITypeDescriptor_out _par_dictTypeDescriptor )
{
  CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName );
  CORBA::StaticAny _sa_dictTypeDescriptor( _marshaller_FoamXServer_ITypeDescriptor, &_par_dictTypeDescriptor.ptr() );
  CORBA::StaticRequest __req( this, "addDictionary" );
  __req.add_in_arg( &_sa_dictName );
  __req.add_out_arg( &_sa_dictTypeDescriptor );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::addDictionary( const char* _par_dictName, FoamXServer::ITypeDescriptor_out _par_dictTypeDescriptor )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addDictionary(_par_dictName, _par_dictTypeDescriptor);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::addDictionary(_par_dictName, _par_dictTypeDescriptor);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::deleteDictionary( const char* _par_dictName )
{
  CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName );
  CORBA::StaticRequest __req( this, "deleteDictionary" );
  __req.add_in_arg( &_sa_dictName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::deleteDictionary( const char* _par_dictName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteDictionary(_par_dictName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::deleteDictionary(_par_dictName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::validate();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IApplication_stub::save()
{
  CORBA::StaticRequest __req( this, "save" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IApplication_stub_clp::save()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IApplication * _myserv = POA_FoamXServer::CaseServer::IApplication::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->save();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IApplication_stub::save();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IGeometricFieldDescriptor
 */

FoamXServer::CaseServer::IGeometricFieldDescriptor::~IGeometricFieldDescriptor()
{
}

void *
FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr
FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0")) {
      _o = new FoamXServer::CaseServer::IGeometricFieldDescriptor_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr
FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IGeometricFieldDescriptor;
}
}
class _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::~_Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IGeometricFieldDescriptor::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor::typecode()
{
  return FoamXServer::CaseServer::_tc_IGeometricFieldDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr &_obj )
{
  FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IGeometricFieldDescriptor
 */

FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::~IGeometricFieldDescriptor_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IGeometricFieldDescriptor *
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IGeometricFieldDescriptor *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::IGeometricFieldDescriptor_stub_clp ()
{
}

FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::IGeometricFieldDescriptor_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::~IGeometricFieldDescriptor_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::typeDescriptor()
{
  FoamXServer::ITypeDescriptor_ptr _res = FoamXServer::ITypeDescriptor::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );

  CORBA::StaticRequest __req( this, "_get_typeDescriptor" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::typeDescriptor()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ITypeDescriptor_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->typeDescriptor();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::typeDescriptor();
}

#endif // MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeDescriptor()
{
  FoamXServer::ITypeDescriptor_ptr _res = FoamXServer::ITypeDescriptor::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );

  CORBA::StaticRequest __req( this, "_get_fieldTypeDescriptor" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::ITypeDescriptor_ptr
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::fieldTypeDescriptor()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::ITypeDescriptor_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->fieldTypeDescriptor();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeDescriptor();
}

#endif // MICO_CONF_NO_POA

FoamXServer::CaseServer::IGeometryDescriptor_ptr FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryDescriptor()
{
  FoamXServer::CaseServer::IGeometryDescriptor_ptr _res = FoamXServer::CaseServer::IGeometryDescriptor::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IGeometryDescriptor, &_res );

  CORBA::StaticRequest __req( this, "_get_geometryDescriptor" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::CaseServer::IGeometryDescriptor_ptr
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::geometryDescriptor()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::CaseServer::IGeometryDescriptor_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->geometryDescriptor();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryDescriptor();
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_fieldTypeName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::fieldTypeName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->fieldTypeName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_fieldTypeName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::fieldTypeName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->fieldTypeName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::fieldTypeName(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryTypeName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_geometryTypeName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::geometryTypeName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->geometryTypeName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryTypeName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryTypeName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_geometryTypeName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::geometryTypeName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->geometryTypeName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::geometryTypeName(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::DimensionSet FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::dimensions()
{
  FoamXServer::DimensionSet _res;
  CORBA::StaticAny __res( _marshaller_FoamXServer_DimensionSet, &_res );

  CORBA::StaticRequest __req( this, "_get_dimensions" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::DimensionSet
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::dimensions()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::DimensionSet __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->dimensions();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::dimensions();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::dimensions( const FoamXServer::DimensionSet& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller_FoamXServer_DimensionSet, &_par__value );
  CORBA::StaticRequest __req( this, "_set_dimensions" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp::dimensions( const FoamXServer::DimensionSet& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricFieldDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->dimensions(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricFieldDescriptor_stub::dimensions(_par__value);
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IPatchDescriptor
 */

FoamXServer::CaseServer::IPatchDescriptor::~IPatchDescriptor()
{
}

void *
FoamXServer::CaseServer::IPatchDescriptor::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IPatchDescriptor_ptr
FoamXServer::CaseServer::IPatchDescriptor::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IPatchDescriptor_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IPatchDescriptor_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0")) {
      _o = new FoamXServer::CaseServer::IPatchDescriptor_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IPatchDescriptor_ptr
FoamXServer::CaseServer::IPatchDescriptor::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IPatchDescriptor;
}
}
class _Marshaller_FoamXServer_CaseServer_IPatchDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IPatchDescriptor_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IPatchDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IPatchDescriptor::~_Marshaller_FoamXServer_CaseServer_IPatchDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IPatchDescriptor::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IPatchDescriptor::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IPatchDescriptor::typecode()
{
  return FoamXServer::CaseServer::_tc_IPatchDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IPatchDescriptor;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IPatchDescriptor_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IPatchDescriptor, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IPatchDescriptor_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IPatchDescriptor, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IPatchDescriptor_ptr &_obj )
{
  FoamXServer::CaseServer::IPatchDescriptor_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IPatchDescriptor, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IPatchDescriptor
 */

FoamXServer::CaseServer::IPatchDescriptor_stub::~IPatchDescriptor_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IPatchDescriptor *
POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IPatchDescriptor *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IPatchDescriptor_stub_clp::IPatchDescriptor_stub_clp ()
{
}

FoamXServer::CaseServer::IPatchDescriptor_stub_clp::IPatchDescriptor_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IPatchDescriptor_stub_clp::~IPatchDescriptor_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchDescriptor_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchDescriptor_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchDescriptor_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchDescriptor_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchDescriptor_stub::displayName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_displayName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::displayName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->displayName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchDescriptor_stub::displayName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchDescriptor_stub::displayName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_displayName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::displayName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->displayName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchDescriptor_stub::displayName(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchDescriptor_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchDescriptor_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchDescriptor_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchDescriptor_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchDescriptor_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IGeometryDescriptor
 */

FoamXServer::CaseServer::IGeometryDescriptor::~IGeometryDescriptor()
{
}

void *
FoamXServer::CaseServer::IGeometryDescriptor::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IGeometryDescriptor_ptr
FoamXServer::CaseServer::IGeometryDescriptor::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IGeometryDescriptor_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IGeometryDescriptor_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0")) {
      _o = new FoamXServer::CaseServer::IGeometryDescriptor_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IGeometryDescriptor_ptr
FoamXServer::CaseServer::IGeometryDescriptor::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IGeometryDescriptor;
}
}
class _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IGeometryDescriptor_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IGeometryDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::~_Marshaller_FoamXServer_CaseServer_IGeometryDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IGeometryDescriptor::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IGeometryDescriptor::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor::typecode()
{
  return FoamXServer::CaseServer::_tc_IGeometryDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometryDescriptor;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IGeometryDescriptor_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometryDescriptor, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IGeometryDescriptor_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometryDescriptor, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IGeometryDescriptor_ptr &_obj )
{
  FoamXServer::CaseServer::IGeometryDescriptor_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IGeometryDescriptor, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IGeometryDescriptor
 */

FoamXServer::CaseServer::IGeometryDescriptor_stub::~IGeometryDescriptor_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IGeometryDescriptor *
POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IGeometryDescriptor *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::IGeometryDescriptor_stub_clp ()
{
}

FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::IGeometryDescriptor_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::~IGeometryDescriptor_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometryDescriptor_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometryDescriptor_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometryDescriptor_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometryDescriptor_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometryDescriptor_stub::displayName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_displayName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::displayName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->displayName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometryDescriptor_stub::displayName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometryDescriptor_stub::displayName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_displayName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::displayName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->displayName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometryDescriptor_stub::displayName(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometryDescriptor_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometryDescriptor_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometryDescriptor_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometryDescriptor_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometryDescriptor * _myserv = POA_FoamXServer::CaseServer::IGeometryDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometryDescriptor_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IPatchPhysicalTypeDescriptor
 */

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::~IPatchPhysicalTypeDescriptor()
{
}

void *
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0")) {
      _o = new FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IPatchPhysicalTypeDescriptor;
}
}
class _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::~_Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor::typecode()
{
  return FoamXServer::CaseServer::_tc_IPatchPhysicalTypeDescriptor;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr &_obj )
{
  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IPatchPhysicalTypeDescriptor
 */

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::~IPatchPhysicalTypeDescriptor_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor *
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::IPatchPhysicalTypeDescriptor_stub_clp ()
{
}

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::IPatchPhysicalTypeDescriptor_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::~IPatchPhysicalTypeDescriptor_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::name( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_name" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::name( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->name(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::name(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::displayName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_displayName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::displayName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->displayName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::displayName();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::displayName( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_displayName" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::displayName( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->displayName(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::displayName(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::description()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_description" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::description()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->description();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::description();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::description( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_description" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::description( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->description(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::description(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchType()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_patchType" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::patchType()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchType();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchType();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchType( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_patchType" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::patchType( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->patchType(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchType(_par__value);
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::parentType()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_parentType" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::parentType()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->parentType();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::parentType();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::parentType( const char* _par__value )
{
  CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value );
  CORBA::StaticRequest __req( this, "_set_parentType" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::parentType( const char* _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->parentType(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::parentType(_par__value);
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringPairList* FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchFieldTypes()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_StringPair );

  CORBA::StaticRequest __req( this, "_get_patchFieldTypes" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringPairList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringPairList*
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::patchFieldTypes()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringPairList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->patchFieldTypes();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchFieldTypes();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchFieldTypes( const FoamXServer::StringPairList& _par__value )
{
  CORBA::StaticAny _sa__value( _marshaller__seq_FoamXServer_StringPair, &_par__value );
  CORBA::StaticRequest __req( this, "_set_patchFieldTypes" );
  __req.add_in_arg( &_sa__value );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp::patchFieldTypes( const FoamXServer::StringPairList& _par__value )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor * _myserv = POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->patchFieldTypes(_par__value);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub::patchFieldTypes(_par__value);
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IGeometricField
 */

FoamXServer::CaseServer::IGeometricField::~IGeometricField()
{
}

void *
FoamXServer::CaseServer::IGeometricField::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseServer/IGeometricField:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseServer::IGeometricField_ptr
FoamXServer::CaseServer::IGeometricField::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseServer::IGeometricField_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseServer/IGeometricField:1.0" )))
      return _duplicate( (FoamXServer::CaseServer::IGeometricField_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseServer/IGeometricField:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseServer/IGeometricField:1.0")) {
      _o = new FoamXServer::CaseServer::IGeometricField_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseServer::IGeometricField_ptr
FoamXServer::CaseServer::IGeometricField::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseServer
{
CORBA::TypeCodeConst _tc_IGeometricField;
}
}
class _Marshaller_FoamXServer_CaseServer_IGeometricField : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseServer::IGeometricField_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseServer_IGeometricField();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseServer_IGeometricField::~_Marshaller_FoamXServer_CaseServer_IGeometricField()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseServer_IGeometricField::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseServer_IGeometricField::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseServer::IGeometricField::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseServer_IGeometricField::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseServer_IGeometricField::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseServer_IGeometricField::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseServer::IGeometricField::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseServer_IGeometricField::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseServer_IGeometricField::typecode()
{
  return FoamXServer::CaseServer::_tc_IGeometricField;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometricField;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseServer::IGeometricField_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometricField, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseServer::IGeometricField_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseServer_IGeometricField, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseServer::IGeometricField_ptr &_obj )
{
  FoamXServer::CaseServer::IGeometricField_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseServer_IGeometricField, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IGeometricField
 */

FoamXServer::CaseServer::IGeometricField_stub::~IGeometricField_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseServer::IGeometricField::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometricField:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseServer::IGeometricField *
POA_FoamXServer::CaseServer::IGeometricField::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseServer/IGeometricField:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseServer::IGeometricField *) p;
  }
  return NULL;
}

FoamXServer::CaseServer::IGeometricField_stub_clp::IGeometricField_stub_clp ()
{
}

FoamXServer::CaseServer::IGeometricField_stub_clp::IGeometricField_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseServer::IGeometricField_stub_clp::~IGeometricField_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CaseServer::IGeometricField_stub::name()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_name" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CaseServer::IGeometricField_stub_clp::name()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricField * _myserv = POA_FoamXServer::CaseServer::IGeometricField::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->name();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricField_stub::name();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricField_stub::getInternalFieldValue( FoamXServer::IDictionaryEntry_out _par_internalFieldValue )
{
  CORBA::StaticAny _sa_internalFieldValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_internalFieldValue.ptr() );
  CORBA::StaticRequest __req( this, "getInternalFieldValue" );
  __req.add_out_arg( &_sa_internalFieldValue );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricField_stub_clp::getInternalFieldValue( FoamXServer::IDictionaryEntry_out _par_internalFieldValue )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricField * _myserv = POA_FoamXServer::CaseServer::IGeometricField::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getInternalFieldValue(_par_internalFieldValue);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricField_stub::getInternalFieldValue(_par_internalFieldValue);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseServer::IGeometricField_stub::getPatchFieldParameters( const char* _par_patchName, FoamXServer::IDictionaryEntry_out _par_patchFieldValue )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_patchFieldValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_patchFieldValue.ptr() );
  CORBA::StaticRequest __req( this, "getPatchFieldParameters" );
  __req.add_in_arg( &_sa_patchName );
  __req.add_out_arg( &_sa_patchFieldValue );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseServer::IGeometricField_stub_clp::getPatchFieldParameters( const char* _par_patchName, FoamXServer::IDictionaryEntry_out _par_patchFieldValue )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricField * _myserv = POA_FoamXServer::CaseServer::IGeometricField::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchFieldParameters(_par_patchName, _par_patchFieldValue);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseServer::IGeometricField_stub::getPatchFieldParameters(_par_patchName, _par_patchFieldValue);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseServer::IGeometricField_stub::modified()
{
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "modified" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseServer::IGeometricField_stub_clp::modified()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseServer::IGeometricField * _myserv = POA_FoamXServer::CaseServer::IGeometricField::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->modified();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseServer::IGeometricField_stub::modified();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class ICasePostServer
 */

FoamXServer::CasePostServer::ICasePostServer::~ICasePostServer()
{
}

void *
FoamXServer::CasePostServer::ICasePostServer::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CasePostServer/ICasePostServer:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CasePostServer::ICasePostServer_ptr
FoamXServer::CasePostServer::ICasePostServer::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CasePostServer::ICasePostServer_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CasePostServer/ICasePostServer:1.0" )))
      return _duplicate( (FoamXServer::CasePostServer::ICasePostServer_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CasePostServer/ICasePostServer:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CasePostServer/ICasePostServer:1.0")) {
      _o = new FoamXServer::CasePostServer::ICasePostServer_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CasePostServer::ICasePostServer_ptr
FoamXServer::CasePostServer::ICasePostServer::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CasePostServer
{
CORBA::TypeCodeConst _tc_ICasePostServer;
}
}
class _Marshaller_FoamXServer_CasePostServer_ICasePostServer : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CasePostServer::ICasePostServer_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CasePostServer_ICasePostServer();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CasePostServer_ICasePostServer::~_Marshaller_FoamXServer_CasePostServer_ICasePostServer()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CasePostServer_ICasePostServer::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CasePostServer_ICasePostServer::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CasePostServer::ICasePostServer::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CasePostServer_ICasePostServer::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CasePostServer_ICasePostServer::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CasePostServer_ICasePostServer::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CasePostServer::ICasePostServer::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CasePostServer_ICasePostServer::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CasePostServer_ICasePostServer::typecode()
{
  return FoamXServer::CasePostServer::_tc_ICasePostServer;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CasePostServer_ICasePostServer;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CasePostServer::ICasePostServer_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CasePostServer_ICasePostServer, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CasePostServer::ICasePostServer_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CasePostServer_ICasePostServer, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CasePostServer::ICasePostServer_ptr &_obj )
{
  FoamXServer::CasePostServer::ICasePostServer_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CasePostServer_ICasePostServer, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class ICasePostServer
 */

FoamXServer::CasePostServer::ICasePostServer_stub::~ICasePostServer_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CasePostServer::ICasePostServer::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CasePostServer/ICasePostServer:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CasePostServer::ICasePostServer *
POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CasePostServer/ICasePostServer:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CasePostServer::ICasePostServer *) p;
  }
  return NULL;
}

FoamXServer::CasePostServer::ICasePostServer_stub_clp::ICasePostServer_stub_clp ()
{
}

FoamXServer::CasePostServer::ICasePostServer_stub_clp::ICasePostServer_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CasePostServer::ICasePostServer_stub_clp::~ICasePostServer_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CasePostServer::ICasePostServer_stub::caseRoot()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_caseRoot" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CasePostServer::ICasePostServer_stub_clp::caseRoot()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->caseRoot();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::caseRoot();
}

#endif // MICO_CONF_NO_POA

char* FoamXServer::CasePostServer::ICasePostServer_stub::caseName()
{
  char* _res = NULL;
  CORBA::StaticAny __res( CORBA::_stc_string, &_res );

  CORBA::StaticRequest __req( this, "_get_caseName" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

char*
FoamXServer::CasePostServer::ICasePostServer_stub_clp::caseName()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      char* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->caseName();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::caseName();
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::CasePostServer::ICasePostServer_stub::nProcs()
{
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "_get_nProcs" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::CasePostServer::ICasePostServer_stub_clp::nProcs()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->nProcs();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::nProcs();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CasePostServer::ICasePostServer_stub::availableTimeSteps()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "_get_availableTimeSteps" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CasePostServer::ICasePostServer_stub_clp::availableTimeSteps()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->availableTimeSteps();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::availableTimeSteps();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::setTime( const char* _par_timeName, CORBA::Long _par_timeIndex )
{
  CORBA::StaticAny _sa_timeName( CORBA::_stc_string, &_par_timeName );
  CORBA::StaticAny _sa_timeIndex( CORBA::_stc_long, &_par_timeIndex );
  CORBA::StaticRequest __req( this, "setTime" );
  __req.add_in_arg( &_sa_timeName );
  __req.add_in_arg( &_sa_timeIndex );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::setTime( const char* _par_timeName, CORBA::Long _par_timeIndex )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->setTime(_par_timeName, _par_timeIndex);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::setTime(_par_timeName, _par_timeIndex);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::setTimeSlave()
{
  CORBA::StaticRequest __req( this, "setTimeSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::setTimeSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->setTimeSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::setTimeSlave();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CasePostServer::ICasePostServer_stub::getPatchNames()
{
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "getPatchNames" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getPatchNames()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getPatchNames();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::getPatchNames();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getPatchNamesSlave()
{
  CORBA::StaticRequest __req( this, "getPatchNamesSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getPatchNamesSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getPatchNamesSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getPatchNamesSlave();
}

#endif // MICO_CONF_NO_POA

FoamXServer::StringList* FoamXServer::CasePostServer::ICasePostServer_stub::getFieldNames( const char* _par_type )
{
  CORBA::StaticAny _sa_type( CORBA::_stc_string, &_par_type );
  CORBA::StaticAny __res( CORBA::_stcseq_string );

  CORBA::StaticRequest __req( this, "getFieldNames" );
  __req.add_in_arg( &_sa_type );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return (FoamXServer::StringList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::StringList*
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getFieldNames( const char* _par_type )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      FoamXServer::StringList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getFieldNames(_par_type);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CasePostServer::ICasePostServer_stub::getFieldNames(_par_type);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getMeshBb( FoamXServer::Point3_out _par_min, FoamXServer::Point3_out _par_max )
{
  CORBA::StaticAny _sa_min( _marshaller__a3_float, _par_min );
  CORBA::StaticAny _sa_max( _marshaller__a3_float, _par_max );
  CORBA::StaticRequest __req( this, "getMeshBb" );
  __req.add_out_arg( &_sa_min );
  __req.add_out_arg( &_sa_max );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getMeshBb( FoamXServer::Point3_out _par_min, FoamXServer::Point3_out _par_max )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getMeshBb(_par_min, _par_max);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getMeshBb(_par_min, _par_max);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getMeshBbSlave()
{
  CORBA::StaticRequest __req( this, "getMeshBbSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getMeshBbSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getMeshBbSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getMeshBbSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getPatchMesh( const char* _par_patchName, CORBA::Double _par_creaseAngle, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_creaseAngle( CORBA::_stc_double, &_par_creaseAngle );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );
  CORBA::StaticRequest __req( this, "getPatchMesh" );
  __req.add_in_arg( &_sa_patchName );
  __req.add_in_arg( &_sa_creaseAngle );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_edges );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_edges = (FoamXServer::LongList*) _sa_edges._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getPatchMesh( const char* _par_patchName, CORBA::Double _par_creaseAngle, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getPatchMesh(_par_patchName, _par_creaseAngle, _par_points, _par_edges);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getPatchMesh(_par_patchName, _par_creaseAngle, _par_points, _par_edges);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getPatchMeshSlave()
{
  CORBA::StaticRequest __req( this, "getPatchMeshSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getPatchMeshSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getPatchMeshSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getPatchMeshSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getCutMesh( const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
  CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );
  CORBA::StaticRequest __req( this, "getCutMesh" );
  __req.add_in_arg( &_sa_basePoint );
  __req.add_in_arg( &_sa_normal );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_edges );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_edges = (FoamXServer::LongList*) _sa_edges._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getCutMesh( const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getCutMesh(_par_basePoint, _par_normal, _par_points, _par_edges);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getCutMesh(_par_basePoint, _par_normal, _par_points, _par_edges);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshSlave()
{
  CORBA::StaticRequest __req( this, "getCutMeshSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getCutMeshSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getCutMeshSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshOutline( const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
  CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );
  CORBA::StaticRequest __req( this, "getCutMeshOutline" );
  __req.add_in_arg( &_sa_basePoint );
  __req.add_in_arg( &_sa_normal );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_edges );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_edges = (FoamXServer::LongList*) _sa_edges._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getCutMeshOutline( const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_edges )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getCutMeshOutline(_par_basePoint, _par_normal, _par_points, _par_edges);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshOutline(_par_basePoint, _par_normal, _par_points, _par_edges);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshOutlineSlave()
{
  CORBA::StaticRequest __req( this, "getCutMeshOutlineSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getCutMeshOutlineSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getCutMeshOutlineSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getCutMeshOutlineSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatch( const char* _par_fieldName, const char* _par_patchName, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
  CORBA::StaticAny _sa_values( CORBA::_stcseq_float );
  CORBA::StaticRequest __req( this, "getTriPatch" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_in_arg( &_sa_patchName );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_triFaces );
  __req.add_out_arg( &_sa_values );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_triFaces = (FoamXServer::LongList*) _sa_triFaces._retn();
  _par_values = (FoamXServer::FloatList*) _sa_values._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getTriPatch( const char* _par_fieldName, const char* _par_patchName, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getTriPatch(_par_fieldName, _par_patchName, _par_points, _par_triFaces, _par_values);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatch(_par_fieldName, _par_patchName, _par_points, _par_triFaces, _par_values);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchSlave()
{
  CORBA::StaticRequest __req( this, "getTriPatchSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getTriPatchSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getTriPatchSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchVec( const char* _par_fieldName, const char* _par_patchName, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
  CORBA::StaticAny _sa_values( CORBA::_stcseq_float );
  CORBA::StaticRequest __req( this, "getTriPatchVec" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_in_arg( &_sa_patchName );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_triFaces );
  __req.add_out_arg( &_sa_values );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_triFaces = (FoamXServer::LongList*) _sa_triFaces._retn();
  _par_values = (FoamXServer::FloatList*) _sa_values._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getTriPatchVec( const char* _par_fieldName, const char* _par_patchName, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getTriPatchVec(_par_fieldName, _par_patchName, _par_points, _par_triFaces, _par_values);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchVec(_par_fieldName, _par_patchName, _par_points, _par_triFaces, _par_values);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchVecSlave()
{
  CORBA::StaticRequest __req( this, "getTriPatchVecSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::getTriPatchVecSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->getTriPatchVecSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::getTriPatchVecSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::cutPlane( const char* _par_fieldName, const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
  CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
  CORBA::StaticAny _sa_values( CORBA::_stcseq_float );
  CORBA::StaticRequest __req( this, "cutPlane" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_in_arg( &_sa_basePoint );
  __req.add_in_arg( &_sa_normal );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_triFaces );
  __req.add_out_arg( &_sa_values );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_triFaces = (FoamXServer::LongList*) _sa_triFaces._retn();
  _par_values = (FoamXServer::FloatList*) _sa_values._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::cutPlane( const char* _par_fieldName, const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->cutPlane(_par_fieldName, _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::cutPlane(_par_fieldName, _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneSlave()
{
  CORBA::StaticRequest __req( this, "cutPlaneSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::cutPlaneSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->cutPlaneSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneVec( const char* _par_fieldName, const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName );
  CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
  CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
  CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
  CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
  CORBA::StaticAny _sa_values( CORBA::_stcseq_float );
  CORBA::StaticRequest __req( this, "cutPlaneVec" );
  __req.add_in_arg( &_sa_fieldName );
  __req.add_in_arg( &_sa_basePoint );
  __req.add_in_arg( &_sa_normal );
  __req.add_out_arg( &_sa_points );
  __req.add_out_arg( &_sa_triFaces );
  __req.add_out_arg( &_sa_values );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  _par_points = (FoamXServer::FloatList*) _sa_points._retn();
  _par_triFaces = (FoamXServer::LongList*) _sa_triFaces._retn();
  _par_values = (FoamXServer::FloatList*) _sa_values._retn();
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::cutPlaneVec( const char* _par_fieldName, const FoamXServer::Point3 _par_basePoint, const FoamXServer::Point3 _par_normal, FoamXServer::FloatList_out _par_points, FoamXServer::LongList_out _par_triFaces, FoamXServer::FloatList_out _par_values )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->cutPlaneVec(_par_fieldName, _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneVec(_par_fieldName, _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneVecSlave()
{
  CORBA::StaticRequest __req( this, "cutPlaneVecSlave" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::cutPlaneVecSlave()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->cutPlaneVecSlave();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::cutPlaneVecSlave();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CasePostServer::ICasePostServer_stub::close()
{
  CORBA::StaticRequest __req( this, "close" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CasePostServer::ICasePostServer_stub_clp::close()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CasePostServer::ICasePostServer * _myserv = POA_FoamXServer::CasePostServer::ICasePostServer::_narrow (_serv);
    if (_myserv) {
      _myserv->close();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CasePostServer::ICasePostServer_stub::close();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class ICaseBrowser
 */

FoamXServer::CaseBrowser::ICaseBrowser::~ICaseBrowser()
{
}

void *
FoamXServer::CaseBrowser::ICaseBrowser::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::CaseBrowser::ICaseBrowser_ptr
FoamXServer::CaseBrowser::ICaseBrowser::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::CaseBrowser::ICaseBrowser_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0" )))
      return _duplicate( (FoamXServer::CaseBrowser::ICaseBrowser_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0")) {
      _o = new FoamXServer::CaseBrowser::ICaseBrowser_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::CaseBrowser::ICaseBrowser_ptr
FoamXServer::CaseBrowser::ICaseBrowser::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace CaseBrowser
{
CORBA::TypeCodeConst _tc_ICaseBrowser;
}
}
class _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::CaseBrowser::ICaseBrowser_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_CaseBrowser_ICaseBrowser();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::~_Marshaller_FoamXServer_CaseBrowser_ICaseBrowser()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::CaseBrowser::ICaseBrowser::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::CaseBrowser::ICaseBrowser::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser::typecode()
{
  return FoamXServer::CaseBrowser::_tc_ICaseBrowser;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseBrowser_ICaseBrowser;

void
operator<<=( CORBA::Any &_a, const FoamXServer::CaseBrowser::ICaseBrowser_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseBrowser_ICaseBrowser, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::CaseBrowser::ICaseBrowser_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_CaseBrowser_ICaseBrowser, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::CaseBrowser::ICaseBrowser_ptr &_obj )
{
  FoamXServer::CaseBrowser::ICaseBrowser_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_CaseBrowser_ICaseBrowser, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class ICaseBrowser
 */

FoamXServer::CaseBrowser::ICaseBrowser_stub::~ICaseBrowser_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::CaseBrowser::ICaseBrowser *
POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::CaseBrowser::ICaseBrowser *) p;
  }
  return NULL;
}

FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::ICaseBrowser_stub_clp ()
{
}

FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::ICaseBrowser_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::~ICaseBrowser_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::CaseServer::IFoamProperties_ptr FoamXServer::CaseBrowser::ICaseBrowser_stub::foamProperties()
{
  FoamXServer::CaseServer::IFoamProperties_ptr _res = FoamXServer::CaseServer::IFoamProperties::_nil();
  CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IFoamProperties, &_res );

  CORBA::StaticRequest __req( this, "_get_foamProperties" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

FoamXServer::CaseServer::IFoamProperties_ptr
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::foamProperties()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      FoamXServer::CaseServer::IFoamProperties_ptr __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->foamProperties();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::foamProperties();
}

#endif // MICO_CONF_NO_POA

FoamXServer::CaseDescriptorList* FoamXServer::CaseBrowser::ICaseBrowser_stub::cases()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_CaseDescriptor );

  CORBA::StaticRequest __req( this, "_get_cases" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::CaseDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::CaseDescriptorList*
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::cases()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      FoamXServer::CaseDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->cases();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::cases();
}

#endif // MICO_CONF_NO_POA

FoamXServer::JobDescriptorList* FoamXServer::CaseBrowser::ICaseBrowser_stub::runningJobs()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_JobDescriptor );

  CORBA::StaticRequest __req( this, "_get_runningJobs" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::JobDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::JobDescriptorList*
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::runningJobs()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      FoamXServer::JobDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->runningJobs();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::runningJobs();
}

#endif // MICO_CONF_NO_POA

FoamXServer::JobDescriptorList* FoamXServer::CaseBrowser::ICaseBrowser_stub::finishedJobs()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_JobDescriptor );

  CORBA::StaticRequest __req( this, "_get_finishedJobs" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::JobDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::JobDescriptorList*
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::finishedJobs()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      FoamXServer::JobDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->finishedJobs();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::finishedJobs();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::getEnv( const char* _par_envName, CORBA::String_out _par_hostName )
{
  CORBA::StaticAny _sa_envName( CORBA::_stc_string, &_par_envName );
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName.ptr() );
  CORBA::StaticRequest __req( this, "getEnv" );
  __req.add_in_arg( &_sa_envName );
  __req.add_out_arg( &_sa_hostName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::getEnv( const char* _par_envName, CORBA::String_out _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getEnv(_par_envName, _par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::getEnv(_par_envName, _par_hostName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::getHostName( CORBA::String_out _par_hostName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName.ptr() );
  CORBA::StaticRequest __req( this, "getHostName" );
  __req.add_out_arg( &_sa_hostName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::getHostName( CORBA::String_out _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getHostName(_par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::getHostName(_par_hostName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::getUserName( CORBA::String_out _par_userName )
{
  CORBA::StaticAny _sa_userName( CORBA::_stc_string, &_par_userName.ptr() );
  CORBA::StaticRequest __req( this, "getUserName" );
  __req.add_out_arg( &_sa_userName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::getUserName( CORBA::String_out _par_userName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->getUserName(_par_userName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::getUserName(_par_userName);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::CaseBrowser::ICaseBrowser_stub::fileModificationDate( const char* _par_fileName )
{
  CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName );
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "fileModificationDate" );
  __req.add_in_arg( &_sa_fileName );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::fileModificationDate( const char* _par_fileName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->fileModificationDate(_par_fileName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::fileModificationDate(_par_fileName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::readFile( const char* _par_fileName, CORBA::String_out _par_contents )
{
  CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName );
  CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents.ptr() );
  CORBA::StaticRequest __req( this, "readFile" );
  __req.add_in_arg( &_sa_fileName );
  __req.add_out_arg( &_sa_contents );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::readFile( const char* _par_fileName, CORBA::String_out _par_contents )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->readFile(_par_fileName, _par_contents);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::readFile(_par_fileName, _par_contents);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::writeFile( const char* _par_fileName, const char* _par_contents )
{
  CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName );
  CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents );
  CORBA::StaticRequest __req( this, "writeFile" );
  __req.add_in_arg( &_sa_fileName );
  __req.add_in_arg( &_sa_contents );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::writeFile( const char* _par_fileName, const char* _par_contents )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->writeFile(_par_fileName, _par_contents);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::writeFile(_par_fileName, _par_contents);
}

#endif // MICO_CONF_NO_POA

CORBA::Long FoamXServer::CaseBrowser::ICaseBrowser_stub::invokeUtility( const char* _par_hostName, const char* _par_utilityName, const FoamXServer::StringList& _par_arguments, const char* _par_logName, CORBA::Boolean _par_backGround )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticAny _sa_utilityName( CORBA::_stc_string, &_par_utilityName );
  CORBA::StaticAny _sa_arguments( CORBA::_stcseq_string, &_par_arguments );
  CORBA::StaticAny _sa_logName( CORBA::_stc_string, &_par_logName );
  CORBA::StaticAny _sa_backGround( CORBA::_stc_boolean, &_par_backGround );
  CORBA::Long _res;
  CORBA::StaticAny __res( CORBA::_stc_long, &_res );

  CORBA::StaticRequest __req( this, "invokeUtility" );
  __req.add_in_arg( &_sa_hostName );
  __req.add_in_arg( &_sa_utilityName );
  __req.add_in_arg( &_sa_arguments );
  __req.add_in_arg( &_sa_logName );
  __req.add_in_arg( &_sa_backGround );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Long
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::invokeUtility( const char* _par_hostName, const char* _par_utilityName, const FoamXServer::StringList& _par_arguments, const char* _par_logName, CORBA::Boolean _par_backGround )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Long __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->invokeUtility(_par_hostName, _par_utilityName, _par_arguments, _par_logName, _par_backGround);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::invokeUtility(_par_hostName, _par_utilityName, _par_arguments, _par_logName, _par_backGround);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::refreshCaseList()
{
  CORBA::StaticRequest __req( this, "refreshCaseList" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::refreshCaseList()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->refreshCaseList();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::refreshCaseList();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::addToCaseList( const char* _par_rootDir )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticRequest __req( this, "addToCaseList" );
  __req.add_in_arg( &_sa_rootDir );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::addToCaseList( const char* _par_rootDir )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addToCaseList(_par_rootDir);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::addToCaseList(_par_rootDir);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::openCase( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticRequest __req( this, "openCase" );
  __req.add_in_arg( &_sa_caseDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::openCase( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->openCase(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::openCase(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::newCase( const char* _par_rootDir, const char* _par_caseName, const char* _par_app )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app );
  CORBA::StaticRequest __req( this, "newCase" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_in_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::newCase( const char* _par_rootDir, const char* _par_caseName, const char* _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->newCase(_par_rootDir, _par_caseName, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::newCase(_par_rootDir, _par_caseName, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::importCase( const char* _par_rootDir, const char* _par_caseName, const char* _par_app )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app );
  CORBA::StaticRequest __req( this, "importCase" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_in_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::importCase( const char* _par_rootDir, const char* _par_caseName, const char* _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->importCase(_par_rootDir, _par_caseName, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::importCase(_par_rootDir, _par_caseName, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::deleteCase( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticRequest __req( this, "deleteCase" );
  __req.add_in_arg( &_sa_caseDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::deleteCase( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->deleteCase(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::deleteCase(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::cloneCase( const FoamXServer::CaseDescriptor& _par_caseDesc, const char* _par_newCaseRootDir, const char* _par_newCaseName, const char* _par_newAppClassName, const char* _par_timeSel )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticAny _sa_newCaseRootDir( CORBA::_stc_string, &_par_newCaseRootDir );
  CORBA::StaticAny _sa_newCaseName( CORBA::_stc_string, &_par_newCaseName );
  CORBA::StaticAny _sa_newAppClassName( CORBA::_stc_string, &_par_newAppClassName );
  CORBA::StaticAny _sa_timeSel( CORBA::_stc_string, &_par_timeSel );
  CORBA::StaticRequest __req( this, "cloneCase" );
  __req.add_in_arg( &_sa_caseDesc );
  __req.add_in_arg( &_sa_newCaseRootDir );
  __req.add_in_arg( &_sa_newCaseName );
  __req.add_in_arg( &_sa_newAppClassName );
  __req.add_in_arg( &_sa_timeSel );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::cloneCase( const FoamXServer::CaseDescriptor& _par_caseDesc, const char* _par_newCaseRootDir, const char* _par_newCaseName, const char* _par_newAppClassName, const char* _par_timeSel )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->cloneCase(_par_caseDesc, _par_newCaseRootDir, _par_newCaseName, _par_newAppClassName, _par_timeSel);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::cloneCase(_par_caseDesc, _par_newCaseRootDir, _par_newCaseName, _par_newAppClassName, _par_timeSel);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseBrowser::ICaseBrowser_stub::getCaseServerReference( const char* _par_rootDir, const char* _par_caseName, FoamXServer::CaseServer::ICaseServer_out _par_caseObj )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_caseObj( _marshaller_FoamXServer_CaseServer_ICaseServer, &_par_caseObj.ptr() );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "getCaseServerReference" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_out_arg( &_sa_caseObj );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::getCaseServerReference( const char* _par_rootDir, const char* _par_caseName, FoamXServer::CaseServer::ICaseServer_out _par_caseObj )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getCaseServerReference(_par_rootDir, _par_caseName, _par_caseObj);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::getCaseServerReference(_par_rootDir, _par_caseName, _par_caseObj);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::openCasePost( const FoamXServer::CaseDescriptor& _par_caseDesc, CORBA::Long _par_nProcs )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticAny _sa_nProcs( CORBA::_stc_long, &_par_nProcs );
  CORBA::StaticRequest __req( this, "openCasePost" );
  __req.add_in_arg( &_sa_caseDesc );
  __req.add_in_arg( &_sa_nProcs );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::openCasePost( const FoamXServer::CaseDescriptor& _par_caseDesc, CORBA::Long _par_nProcs )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->openCasePost(_par_caseDesc, _par_nProcs);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::openCasePost(_par_caseDesc, _par_nProcs);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseBrowser::ICaseBrowser_stub::getCasePostServerReference( const char* _par_rootDir, const char* _par_caseName, CORBA::Long _par_nProcs, FoamXServer::CasePostServer::ICasePostServer_out _par_caseObj )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_nProcs( CORBA::_stc_long, &_par_nProcs );
  CORBA::StaticAny _sa_caseObj( _marshaller_FoamXServer_CasePostServer_ICasePostServer, &_par_caseObj.ptr() );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "getCasePostServerReference" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_in_arg( &_sa_nProcs );
  __req.add_out_arg( &_sa_caseObj );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::getCasePostServerReference( const char* _par_rootDir, const char* _par_caseName, CORBA::Long _par_nProcs, FoamXServer::CasePostServer::ICasePostServer_out _par_caseObj )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getCasePostServerReference(_par_rootDir, _par_caseName, _par_nProcs, _par_caseObj);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::getCasePostServerReference(_par_rootDir, _par_caseName, _par_nProcs, _par_caseObj);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseBrowser::ICaseBrowser_stub::caseLocked( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "caseLocked" );
  __req.add_in_arg( &_sa_caseDesc );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::caseLocked( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->caseLocked(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::caseLocked(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::unlockCase( const char* _par_rootDir, const char* _par_caseName )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticRequest __req( this, "unlockCase" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::unlockCase( const char* _par_rootDir, const char* _par_caseName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->unlockCase(_par_rootDir, _par_caseName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::unlockCase(_par_rootDir, _par_caseName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::unlockCaseDescriptor( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticRequest __req( this, "unlockCaseDescriptor" );
  __req.add_in_arg( &_sa_caseDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::unlockCaseDescriptor( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->unlockCaseDescriptor(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::unlockCaseDescriptor(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::addCase( const char* _par_rootDir, const char* _par_rawRootDir, const char* _par_caseName, const char* _par_app )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app );
  CORBA::StaticRequest __req( this, "addCase" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_rawRootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_in_arg( &_sa_app );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::addCase( const char* _par_rootDir, const char* _par_rawRootDir, const char* _par_caseName, const char* _par_app )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->addCase(_par_rootDir, _par_rawRootDir, _par_caseName, _par_app);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::addCase(_par_rootDir, _par_rawRootDir, _par_caseName, _par_app);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::caseOpen( const char* _par_rootDir, const char* _par_caseName )
{
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticRequest __req( this, "caseOpen" );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::caseOpen( const char* _par_rootDir, const char* _par_caseName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->caseOpen(_par_rootDir, _par_caseName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::caseOpen(_par_rootDir, _par_caseName);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::CaseBrowser::ICaseBrowser_stub::isCaseInError( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "isCaseInError" );
  __req.add_in_arg( &_sa_caseDesc );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::isCaseInError( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->isCaseInError(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::CaseBrowser::ICaseBrowser_stub::isCaseInError(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::caseIsInError( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
  CORBA::StaticRequest __req( this, "caseIsInError" );
  __req.add_in_arg( &_sa_caseDesc );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::caseIsInError( const FoamXServer::CaseDescriptor& _par_caseDesc )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->caseIsInError(_par_caseDesc);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::caseIsInError(_par_caseDesc);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::refreshJobsLists()
{
  CORBA::StaticRequest __req( this, "refreshJobsLists" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::refreshJobsLists()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->refreshJobsLists();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::refreshJobsLists();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::checkRunningJobs()
{
  CORBA::StaticRequest __req( this, "checkRunningJobs" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::checkRunningJobs()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->checkRunningJobs();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::checkRunningJobs();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeRunningJobs()
{
  CORBA::StaticRequest __req( this, "purgeRunningJobs" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::purgeRunningJobs()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->purgeRunningJobs();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeRunningJobs();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeFinishedJob( const FoamXServer::JobID& _par_jobID )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticRequest __req( this, "purgeFinishedJob" );
  __req.add_in_arg( &_sa_jobID );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::purgeFinishedJob( const FoamXServer::JobID& _par_jobID )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->purgeFinishedJob(_par_jobID);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeFinishedJob(_par_jobID);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeFinishedJobs( CORBA::Long _par_nDays )
{
  CORBA::StaticAny _sa_nDays( CORBA::_stc_long, &_par_nDays );
  CORBA::StaticRequest __req( this, "purgeFinishedJobs" );
  __req.add_in_arg( &_sa_nDays );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::purgeFinishedJobs( CORBA::Long _par_nDays )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->purgeFinishedJobs(_par_nDays);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::purgeFinishedJobs(_par_nDays);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::kill( const FoamXServer::JobID& _par_jobID )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticRequest __req( this, "kill" );
  __req.add_in_arg( &_sa_jobID );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::kill( const FoamXServer::JobID& _par_jobID )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->kill(_par_jobID);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::kill(_par_jobID);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::suspend( const FoamXServer::JobID& _par_jobID )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticRequest __req( this, "suspend" );
  __req.add_in_arg( &_sa_jobID );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::suspend( const FoamXServer::JobID& _par_jobID )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->suspend(_par_jobID);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::suspend(_par_jobID);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::cont( const FoamXServer::JobID& _par_jobID )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticRequest __req( this, "cont" );
  __req.add_in_arg( &_sa_jobID );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::cont( const FoamXServer::JobID& _par_jobID )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->cont(_par_jobID);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::cont(_par_jobID);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::end( const FoamXServer::JobID& _par_jobID, const char* _par_rootDir, const char* _par_caseName, CORBA::Boolean _par_now )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir );
  CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName );
  CORBA::StaticAny _sa_now( CORBA::_stc_boolean, &_par_now );
  CORBA::StaticRequest __req( this, "end" );
  __req.add_in_arg( &_sa_jobID );
  __req.add_in_arg( &_sa_rootDir );
  __req.add_in_arg( &_sa_caseName );
  __req.add_in_arg( &_sa_now );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::end( const FoamXServer::JobID& _par_jobID, const char* _par_rootDir, const char* _par_caseName, CORBA::Boolean _par_now )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->end(_par_jobID, _par_rootDir, _par_caseName, _par_now);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::end(_par_jobID, _par_rootDir, _par_caseName, _par_now);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::setStatus( const FoamXServer::JobID& _par_jobID, FoamXServer::JobStatus _par_jobStatus )
{
  CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
  CORBA::StaticAny _sa_jobStatus( _marshaller_FoamXServer_JobStatus, &_par_jobStatus );
  CORBA::StaticRequest __req( this, "setStatus" );
  __req.add_in_arg( &_sa_jobID );
  __req.add_in_arg( &_sa_jobStatus );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::setStatus( const FoamXServer::JobID& _par_jobID, FoamXServer::JobStatus _par_jobStatus )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->setStatus(_par_jobID, _par_jobStatus);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::setStatus(_par_jobID, _par_jobStatus);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::validate();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::save()
{
  CORBA::StaticRequest __req( this, "save" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::save()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->save();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::save();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::CaseBrowser::ICaseBrowser_stub::close()
{
  CORBA::StaticRequest __req( this, "close" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::CaseBrowser::ICaseBrowser_stub_clp::close()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::CaseBrowser::ICaseBrowser * _myserv = POA_FoamXServer::CaseBrowser::ICaseBrowser::_narrow (_serv);
    if (_myserv) {
      _myserv->close();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::CaseBrowser::ICaseBrowser_stub::close();
}

#endif // MICO_CONF_NO_POA


/*
 * Base interface for class IHostBrowser
 */

FoamXServer::HostBrowser::IHostBrowser::~IHostBrowser()
{
}

void *
FoamXServer::HostBrowser::IHostBrowser::_narrow_helper( const char *_repoid )
{
  if( strcmp( _repoid, "IDL:FoamXServer/HostBrowser/IHostBrowser:1.0" ) == 0 )
    return (void *)this;
  return NULL;
}

FoamXServer::HostBrowser::IHostBrowser_ptr
FoamXServer::HostBrowser::IHostBrowser::_narrow( CORBA::Object_ptr _obj )
{
  FoamXServer::HostBrowser::IHostBrowser_ptr _o;
  if( !CORBA::is_nil( _obj ) ) {
    void *_p;
    if( (_p = _obj->_narrow_helper( "IDL:FoamXServer/HostBrowser/IHostBrowser:1.0" )))
      return _duplicate( (FoamXServer::HostBrowser::IHostBrowser_ptr) _p );
    if (!strcmp (_obj->_repoid(), "IDL:FoamXServer/HostBrowser/IHostBrowser:1.0") || _obj->_is_a_remote ("IDL:FoamXServer/HostBrowser/IHostBrowser:1.0")) {
      _o = new FoamXServer::HostBrowser::IHostBrowser_stub;
      _o->CORBA::Object::operator=( *_obj );
      return _o;
    }
  }
  return _nil();
}

FoamXServer::HostBrowser::IHostBrowser_ptr
FoamXServer::HostBrowser::IHostBrowser::_narrow( CORBA::AbstractBase_ptr _obj )
{
  return _narrow (_obj->_to_object());
}

namespace FoamXServer
{
namespace HostBrowser
{
CORBA::TypeCodeConst _tc_IHostBrowser;
}
}
class _Marshaller_FoamXServer_HostBrowser_IHostBrowser : public ::CORBA::StaticTypeInfo {
    typedef FoamXServer::HostBrowser::IHostBrowser_ptr _MICO_T;
  public:
    ~_Marshaller_FoamXServer_HostBrowser_IHostBrowser();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    void release (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller_FoamXServer_HostBrowser_IHostBrowser::~_Marshaller_FoamXServer_HostBrowser_IHostBrowser()
{
}

::CORBA::StaticValueType _Marshaller_FoamXServer_HostBrowser_IHostBrowser::create() const
{
  return (StaticValueType) new _MICO_T( 0 );
}

void _Marshaller_FoamXServer_HostBrowser_IHostBrowser::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = ::FoamXServer::HostBrowser::IHostBrowser::_duplicate( *(_MICO_T*) s );
}

void _Marshaller_FoamXServer_HostBrowser_IHostBrowser::free( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
  delete (_MICO_T*) v;
}

void _Marshaller_FoamXServer_HostBrowser_IHostBrowser::release( StaticValueType v ) const
{
  ::CORBA::release( *(_MICO_T *) v );
}

::CORBA::Boolean _Marshaller_FoamXServer_HostBrowser_IHostBrowser::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj;
  if (!::CORBA::_stc_Object->demarshal(dc, &obj))
    return FALSE;
  *(_MICO_T *) v = ::FoamXServer::HostBrowser::IHostBrowser::_narrow( obj );
  ::CORBA::Boolean ret = ::CORBA::is_nil (obj) || !::CORBA::is_nil (*(_MICO_T *)v);
  ::CORBA::release (obj);
  return ret;
}

void _Marshaller_FoamXServer_HostBrowser_IHostBrowser::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::Object_ptr obj = *(_MICO_T *) v;
  ::CORBA::_stc_Object->marshal( ec, &obj );
}

::CORBA::TypeCode_ptr _Marshaller_FoamXServer_HostBrowser_IHostBrowser::typecode()
{
  return FoamXServer::HostBrowser::_tc_IHostBrowser;
}

::CORBA::StaticTypeInfo *_marshaller_FoamXServer_HostBrowser_IHostBrowser;

void
operator<<=( CORBA::Any &_a, const FoamXServer::HostBrowser::IHostBrowser_ptr _obj )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_HostBrowser_IHostBrowser, &_obj);
  _a.from_static_any (_sa);
}

void
operator<<=( CORBA::Any &_a, FoamXServer::HostBrowser::IHostBrowser_ptr* _obj_ptr )
{
  CORBA::StaticAny _sa (_marshaller_FoamXServer_HostBrowser_IHostBrowser, _obj_ptr);
  _a.from_static_any (_sa);
  CORBA::release (*_obj_ptr);
}

CORBA::Boolean
operator>>=( const CORBA::Any &_a, FoamXServer::HostBrowser::IHostBrowser_ptr &_obj )
{
  FoamXServer::HostBrowser::IHostBrowser_ptr *p;
  if (_a.to_static_any (_marshaller_FoamXServer_HostBrowser_IHostBrowser, (void *&)p)) {
    _obj = *p;
    return TRUE;
  }
  return FALSE;
}


/*
 * Stub interface for class IHostBrowser
 */

FoamXServer::HostBrowser::IHostBrowser_stub::~IHostBrowser_stub()
{
}

#ifndef MICO_CONF_NO_POA

void *
POA_FoamXServer::HostBrowser::IHostBrowser::_narrow_helper (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/HostBrowser/IHostBrowser:1.0") == 0) {
    return (void *) this;
  }
  return NULL;
}

POA_FoamXServer::HostBrowser::IHostBrowser *
POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (PortableServer::Servant serv) 
{
  void * p;
  if ((p = serv->_narrow_helper ("IDL:FoamXServer/HostBrowser/IHostBrowser:1.0")) != NULL) {
    serv->_add_ref ();
    return (POA_FoamXServer::HostBrowser::IHostBrowser *) p;
  }
  return NULL;
}

FoamXServer::HostBrowser::IHostBrowser_stub_clp::IHostBrowser_stub_clp ()
{
}

FoamXServer::HostBrowser::IHostBrowser_stub_clp::IHostBrowser_stub_clp (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
  : CORBA::Object(*obj), PortableServer::StubBase(poa)
{
}

FoamXServer::HostBrowser::IHostBrowser_stub_clp::~IHostBrowser_stub_clp ()
{
}

#endif // MICO_CONF_NO_POA

FoamXServer::HostDescriptorList* FoamXServer::HostBrowser::IHostBrowser_stub::hosts()
{
  CORBA::StaticAny __res( _marshaller__seq_FoamXServer_HostDescriptor );

  CORBA::StaticRequest __req( this, "_get_hosts" );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    0);
  return (FoamXServer::HostDescriptorList*) __res._retn();
}


#ifndef MICO_CONF_NO_POA

FoamXServer::HostDescriptorList*
FoamXServer::HostBrowser::IHostBrowser_stub_clp::hosts()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      FoamXServer::HostDescriptorList* __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->hosts();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::HostBrowser::IHostBrowser_stub::hosts();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::refreshHostList()
{
  CORBA::StaticRequest __req( this, "refreshHostList" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::refreshHostList()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->refreshHostList();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::refreshHostList();
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::HostBrowser::IHostBrowser_stub::isHostAlive( const char* _par_hostName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "isHostAlive" );
  __req.add_in_arg( &_sa_hostName );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::HostBrowser::IHostBrowser_stub_clp::isHostAlive( const char* _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->isHostAlive(_par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::HostBrowser::IHostBrowser_stub::isHostAlive(_par_hostName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::hostIsAlive( const char* _par_hostName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticRequest __req( this, "hostIsAlive" );
  __req.add_in_arg( &_sa_hostName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::hostIsAlive( const char* _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->hostIsAlive(_par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::hostIsAlive(_par_hostName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::hostIsDead( const char* _par_hostName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticRequest __req( this, "hostIsDead" );
  __req.add_in_arg( &_sa_hostName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::hostIsDead( const char* _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->hostIsDead(_par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::hostIsDead(_par_hostName);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::openCaseBrowser( const char* _par_hostName )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticRequest __req( this, "openCaseBrowser" );
  __req.add_in_arg( &_sa_hostName );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::openCaseBrowser( const char* _par_hostName )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->openCaseBrowser(_par_hostName);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::openCaseBrowser(_par_hostName);
}

#endif // MICO_CONF_NO_POA

CORBA::Boolean FoamXServer::HostBrowser::IHostBrowser_stub::getCaseBrowserReference( const char* _par_hostName, FoamXServer::CaseBrowser::ICaseBrowser_out _par_browserObj )
{
  CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );
  CORBA::StaticAny _sa_browserObj( _marshaller_FoamXServer_CaseBrowser_ICaseBrowser, &_par_browserObj.ptr() );
  CORBA::Boolean _res;
  CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );

  CORBA::StaticRequest __req( this, "getCaseBrowserReference" );
  __req.add_in_arg( &_sa_hostName );
  __req.add_out_arg( &_sa_browserObj );
  __req.set_result( &__res );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXSYSError, "IDL:FoamXServer/FoamXSYSError:1.0",
    0);
  return _res;
}


#ifndef MICO_CONF_NO_POA

CORBA::Boolean
FoamXServer::HostBrowser::IHostBrowser_stub_clp::getCaseBrowserReference( const char* _par_hostName, FoamXServer::CaseBrowser::ICaseBrowser_out _par_browserObj )
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      CORBA::Boolean __res;

      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        __res = _myserv->getCaseBrowserReference(_par_hostName, _par_browserObj);
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return __res;
    }
    _postinvoke ();
  }

  return FoamXServer::HostBrowser::IHostBrowser_stub::getCaseBrowserReference(_par_hostName, _par_browserObj);
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::validate()
{
  CORBA::StaticRequest __req( this, "validate" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::validate()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->validate();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::validate();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::save()
{
  CORBA::StaticRequest __req( this, "save" );

  __req.invoke();

  mico_sii_throw( &__req, 
    _marshaller_FoamXServer_FoamXError, "IDL:FoamXServer/FoamXError:1.0",
    _marshaller_FoamXServer_FoamXIOError, "IDL:FoamXServer/FoamXIOError:1.0",
    _marshaller_FoamXServer_ValidationError, "IDL:FoamXServer/ValidationError:1.0",
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::save()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      #ifdef HAVE_EXCEPTIONS
      try {
      #endif
        _myserv->save();
      #ifdef HAVE_EXCEPTIONS
      }
      catch (...) {
        _myserv->_remove_ref();
        _postinvoke();
        throw;
      }
      #endif

      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::save();
}

#endif // MICO_CONF_NO_POA

void FoamXServer::HostBrowser::IHostBrowser_stub::close()
{
  CORBA::StaticRequest __req( this, "close" );

  __req.oneway();

  mico_sii_throw( &__req, 
    0);
}


#ifndef MICO_CONF_NO_POA

void
FoamXServer::HostBrowser::IHostBrowser_stub_clp::close()
{
  PortableServer::Servant _serv = _preinvoke ();
  if (_serv) {
    POA_FoamXServer::HostBrowser::IHostBrowser * _myserv = POA_FoamXServer::HostBrowser::IHostBrowser::_narrow (_serv);
    if (_myserv) {
      _myserv->close();
      _myserv->_remove_ref();
      _postinvoke ();
      return;
    }
    _postinvoke ();
  }

  FoamXServer::HostBrowser::IHostBrowser_stub::close();
}

#endif // MICO_CONF_NO_POA

class _Marshaller__seq_FoamXServer_FoamXAny : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_FoamXAny();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_FoamXAny::~_Marshaller__seq_FoamXServer_FoamXAny()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_FoamXAny::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_FoamXAny::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_FoamXAny::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_FoamXAny::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_FoamXAny->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_FoamXAny::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_FoamXAny->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_FoamXAny::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000034020000010000000f00000024020000010000001d00"
    "000049444c3a466f616d585365727665722f466f616d58416e793a312e30"
    "0000000009000000466f616d58416e790000000002000000050000007479"
    "70650000000011000000c3010000010000001e00000049444c3a466f616d"
    "585365727665722f466f616d58547970653a312e300000000a000000466f"
    "616d5854797065000000150000000f000000547970655f556e646566696e"
    "656400000d000000547970655f426f6f6c65616e000000000b0000005479"
    "70655f4c6162656c00000c000000547970655f5363616c6172000a000000"
    "547970655f436861720000000a000000547970655f576f72640000000c00"
    "0000547970655f537472696e67000d000000547970655f526f6f74446972"
    "0000000011000000547970655f526f6f74416e6443617365000000000e00"
    "0000547970655f436173654e616d650000000e000000547970655f486f73"
    "744e616d650000000a000000547970655f46696c650000000f0000005479"
    "70655f4469726563746f727900000a000000547970655f54696d65000000"
    "12000000547970655f44696d656e73696f6e5365740000000f0000005479"
    "70655f46697865644c69737400000a000000547970655f4c697374000000"
    "10000000547970655f44696374696f6e617279000f000000547970655f53"
    "656c656374696f6e00000e000000547970655f436f6d706f756e64000000"
    "0b000000547970655f4669656c6400000600000076616c75650000000b00"
    "000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_FoamXAny::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_FoamXAny;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_FoamXAny, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_FoamXAny, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_FoamXAny, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_ITypeDescriptor : public ::CORBA::StaticTypeInfo {
    typedef IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_ITypeDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_ITypeDescriptor::~_Marshaller__seq_FoamXServer_ITypeDescriptor()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_ITypeDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_ITypeDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_ITypeDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_ITypeDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_ITypeDescriptor->demarshal( dc, &(*(_MICO_T*)v)[i]._for_demarshal() ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_ITypeDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_ITypeDescriptor->marshal( ec, &(*(_MICO_T*)v)[i].inout() );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_ITypeDescriptor::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000050000000010000000e00000040000000010000002400"
    "000049444c3a466f616d585365727665722f495479706544657363726970"
    "746f723a312e300010000000495479706544657363726970746f72000000"
    "0000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_ITypeDescriptor::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_ITypeDescriptor;

void operator<<=( CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_ITypeDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_ITypeDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_ITypeDescriptor, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_IDictionaryEntry : public ::CORBA::StaticTypeInfo {
    typedef IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_IDictionaryEntry();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_IDictionaryEntry::~_Marshaller__seq_FoamXServer_IDictionaryEntry()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_IDictionaryEntry::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_IDictionaryEntry::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_IDictionaryEntry::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_IDictionaryEntry::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_IDictionaryEntry->demarshal( dc, &(*(_MICO_T*)v)[i]._for_demarshal() ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_IDictionaryEntry::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_IDictionaryEntry->marshal( ec, &(*(_MICO_T*)v)[i].inout() );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_IDictionaryEntry::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000058000000010000000e00000045000000010000002500"
    "000049444c3a466f616d585365727665722f4944696374696f6e61727945"
    "6e7472793a312e3000000000110000004944696374696f6e617279456e74"
    "72790000000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_IDictionaryEntry::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_IDictionaryEntry;

void operator<<=( CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_IDictionaryEntry, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_IDictionaryEntry, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_IDictionaryEntry, (void *&)_s);
}


class _Marshaller__a3_float : public ::CORBA::StaticTypeInfo {
    typedef CORBA::Float _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__a3_float();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__a3_float::~_Marshaller__a3_float()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__a3_float::create() const
{
  return (StaticValueType) new _MICO_T[ 3 ];
}

void _Marshaller__a3_float::assign( StaticValueType d, const StaticValueType s ) const
{
  for( int i = 0; i < 3; i++ )
    ((CORBA::Float *) d)[ i ] = ((CORBA::Float *) s)[ i ];
}

void _Marshaller__a3_float::free( StaticValueType v ) const
{
  delete[] (_MICO_T *) v;
}

::CORBA::Boolean _Marshaller__a3_float::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  if( !dc.arr_begin() )
    return FALSE;
  if (!dc.get_floats (&((_MICO_T *)v)[0], 3))
    return FALSE;
  return dc.arr_end();
}

void _Marshaller__a3_float::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ec.arr_begin();
  ec.put_floats (&((_MICO_T *)v)[0], 3);
  ec.arr_end();
}

::CORBA::TypeCode_ptr _Marshaller__a3_float::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "01000000140000000c000000010000000600000003000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__a3_float::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__a3_float;

class _Marshaller__seq_FoamXServer_StringPair : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_StringPair();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_StringPair::~_Marshaller__seq_FoamXServer_StringPair()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_StringPair::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_StringPair::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_StringPair::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_StringPair::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_StringPair->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_StringPair::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_StringPair->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_StringPair::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000074000000010000000f00000064000000010000001f00"
    "000049444c3a466f616d585365727665722f537472696e67506169723a31"
    "2e3000000b000000537472696e6750616972000002000000050000006e61"
    "6d650000000012000000000000000600000076616c756500000012000000"
    "0000000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_StringPair::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_StringPair;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_StringPair, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_StringPair, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_StringPair, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_HostDescriptor : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_HostDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_HostDescriptor::~_Marshaller__seq_FoamXServer_HostDescriptor()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_HostDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_HostDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_HostDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_HostDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_HostDescriptor->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_HostDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_HostDescriptor->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_HostDescriptor::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000078000000010000000f00000068000000010000002300"
    "000049444c3a466f616d585365727665722f486f73744465736372697074"
    "6f723a312e3000000f000000486f737444657363726970746f7200000200"
    "0000050000006e616d6500000000120000000000000006000000616c6976"
    "650000000800000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_HostDescriptor::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_HostDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_HostDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_HostDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_HostDescriptor, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_ApplicationDescriptor : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_ApplicationDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_ApplicationDescriptor::~_Marshaller__seq_FoamXServer_ApplicationDescriptor()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_ApplicationDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_ApplicationDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_ApplicationDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_ApplicationDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_ApplicationDescriptor->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_ApplicationDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_ApplicationDescriptor->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_ApplicationDescriptor::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "0100000013000000b8000000010000000f000000a8000000010000002a00"
    "000049444c3a466f616d585365727665722f4170706c69636174696f6e44"
    "657363726970746f723a312e30000000160000004170706c69636174696f"
    "6e44657363726970746f7200000004000000050000006e616d6500000000"
    "12000000000000000900000063617465676f727900000000120000000000"
    "000005000000706174680000000012000000000000000c00000073797374"
    "656d436c617373000800000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_ApplicationDescriptor::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_ApplicationDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_ApplicationDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_ApplicationDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_ApplicationDescriptor, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_CaseDescriptor : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_CaseDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_CaseDescriptor::~_Marshaller__seq_FoamXServer_CaseDescriptor()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_CaseDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_CaseDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_CaseDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_CaseDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_CaseDescriptor->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_CaseDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_CaseDescriptor->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_CaseDescriptor::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "0100000013000000e8000000010000000f000000d8000000010000002300"
    "000049444c3a466f616d585365727665722f436173654465736372697074"
    "6f723a312e3000000f0000004361736544657363726970746f7200000800"
    "000008000000726f6f744469720012000000000000000b00000072617752"
    "6f6f744469720000120000000000000009000000636173654e616d650000"
    "000012000000000000000400000061707000120000000000000007000000"
    "6e50726f6373000003000000080000006d616e6167656400080000000700"
    "00006c6f636b6564000008000000060000006572726f7200000008000000"
    "00000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_CaseDescriptor::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_CaseDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_CaseDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_CaseDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_CaseDescriptor, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_JobID : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_JobID();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_JobID::~_Marshaller__seq_FoamXServer_JobID()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_JobID::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_JobID::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_JobID::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_JobID::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_JobID->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_JobID::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_JobID->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_JobID::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "010000001300000070000000010000000f00000060000000010000001a00"
    "000049444c3a466f616d585365727665722f4a6f6249443a312e30000000"
    "060000004a6f6249440000000200000009000000686f73744e616d650000"
    "000012000000000000000a00000070726f63657373494400000003000000"
    "00000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_JobID::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_JobID;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_JobID, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_JobID, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_JobID, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_JobDescriptor : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_JobDescriptor();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_JobDescriptor::~_Marshaller__seq_FoamXServer_JobDescriptor()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_JobDescriptor::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_JobDescriptor::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_JobDescriptor::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_JobDescriptor::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_JobDescriptor->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_JobDescriptor::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_JobDescriptor->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_JobDescriptor::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "0100000013000000ac030000010000000f0000009c030000010000002200"
    "000049444c3a466f616d585365727665722f4a6f6244657363726970746f"
    "723a312e300000000e0000004a6f6244657363726970746f720000001300"
    "0000060000006a6f6249440000000f00000060000000010000001a000000"
    "49444c3a466f616d585365727665722f4a6f6249443a312e300000000600"
    "00004a6f6249440000000200000009000000686f73744e616d6500000000"
    "12000000000000000a00000070726f636573734944000000030000000500"
    "000070706964000000000300000005000000706769640000000003000000"
    "0a00000073746172744461746500000012000000000000000a0000007374"
    "61727454696d65000000120000000000000009000000757365724e616d65"
    "0000000012000000000000000c000000666f616d56657273696f6e001200"
    "00000000000005000000636f646500000000120000000000000008000000"
    "6172674c6973740012000000000000000b00000063757272656e74446972"
    "0000120000000000000008000000726f6f74446972001200000000000000"
    "09000000636173654e616d65000000001200000000000000070000006e50"
    "726f637300000300000007000000736c61766573000015000000b0000000"
    "010000001e00000049444c3a466f616d585365727665722f4a6f6249444c"
    "6973743a312e300000000a0000004a6f6249444c69737400000013000000"
    "70000000010000000f00000060000000010000001a00000049444c3a466f"
    "616d585365727665722f4a6f6249443a312e30000000060000004a6f6249"
    "440000000200000009000000686f73744e616d6500000000120000000000"
    "00000a00000070726f63657373494400000003000000000000000e000000"
    "6e436f756e74656450726f6373000000030000000800000063707554696d"
    "65000700000008000000656e644461746500120000000000000008000000"
    "656e6454696d650012000000000000000700000073746174757300001100"
    "0000c0000000010000001e00000049444c3a466f616d585365727665722f"
    "4a6f625374617475733a312e300000000a0000004a6f6253746174757300"
    "0000070000000e0000004a4f425f554e444546494e45440000000e000000"
    "4a4f425f4c41554e4348494e470000000c0000004a4f425f52554e4e494e"
    "47000d0000004a4f425f53544f5050494e47000000000e0000004a4f425f"
    "53555350454e4445440000000d0000004a4f425f46494e49534845440000"
    "00000c0000004a4f425f41424f525445440000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_JobDescriptor::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_JobDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_JobDescriptor, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_JobDescriptor, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_JobDescriptor, (void *&)_s);
}


class _Marshaller__seq_FoamXServer_DimensionSet : public ::CORBA::StaticTypeInfo {
    typedef SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> _MICO_T;
    static ::CORBA::TypeCode_ptr _tc;
  public:
    ~_Marshaller__seq_FoamXServer_DimensionSet();
    StaticValueType create () const;
    void assign (StaticValueType dst, const StaticValueType src) const;
    void free (StaticValueType) const;
    ::CORBA::Boolean demarshal (::CORBA::DataDecoder&, StaticValueType) const;
    void marshal (::CORBA::DataEncoder &, StaticValueType) const;
    ::CORBA::TypeCode_ptr typecode ();
};


_Marshaller__seq_FoamXServer_DimensionSet::~_Marshaller__seq_FoamXServer_DimensionSet()
{
  if (_tc)
    delete _tc;
}

::CORBA::StaticValueType _Marshaller__seq_FoamXServer_DimensionSet::create() const
{
  return (StaticValueType) new _MICO_T;
}

void _Marshaller__seq_FoamXServer_DimensionSet::assign( StaticValueType d, const StaticValueType s ) const
{
  *(_MICO_T*) d = *(_MICO_T*) s;
}

void _Marshaller__seq_FoamXServer_DimensionSet::free( StaticValueType v ) const
{
  delete (_MICO_T*) v;
}

::CORBA::Boolean _Marshaller__seq_FoamXServer_DimensionSet::demarshal( ::CORBA::DataDecoder &dc, StaticValueType v ) const
{
  ::CORBA::ULong len;
  if( !dc.seq_begin( len ) )
    return FALSE;
  ((_MICO_T *) v)->length( len );
  for( ::CORBA::ULong i = 0; i < len; i++ ) {
    if( !_marshaller_FoamXServer_DimensionSet->demarshal( dc, &(*(_MICO_T*)v)[i] ) )
      return FALSE;
  }
  return dc.seq_end();
}

void _Marshaller__seq_FoamXServer_DimensionSet::marshal( ::CORBA::DataEncoder &ec, StaticValueType v ) const
{
  ::CORBA::ULong len = ((_MICO_T *) v)->length();
  ec.seq_begin( len );
  for( ::CORBA::ULong i = 0; i < len; i++ )
    _marshaller_FoamXServer_DimensionSet->marshal( ec, &(*(_MICO_T*)v)[i] );
  ec.seq_end();
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_DimensionSet::typecode()
{
  if (!_tc)
    _tc = (new ::CORBA::TypeCode (
    "0100000013000000d4000000010000000f000000c4000000010000002100"
    "000049444c3a466f616d585365727665722f44696d656e73696f6e536574"
    "3a312e30000000000d00000044696d656e73696f6e536574000000000700"
    "0000050000006d6173730000000007000000070000006c656e6774680000"
    "070000000500000074696d6500000000070000000c00000074656d706572"
    "61747572650007000000060000006d6f6c65730000000700000008000000"
    "63757272656e740007000000120000006c756d696e6f7573496e74656e73"
    "6974790000000700000000000000"))->mk_constant();
  return _tc;
}

::CORBA::TypeCode_ptr _Marshaller__seq_FoamXServer_DimensionSet::_tc = 0;
::CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_DimensionSet;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_DimensionSet, &_s);
  _a.from_static_any (_sa);
}

void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> *_s )
{
  _a <<= *_s;
  delete _s;
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> &_s )
{
  CORBA::StaticAny _sa (_marshaller__seq_FoamXServer_DimensionSet, &_s);
  return _a.to_static_any (_sa);
}

CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> *&_s )
{
  return _a.to_static_any (_marshaller__seq_FoamXServer_DimensionSet, (void *&)_s);
}


struct __tc_init_FOAMX {
  __tc_init_FOAMX()
  {
    FoamXServer::_tc_FoamXType = 
    "0100000011000000c3010000010000001e00000049444c3a466f616d5853"
    "65727665722f466f616d58547970653a312e300000000a000000466f616d"
    "5854797065000000150000000f000000547970655f556e646566696e6564"
    "00000d000000547970655f426f6f6c65616e000000000b00000054797065"
    "5f4c6162656c00000c000000547970655f5363616c6172000a0000005479"
    "70655f436861720000000a000000547970655f576f72640000000c000000"
    "547970655f537472696e67000d000000547970655f526f6f744469720000"
    "000011000000547970655f526f6f74416e6443617365000000000e000000"
    "547970655f436173654e616d650000000e000000547970655f486f73744e"
    "616d650000000a000000547970655f46696c650000000f00000054797065"
    "5f4469726563746f727900000a000000547970655f54696d650000001200"
    "0000547970655f44696d656e73696f6e5365740000000f00000054797065"
    "5f46697865644c69737400000a000000547970655f4c6973740000001000"
    "0000547970655f44696374696f6e617279000f000000547970655f53656c"
    "656374696f6e00000e000000547970655f436f6d706f756e640000000b00"
    "0000547970655f4669656c6400";
    _marshaller_FoamXServer_FoamXType = new _Marshaller_FoamXServer_FoamXType;
    FoamXServer::_tc_FoamXAny = 
    "010000000f00000024020000010000001d00000049444c3a466f616d5853"
    "65727665722f466f616d58416e793a312e300000000009000000466f616d"
    "58416e79000000000200000005000000747970650000000011000000c301"
    "0000010000001e00000049444c3a466f616d585365727665722f466f616d"
    "58547970653a312e300000000a000000466f616d58547970650000001500"
    "00000f000000547970655f556e646566696e656400000d00000054797065"
    "5f426f6f6c65616e000000000b000000547970655f4c6162656c00000c00"
    "0000547970655f5363616c6172000a000000547970655f43686172000000"
    "0a000000547970655f576f72640000000c000000547970655f537472696e"
    "67000d000000547970655f526f6f74446972000000001100000054797065"
    "5f526f6f74416e6443617365000000000e000000547970655f436173654e"
    "616d650000000e000000547970655f486f73744e616d650000000a000000"
    "547970655f46696c650000000f000000547970655f4469726563746f7279"
    "00000a000000547970655f54696d6500000012000000547970655f44696d"
    "656e73696f6e5365740000000f000000547970655f46697865644c697374"
    "00000a000000547970655f4c69737400000010000000547970655f446963"
    "74696f6e617279000f000000547970655f53656c656374696f6e00000e00"
    "0000547970655f436f6d706f756e640000000b000000547970655f466965"
    "6c6400000600000076616c75650000000b000000";
    _marshaller_FoamXServer_FoamXAny = new _Marshaller_FoamXServer_FoamXAny;
    FoamXServer::_tc_FoamXAnyList = 
    "01000000150000007c020000010000002100000049444c3a466f616d5853"
    "65727665722f466f616d58416e794c6973743a312e30000000000d000000"
    "466f616d58416e794c697374000000001300000034020000010000000f00"
    "000024020000010000001d00000049444c3a466f616d585365727665722f"
    "466f616d58416e793a312e300000000009000000466f616d58416e790000"
    "00000200000005000000747970650000000011000000c301000001000000"
    "1e00000049444c3a466f616d585365727665722f466f616d58547970653a"
    "312e300000000a000000466f616d5854797065000000150000000f000000"
    "547970655f556e646566696e656400000d000000547970655f426f6f6c65"
    "616e000000000b000000547970655f4c6162656c00000c00000054797065"
    "5f5363616c6172000a000000547970655f436861720000000a0000005479"
    "70655f576f72640000000c000000547970655f537472696e67000d000000"
    "547970655f526f6f744469720000000011000000547970655f526f6f7441"
    "6e6443617365000000000e000000547970655f436173654e616d65000000"
    "0e000000547970655f486f73744e616d650000000a000000547970655f46"
    "696c650000000f000000547970655f4469726563746f727900000a000000"
    "547970655f54696d6500000012000000547970655f44696d656e73696f6e"
    "5365740000000f000000547970655f46697865644c69737400000a000000"
    "547970655f4c69737400000010000000547970655f44696374696f6e6172"
    "79000f000000547970655f53656c656374696f6e00000e00000054797065"
    "5f436f6d706f756e640000000b000000547970655f4669656c6400000600"
    "000076616c75650000000b00000000000000";
    FoamXServer::_tc_StringList = 
    "010000001500000050000000010000001f00000049444c3a466f616d5853"
    "65727665722f537472696e674c6973743a312e3000000b00000053747269"
    "6e674c697374000013000000100000000100000012000000000000000000"
    "0000";
    FoamXServer::_tc_TypeDescriptorList = 
    "0100000015000000a0000000010000002700000049444c3a466f616d5853"
    "65727665722f5479706544657363726970746f724c6973743a312e300000"
    "130000005479706544657363726970746f724c6973740000130000005000"
    "0000010000000e00000040000000010000002400000049444c3a466f616d"
    "585365727665722f495479706544657363726970746f723a312e30001000"
    "0000495479706544657363726970746f720000000000";
    FoamXServer::_tc_ErrorCode = 
    "0100000011000000dd000000010000001e00000049444c3a466f616d5853"
    "65727665722f4572726f72436f64653a312e300000000a0000004572726f"
    "72436f64650000000900000005000000535f4f4b0000000007000000455f"
    "4641494c000007000000455f464f414d00000e000000455f494e56414c49"
    "445f4152470000000e000000455f494e56414c49445f5054520000000e00"
    "0000455f494e56414c49445f52454600000016000000455f494e4445585f"
    "4f55545f4f465f424f554e44530000000f000000455f554e4b4e4f574e5f"
    "4e414d4500000d000000455f554e455850454354454400";
    _marshaller_FoamXServer_ErrorCode = new _Marshaller_FoamXServer_ErrorCode;
    FoamXServer::_tc_FoamXError = 
    "010000001600000090010000010000001f00000049444c3a466f616d5853"
    "65727665722f466f616d584572726f723a312e3000000b000000466f616d"
    "584572726f720000050000000a0000006572726f72436f64650000001100"
    "0000dd000000010000001e00000049444c3a466f616d585365727665722f"
    "4572726f72436f64653a312e300000000a0000004572726f72436f646500"
    "00000900000005000000535f4f4b0000000007000000455f4641494c0000"
    "07000000455f464f414d00000e000000455f494e56414c49445f41524700"
    "00000e000000455f494e56414c49445f5054520000000e000000455f494e"
    "56414c49445f52454600000016000000455f494e4445585f4f55545f4f46"
    "5f424f554e44530000000f000000455f554e4b4e4f574e5f4e414d450000"
    "0d000000455f554e4558504543544544000000000d0000006572726f724d"
    "6573736167650000000012000000000000000b0000006d6574686f644e61"
    "6d65000012000000000000000900000066696c654e616d65000000001200"
    "000000000000070000006c696e654e6f000003000000";
    _marshaller_FoamXServer_FoamXError = new _Marshaller_FoamXServer_FoamXError;
    FoamXServer::_tc_FoamXIOError = 
    "0100000016000000e4000000010000002100000049444c3a466f616d5853"
    "65727665722f466f616d58494f4572726f723a312e30000000000d000000"
    "466f616d58494f4572726f7200000000070000000d0000006572726f724d"
    "6573736167650000000012000000000000000b000000696f46696c654e61"
    "6d65000012000000000000000e000000696f53746172744c696e654e6f00"
    "0000030000000c000000696f456e644c696e654e6f00030000000b000000"
    "6d6574686f644e616d65000012000000000000000900000066696c654e61"
    "6d65000000001200000000000000070000006c696e654e6f000003000000"
    ;
    _marshaller_FoamXServer_FoamXIOError = new _Marshaller_FoamXServer_FoamXIOError;
    FoamXServer::_tc_ValidationError = 
    "010000001600000070010000010000002400000049444c3a466f616d5853"
    "65727665722f56616c69646174696f6e4572726f723a312e300010000000"
    "56616c69646174696f6e4572726f7200030000000a0000006572726f7243"
    "6f646500000011000000dd000000010000001e00000049444c3a466f616d"
    "585365727665722f4572726f72436f64653a312e300000000a0000004572"
    "726f72436f64650000000900000005000000535f4f4b0000000007000000"
    "455f4641494c000007000000455f464f414d00000e000000455f494e5641"
    "4c49445f4152470000000e000000455f494e56414c49445f505452000000"
    "0e000000455f494e56414c49445f52454600000016000000455f494e4445"
    "585f4f55545f4f465f424f554e44530000000f000000455f554e4b4e4f57"
    "4e5f4e414d4500000d000000455f554e4558504543544544000000000d00"
    "00006572726f724d65737361676500000000120000000000000009000000"
    "6974656d50617468000000001200000000000000";
    _marshaller_FoamXServer_ValidationError = new _Marshaller_FoamXServer_ValidationError;
    FoamXServer::_tc_ITypeDescriptor = 
    "010000000e00000040000000010000002400000049444c3a466f616d5853"
    "65727665722f495479706544657363726970746f723a312e300010000000"
    "495479706544657363726970746f7200";
    _marshaller_FoamXServer_ITypeDescriptor = new _Marshaller_FoamXServer_ITypeDescriptor;
    FoamXServer::_tc_DictionaryEntryList = 
    "0100000015000000a8000000010000002800000049444c3a466f616d5853"
    "65727665722f44696374696f6e617279456e7472794c6973743a312e3000"
    "1400000044696374696f6e617279456e7472794c69737400130000005800"
    "0000010000000e00000045000000010000002500000049444c3a466f616d"
    "585365727665722f4944696374696f6e617279456e7472793a312e300000"
    "0000110000004944696374696f6e617279456e7472790000000000000000"
    ;
    FoamXServer::_tc_IDictionaryEntry = 
    "010000000e00000045000000010000002500000049444c3a466f616d5853"
    "65727665722f4944696374696f6e617279456e7472793a312e3000000000"
    "110000004944696374696f6e617279456e74727900";
    _marshaller_FoamXServer_IDictionaryEntry = new _Marshaller_FoamXServer_IDictionaryEntry;
    FoamXServer::_tc_JobStatus = 
    "0100000011000000c0000000010000001e00000049444c3a466f616d5853"
    "65727665722f4a6f625374617475733a312e300000000a0000004a6f6253"
    "7461747573000000070000000e0000004a4f425f554e444546494e454400"
    "00000e0000004a4f425f4c41554e4348494e470000000c0000004a4f425f"
    "52554e4e494e47000d0000004a4f425f53544f5050494e47000000000e00"
    "00004a4f425f53555350454e4445440000000d0000004a4f425f46494e49"
    "53484544000000000c0000004a4f425f41424f5254454400";
    _marshaller_FoamXServer_JobStatus = new _Marshaller_FoamXServer_JobStatus;
    FoamXServer::_tc_MessageType = 
    "01000000110000006c000000010000002000000049444c3a466f616d5853"
    "65727665722f4d657373616765547970653a312e30000c0000004d657373"
    "6167655479706500030000000d0000004d5f444941474e4f535449430000"
    "00000a0000004d5f5741524e494e47000000080000004d5f4552524f5200"
    ;
    _marshaller_FoamXServer_MessageType = new _Marshaller_FoamXServer_MessageType;
    FoamXServer::_tc_DoubleList = 
    "01000000150000004c000000010000001f00000049444c3a466f616d5853"
    "65727665722f446f75626c654c6973743a312e3000000b000000446f7562"
    "6c654c6973740000130000000c000000010000000700000000000000";
    FoamXServer::_tc_FloatList = 
    "01000000150000004c000000010000001e00000049444c3a466f616d5853"
    "65727665722f466c6f61744c6973743a312e300000000a000000466c6f61"
    "744c697374000000130000000c000000010000000600000000000000";
    FoamXServer::_tc_LongList = 
    "01000000150000004c000000010000001d00000049444c3a466f616d5853"
    "65727665722f4c6f6e674c6973743a312e3000000000090000004c6f6e67"
    "4c69737400000000130000000c000000010000000300000000000000";
    FoamXServer::_tc_Point3 = 
    "010000001500000044000000010000001b00000049444c3a466f616d5853"
    "65727665722f506f696e74333a312e30000007000000506f696e74330000"
    "140000000c000000010000000600000003000000";
    FoamXServer::_tc_StringPair = 
    "010000000f00000064000000010000001f00000049444c3a466f616d5853"
    "65727665722f537472696e67506169723a312e3000000b00000053747269"
    "6e6750616972000002000000050000006e616d6500000000120000000000"
    "00000600000076616c75650000001200000000000000";
    _marshaller_FoamXServer_StringPair = new _Marshaller_FoamXServer_StringPair;
    FoamXServer::_tc_StringPairList = 
    "0100000015000000bc000000010000002300000049444c3a466f616d5853"
    "65727665722f537472696e67506169724c6973743a312e3000000f000000"
    "537472696e67506169724c69737400001300000074000000010000000f00"
    "000064000000010000001f00000049444c3a466f616d585365727665722f"
    "537472696e67506169723a312e3000000b000000537472696e6750616972"
    "000002000000050000006e616d6500000000120000000000000006000000"
    "76616c7565000000120000000000000000000000";
    FoamXServer::_tc_HostDescriptor = 
    "010000000f00000068000000010000002300000049444c3a466f616d5853"
    "65727665722f486f737444657363726970746f723a312e3000000f000000"
    "486f737444657363726970746f72000002000000050000006e616d650000"
    "0000120000000000000006000000616c69766500000008000000";
    _marshaller_FoamXServer_HostDescriptor = new _Marshaller_FoamXServer_HostDescriptor;
    FoamXServer::_tc_HostDescriptorList = 
    "0100000015000000c8000000010000002700000049444c3a466f616d5853"
    "65727665722f486f737444657363726970746f724c6973743a312e300000"
    "13000000486f737444657363726970746f724c6973740000130000007800"
    "0000010000000f00000068000000010000002300000049444c3a466f616d"
    "585365727665722f486f737444657363726970746f723a312e3000000f00"
    "0000486f737444657363726970746f72000002000000050000006e616d65"
    "00000000120000000000000006000000616c697665000000080000000000"
    "0000";
    FoamXServer::_tc_ApplicationDescriptor = 
    "010000000f000000a8000000010000002a00000049444c3a466f616d5853"
    "65727665722f4170706c69636174696f6e44657363726970746f723a312e"
    "30000000160000004170706c69636174696f6e44657363726970746f7200"
    "000004000000050000006e616d6500000000120000000000000009000000"
    "63617465676f727900000000120000000000000005000000706174680000"
    "000012000000000000000c00000073797374656d436c6173730008000000"
    ;
    _marshaller_FoamXServer_ApplicationDescriptor = new _Marshaller_FoamXServer_ApplicationDescriptor;
    FoamXServer::_tc_ApplicationDescriptorList = 
    "010000001500000018010000010000002e00000049444c3a466f616d5853"
    "65727665722f4170706c69636174696f6e44657363726970746f724c6973"
    "743a312e300000001a0000004170706c69636174696f6e44657363726970"
    "746f724c69737400000013000000b8000000010000000f000000a8000000"
    "010000002a00000049444c3a466f616d585365727665722f4170706c6963"
    "6174696f6e44657363726970746f723a312e30000000160000004170706c"
    "69636174696f6e44657363726970746f7200000004000000050000006e61"
    "6d650000000012000000000000000900000063617465676f727900000000"
    "120000000000000005000000706174680000000012000000000000000c00"
    "000073797374656d436c617373000800000000000000";
    FoamXServer::_tc_CaseDescriptor = 
    "010000000f000000d8000000010000002300000049444c3a466f616d5853"
    "65727665722f4361736544657363726970746f723a312e3000000f000000"
    "4361736544657363726970746f7200000800000008000000726f6f744469"
    "720012000000000000000b000000726177526f6f74446972000012000000"
    "0000000009000000636173654e616d650000000012000000000000000400"
    "0000617070001200000000000000070000006e50726f6373000003000000"
    "080000006d616e616765640008000000070000006c6f636b656400000800"
    "0000060000006572726f7200000008000000";
    _marshaller_FoamXServer_CaseDescriptor = new _Marshaller_FoamXServer_CaseDescriptor;
    FoamXServer::_tc_CaseDescriptorList = 
    "010000001500000038010000010000002700000049444c3a466f616d5853"
    "65727665722f4361736544657363726970746f724c6973743a312e300000"
    "130000004361736544657363726970746f724c697374000013000000e800"
    "0000010000000f000000d8000000010000002300000049444c3a466f616d"
    "585365727665722f4361736544657363726970746f723a312e3000000f00"
    "00004361736544657363726970746f7200000800000008000000726f6f74"
    "4469720012000000000000000b000000726177526f6f7444697200001200"
    "00000000000009000000636173654e616d65000000001200000000000000"
    "04000000617070001200000000000000070000006e50726f637300000300"
    "0000080000006d616e616765640008000000070000006c6f636b65640000"
    "08000000060000006572726f720000000800000000000000";
    FoamXServer::_tc_JobID = 
    "010000000f00000060000000010000001a00000049444c3a466f616d5853"
    "65727665722f4a6f6249443a312e30000000060000004a6f624944000000"
    "0200000009000000686f73744e616d650000000012000000000000000a00"
    "000070726f63657373494400000003000000";
    _marshaller_FoamXServer_JobID = new _Marshaller_FoamXServer_JobID;
    FoamXServer::_tc_JobIDList = 
    "0100000015000000b0000000010000001e00000049444c3a466f616d5853"
    "65727665722f4a6f6249444c6973743a312e300000000a0000004a6f6249"
    "444c6973740000001300000070000000010000000f000000600000000100"
    "00001a00000049444c3a466f616d585365727665722f4a6f6249443a312e"
    "30000000060000004a6f6249440000000200000009000000686f73744e61"
    "6d650000000012000000000000000a00000070726f636573734944000000"
    "0300000000000000";
    FoamXServer::_tc_JobDescriptor = 
    "010000000f0000009c030000010000002200000049444c3a466f616d5853"
    "65727665722f4a6f6244657363726970746f723a312e300000000e000000"
    "4a6f6244657363726970746f7200000013000000060000006a6f62494400"
    "00000f00000060000000010000001a00000049444c3a466f616d58536572"
    "7665722f4a6f6249443a312e30000000060000004a6f6249440000000200"
    "000009000000686f73744e616d650000000012000000000000000a000000"
    "70726f636573734944000000030000000500000070706964000000000300"
    "0000050000007067696400000000030000000a0000007374617274446174"
    "6500000012000000000000000a000000737461727454696d650000001200"
    "00000000000009000000757365724e616d65000000001200000000000000"
    "0c000000666f616d56657273696f6e00120000000000000005000000636f"
    "6465000000001200000000000000080000006172674c6973740012000000"
    "000000000b00000063757272656e74446972000012000000000000000800"
    "0000726f6f7444697200120000000000000009000000636173654e616d65"
    "000000001200000000000000070000006e50726f63730000030000000700"
    "0000736c61766573000015000000b0000000010000001e00000049444c3a"
    "466f616d585365727665722f4a6f6249444c6973743a312e300000000a00"
    "00004a6f6249444c6973740000001300000070000000010000000f000000"
    "60000000010000001a00000049444c3a466f616d585365727665722f4a6f"
    "6249443a312e30000000060000004a6f6249440000000200000009000000"
    "686f73744e616d650000000012000000000000000a00000070726f636573"
    "73494400000003000000000000000e0000006e436f756e74656450726f63"
    "73000000030000000800000063707554696d65000700000008000000656e"
    "644461746500120000000000000008000000656e6454696d650012000000"
    "0000000007000000737461747573000011000000c0000000010000001e00"
    "000049444c3a466f616d585365727665722f4a6f625374617475733a312e"
    "300000000a0000004a6f62537461747573000000070000000e0000004a4f"
    "425f554e444546494e45440000000e0000004a4f425f4c41554e4348494e"
    "470000000c0000004a4f425f52554e4e494e47000d0000004a4f425f5354"
    "4f5050494e47000000000e0000004a4f425f53555350454e444544000000"
    "0d0000004a4f425f46494e4953484544000000000c0000004a4f425f4142"
    "4f5254454400";
    _marshaller_FoamXServer_JobDescriptor = new _Marshaller_FoamXServer_JobDescriptor;
    FoamXServer::_tc_JobDescriptorList = 
    "0100000015000000fc030000010000002600000049444c3a466f616d5853"
    "65727665722f4a6f6244657363726970746f724c6973743a312e30000000"
    "120000004a6f6244657363726970746f724c69737400000013000000ac03"
    "0000010000000f0000009c030000010000002200000049444c3a466f616d"
    "585365727665722f4a6f6244657363726970746f723a312e300000000e00"
    "00004a6f6244657363726970746f7200000013000000060000006a6f6249"
    "440000000f00000060000000010000001a00000049444c3a466f616d5853"
    "65727665722f4a6f6249443a312e30000000060000004a6f624944000000"
    "0200000009000000686f73744e616d650000000012000000000000000a00"
    "000070726f63657373494400000003000000050000007070696400000000"
    "03000000050000007067696400000000030000000a000000737461727444"
    "61746500000012000000000000000a000000737461727454696d65000000"
    "120000000000000009000000757365724e616d6500000000120000000000"
    "00000c000000666f616d56657273696f6e00120000000000000005000000"
    "636f6465000000001200000000000000080000006172674c697374001200"
    "0000000000000b00000063757272656e7444697200001200000000000000"
    "08000000726f6f7444697200120000000000000009000000636173654e61"
    "6d65000000001200000000000000070000006e50726f6373000003000000"
    "07000000736c61766573000015000000b0000000010000001e0000004944"
    "4c3a466f616d585365727665722f4a6f6249444c6973743a312e30000000"
    "0a0000004a6f6249444c6973740000001300000070000000010000000f00"
    "000060000000010000001a00000049444c3a466f616d585365727665722f"
    "4a6f6249443a312e30000000060000004a6f624944000000020000000900"
    "0000686f73744e616d650000000012000000000000000a00000070726f63"
    "657373494400000003000000000000000e0000006e436f756e7465645072"
    "6f6373000000030000000800000063707554696d65000700000008000000"
    "656e644461746500120000000000000008000000656e6454696d65001200"
    "00000000000007000000737461747573000011000000c000000001000000"
    "1e00000049444c3a466f616d585365727665722f4a6f625374617475733a"
    "312e300000000a0000004a6f62537461747573000000070000000e000000"
    "4a4f425f554e444546494e45440000000e0000004a4f425f4c41554e4348"
    "494e470000000c0000004a4f425f52554e4e494e47000d0000004a4f425f"
    "53544f5050494e47000000000e0000004a4f425f53555350454e44454400"
    "00000d0000004a4f425f46494e4953484544000000000c0000004a4f425f"
    "41424f525445440000000000";
    FoamXServer::_tc_DimensionSet = 
    "010000000f000000c4000000010000002100000049444c3a466f616d5853"
    "65727665722f44696d656e73696f6e5365743a312e30000000000d000000"
    "44696d656e73696f6e5365740000000007000000050000006d6173730000"
    "000007000000070000006c656e6774680000070000000500000074696d65"
    "00000000070000000c00000074656d706572617475726500070000000600"
    "00006d6f6c6573000000070000000800000063757272656e740007000000"
    "120000006c756d696e6f7573496e74656e7369747900000007000000";
    _marshaller_FoamXServer_DimensionSet = new _Marshaller_FoamXServer_DimensionSet;
    FoamXServer::_tc_DimensionSetList = 
    "010000001500000024010000010000002500000049444c3a466f616d5853"
    "65727665722f44696d656e73696f6e5365744c6973743a312e3000000000"
    "1100000044696d656e73696f6e5365744c6973740000000013000000d400"
    "0000010000000f000000c4000000010000002100000049444c3a466f616d"
    "585365727665722f44696d656e73696f6e5365743a312e30000000000d00"
    "000044696d656e73696f6e5365740000000007000000050000006d617373"
    "0000000007000000070000006c656e677468000007000000050000007469"
    "6d6500000000070000000c00000074656d70657261747572650007000000"
    "060000006d6f6c6573000000070000000800000063757272656e74000700"
    "0000120000006c756d696e6f7573496e74656e7369747900000007000000"
    "00000000";
    FoamXServer::_tc_FoamXSYSError = 
    "0100000016000000b0010000010000002200000049444c3a466f616d5853"
    "65727665722f466f616d585359534572726f723a312e300000000e000000"
    "466f616d585359534572726f72000000060000000a0000006572726f7243"
    "6f646500000011000000dd000000010000001e00000049444c3a466f616d"
    "585365727665722f4572726f72436f64653a312e300000000a0000004572"
    "726f72436f64650000000900000005000000535f4f4b0000000007000000"
    "455f4641494c000007000000455f464f414d00000e000000455f494e5641"
    "4c49445f4152470000000e000000455f494e56414c49445f505452000000"
    "0e000000455f494e56414c49445f52454600000016000000455f494e4445"
    "585f4f55545f4f465f424f554e44530000000f000000455f554e4b4e4f57"
    "4e5f4e414d4500000d000000455f554e4558504543544544000000000d00"
    "00006572726f724d65737361676500000000120000000000000009000000"
    "686f73744e616d650000000012000000000000000b0000006d6574686f64"
    "4e616d65000012000000000000000900000066696c654e616d6500000000"
    "1200000000000000070000006c696e654e6f000003000000";
    _marshaller_FoamXServer_FoamXSYSError = new _Marshaller_FoamXServer_FoamXSYSError;
    FoamXServer::CaseServer::_tc_ICaseServer = 
    "010000000e00000044000000010000002b00000049444c3a466f616d5853"
    "65727665722f436173655365727665722f49436173655365727665723a31"
    "2e3000000c000000494361736553657276657200";
    _marshaller_FoamXServer_CaseServer_ICaseServer = new _Marshaller_FoamXServer_CaseServer_ICaseServer;
    FoamXServer::CaseServer::_tc_IFoamProperties = 
    "010000000e0000004c000000010000002f00000049444c3a466f616d5853"
    "65727665722f436173655365727665722f49466f616d50726f7065727469"
    "65733a312e3000001000000049466f616d50726f7065727469657300";
    _marshaller_FoamXServer_CaseServer_IFoamProperties = new _Marshaller_FoamXServer_CaseServer_IFoamProperties;
    FoamXServer::CaseServer::_tc_IApplication = 
    "010000000e00000045000000010000002c00000049444c3a466f616d5853"
    "65727665722f436173655365727665722f494170706c69636174696f6e3a"
    "312e30000d000000494170706c69636174696f6e00";
    _marshaller_FoamXServer_CaseServer_IApplication = new _Marshaller_FoamXServer_CaseServer_IApplication;
    FoamXServer::CaseServer::_tc_IGeometricFieldDescriptor = 
    "010000000e00000062000000010000003900000049444c3a466f616d5853"
    "65727665722f436173655365727665722f4947656f6d6574726963466965"
    "6c6444657363726970746f723a312e30000000001a0000004947656f6d65"
    "747269634669656c6444657363726970746f7200";
    _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor = new _Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor;
    FoamXServer::CaseServer::_tc_IPatchDescriptor = 
    "010000000e0000004d000000010000003000000049444c3a466f616d5853"
    "65727665722f436173655365727665722f49506174636844657363726970"
    "746f723a312e30001100000049506174636844657363726970746f7200";
    _marshaller_FoamXServer_CaseServer_IPatchDescriptor = new _Marshaller_FoamXServer_CaseServer_IPatchDescriptor;
    FoamXServer::CaseServer::_tc_IGeometryDescriptor = 
    "010000000e00000054000000010000003300000049444c3a466f616d5853"
    "65727665722f436173655365727665722f4947656f6d6574727944657363"
    "726970746f723a312e300000140000004947656f6d657472794465736372"
    "6970746f7200";
    _marshaller_FoamXServer_CaseServer_IGeometryDescriptor = new _Marshaller_FoamXServer_CaseServer_IGeometryDescriptor;
    FoamXServer::CaseServer::_tc_IPatchPhysicalTypeDescriptor = 
    "010000000e00000065000000010000003c00000049444c3a466f616d5853"
    "65727665722f436173655365727665722f49506174636850687973696361"
    "6c5479706544657363726970746f723a312e30001d000000495061746368"
    "506879736963616c5479706544657363726970746f7200";
    _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor = new _Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor;
    FoamXServer::CaseServer::_tc_IGeometricField = 
    "010000000e0000004c000000010000002f00000049444c3a466f616d5853"
    "65727665722f436173655365727665722f4947656f6d6574726963466965"
    "6c643a312e300000100000004947656f6d65747269634669656c6400";
    _marshaller_FoamXServer_CaseServer_IGeometricField = new _Marshaller_FoamXServer_CaseServer_IGeometricField;
    FoamXServer::CasePostServer::_tc_ICasePostServer = 
    "010000000e00000050000000010000003300000049444c3a466f616d5853"
    "65727665722f43617365506f73745365727665722f4943617365506f7374"
    "5365727665723a312e300000100000004943617365506f73745365727665"
    "7200";
    _marshaller_FoamXServer_CasePostServer_ICasePostServer = new _Marshaller_FoamXServer_CasePostServer_ICasePostServer;
    FoamXServer::CaseBrowser::_tc_ICaseBrowser = 
    "010000000e00000049000000010000002d00000049444c3a466f616d5853"
    "65727665722f4361736542726f777365722f494361736542726f77736572"
    "3a312e30000000000d000000494361736542726f7773657200";
    _marshaller_FoamXServer_CaseBrowser_ICaseBrowser = new _Marshaller_FoamXServer_CaseBrowser_ICaseBrowser;
    FoamXServer::HostBrowser::_tc_IHostBrowser = 
    "010000000e00000049000000010000002d00000049444c3a466f616d5853"
    "65727665722f486f737442726f777365722f49486f737442726f77736572"
    "3a312e30000000000d00000049486f737442726f7773657200";
    _marshaller_FoamXServer_HostBrowser_IHostBrowser = new _Marshaller_FoamXServer_HostBrowser_IHostBrowser;
    _marshaller__seq_FoamXServer_FoamXAny = new _Marshaller__seq_FoamXServer_FoamXAny;
    _marshaller__seq_FoamXServer_ITypeDescriptor = new _Marshaller__seq_FoamXServer_ITypeDescriptor;
    _marshaller__seq_FoamXServer_IDictionaryEntry = new _Marshaller__seq_FoamXServer_IDictionaryEntry;
    _marshaller__a3_float = new _Marshaller__a3_float;
    _marshaller__seq_FoamXServer_StringPair = new _Marshaller__seq_FoamXServer_StringPair;
    _marshaller__seq_FoamXServer_HostDescriptor = new _Marshaller__seq_FoamXServer_HostDescriptor;
    _marshaller__seq_FoamXServer_ApplicationDescriptor = new _Marshaller__seq_FoamXServer_ApplicationDescriptor;
    _marshaller__seq_FoamXServer_CaseDescriptor = new _Marshaller__seq_FoamXServer_CaseDescriptor;
    _marshaller__seq_FoamXServer_JobID = new _Marshaller__seq_FoamXServer_JobID;
    _marshaller__seq_FoamXServer_JobDescriptor = new _Marshaller__seq_FoamXServer_JobDescriptor;
    _marshaller__seq_FoamXServer_DimensionSet = new _Marshaller__seq_FoamXServer_DimensionSet;
  }

  ~__tc_init_FOAMX()
  {
    delete static_cast<_Marshaller_FoamXServer_FoamXType*>(_marshaller_FoamXServer_FoamXType);
    delete static_cast<_Marshaller_FoamXServer_FoamXAny*>(_marshaller_FoamXServer_FoamXAny);
    delete static_cast<_Marshaller_FoamXServer_ErrorCode*>(_marshaller_FoamXServer_ErrorCode);
    delete static_cast<_Marshaller_FoamXServer_FoamXError*>(_marshaller_FoamXServer_FoamXError);
    delete static_cast<_Marshaller_FoamXServer_FoamXIOError*>(_marshaller_FoamXServer_FoamXIOError);
    delete static_cast<_Marshaller_FoamXServer_ValidationError*>(_marshaller_FoamXServer_ValidationError);
    delete static_cast<_Marshaller_FoamXServer_ITypeDescriptor*>(_marshaller_FoamXServer_ITypeDescriptor);
    delete static_cast<_Marshaller_FoamXServer_IDictionaryEntry*>(_marshaller_FoamXServer_IDictionaryEntry);
    delete static_cast<_Marshaller_FoamXServer_JobStatus*>(_marshaller_FoamXServer_JobStatus);
    delete static_cast<_Marshaller_FoamXServer_MessageType*>(_marshaller_FoamXServer_MessageType);
    delete static_cast<_Marshaller_FoamXServer_StringPair*>(_marshaller_FoamXServer_StringPair);
    delete static_cast<_Marshaller_FoamXServer_HostDescriptor*>(_marshaller_FoamXServer_HostDescriptor);
    delete static_cast<_Marshaller_FoamXServer_ApplicationDescriptor*>(_marshaller_FoamXServer_ApplicationDescriptor);
    delete static_cast<_Marshaller_FoamXServer_CaseDescriptor*>(_marshaller_FoamXServer_CaseDescriptor);
    delete static_cast<_Marshaller_FoamXServer_JobID*>(_marshaller_FoamXServer_JobID);
    delete static_cast<_Marshaller_FoamXServer_JobDescriptor*>(_marshaller_FoamXServer_JobDescriptor);
    delete static_cast<_Marshaller_FoamXServer_DimensionSet*>(_marshaller_FoamXServer_DimensionSet);
    delete static_cast<_Marshaller_FoamXServer_FoamXSYSError*>(_marshaller_FoamXServer_FoamXSYSError);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_ICaseServer*>(_marshaller_FoamXServer_CaseServer_ICaseServer);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IFoamProperties*>(_marshaller_FoamXServer_CaseServer_IFoamProperties);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IApplication*>(_marshaller_FoamXServer_CaseServer_IApplication);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor*>(_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IPatchDescriptor*>(_marshaller_FoamXServer_CaseServer_IPatchDescriptor);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IGeometryDescriptor*>(_marshaller_FoamXServer_CaseServer_IGeometryDescriptor);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor*>(_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor);
    delete static_cast<_Marshaller_FoamXServer_CaseServer_IGeometricField*>(_marshaller_FoamXServer_CaseServer_IGeometricField);
    delete static_cast<_Marshaller_FoamXServer_CasePostServer_ICasePostServer*>(_marshaller_FoamXServer_CasePostServer_ICasePostServer);
    delete static_cast<_Marshaller_FoamXServer_CaseBrowser_ICaseBrowser*>(_marshaller_FoamXServer_CaseBrowser_ICaseBrowser);
    delete static_cast<_Marshaller_FoamXServer_HostBrowser_IHostBrowser*>(_marshaller_FoamXServer_HostBrowser_IHostBrowser);
    delete static_cast<_Marshaller__seq_FoamXServer_FoamXAny*>(_marshaller__seq_FoamXServer_FoamXAny);
    delete static_cast<_Marshaller__seq_FoamXServer_ITypeDescriptor*>(_marshaller__seq_FoamXServer_ITypeDescriptor);
    delete static_cast<_Marshaller__seq_FoamXServer_IDictionaryEntry*>(_marshaller__seq_FoamXServer_IDictionaryEntry);
    delete static_cast<_Marshaller__a3_float*>(_marshaller__a3_float);
    delete static_cast<_Marshaller__seq_FoamXServer_StringPair*>(_marshaller__seq_FoamXServer_StringPair);
    delete static_cast<_Marshaller__seq_FoamXServer_HostDescriptor*>(_marshaller__seq_FoamXServer_HostDescriptor);
    delete static_cast<_Marshaller__seq_FoamXServer_ApplicationDescriptor*>(_marshaller__seq_FoamXServer_ApplicationDescriptor);
    delete static_cast<_Marshaller__seq_FoamXServer_CaseDescriptor*>(_marshaller__seq_FoamXServer_CaseDescriptor);
    delete static_cast<_Marshaller__seq_FoamXServer_JobID*>(_marshaller__seq_FoamXServer_JobID);
    delete static_cast<_Marshaller__seq_FoamXServer_JobDescriptor*>(_marshaller__seq_FoamXServer_JobDescriptor);
    delete static_cast<_Marshaller__seq_FoamXServer_DimensionSet*>(_marshaller__seq_FoamXServer_DimensionSet);
  }
};

static __tc_init_FOAMX __init_FOAMX;

//--------------------------------------------------------
//  Implementation of skeletons
//--------------------------------------------------------

// PortableServer Skeleton Class for interface FoamXServer::ITypeDescriptor
POA_FoamXServer::ITypeDescriptor::~ITypeDescriptor()
{
}

::FoamXServer::ITypeDescriptor_ptr
POA_FoamXServer::ITypeDescriptor::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::ITypeDescriptor::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::ITypeDescriptor::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/ITypeDescriptor:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::ITypeDescriptor::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/ITypeDescriptor:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::ITypeDescriptor::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/ITypeDescriptor:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::ITypeDescriptor::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::ITypeDescriptor_stub_clp (poa, obj);
}

bool
POA_FoamXServer::ITypeDescriptor::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 71)) {
    case 1:
      if( strcmp( __req->op_name(), "_set_optional" ) == 0 ) {
        CORBA::Boolean _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        optional( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "_get_type" ) == 0 ) {
        ::FoamXServer::FoamXType _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXType, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = type();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_maxValue" ) == 0 ) {
        ::FoamXServer::FoamXAny _par__value;
        CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        maxValue( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "_set_path" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        path( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "addSubType" ) == 0 ) {
        ::FoamXServer::FoamXType _par_type;
        CORBA::StaticAny _sa_type( _marshaller_FoamXServer_FoamXType, &_par_type );
        ::FoamXServer::ITypeDescriptor_ptr _par_subEntry;
        CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_ITypeDescriptor, &_par_subEntry );

        __req->add_in_arg( &_sa_type );
        __req->add_out_arg( &_sa_subEntry );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addSubType( _par_type, _par_subEntry );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_subEntry );
        return true;
      }
      break;
    case 10:
      if( strcmp( __req->op_name(), "_get_isCompoundType" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = isCompoundType();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_lookupDict" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = lookupDict();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 16:
      if( strcmp( __req->op_name(), "_set_lookupDict" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        lookupDict( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 17:
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_elementType" ) == 0 ) {
        ::FoamXServer::ITypeDescriptor_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = elementType();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 19:
      if( strcmp( __req->op_name(), "_get_editable" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = editable();
        __req->write_results();
        return true;
      }
      break;
    case 23:
      if( strcmp( __req->op_name(), "_get_helpURL" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = helpURL();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 25:
      if( strcmp( __req->op_name(), "getDefaultValue" ) == 0 ) {
        ::FoamXServer::IDictionaryEntry_ptr _par_defaultValue;
        CORBA::StaticAny _sa_defaultValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_defaultValue );

        __req->add_out_arg( &_sa_defaultValue );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getDefaultValue( _par_defaultValue );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_defaultValue );
        return true;
      }
      break;
    case 26:
      if( strcmp( __req->op_name(), "_get_minValue" ) == 0 ) {
        ::FoamXServer::FoamXAny* _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = minValue();
        __res.value( _marshaller_FoamXServer_FoamXAny, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 27:
      if( strcmp( __req->op_name(), "_get_path" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = path();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 28:
      if( strcmp( __req->op_name(), "_get_visible" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = visible();
        __req->write_results();
        return true;
      }
      break;
    case 31:
      if( strcmp( __req->op_name(), "_set_valueList" ) == 0 ) {
        ::FoamXServer::FoamXAnyList _par__value;
        CORBA::StaticAny _sa__value( _marshaller__seq_FoamXServer_FoamXAny, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        valueList( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 32:
      if( strcmp( __req->op_name(), "_get_dictionaryPath" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = dictionaryPath();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 33:
      if( strcmp( __req->op_name(), "_get_category" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = category();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_set_comment" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        comment( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 34:
      if( strcmp( __req->op_name(), "_get_subTypes" ) == 0 ) {
        ::FoamXServer::TypeDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ITypeDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = subTypes();
        __res.value( _marshaller__seq_FoamXServer_ITypeDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      if( strcmp( __req->op_name(), "_set_helpURL" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        helpURL( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 36:
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 38:
      if( strcmp( __req->op_name(), "_get_displayName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = displayName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_optional" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = optional();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_maxValue" ) == 0 ) {
        ::FoamXServer::FoamXAny* _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = maxValue();
        __res.value( _marshaller_FoamXServer_FoamXAny, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 39:
      if( strcmp( __req->op_name(), "_set_elementLabels" ) == 0 ) {
        ::FoamXServer::StringList _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stcseq_string, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        elementLabels( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 40:
      if( strcmp( __req->op_name(), "_set_type" ) == 0 ) {
        ::FoamXServer::FoamXType _par__value;
        CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXType, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        type( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 41:
      if( strcmp( __req->op_name(), "_set_category" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        category( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 43:
      if( strcmp( __req->op_name(), "_get_isPrimitiveType" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = isPrimitiveType();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_comment" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = comment();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 44:
      if( strcmp( __req->op_name(), "_get_numElements" ) == 0 ) {
        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = numElements();
        __req->write_results();
        return true;
      }
      break;
    case 46:
      if( strcmp( __req->op_name(), "_set_displayName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        displayName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 50:
      if( strcmp( __req->op_name(), "_get_elementLabels" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = elementLabels();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 51:
      if( strcmp( __req->op_name(), "_get_iconURL" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = iconURL();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 52:
      if( strcmp( __req->op_name(), "hasDefaultValue" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = hasDefaultValue();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_numElements" ) == 0 ) {
        CORBA::Long _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_long, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        numElements( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 60:
      if( strcmp( __req->op_name(), "_set_minValue" ) == 0 ) {
        ::FoamXServer::FoamXAny _par__value;
        CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        minValue( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 61:
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_editable" ) == 0 ) {
        CORBA::Boolean _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        editable( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 64:
      if( strcmp( __req->op_name(), "_set_visible" ) == 0 ) {
        CORBA::Boolean _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        visible( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 65:
      if( strcmp( __req->op_name(), "_set_iconURL" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        iconURL( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 66:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_valueList" ) == 0 ) {
        ::FoamXServer::FoamXAnyList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_FoamXAny );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = valueList();
        __res.value( _marshaller__seq_FoamXServer_FoamXAny, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      if( strcmp( __req->op_name(), "_set_dictionaryPath" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        dictionaryPath( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 69:
      if( strcmp( __req->op_name(), "removeSubType" ) == 0 ) {
        ::FoamXServer::ITypeDescriptor_var _par_subEntry;
        CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_ITypeDescriptor, &_par_subEntry._for_demarshal() );

        __req->add_in_arg( &_sa_subEntry );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          removeSubType( _par_subEntry.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::ITypeDescriptor::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::IDictionaryEntry
POA_FoamXServer::IDictionaryEntry::~IDictionaryEntry()
{
}

::FoamXServer::IDictionaryEntry_ptr
POA_FoamXServer::IDictionaryEntry::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::IDictionaryEntry::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::IDictionaryEntry::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/IDictionaryEntry:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::IDictionaryEntry::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/IDictionaryEntry:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::IDictionaryEntry::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/IDictionaryEntry:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::IDictionaryEntry::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::IDictionaryEntry_stub_clp (poa, obj);
}

bool
POA_FoamXServer::IDictionaryEntry::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 23)) {
    case 0:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 1:
      if( strcmp( __req->op_name(), "modified" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = modified();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 3:
      if( strcmp( __req->op_name(), "save" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          save();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_typeDescriptor" ) == 0 ) {
        ::FoamXServer::ITypeDescriptor_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = typeDescriptor();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "_get_selection" ) == 0 ) {
        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = selection();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_selection" ) == 0 ) {
        CORBA::Long _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_long, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        selection( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 5:
      if( strcmp( __req->op_name(), "nSubElements" ) == 0 ) {
        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = nSubElements();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "_get_value" ) == 0 ) {
        ::FoamXServer::FoamXAny* _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_FoamXAny );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = value();
        __res.value( _marshaller_FoamXServer_FoamXAny, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "_get_subElements" ) == 0 ) {
        ::FoamXServer::DictionaryEntryList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_IDictionaryEntry );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = subElements();
        __res.value( _marshaller__seq_FoamXServer_IDictionaryEntry, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 18:
      if( strcmp( __req->op_name(), "_set_value" ) == 0 ) {
        ::FoamXServer::FoamXAny _par__value;
        CORBA::StaticAny _sa__value( _marshaller_FoamXServer_FoamXAny, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        value( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 19:
      if( strcmp( __req->op_name(), "packedList" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = packedList();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 21:
      if( strcmp( __req->op_name(), "addElement" ) == 0 ) {
        ::FoamXServer::IDictionaryEntry_ptr _par_subEntry;
        CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_IDictionaryEntry, &_par_subEntry );

        __req->add_out_arg( &_sa_subEntry );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addElement( _par_subEntry );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_subEntry );
        return true;
      }
      if( strcmp( __req->op_name(), "removeElement" ) == 0 ) {
        ::FoamXServer::IDictionaryEntry_var _par_subEntry;
        CORBA::StaticAny _sa_subEntry( _marshaller_FoamXServer_IDictionaryEntry, &_par_subEntry._for_demarshal() );

        __req->add_in_arg( &_sa_subEntry );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          removeElement( _par_subEntry.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 22:
      if( strcmp( __req->op_name(), "setValue" ) == 0 ) {
        ::FoamXServer::FoamXAny _par_value;
        CORBA::StaticAny _sa_value( _marshaller_FoamXServer_FoamXAny, &_par_value );

        __req->add_in_arg( &_sa_value );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          setValue( _par_value );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::IDictionaryEntry::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::ICaseServer
POA_FoamXServer::CaseServer::ICaseServer::~ICaseServer()
{
}

::FoamXServer::CaseServer::ICaseServer_ptr
POA_FoamXServer::CaseServer::ICaseServer::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::ICaseServer::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::ICaseServer::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/ICaseServer:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::ICaseServer::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/ICaseServer:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::ICaseServer::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/ICaseServer:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::ICaseServer::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::ICaseServer_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::ICaseServer::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 43)) {
    case 1:
      if( strcmp( __req->op_name(), "readMesh" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          readMesh();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 3:
      if( strcmp( __req->op_name(), "getTime" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getTime();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 5:
      if( strcmp( __req->op_name(), "_get_caseName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = caseName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "_get_foamProperties" ) == 0 ) {
        ::FoamXServer::CaseServer::IFoamProperties_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IFoamProperties, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = foamProperties();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 7:
      if( strcmp( __req->op_name(), "getFieldValues" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::CaseServer::IGeometricField_ptr _par_fieldValues;
        CORBA::StaticAny _sa_fieldValues( _marshaller_FoamXServer_CaseServer_IGeometricField, &_par_fieldValues );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_out_arg( &_sa_fieldValues );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getFieldValues( _par_fieldName.inout(), _par_fieldValues );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_fieldValues );
        return true;
      }
      if( strcmp( __req->op_name(), "save" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          save();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_meshDefined" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = meshDefined();
        __req->write_results();
        return true;
      }
      break;
    case 10:
      if( strcmp( __req->op_name(), "setPatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        CORBA::String_var _par_patchPhysicalType;
        CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType._for_demarshal() );

        __req->add_in_arg( &_sa_patchName );
        __req->add_in_arg( &_sa_patchPhysicalType );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          setPatchPhysicalType( _par_patchName.inout(), _par_patchPhysicalType.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_patchNames" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchNames();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "_get_caseRoot" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = caseRoot();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 12:
      if( strcmp( __req->op_name(), "getDictionary" ) == 0 ) {
        CORBA::String_var _par_dictionaryName;
        CORBA::StaticAny _sa_dictionaryName( CORBA::_stc_string, &_par_dictionaryName._for_demarshal() );
        CORBA::Boolean _par_forceRead;
        CORBA::StaticAny _sa_forceRead( CORBA::_stc_boolean, &_par_forceRead );
        ::FoamXServer::IDictionaryEntry_ptr _par_dictRoot;
        CORBA::StaticAny _sa_dictRoot( _marshaller_FoamXServer_IDictionaryEntry, &_par_dictRoot );

        __req->add_in_arg( &_sa_dictionaryName );
        __req->add_in_arg( &_sa_forceRead );
        __req->add_out_arg( &_sa_dictRoot );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getDictionary( _par_dictionaryName.inout(), _par_forceRead, _par_dictRoot );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_dictRoot );
        return true;
      }
      if( strcmp( __req->op_name(), "killCase" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          killCase();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 16:
      if( strcmp( __req->op_name(), "_get_availableTimeSteps" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = availableTimeSteps();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 19:
      if( strcmp( __req->op_name(), "setTime" ) == 0 ) {
        CORBA::String_var _par_timeName;
        CORBA::StaticAny _sa_timeName( CORBA::_stc_string, &_par_timeName._for_demarshal() );
        CORBA::Long _par_timeIndex;
        CORBA::StaticAny _sa_timeIndex( CORBA::_stc_long, &_par_timeIndex );

        __req->add_in_arg( &_sa_timeName );
        __req->add_in_arg( &_sa_timeIndex );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          setTime( _par_timeName.inout(), _par_timeIndex );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 24:
      if( strcmp( __req->op_name(), "importMesh" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );

        __req->add_in_arg( &_sa_hostName );
        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          importMesh( _par_hostName.inout(), _par_rootDir.inout(), _par_caseName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 26:
      if( strcmp( __req->op_name(), "readFile" ) == 0 ) {
        CORBA::String_var _par_name;
        CORBA::StaticAny _sa_name( CORBA::_stc_string, &_par_name._for_demarshal() );
        char* _par_contents;
        CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents );

        __req->add_in_arg( &_sa_name );
        __req->add_out_arg( &_sa_contents );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          readFile( _par_name.inout(), _par_contents );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_contents );
        return true;
      }
      break;
    case 27:
      if( strcmp( __req->op_name(), "runCase" ) == 0 ) {
        CORBA::String_var _par_arguments;
        CORBA::StaticAny _sa_arguments( CORBA::_stc_string, &_par_arguments._for_demarshal() );

        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->add_in_arg( &_sa_arguments );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = runCase( _par_arguments.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 28:
      if( strcmp( __req->op_name(), "writeFile" ) == 0 ) {
        CORBA::String_var _par_name;
        CORBA::StaticAny _sa_name( CORBA::_stc_string, &_par_name._for_demarshal() );
        CORBA::String_var _par_contents;
        CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents._for_demarshal() );

        __req->add_in_arg( &_sa_name );
        __req->add_in_arg( &_sa_contents );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          writeFile( _par_name.inout(), _par_contents.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 29:
      if( strcmp( __req->op_name(), "getPatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        char* _par_patchPhysicalType;
        CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType );

        __req->add_in_arg( &_sa_patchName );
        __req->add_out_arg( &_sa_patchPhysicalType );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchPhysicalType( _par_patchName.inout(), _par_patchPhysicalType );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_patchPhysicalType );
        return true;
      }
      if( strcmp( __req->op_name(), "_set_managed" ) == 0 ) {
        CORBA::Boolean _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_boolean, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        managed( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 30:
      if( strcmp( __req->op_name(), "fileModificationDate" ) == 0 ) {
        CORBA::String_var _par_fileName;
        CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName._for_demarshal() );

        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->add_in_arg( &_sa_fileName );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = fileModificationDate( _par_fileName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 31:
      if( strcmp( __req->op_name(), "modified" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = modified();
        __req->write_results();
        return true;
      }
      break;
    case 33:
      if( strcmp( __req->op_name(), "_get_application" ) == 0 ) {
        ::FoamXServer::CaseServer::IApplication_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IApplication, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = application();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 34:
      if( strcmp( __req->op_name(), "_get_managed" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = managed();
        __req->write_results();
        return true;
      }
      break;
    case 35:
      if( strcmp( __req->op_name(), "deletePatch" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );

        __req->add_in_arg( &_sa_patchName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deletePatch( _par_patchName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "close" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        close();
        __req->write_results();
        return true;
      }
      break;
    case 42:
      if( strcmp( __req->op_name(), "addPatch" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        CORBA::String_var _par_patchPhysicalType;
        CORBA::StaticAny _sa_patchPhysicalType( CORBA::_stc_string, &_par_patchPhysicalType._for_demarshal() );

        __req->add_in_arg( &_sa_patchName );
        __req->add_in_arg( &_sa_patchPhysicalType );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addPatch( _par_patchName.inout(), _par_patchPhysicalType.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "deleteAllPatches" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteAllPatches();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::ICaseServer::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IFoamProperties
POA_FoamXServer::CaseServer::IFoamProperties::~IFoamProperties()
{
}

::FoamXServer::CaseServer::IFoamProperties_ptr
POA_FoamXServer::CaseServer::IFoamProperties::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IFoamProperties::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IFoamProperties::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IFoamProperties:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IFoamProperties::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IFoamProperties:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IFoamProperties::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IFoamProperties:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IFoamProperties::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IFoamProperties_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IFoamProperties::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 41)) {
    case 6:
      if( strcmp( __req->op_name(), "getFoamControlDict" ) == 0 ) {
        ::FoamXServer::IDictionaryEntry_ptr _par_controlDict;
        CORBA::StaticAny _sa_controlDict( _marshaller_FoamXServer_IDictionaryEntry, &_par_controlDict );

        __req->add_out_arg( &_sa_controlDict );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getFoamControlDict( _par_controlDict );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_controlDict );
        return true;
      }
      break;
    case 7:
      if( strcmp( __req->op_name(), "getFoamType" ) == 0 ) {
        CORBA::String_var _par_foamTypeName;
        CORBA::StaticAny _sa_foamTypeName( CORBA::_stc_string, &_par_foamTypeName._for_demarshal() );
        ::FoamXServer::ITypeDescriptor_ptr _par_typeDesc;
        CORBA::StaticAny _sa_typeDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_typeDesc );

        __req->add_in_arg( &_sa_foamTypeName );
        __req->add_out_arg( &_sa_typeDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getFoamType( _par_foamTypeName.inout(), _par_typeDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_typeDesc );
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "saveUserProperties" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          saveUserProperties();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_utilities" ) == 0 ) {
        ::FoamXServer::ApplicationDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ApplicationDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = utilities();
        __res.value( _marshaller__seq_FoamXServer_ApplicationDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "findPatchFieldType" ) == 0 ) {
        CORBA::String_var _par_patchFieldTypeName;
        CORBA::StaticAny _sa_patchFieldTypeName( CORBA::_stc_string, &_par_patchFieldTypeName._for_demarshal() );
        ::FoamXServer::ITypeDescriptor_ptr _par_patchFieldDesc;
        CORBA::StaticAny _sa_patchFieldDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_patchFieldDesc );

        __req->add_in_arg( &_sa_patchFieldTypeName );
        __req->add_out_arg( &_sa_patchFieldDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          findPatchFieldType( _par_patchFieldTypeName.inout(), _par_patchFieldDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchFieldDesc );
        return true;
      }
      break;
    case 12:
      if( strcmp( __req->op_name(), "getUtilityControlDict" ) == 0 ) {
        CORBA::String_var _par_utilityName;
        CORBA::StaticAny _sa_utilityName( CORBA::_stc_string, &_par_utilityName._for_demarshal() );
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        ::FoamXServer::IDictionaryEntry_ptr _par_controlDict;
        CORBA::StaticAny _sa_controlDict( _marshaller_FoamXServer_IDictionaryEntry, &_par_controlDict );

        __req->add_in_arg( &_sa_utilityName );
        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_out_arg( &_sa_controlDict );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getUtilityControlDict( _par_utilityName.inout(), _par_rootDir.inout(), _par_caseName.inout(), _par_controlDict );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_controlDict );
        return true;
      }
      break;
    case 13:
      if( strcmp( __req->op_name(), "findPatchType" ) == 0 ) {
        CORBA::String_var _par_patchTypeName;
        CORBA::StaticAny _sa_patchTypeName( CORBA::_stc_string, &_par_patchTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IPatchDescriptor_ptr _par_patchDesc;
        CORBA::StaticAny _sa_patchDesc( _marshaller_FoamXServer_CaseServer_IPatchDescriptor, &_par_patchDesc );

        __req->add_in_arg( &_sa_patchTypeName );
        __req->add_out_arg( &_sa_patchDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          findPatchType( _par_patchTypeName.inout(), _par_patchDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchDesc );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_patchTypes" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchTypes();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 17:
      if( strcmp( __req->op_name(), "cloneApplication" ) == 0 ) {
        CORBA::String_var _par_appNameSrc;
        CORBA::StaticAny _sa_appNameSrc( CORBA::_stc_string, &_par_appNameSrc._for_demarshal() );
        CORBA::String_var _par_appNameDest;
        CORBA::StaticAny _sa_appNameDest( CORBA::_stc_string, &_par_appNameDest._for_demarshal() );
        CORBA::String_var _par_appDestPath;
        CORBA::StaticAny _sa_appDestPath( CORBA::_stc_string, &_par_appDestPath._for_demarshal() );
        ::FoamXServer::CaseServer::IApplication_ptr _par_app;
        CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app );

        __req->add_in_arg( &_sa_appNameSrc );
        __req->add_in_arg( &_sa_appNameDest );
        __req->add_in_arg( &_sa_appDestPath );
        __req->add_out_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          cloneApplication( _par_appNameSrc.inout(), _par_appNameDest.inout(), _par_appDestPath.inout(), _par_app );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_app );
        return true;
      }
      break;
    case 18:
      if( strcmp( __req->op_name(), "addRootDirectory" ) == 0 ) {
        CORBA::String_var _par_rawRootDir;
        CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir._for_demarshal() );

        __req->add_in_arg( &_sa_rawRootDir );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addRootDirectory( _par_rawRootDir.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_applicationes" ) == 0 ) {
        ::FoamXServer::ApplicationDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_ApplicationDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = applicationes();
        __res.value( _marshaller__seq_FoamXServer_ApplicationDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 19:
      if( strcmp( __req->op_name(), "getGeometryType" ) == 0 ) {
        CORBA::String_var _par_geometryTypeName;
        CORBA::StaticAny _sa_geometryTypeName( CORBA::_stc_string, &_par_geometryTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IGeometryDescriptor_ptr _par_geometryDesc;
        CORBA::StaticAny _sa_geometryDesc( _marshaller_FoamXServer_CaseServer_IGeometryDescriptor, &_par_geometryDesc );

        __req->add_in_arg( &_sa_geometryTypeName );
        __req->add_out_arg( &_sa_geometryDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getGeometryType( _par_geometryTypeName.inout(), _par_geometryDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_geometryDesc );
        return true;
      }
      break;
    case 21:
      if( strcmp( __req->op_name(), "_get_availableModules" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = availableModules();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      if( strcmp( __req->op_name(), "_get_geometryTypes" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = geometryTypes();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 22:
      if( strcmp( __req->op_name(), "_get_rawRootDirectories" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = rawRootDirectories();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 26:
      if( strcmp( __req->op_name(), "getPatchFieldType" ) == 0 ) {
        CORBA::String_var _par_patchFieldTypeName;
        CORBA::StaticAny _sa_patchFieldTypeName( CORBA::_stc_string, &_par_patchFieldTypeName._for_demarshal() );
        ::FoamXServer::ITypeDescriptor_ptr _par_patchFieldDesc;
        CORBA::StaticAny _sa_patchFieldDesc( _marshaller_FoamXServer_ITypeDescriptor, &_par_patchFieldDesc );

        __req->add_in_arg( &_sa_patchFieldTypeName );
        __req->add_out_arg( &_sa_patchFieldDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchFieldType( _par_patchFieldTypeName.inout(), _par_patchFieldDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchFieldDesc );
        return true;
      }
      break;
    case 27:
      if( strcmp( __req->op_name(), "saveSystemProperties" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          saveSystemProperties();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 29:
      if( strcmp( __req->op_name(), "addApplication" ) == 0 ) {
        CORBA::String_var _par_appName;
        CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName._for_demarshal() );
        ::FoamXServer::CaseServer::IApplication_ptr _par_app;
        CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app );

        __req->add_in_arg( &_sa_appName );
        __req->add_out_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addApplication( _par_appName.inout(), _par_app );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_app );
        return true;
      }
      if( strcmp( __req->op_name(), "deleteApplication" ) == 0 ) {
        CORBA::String_var _par_appName;
        CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName._for_demarshal() );

        __req->add_in_arg( &_sa_appName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteApplication( _par_appName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 34:
      if( strcmp( __req->op_name(), "getApplication" ) == 0 ) {
        CORBA::String_var _par_appName;
        CORBA::StaticAny _sa_appName( CORBA::_stc_string, &_par_appName._for_demarshal() );
        ::FoamXServer::CaseServer::IApplication_ptr _par_app;
        CORBA::StaticAny _sa_app( _marshaller_FoamXServer_CaseServer_IApplication, &_par_app );

        __req->add_in_arg( &_sa_appName );
        __req->add_out_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getApplication( _par_appName.inout(), _par_app );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_app );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_rootDirectories" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = rootDirectories();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 36:
      if( strcmp( __req->op_name(), "deleteRootDirectory" ) == 0 ) {
        CORBA::String_var _par_rawRootDir;
        CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir._for_demarshal() );

        __req->add_in_arg( &_sa_rawRootDir );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteRootDirectory( _par_rawRootDir.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 37:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 39:
      if( strcmp( __req->op_name(), "_get_foamTypes" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = foamTypes();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 40:
      if( strcmp( __req->op_name(), "getPatchType" ) == 0 ) {
        CORBA::String_var _par_patchTypeName;
        CORBA::StaticAny _sa_patchTypeName( CORBA::_stc_string, &_par_patchTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IPatchDescriptor_ptr _par_patchDesc;
        CORBA::StaticAny _sa_patchDesc( _marshaller_FoamXServer_CaseServer_IPatchDescriptor, &_par_patchDesc );

        __req->add_in_arg( &_sa_patchTypeName );
        __req->add_out_arg( &_sa_patchDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchType( _par_patchTypeName.inout(), _par_patchDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchDesc );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_patchFieldTypes" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchFieldTypes();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IFoamProperties::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IApplication
POA_FoamXServer::CaseServer::IApplication::~IApplication()
{
}

::FoamXServer::CaseServer::IApplication_ptr
POA_FoamXServer::CaseServer::IApplication::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IApplication::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IApplication::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IApplication:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IApplication::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IApplication:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IApplication::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IApplication:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IApplication::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IApplication_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IApplication::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 37)) {
    case 0:
      if( strcmp( __req->op_name(), "deleteField" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );

        __req->add_in_arg( &_sa_fieldName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteField( _par_fieldName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_fields" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = fields();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 1:
      if( strcmp( __req->op_name(), "getPatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchPhysicalTypeName;
        CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _par_patchPhysicalTypeDescriptor;
        CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor );

        __req->add_in_arg( &_sa_patchPhysicalTypeName );
        __req->add_out_arg( &_sa_patchPhysicalTypeDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchPhysicalType( _par_patchPhysicalTypeName.inout(), _par_patchPhysicalTypeDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchPhysicalTypeDescriptor );
        return true;
      }
      break;
    case 2:
      if( strcmp( __req->op_name(), "addPatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchPhysicalTypeName;
        CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _par_patchPhysicalTypeDescriptor;
        CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor );

        __req->add_in_arg( &_sa_patchPhysicalTypeName );
        __req->add_out_arg( &_sa_patchPhysicalTypeDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addPatchPhysicalType( _par_patchPhysicalTypeName.inout(), _par_patchPhysicalTypeDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchPhysicalTypeDescriptor );
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "_get_dictionaries" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = dictionaries();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 5:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "_set_modules" ) == 0 ) {
        ::FoamXServer::StringList _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stcseq_string, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        modules( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "deleteDictionary" ) == 0 ) {
        CORBA::String_var _par_dictName;
        CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName._for_demarshal() );

        __req->add_in_arg( &_sa_dictName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteDictionary( _par_dictName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "getField" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _par_fieldDescriptor;
        CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_out_arg( &_sa_fieldDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getField( _par_fieldName.inout(), _par_fieldDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_fieldDescriptor );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_patchPhysicalTypes" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchPhysicalTypes();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 12:
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "_set_category" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        category( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 16:
      if( strcmp( __req->op_name(), "addField" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _par_fieldDescriptor;
        CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_out_arg( &_sa_fieldDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addField( _par_fieldName.inout(), _par_fieldDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_fieldDescriptor );
        return true;
      }
      break;
    case 17:
      if( strcmp( __req->op_name(), "_get_modules" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = modules();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 20:
      if( strcmp( __req->op_name(), "getDictionary" ) == 0 ) {
        CORBA::String_var _par_dictName;
        CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName._for_demarshal() );
        ::FoamXServer::ITypeDescriptor_ptr _par_dictTypeDescriptor;
        CORBA::StaticAny _sa_dictTypeDescriptor( _marshaller_FoamXServer_ITypeDescriptor, &_par_dictTypeDescriptor );

        __req->add_in_arg( &_sa_dictName );
        __req->add_out_arg( &_sa_dictTypeDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getDictionary( _par_dictName.inout(), _par_dictTypeDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_dictTypeDescriptor );
        return true;
      }
      break;
    case 25:
      if( strcmp( __req->op_name(), "_get_category" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = category();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 26:
      if( strcmp( __req->op_name(), "save" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          save();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 28:
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 32:
      if( strcmp( __req->op_name(), "findField" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _par_fieldDescriptor;
        CORBA::StaticAny _sa_fieldDescriptor( _marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor, &_par_fieldDescriptor );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_out_arg( &_sa_fieldDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          findField( _par_fieldName.inout(), _par_fieldDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_fieldDescriptor );
        return true;
      }
      break;
    case 34:
      if( strcmp( __req->op_name(), "findPatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchPhysicalTypeName;
        CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName._for_demarshal() );
        ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _par_patchPhysicalTypeDescriptor;
        CORBA::StaticAny _sa_patchPhysicalTypeDescriptor( _marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor, &_par_patchPhysicalTypeDescriptor );

        __req->add_in_arg( &_sa_patchPhysicalTypeName );
        __req->add_out_arg( &_sa_patchPhysicalTypeDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          findPatchPhysicalType( _par_patchPhysicalTypeName.inout(), _par_patchPhysicalTypeDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchPhysicalTypeDescriptor );
        return true;
      }
      break;
    case 35:
      if( strcmp( __req->op_name(), "addDictionary" ) == 0 ) {
        CORBA::String_var _par_dictName;
        CORBA::StaticAny _sa_dictName( CORBA::_stc_string, &_par_dictName._for_demarshal() );
        ::FoamXServer::ITypeDescriptor_ptr _par_dictTypeDescriptor;
        CORBA::StaticAny _sa_dictTypeDescriptor( _marshaller_FoamXServer_ITypeDescriptor, &_par_dictTypeDescriptor );

        __req->add_in_arg( &_sa_dictName );
        __req->add_out_arg( &_sa_dictTypeDescriptor );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addDictionary( _par_dictName.inout(), _par_dictTypeDescriptor );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_dictTypeDescriptor );
        return true;
      }
      break;
    case 36:
      if( strcmp( __req->op_name(), "deletePatchPhysicalType" ) == 0 ) {
        CORBA::String_var _par_patchPhysicalTypeName;
        CORBA::StaticAny _sa_patchPhysicalTypeName( CORBA::_stc_string, &_par_patchPhysicalTypeName._for_demarshal() );

        __req->add_in_arg( &_sa_patchPhysicalTypeName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deletePatchPhysicalType( _par_patchPhysicalTypeName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_systemClass" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = systemClass();
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IApplication::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IGeometricFieldDescriptor
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::~IGeometricFieldDescriptor()
{
}

::FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IGeometricFieldDescriptor::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IGeometricFieldDescriptor:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IGeometricFieldDescriptor_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 19)) {
    case 1:
      if( strcmp( __req->op_name(), "_get_geometryDescriptor" ) == 0 ) {
        ::FoamXServer::CaseServer::IGeometryDescriptor_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IGeometryDescriptor, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = geometryDescriptor();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 2:
      if( strcmp( __req->op_name(), "_get_dimensions" ) == 0 ) {
        ::FoamXServer::DimensionSet _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_DimensionSet, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = dimensions();
        __req->write_results();
        return true;
      }
      break;
    case 3:
      if( strcmp( __req->op_name(), "_get_typeDescriptor" ) == 0 ) {
        ::FoamXServer::ITypeDescriptor_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = typeDescriptor();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 5:
      if( strcmp( __req->op_name(), "_set_fieldTypeName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        fieldTypeName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 8:
      if( strcmp( __req->op_name(), "_get_fieldTypeDescriptor" ) == 0 ) {
        ::FoamXServer::ITypeDescriptor_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_ITypeDescriptor, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = fieldTypeDescriptor();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 10:
      if( strcmp( __req->op_name(), "_set_geometryTypeName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        geometryTypeName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 13:
      if( strcmp( __req->op_name(), "_get_fieldTypeName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = fieldTypeName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_geometryTypeName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = geometryTypeName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 17:
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 18:
      if( strcmp( __req->op_name(), "_set_dimensions" ) == 0 ) {
        ::FoamXServer::DimensionSet _par__value;
        CORBA::StaticAny _sa__value( _marshaller_FoamXServer_DimensionSet, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        dimensions( _par__value );
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IGeometricFieldDescriptor::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IPatchDescriptor
POA_FoamXServer::CaseServer::IPatchDescriptor::~IPatchDescriptor()
{
}

::FoamXServer::CaseServer::IPatchDescriptor_ptr
POA_FoamXServer::CaseServer::IPatchDescriptor::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IPatchDescriptor::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IPatchDescriptor::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IPatchDescriptor::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IPatchDescriptor::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IPatchDescriptor:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IPatchDescriptor::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IPatchDescriptor_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IPatchDescriptor::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 11)) {
    case 4:
      if( strcmp( __req->op_name(), "_get_displayName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = displayName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 7:
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_set_displayName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        displayName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IPatchDescriptor::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IGeometryDescriptor
POA_FoamXServer::CaseServer::IGeometryDescriptor::~IGeometryDescriptor()
{
}

::FoamXServer::CaseServer::IGeometryDescriptor_ptr
POA_FoamXServer::CaseServer::IGeometryDescriptor::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IGeometryDescriptor::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IGeometryDescriptor::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IGeometryDescriptor::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IGeometryDescriptor::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IGeometryDescriptor:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IGeometryDescriptor::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IGeometryDescriptor_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IGeometryDescriptor::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 11)) {
    case 4:
      if( strcmp( __req->op_name(), "_get_displayName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = displayName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 7:
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_set_displayName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        displayName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IGeometryDescriptor::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::~IPatchPhysicalTypeDescriptor()
{
}

::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IPatchPhysicalTypeDescriptor:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 19)) {
    case 0:
      if( strcmp( __req->op_name(), "_get_displayName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = displayName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 2:
      if( strcmp( __req->op_name(), "_get_patchType" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchType();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_set_patchType" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        patchType( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "_set_patchFieldTypes" ) == 0 ) {
        ::FoamXServer::StringPairList _par__value;
        CORBA::StaticAny _sa__value( _marshaller__seq_FoamXServer_StringPair, &_par__value );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        patchFieldTypes( _par__value );
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "_set_displayName" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        displayName( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 10:
      if( strcmp( __req->op_name(), "_get_patchFieldTypes" ) == 0 ) {
        ::FoamXServer::StringPairList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_StringPair );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = patchFieldTypes();
        __res.value( _marshaller__seq_FoamXServer_StringPair, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      if( strcmp( __req->op_name(), "_get_description" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = description();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 15:
      if( strcmp( __req->op_name(), "_set_parentType" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        parentType( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 17:
      if( strcmp( __req->op_name(), "_set_name" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        name( _par__value.inout() );
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_set_description" ) == 0 ) {
        CORBA::String_var _par__value;
        CORBA::StaticAny _sa__value( CORBA::_stc_string, &_par__value._for_demarshal() );

        __req->add_in_arg( &_sa__value );

        if( !__req->read_args() )
          return true;

        description( _par__value.inout() );
        __req->write_results();
        return true;
      }
      break;
    case 18:
      if( strcmp( __req->op_name(), "_get_parentType" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = parentType();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseServer::IGeometricField
POA_FoamXServer::CaseServer::IGeometricField::~IGeometricField()
{
}

::FoamXServer::CaseServer::IGeometricField_ptr
POA_FoamXServer::CaseServer::IGeometricField::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseServer::IGeometricField::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseServer::IGeometricField::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseServer/IGeometricField:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseServer::IGeometricField::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseServer/IGeometricField:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseServer::IGeometricField::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseServer/IGeometricField:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseServer::IGeometricField::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseServer::IGeometricField_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseServer::IGeometricField::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 7)) {
    case 2:
      if( strcmp( __req->op_name(), "getPatchFieldParameters" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        ::FoamXServer::IDictionaryEntry_ptr _par_patchFieldValue;
        CORBA::StaticAny _sa_patchFieldValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_patchFieldValue );

        __req->add_in_arg( &_sa_patchName );
        __req->add_out_arg( &_sa_patchFieldValue );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchFieldParameters( _par_patchName.inout(), _par_patchFieldValue );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_patchFieldValue );
        return true;
      }
      break;
    case 5:
      if( strcmp( __req->op_name(), "getInternalFieldValue" ) == 0 ) {
        ::FoamXServer::IDictionaryEntry_ptr _par_internalFieldValue;
        CORBA::StaticAny _sa_internalFieldValue( _marshaller_FoamXServer_IDictionaryEntry, &_par_internalFieldValue );

        __req->add_out_arg( &_sa_internalFieldValue );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getInternalFieldValue( _par_internalFieldValue );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_internalFieldValue );
        return true;
      }
      break;
    case 6:
      if( strcmp( __req->op_name(), "modified" ) == 0 ) {
        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = modified();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_name" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = name();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseServer::IGeometricField::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CasePostServer::ICasePostServer
POA_FoamXServer::CasePostServer::ICasePostServer::~ICasePostServer()
{
}

::FoamXServer::CasePostServer::ICasePostServer_ptr
POA_FoamXServer::CasePostServer::ICasePostServer::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CasePostServer::ICasePostServer::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CasePostServer::ICasePostServer::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CasePostServer/ICasePostServer:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CasePostServer::ICasePostServer::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CasePostServer/ICasePostServer:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CasePostServer::ICasePostServer::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CasePostServer/ICasePostServer:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CasePostServer::ICasePostServer::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CasePostServer::ICasePostServer_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CasePostServer::ICasePostServer::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 41)) {
    case 0:
      if( strcmp( __req->op_name(), "getCutMesh" ) == 0 ) {
        ::FoamXServer::Point3_slice* _par_basePoint = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
        ::FoamXServer::Point3_slice* _par_normal = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_edges;
        CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );

        __req->add_in_arg( &_sa_basePoint );
        __req->add_in_arg( &_sa_normal );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_edges );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getCutMesh( _par_basePoint, _par_normal, _par_points, _par_edges );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_edges.value( CORBA::_stcseq_long, _par_edges );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        ::FoamXServer::Point3_free( _par_basePoint );
        ::FoamXServer::Point3_free( _par_normal );
        delete _par_points;
        delete _par_edges;
        return true;
      }
      break;
    case 1:
      if( strcmp( __req->op_name(), "_get_nProcs" ) == 0 ) {
        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = nProcs();
        __req->write_results();
        return true;
      }
      break;
    case 2:
      if( strcmp( __req->op_name(), "getCutMeshOutlineSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getCutMeshOutlineSlave();
        __req->write_results();
        return true;
      }
      break;
    case 3:
      if( strcmp( __req->op_name(), "getPatchNamesSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getPatchNamesSlave();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "getTriPatchVec" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_triFaces;
        CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
        ::FoamXServer::FloatList* _par_values;
        CORBA::StaticAny _sa_values( CORBA::_stcseq_float );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_in_arg( &_sa_patchName );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_triFaces );
        __req->add_out_arg( &_sa_values );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getTriPatchVec( _par_fieldName.inout(), _par_patchName.inout(), _par_points, _par_triFaces, _par_values );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_triFaces.value( CORBA::_stcseq_long, _par_triFaces );
          _sa_values.value( CORBA::_stcseq_float, _par_values );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        delete _par_points;
        delete _par_triFaces;
        delete _par_values;
        return true;
      }
      if( strcmp( __req->op_name(), "cutPlaneVecSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        cutPlaneVecSlave();
        __req->write_results();
        return true;
      }
      break;
    case 8:
      if( strcmp( __req->op_name(), "getCutMeshSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getCutMeshSlave();
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "getPatchMeshSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getPatchMeshSlave();
        __req->write_results();
        return true;
      }
      break;
    case 10:
      if( strcmp( __req->op_name(), "getMeshBb" ) == 0 ) {
        ::FoamXServer::Point3_slice* _par_min = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_min( _marshaller__a3_float, _par_min );
        ::FoamXServer::Point3_slice* _par_max = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_max( _marshaller__a3_float, _par_max );

        __req->add_out_arg( &_sa_min );
        __req->add_out_arg( &_sa_max );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getMeshBb( _par_min, _par_max );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        ::FoamXServer::Point3_free( _par_min );
        ::FoamXServer::Point3_free( _par_max );
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "getMeshBbSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getMeshBbSlave();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "getPatchMesh" ) == 0 ) {
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        CORBA::Double _par_creaseAngle;
        CORBA::StaticAny _sa_creaseAngle( CORBA::_stc_double, &_par_creaseAngle );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_edges;
        CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );

        __req->add_in_arg( &_sa_patchName );
        __req->add_in_arg( &_sa_creaseAngle );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_edges );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getPatchMesh( _par_patchName.inout(), _par_creaseAngle, _par_points, _par_edges );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_edges.value( CORBA::_stcseq_long, _par_edges );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        delete _par_points;
        delete _par_edges;
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "close" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        close();
        __req->write_results();
        return true;
      }
      break;
    case 15:
      if( strcmp( __req->op_name(), "setTime" ) == 0 ) {
        CORBA::String_var _par_timeName;
        CORBA::StaticAny _sa_timeName( CORBA::_stc_string, &_par_timeName._for_demarshal() );
        CORBA::Long _par_timeIndex;
        CORBA::StaticAny _sa_timeIndex( CORBA::_stc_long, &_par_timeIndex );

        __req->add_in_arg( &_sa_timeName );
        __req->add_in_arg( &_sa_timeIndex );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          setTime( _par_timeName.inout(), _par_timeIndex );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "cutPlaneVec" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::Point3_slice* _par_basePoint = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
        ::FoamXServer::Point3_slice* _par_normal = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_triFaces;
        CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
        ::FoamXServer::FloatList* _par_values;
        CORBA::StaticAny _sa_values( CORBA::_stcseq_float );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_in_arg( &_sa_basePoint );
        __req->add_in_arg( &_sa_normal );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_triFaces );
        __req->add_out_arg( &_sa_values );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          cutPlaneVec( _par_fieldName.inout(), _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_triFaces.value( CORBA::_stcseq_long, _par_triFaces );
          _sa_values.value( CORBA::_stcseq_float, _par_values );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        ::FoamXServer::Point3_free( _par_basePoint );
        ::FoamXServer::Point3_free( _par_normal );
        delete _par_points;
        delete _par_triFaces;
        delete _par_values;
        return true;
      }
      if( strcmp( __req->op_name(), "_get_caseRoot" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = caseRoot();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    case 19:
      if( strcmp( __req->op_name(), "getTriPatch" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        CORBA::String_var _par_patchName;
        CORBA::StaticAny _sa_patchName( CORBA::_stc_string, &_par_patchName._for_demarshal() );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_triFaces;
        CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
        ::FoamXServer::FloatList* _par_values;
        CORBA::StaticAny _sa_values( CORBA::_stcseq_float );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_in_arg( &_sa_patchName );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_triFaces );
        __req->add_out_arg( &_sa_values );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getTriPatch( _par_fieldName.inout(), _par_patchName.inout(), _par_points, _par_triFaces, _par_values );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_triFaces.value( CORBA::_stcseq_long, _par_triFaces );
          _sa_values.value( CORBA::_stcseq_float, _par_values );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        delete _par_points;
        delete _par_triFaces;
        delete _par_values;
        return true;
      }
      break;
    case 22:
      if( strcmp( __req->op_name(), "getTriPatchVecSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getTriPatchVecSlave();
        __req->write_results();
        return true;
      }
      break;
    case 24:
      if( strcmp( __req->op_name(), "cutPlaneSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        cutPlaneSlave();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_availableTimeSteps" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = availableTimeSteps();
        __res.value( CORBA::_stcseq_string, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 30:
      if( strcmp( __req->op_name(), "setTimeSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        setTimeSlave();
        __req->write_results();
        return true;
      }
      break;
    case 33:
      if( strcmp( __req->op_name(), "getTriPatchSlave" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        getTriPatchSlave();
        __req->write_results();
        return true;
      }
      break;
    case 35:
      if( strcmp( __req->op_name(), "getPatchNames" ) == 0 ) {
        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getPatchNames();
          __res.value( CORBA::_stcseq_string, _res );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        delete _res;
        return true;
      }
      if( strcmp( __req->op_name(), "getFieldNames" ) == 0 ) {
        CORBA::String_var _par_type;
        CORBA::StaticAny _sa_type( CORBA::_stc_string, &_par_type._for_demarshal() );

        ::FoamXServer::StringList* _res;
        CORBA::StaticAny __res( CORBA::_stcseq_string );
        __req->add_in_arg( &_sa_type );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getFieldNames( _par_type.inout() );
          __res.value( CORBA::_stcseq_string, _res );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        delete _res;
        return true;
      }
      if( strcmp( __req->op_name(), "cutPlane" ) == 0 ) {
        CORBA::String_var _par_fieldName;
        CORBA::StaticAny _sa_fieldName( CORBA::_stc_string, &_par_fieldName._for_demarshal() );
        ::FoamXServer::Point3_slice* _par_basePoint = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
        ::FoamXServer::Point3_slice* _par_normal = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_triFaces;
        CORBA::StaticAny _sa_triFaces( CORBA::_stcseq_long );
        ::FoamXServer::FloatList* _par_values;
        CORBA::StaticAny _sa_values( CORBA::_stcseq_float );

        __req->add_in_arg( &_sa_fieldName );
        __req->add_in_arg( &_sa_basePoint );
        __req->add_in_arg( &_sa_normal );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_triFaces );
        __req->add_out_arg( &_sa_values );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          cutPlane( _par_fieldName.inout(), _par_basePoint, _par_normal, _par_points, _par_triFaces, _par_values );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_triFaces.value( CORBA::_stcseq_long, _par_triFaces );
          _sa_values.value( CORBA::_stcseq_float, _par_values );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        ::FoamXServer::Point3_free( _par_basePoint );
        ::FoamXServer::Point3_free( _par_normal );
        delete _par_points;
        delete _par_triFaces;
        delete _par_values;
        return true;
      }
      break;
    case 38:
      if( strcmp( __req->op_name(), "getCutMeshOutline" ) == 0 ) {
        ::FoamXServer::Point3_slice* _par_basePoint = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_basePoint( _marshaller__a3_float, _par_basePoint );
        ::FoamXServer::Point3_slice* _par_normal = ::FoamXServer::Point3_alloc();
        CORBA::StaticAny _sa_normal( _marshaller__a3_float, _par_normal );
        ::FoamXServer::FloatList* _par_points;
        CORBA::StaticAny _sa_points( CORBA::_stcseq_float );
        ::FoamXServer::LongList* _par_edges;
        CORBA::StaticAny _sa_edges( CORBA::_stcseq_long );

        __req->add_in_arg( &_sa_basePoint );
        __req->add_in_arg( &_sa_normal );
        __req->add_out_arg( &_sa_points );
        __req->add_out_arg( &_sa_edges );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getCutMeshOutline( _par_basePoint, _par_normal, _par_points, _par_edges );
          _sa_points.value( CORBA::_stcseq_float, _par_points );
          _sa_edges.value( CORBA::_stcseq_long, _par_edges );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        ::FoamXServer::Point3_free( _par_basePoint );
        ::FoamXServer::Point3_free( _par_normal );
        delete _par_points;
        delete _par_edges;
        return true;
      }
      break;
    case 40:
      if( strcmp( __req->op_name(), "_get_caseName" ) == 0 ) {
        char* _res;
        CORBA::StaticAny __res( CORBA::_stc_string, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = caseName();
        __req->write_results();
        CORBA::string_free( _res );
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CasePostServer::ICasePostServer::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::CaseBrowser::ICaseBrowser
POA_FoamXServer::CaseBrowser::ICaseBrowser::~ICaseBrowser()
{
}

::FoamXServer::CaseBrowser::ICaseBrowser_ptr
POA_FoamXServer::CaseBrowser::ICaseBrowser::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::CaseBrowser::ICaseBrowser::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::CaseBrowser::ICaseBrowser::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::CaseBrowser::ICaseBrowser::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::CaseBrowser::ICaseBrowser::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/CaseBrowser/ICaseBrowser:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::CaseBrowser::ICaseBrowser::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::CaseBrowser::ICaseBrowser_stub_clp (poa, obj);
}

bool
POA_FoamXServer::CaseBrowser::ICaseBrowser::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 61)) {
    case 1:
      if( strcmp( __req->op_name(), "refreshJobsLists" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          refreshJobsLists();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 2:
      if( strcmp( __req->op_name(), "addToCaseList" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addToCaseList( _par_rootDir.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 4:
      if( strcmp( __req->op_name(), "openCasePost" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
        CORBA::Long _par_nProcs;
        CORBA::StaticAny _sa_nProcs( CORBA::_stc_long, &_par_nProcs );

        __req->add_in_arg( &_sa_caseDesc );
        __req->add_in_arg( &_sa_nProcs );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          openCasePost( _par_caseDesc, _par_nProcs );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 8:
      if( strcmp( __req->op_name(), "getUserName" ) == 0 ) {
        char* _par_userName;
        CORBA::StaticAny _sa_userName( CORBA::_stc_string, &_par_userName );

        __req->add_out_arg( &_sa_userName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getUserName( _par_userName );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_userName );
        return true;
      }
      if( strcmp( __req->op_name(), "readFile" ) == 0 ) {
        CORBA::String_var _par_fileName;
        CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName._for_demarshal() );
        char* _par_contents;
        CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents );

        __req->add_in_arg( &_sa_fileName );
        __req->add_out_arg( &_sa_contents );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          readFile( _par_fileName.inout(), _par_contents );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_contents );
        return true;
      }
      if( strcmp( __req->op_name(), "unlockCaseDescriptor" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        __req->add_in_arg( &_sa_caseDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          unlockCaseDescriptor( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 11:
      if( strcmp( __req->op_name(), "newCase" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        CORBA::String_var _par_app;
        CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_in_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          newCase( _par_rootDir.inout(), _par_caseName.inout(), _par_app.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "caseLocked" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_caseDesc );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = caseLocked( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 12:
      if( strcmp( __req->op_name(), "cont" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );

        __req->add_in_arg( &_sa_jobID );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          cont( _par_jobID );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "invokeUtility" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );
        CORBA::String_var _par_utilityName;
        CORBA::StaticAny _sa_utilityName( CORBA::_stc_string, &_par_utilityName._for_demarshal() );
        ::FoamXServer::StringList _par_arguments;
        CORBA::StaticAny _sa_arguments( CORBA::_stcseq_string, &_par_arguments );
        CORBA::String_var _par_logName;
        CORBA::StaticAny _sa_logName( CORBA::_stc_string, &_par_logName._for_demarshal() );
        CORBA::Boolean _par_backGround;
        CORBA::StaticAny _sa_backGround( CORBA::_stc_boolean, &_par_backGround );

        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->add_in_arg( &_sa_hostName );
        __req->add_in_arg( &_sa_utilityName );
        __req->add_in_arg( &_sa_arguments );
        __req->add_in_arg( &_sa_logName );
        __req->add_in_arg( &_sa_backGround );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = invokeUtility( _par_hostName.inout(), _par_utilityName.inout(), _par_arguments, _par_logName.inout(), _par_backGround );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "openCase" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        __req->add_in_arg( &_sa_caseDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          openCase( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "suspend" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );

        __req->add_in_arg( &_sa_jobID );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          suspend( _par_jobID );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 15:
      if( strcmp( __req->op_name(), "deleteCase" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        __req->add_in_arg( &_sa_caseDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          deleteCase( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 16:
      if( strcmp( __req->op_name(), "getEnv" ) == 0 ) {
        CORBA::String_var _par_envName;
        CORBA::StaticAny _sa_envName( CORBA::_stc_string, &_par_envName._for_demarshal() );
        char* _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );

        __req->add_in_arg( &_sa_envName );
        __req->add_out_arg( &_sa_hostName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getEnv( _par_envName.inout(), _par_hostName );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_hostName );
        return true;
      }
      if( strcmp( __req->op_name(), "getCaseServerReference" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        ::FoamXServer::CaseServer::ICaseServer_ptr _par_caseObj;
        CORBA::StaticAny _sa_caseObj( _marshaller_FoamXServer_CaseServer_ICaseServer, &_par_caseObj );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_out_arg( &_sa_caseObj );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getCaseServerReference( _par_rootDir.inout(), _par_caseName.inout(), _par_caseObj );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_caseObj );
        return true;
      }
      break;
    case 18:
      if( strcmp( __req->op_name(), "close" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        close();
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_cases" ) == 0 ) {
        ::FoamXServer::CaseDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_CaseDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = cases();
        __res.value( _marshaller__seq_FoamXServer_CaseDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 22:
      if( strcmp( __req->op_name(), "end" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        CORBA::Boolean _par_now;
        CORBA::StaticAny _sa_now( CORBA::_stc_boolean, &_par_now );

        __req->add_in_arg( &_sa_jobID );
        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_in_arg( &_sa_now );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          end( _par_jobID, _par_rootDir.inout(), _par_caseName.inout(), _par_now );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 23:
      if( strcmp( __req->op_name(), "getCasePostServerReference" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        CORBA::Long _par_nProcs;
        CORBA::StaticAny _sa_nProcs( CORBA::_stc_long, &_par_nProcs );
        ::FoamXServer::CasePostServer::ICasePostServer_ptr _par_caseObj;
        CORBA::StaticAny _sa_caseObj( _marshaller_FoamXServer_CasePostServer_ICasePostServer, &_par_caseObj );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_in_arg( &_sa_nProcs );
        __req->add_out_arg( &_sa_caseObj );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getCasePostServerReference( _par_rootDir.inout(), _par_caseName.inout(), _par_nProcs, _par_caseObj );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_caseObj );
        return true;
      }
      break;
    case 25:
      if( strcmp( __req->op_name(), "purgeFinishedJob" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );

        __req->add_in_arg( &_sa_jobID );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          purgeFinishedJob( _par_jobID );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 26:
      if( strcmp( __req->op_name(), "_get_foamProperties" ) == 0 ) {
        ::FoamXServer::CaseServer::IFoamProperties_ptr _res;
        CORBA::StaticAny __res( _marshaller_FoamXServer_CaseServer_IFoamProperties, &_res );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = foamProperties();
        __req->write_results();
        CORBA::release( _res );
        return true;
      }
      break;
    case 27:
      if( strcmp( __req->op_name(), "refreshCaseList" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          refreshCaseList();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 28:
      if( strcmp( __req->op_name(), "unlockCase" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          unlockCase( _par_rootDir.inout(), _par_caseName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 31:
      if( strcmp( __req->op_name(), "purgeFinishedJobs" ) == 0 ) {
        CORBA::Long _par_nDays;
        CORBA::StaticAny _sa_nDays( CORBA::_stc_long, &_par_nDays );

        __req->add_in_arg( &_sa_nDays );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          purgeFinishedJobs( _par_nDays );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 33:
      if( strcmp( __req->op_name(), "kill" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );

        __req->add_in_arg( &_sa_jobID );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          kill( _par_jobID );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 36:
      if( strcmp( __req->op_name(), "isCaseInError" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_caseDesc );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = isCaseInError( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 37:
      if( strcmp( __req->op_name(), "cloneCase" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );
        CORBA::String_var _par_newCaseRootDir;
        CORBA::StaticAny _sa_newCaseRootDir( CORBA::_stc_string, &_par_newCaseRootDir._for_demarshal() );
        CORBA::String_var _par_newCaseName;
        CORBA::StaticAny _sa_newCaseName( CORBA::_stc_string, &_par_newCaseName._for_demarshal() );
        CORBA::String_var _par_newAppClassName;
        CORBA::StaticAny _sa_newAppClassName( CORBA::_stc_string, &_par_newAppClassName._for_demarshal() );
        CORBA::String_var _par_timeSel;
        CORBA::StaticAny _sa_timeSel( CORBA::_stc_string, &_par_timeSel._for_demarshal() );

        __req->add_in_arg( &_sa_caseDesc );
        __req->add_in_arg( &_sa_newCaseRootDir );
        __req->add_in_arg( &_sa_newCaseName );
        __req->add_in_arg( &_sa_newAppClassName );
        __req->add_in_arg( &_sa_timeSel );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          cloneCase( _par_caseDesc, _par_newCaseRootDir.inout(), _par_newCaseName.inout(), _par_newAppClassName.inout(), _par_timeSel.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 39:
      if( strcmp( __req->op_name(), "fileModificationDate" ) == 0 ) {
        CORBA::String_var _par_fileName;
        CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName._for_demarshal() );

        CORBA::Long _res;
        CORBA::StaticAny __res( CORBA::_stc_long, &_res );
        __req->add_in_arg( &_sa_fileName );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = fileModificationDate( _par_fileName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 40:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "save" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          save();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 44:
      if( strcmp( __req->op_name(), "purgeRunningJobs" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          purgeRunningJobs();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 46:
      if( strcmp( __req->op_name(), "importCase" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        CORBA::String_var _par_app;
        CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_in_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          importCase( _par_rootDir.inout(), _par_caseName.inout(), _par_app.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 49:
      if( strcmp( __req->op_name(), "addCase" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_rawRootDir;
        CORBA::StaticAny _sa_rawRootDir( CORBA::_stc_string, &_par_rawRootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );
        CORBA::String_var _par_app;
        CORBA::StaticAny _sa_app( CORBA::_stc_string, &_par_app._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_rawRootDir );
        __req->add_in_arg( &_sa_caseName );
        __req->add_in_arg( &_sa_app );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          addCase( _par_rootDir.inout(), _par_rawRootDir.inout(), _par_caseName.inout(), _par_app.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "caseOpen" ) == 0 ) {
        CORBA::String_var _par_rootDir;
        CORBA::StaticAny _sa_rootDir( CORBA::_stc_string, &_par_rootDir._for_demarshal() );
        CORBA::String_var _par_caseName;
        CORBA::StaticAny _sa_caseName( CORBA::_stc_string, &_par_caseName._for_demarshal() );

        __req->add_in_arg( &_sa_rootDir );
        __req->add_in_arg( &_sa_caseName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          caseOpen( _par_rootDir.inout(), _par_caseName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 50:
      if( strcmp( __req->op_name(), "checkRunningJobs" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          checkRunningJobs();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 54:
      if( strcmp( __req->op_name(), "caseIsInError" ) == 0 ) {
        ::FoamXServer::CaseDescriptor _par_caseDesc;
        CORBA::StaticAny _sa_caseDesc( _marshaller_FoamXServer_CaseDescriptor, &_par_caseDesc );

        __req->add_in_arg( &_sa_caseDesc );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          caseIsInError( _par_caseDesc );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_finishedJobs" ) == 0 ) {
        ::FoamXServer::JobDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_JobDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = finishedJobs();
        __res.value( _marshaller__seq_FoamXServer_JobDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 56:
      if( strcmp( __req->op_name(), "setStatus" ) == 0 ) {
        ::FoamXServer::JobID _par_jobID;
        CORBA::StaticAny _sa_jobID( _marshaller_FoamXServer_JobID, &_par_jobID );
        ::FoamXServer::JobStatus _par_jobStatus;
        CORBA::StaticAny _sa_jobStatus( _marshaller_FoamXServer_JobStatus, &_par_jobStatus );

        __req->add_in_arg( &_sa_jobID );
        __req->add_in_arg( &_sa_jobStatus );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          setStatus( _par_jobID, _par_jobStatus );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 57:
      if( strcmp( __req->op_name(), "getHostName" ) == 0 ) {
        char* _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName );

        __req->add_out_arg( &_sa_hostName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          getHostName( _par_hostName );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::string_free( _par_hostName );
        return true;
      }
      if( strcmp( __req->op_name(), "writeFile" ) == 0 ) {
        CORBA::String_var _par_fileName;
        CORBA::StaticAny _sa_fileName( CORBA::_stc_string, &_par_fileName._for_demarshal() );
        CORBA::String_var _par_contents;
        CORBA::StaticAny _sa_contents( CORBA::_stc_string, &_par_contents._for_demarshal() );

        __req->add_in_arg( &_sa_fileName );
        __req->add_in_arg( &_sa_contents );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          writeFile( _par_fileName.inout(), _par_contents.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_runningJobs" ) == 0 ) {
        ::FoamXServer::JobDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_JobDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = runningJobs();
        __res.value( _marshaller__seq_FoamXServer_JobDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::CaseBrowser::ICaseBrowser::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}


// PortableServer Skeleton Class for interface FoamXServer::HostBrowser::IHostBrowser
POA_FoamXServer::HostBrowser::IHostBrowser::~IHostBrowser()
{
}

::FoamXServer::HostBrowser::IHostBrowser_ptr
POA_FoamXServer::HostBrowser::IHostBrowser::_this ()
{
  CORBA::Object_var obj = PortableServer::ServantBase::_this();
  return ::FoamXServer::HostBrowser::IHostBrowser::_narrow (obj);
}

CORBA::Boolean
POA_FoamXServer::HostBrowser::IHostBrowser::_is_a (const char * repoid)
{
  if (strcmp (repoid, "IDL:FoamXServer/HostBrowser/IHostBrowser:1.0") == 0) {
    return TRUE;
  }
  return FALSE;
}

CORBA::InterfaceDef_ptr
POA_FoamXServer::HostBrowser::IHostBrowser::_get_interface ()
{
  CORBA::InterfaceDef_ptr ifd = PortableServer::ServantBase::_get_interface ("IDL:FoamXServer/HostBrowser/IHostBrowser:1.0");

  if (CORBA::is_nil (ifd)) {
    mico_throw (CORBA::OBJ_ADAPTER (0, CORBA::COMPLETED_NO));
  }

  return ifd;
}

CORBA::RepositoryId
POA_FoamXServer::HostBrowser::IHostBrowser::_primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr)
{
  return CORBA::string_dup ("IDL:FoamXServer/HostBrowser/IHostBrowser:1.0");
}

CORBA::Object_ptr
POA_FoamXServer::HostBrowser::IHostBrowser::_make_stub (PortableServer::POA_ptr poa, CORBA::Object_ptr obj)
{
  return new ::FoamXServer::HostBrowser::IHostBrowser_stub_clp (poa, obj);
}

bool
POA_FoamXServer::HostBrowser::IHostBrowser::dispatch (CORBA::StaticServerRequest_ptr __req)
{
  #ifdef HAVE_EXCEPTIONS
  try {
  #endif
    switch (mico_string_hash (__req->op_name(), 17)) {
    case 3:
      if( strcmp( __req->op_name(), "close" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        close();
        __req->write_results();
        return true;
      }
      break;
    case 7:
      if( strcmp( __req->op_name(), "isHostAlive" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_hostName );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = isHostAlive( _par_hostName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 8:
      if( strcmp( __req->op_name(), "validate" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          validate();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 9:
      if( strcmp( __req->op_name(), "hostIsDead" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );

        __req->add_in_arg( &_sa_hostName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          hostIsDead( _par_hostName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 12:
      if( strcmp( __req->op_name(), "getCaseBrowserReference" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );
        ::FoamXServer::CaseBrowser::ICaseBrowser_ptr _par_browserObj;
        CORBA::StaticAny _sa_browserObj( _marshaller_FoamXServer_CaseBrowser_ICaseBrowser, &_par_browserObj );

        CORBA::Boolean _res;
        CORBA::StaticAny __res( CORBA::_stc_boolean, &_res );
        __req->add_in_arg( &_sa_hostName );
        __req->add_out_arg( &_sa_browserObj );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          _res = getCaseBrowserReference( _par_hostName.inout(), _par_browserObj );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        CORBA::release( _par_browserObj );
        return true;
      }
      break;
    case 13:
      if( strcmp( __req->op_name(), "hostIsAlive" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );

        __req->add_in_arg( &_sa_hostName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          hostIsAlive( _par_hostName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "openCaseBrowser" ) == 0 ) {
        CORBA::String_var _par_hostName;
        CORBA::StaticAny _sa_hostName( CORBA::_stc_string, &_par_hostName._for_demarshal() );

        __req->add_in_arg( &_sa_hostName );

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          openCaseBrowser( _par_hostName.inout() );
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXSYSError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    case 14:
      if( strcmp( __req->op_name(), "refreshHostList" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          refreshHostList();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      if( strcmp( __req->op_name(), "_get_hosts" ) == 0 ) {
        ::FoamXServer::HostDescriptorList* _res;
        CORBA::StaticAny __res( _marshaller__seq_FoamXServer_HostDescriptor );
        __req->set_result( &__res );

        if( !__req->read_args() )
          return true;

        _res = hosts();
        __res.value( _marshaller__seq_FoamXServer_HostDescriptor, _res );
        __req->write_results();
        delete _res;
        return true;
      }
      break;
    case 16:
      if( strcmp( __req->op_name(), "save" ) == 0 ) {

        if( !__req->read_args() )
          return true;

        #ifdef HAVE_EXCEPTIONS
        try {
        #endif
          save();
        #ifdef HAVE_EXCEPTIONS
        } catch( ::FoamXServer::FoamXError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::FoamXIOError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        } catch( ::FoamXServer::ValidationError_catch &_ex ) {
          __req->set_exception( _ex->_clone() );
          __req->write_results();
          return true;
        }
        #endif
        __req->write_results();
        return true;
      }
      break;
    }
  #ifdef HAVE_EXCEPTIONS
  } catch( CORBA::SystemException_catch &_ex ) {
    __req->set_exception( _ex->_clone() );
    __req->write_results();
    return true;
  } catch( ... ) {
    CORBA::UNKNOWN _ex (CORBA::OMGVMCID | 1, CORBA::COMPLETED_MAYBE);
    __req->set_exception (_ex->_clone());
    __req->write_results ();
    return true;
  }
  #endif

  return false;
}

void
POA_FoamXServer::HostBrowser::IHostBrowser::invoke (CORBA::StaticServerRequest_ptr __req)
{
  if (dispatch (__req)) {
      return;
  }

  CORBA::Exception * ex = 
    new CORBA::BAD_OPERATION (0, CORBA::COMPLETED_NO);
  __req->set_exception (ex);
  __req->write_results();
}

