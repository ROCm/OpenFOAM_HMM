/*
 *  MICO --- an Open Source CORBA implementation
 *  Copyright (c) 1997-2006 by The Mico Team
 *
 *  This file was automatically generated. DO NOT EDIT!
 */

#include <CORBA.h>
#include <mico/throw.h>

#ifndef __FOAMX_H__
#define __FOAMX_H__






namespace FoamXServer
{

class ITypeDescriptor;
typedef ITypeDescriptor *ITypeDescriptor_ptr;
typedef ITypeDescriptor_ptr ITypeDescriptorRef;
typedef ObjVar< ITypeDescriptor > ITypeDescriptor_var;
typedef ObjOut< ITypeDescriptor > ITypeDescriptor_out;

class IDictionaryEntry;
typedef IDictionaryEntry *IDictionaryEntry_ptr;
typedef IDictionaryEntry_ptr IDictionaryEntryRef;
typedef ObjVar< IDictionaryEntry > IDictionaryEntry_var;
typedef ObjOut< IDictionaryEntry > IDictionaryEntry_out;



namespace CaseServer
{

class IApplication;
typedef IApplication *IApplication_ptr;
typedef IApplication_ptr IApplicationRef;
typedef ObjVar< IApplication > IApplication_var;
typedef ObjOut< IApplication > IApplication_out;

class IFoamProperties;
typedef IFoamProperties *IFoamProperties_ptr;
typedef IFoamProperties_ptr IFoamPropertiesRef;
typedef ObjVar< IFoamProperties > IFoamProperties_var;
typedef ObjOut< IFoamProperties > IFoamProperties_out;

class IGeometricField;
typedef IGeometricField *IGeometricField_ptr;
typedef IGeometricField_ptr IGeometricFieldRef;
typedef ObjVar< IGeometricField > IGeometricField_var;
typedef ObjOut< IGeometricField > IGeometricField_out;

class ICaseServer;
typedef ICaseServer *ICaseServer_ptr;
typedef ICaseServer_ptr ICaseServerRef;
typedef ObjVar< ICaseServer > ICaseServer_var;
typedef ObjOut< ICaseServer > ICaseServer_out;

class IGeometryDescriptor;
typedef IGeometryDescriptor *IGeometryDescriptor_ptr;
typedef IGeometryDescriptor_ptr IGeometryDescriptorRef;
typedef ObjVar< IGeometryDescriptor > IGeometryDescriptor_var;
typedef ObjOut< IGeometryDescriptor > IGeometryDescriptor_out;

class IPatchDescriptor;
typedef IPatchDescriptor *IPatchDescriptor_ptr;
typedef IPatchDescriptor_ptr IPatchDescriptorRef;
typedef ObjVar< IPatchDescriptor > IPatchDescriptor_var;
typedef ObjOut< IPatchDescriptor > IPatchDescriptor_out;

class IGeometricFieldDescriptor;
typedef IGeometricFieldDescriptor *IGeometricFieldDescriptor_ptr;
typedef IGeometricFieldDescriptor_ptr IGeometricFieldDescriptorRef;
typedef ObjVar< IGeometricFieldDescriptor > IGeometricFieldDescriptor_var;
typedef ObjOut< IGeometricFieldDescriptor > IGeometricFieldDescriptor_out;

class IPatchPhysicalTypeDescriptor;
typedef IPatchPhysicalTypeDescriptor *IPatchPhysicalTypeDescriptor_ptr;
typedef IPatchPhysicalTypeDescriptor_ptr IPatchPhysicalTypeDescriptorRef;
typedef ObjVar< IPatchPhysicalTypeDescriptor > IPatchPhysicalTypeDescriptor_var;
typedef ObjOut< IPatchPhysicalTypeDescriptor > IPatchPhysicalTypeDescriptor_out;

}


namespace CasePostServer
{

class ICasePostServer;
typedef ICasePostServer *ICasePostServer_ptr;
typedef ICasePostServer_ptr ICasePostServerRef;
typedef ObjVar< ICasePostServer > ICasePostServer_var;
typedef ObjOut< ICasePostServer > ICasePostServer_out;

}


namespace CaseBrowser
{

class ICaseBrowser;
typedef ICaseBrowser *ICaseBrowser_ptr;
typedef ICaseBrowser_ptr ICaseBrowserRef;
typedef ObjVar< ICaseBrowser > ICaseBrowser_var;
typedef ObjOut< ICaseBrowser > ICaseBrowser_out;

}


namespace HostBrowser
{

class IHostBrowser;
typedef IHostBrowser *IHostBrowser_ptr;
typedef IHostBrowser_ptr IHostBrowserRef;
typedef ObjVar< IHostBrowser > IHostBrowser_var;
typedef ObjOut< IHostBrowser > IHostBrowser_out;

}
}






namespace FoamXServer
{

enum FoamXType {
  Type_Undefined = 0,
  Type_Boolean,
  Type_Label,
  Type_Scalar,
  Type_Char,
  Type_Word,
  Type_String,
  Type_RootDir,
  Type_RootAndCase,
  Type_CaseName,
  Type_HostName,
  Type_File,
  Type_Directory,
  Type_Time,
  Type_DimensionSet,
  Type_FixedList,
  Type_List,
  Type_Dictionary,
  Type_Selection,
  Type_Compound,
  Type_Field
};

typedef FoamXType& FoamXType_out;

extern CORBA::TypeCodeConst _tc_FoamXType;

struct FoamXAny;
typedef TVarVar< FoamXAny > FoamXAny_var;
typedef TVarOut< FoamXAny > FoamXAny_out;


struct FoamXAny {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef FoamXAny_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  FoamXAny();
  ~FoamXAny();
  FoamXAny( const FoamXAny& s );
  FoamXAny& operator=( const FoamXAny& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  FoamXType type;
  CORBA::Any value;
};

extern CORBA::TypeCodeConst _tc_FoamXAny;

typedef SequenceTmpl< FoamXAny,MICO_TID_DEF> FoamXAnyList;
typedef TSeqVar< SequenceTmpl< FoamXAny,MICO_TID_DEF> > FoamXAnyList_var;
typedef TSeqOut< SequenceTmpl< FoamXAny,MICO_TID_DEF> > FoamXAnyList_out;

extern CORBA::TypeCodeConst _tc_FoamXAnyList;

typedef StringSequenceTmpl<CORBA::String_var> StringList;
typedef TSeqVar< StringSequenceTmpl<CORBA::String_var> > StringList_var;
typedef TSeqOut< StringSequenceTmpl<CORBA::String_var> > StringList_out;

extern CORBA::TypeCodeConst _tc_StringList;

typedef IfaceSequenceTmpl< ITypeDescriptor_var,ITypeDescriptor_ptr> TypeDescriptorList;
typedef TSeqVar< IfaceSequenceTmpl< ITypeDescriptor_var,ITypeDescriptor_ptr> > TypeDescriptorList_var;
typedef TSeqOut< IfaceSequenceTmpl< ITypeDescriptor_var,ITypeDescriptor_ptr> > TypeDescriptorList_out;

extern CORBA::TypeCodeConst _tc_TypeDescriptorList;

enum ErrorCode {
  S_OK = 0,
  E_FAIL,
  E_FOAM,
  E_INVALID_ARG,
  E_INVALID_PTR,
  E_INVALID_REF,
  E_INDEX_OUT_OF_BOUNDS,
  E_UNKNOWN_NAME,
  E_UNEXPECTED
};

typedef ErrorCode& ErrorCode_out;

extern CORBA::TypeCodeConst _tc_ErrorCode;

struct FoamXError : public CORBA::UserException {
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  FoamXError();
  ~FoamXError();
  FoamXError( const FoamXError& s );
  FoamXError& operator=( const FoamXError& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  #ifndef HAVE_EXPLICIT_STRUCT_OPS
  FoamXError();
  #endif //HAVE_EXPLICIT_STRUCT_OPS
  FoamXError( ::FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2, const char* _m3, CORBA::Long _m4 );

  #ifdef HAVE_STD_EH
  FoamXError *operator->() { return this; }
  FoamXError& operator*() { return *this; }
  operator FoamXError*() { return this; }
  #endif // HAVE_STD_EH

  void _throwit() const;
  const char *_repoid() const;
  void _encode( CORBA::DataEncoder &en ) const;
  void _encode_any( CORBA::Any &a ) const;
  CORBA::Exception *_clone() const;
  static FoamXError *_downcast( CORBA::Exception *ex );
  static const FoamXError *_downcast( const CORBA::Exception *ex );
  ErrorCode errorCode;
  CORBA::String_var errorMessage;
  CORBA::String_var methodName;
  CORBA::String_var fileName;
  CORBA::Long lineNo;
};

#ifdef HAVE_STD_EH
typedef FoamXError FoamXError_catch;
#else
typedef ExceptVar< FoamXError > FoamXError_var;
typedef TVarOut< FoamXError > FoamXError_out;
typedef FoamXError_var FoamXError_catch;
#endif // HAVE_STD_EH

extern CORBA::TypeCodeConst _tc_FoamXError;

struct FoamXIOError : public CORBA::UserException {
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  FoamXIOError();
  ~FoamXIOError();
  FoamXIOError( const FoamXIOError& s );
  FoamXIOError& operator=( const FoamXIOError& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  #ifndef HAVE_EXPLICIT_STRUCT_OPS
  FoamXIOError();
  #endif //HAVE_EXPLICIT_STRUCT_OPS
  FoamXIOError( const char* _m0, const char* _m1, CORBA::Long _m2, CORBA::Long _m3, const char* _m4, const char* _m5, CORBA::Long _m6 );

  #ifdef HAVE_STD_EH
  FoamXIOError *operator->() { return this; }
  FoamXIOError& operator*() { return *this; }
  operator FoamXIOError*() { return this; }
  #endif // HAVE_STD_EH

  void _throwit() const;
  const char *_repoid() const;
  void _encode( CORBA::DataEncoder &en ) const;
  void _encode_any( CORBA::Any &a ) const;
  CORBA::Exception *_clone() const;
  static FoamXIOError *_downcast( CORBA::Exception *ex );
  static const FoamXIOError *_downcast( const CORBA::Exception *ex );
  CORBA::String_var errorMessage;
  CORBA::String_var ioFileName;
  CORBA::Long ioStartLineNo;
  CORBA::Long ioEndLineNo;
  CORBA::String_var methodName;
  CORBA::String_var fileName;
  CORBA::Long lineNo;
};

#ifdef HAVE_STD_EH
typedef FoamXIOError FoamXIOError_catch;
#else
typedef ExceptVar< FoamXIOError > FoamXIOError_var;
typedef TVarOut< FoamXIOError > FoamXIOError_out;
typedef FoamXIOError_var FoamXIOError_catch;
#endif // HAVE_STD_EH

extern CORBA::TypeCodeConst _tc_FoamXIOError;

struct ValidationError : public CORBA::UserException {
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  ValidationError();
  ~ValidationError();
  ValidationError( const ValidationError& s );
  ValidationError& operator=( const ValidationError& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  #ifndef HAVE_EXPLICIT_STRUCT_OPS
  ValidationError();
  #endif //HAVE_EXPLICIT_STRUCT_OPS
  ValidationError( ::FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2 );

  #ifdef HAVE_STD_EH
  ValidationError *operator->() { return this; }
  ValidationError& operator*() { return *this; }
  operator ValidationError*() { return this; }
  #endif // HAVE_STD_EH

  void _throwit() const;
  const char *_repoid() const;
  void _encode( CORBA::DataEncoder &en ) const;
  void _encode_any( CORBA::Any &a ) const;
  CORBA::Exception *_clone() const;
  static ValidationError *_downcast( CORBA::Exception *ex );
  static const ValidationError *_downcast( const CORBA::Exception *ex );
  ErrorCode errorCode;
  CORBA::String_var errorMessage;
  CORBA::String_var itemPath;
};

#ifdef HAVE_STD_EH
typedef ValidationError ValidationError_catch;
#else
typedef ExceptVar< ValidationError > ValidationError_var;
typedef TVarOut< ValidationError > ValidationError_out;
typedef ValidationError_var ValidationError_catch;
#endif // HAVE_STD_EH

extern CORBA::TypeCodeConst _tc_ValidationError;


/*
 * Base class and common definitions for interface ITypeDescriptor
 */

class ITypeDescriptor : 
  virtual public CORBA::Object
{
  public:
    virtual ~ITypeDescriptor();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef ITypeDescriptor_ptr _ptr_type;
    typedef ITypeDescriptor_var _var_type;
    #endif

    static ITypeDescriptor_ptr _narrow( CORBA::Object_ptr obj );
    static ITypeDescriptor_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static ITypeDescriptor_ptr _duplicate( ITypeDescriptor_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static ITypeDescriptor_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::FoamXType type() = 0;
    virtual void type( ::FoamXServer::FoamXType value ) = 0;
    virtual CORBA::Boolean isPrimitiveType() = 0;
    virtual CORBA::Boolean isCompoundType() = 0;
    virtual char* path() = 0;
    virtual void path( const char* value ) = 0;
    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* comment() = 0;
    virtual void comment( const char* value ) = 0;
    virtual char* category() = 0;
    virtual void category( const char* value ) = 0;
    virtual char* helpURL() = 0;
    virtual void helpURL( const char* value ) = 0;
    virtual char* iconURL() = 0;
    virtual void iconURL( const char* value ) = 0;
    virtual CORBA::Boolean optional() = 0;
    virtual void optional( CORBA::Boolean value ) = 0;
    virtual CORBA::Boolean visible() = 0;
    virtual void visible( CORBA::Boolean value ) = 0;
    virtual CORBA::Boolean editable() = 0;
    virtual void editable( CORBA::Boolean value ) = 0;
    virtual ::FoamXServer::FoamXAny* minValue() = 0;
    virtual void minValue( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual ::FoamXServer::FoamXAny* maxValue() = 0;
    virtual void maxValue( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual char* lookupDict() = 0;
    virtual void lookupDict( const char* value ) = 0;
    virtual ::FoamXServer::FoamXAnyList* valueList() = 0;
    virtual void valueList( const ::FoamXServer::FoamXAnyList& value ) = 0;
    virtual char* dictionaryPath() = 0;
    virtual void dictionaryPath( const char* value ) = 0;
    virtual CORBA::Long numElements() = 0;
    virtual void numElements( CORBA::Long value ) = 0;
    virtual ::FoamXServer::StringList* elementLabels() = 0;
    virtual void elementLabels( const ::FoamXServer::StringList& value ) = 0;
    virtual ::FoamXServer::TypeDescriptorList* subTypes() = 0;
    virtual ::FoamXServer::ITypeDescriptor_ptr elementType() = 0;

    virtual CORBA::Boolean hasDefaultValue() = 0;
    virtual void getDefaultValue( IDictionaryEntry_out defaultValue ) = 0;
    virtual void addSubType( FoamXType type, ITypeDescriptor_out subEntry ) = 0;
    virtual void removeSubType( ITypeDescriptor_ptr subEntry ) = 0;
    virtual void validate() = 0;

  protected:
    ITypeDescriptor() {};
  private:
    ITypeDescriptor( const ITypeDescriptor& );
    void operator=( const ITypeDescriptor& );
};

extern CORBA::TypeCodeConst _tc_ITypeDescriptor;

// Stub for interface ITypeDescriptor
class ITypeDescriptor_stub:
  virtual public ITypeDescriptor
{
  public:
    virtual ~ITypeDescriptor_stub();
    ::FoamXServer::FoamXType type();
    void type( ::FoamXServer::FoamXType value );
    CORBA::Boolean isPrimitiveType();
    CORBA::Boolean isCompoundType();
    char* path();
    void path( const char* value );
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );
    char* comment();
    void comment( const char* value );
    char* category();
    void category( const char* value );
    char* helpURL();
    void helpURL( const char* value );
    char* iconURL();
    void iconURL( const char* value );
    CORBA::Boolean optional();
    void optional( CORBA::Boolean value );
    CORBA::Boolean visible();
    void visible( CORBA::Boolean value );
    CORBA::Boolean editable();
    void editable( CORBA::Boolean value );
    ::FoamXServer::FoamXAny* minValue();
    void minValue( const ::FoamXServer::FoamXAny& value );
    ::FoamXServer::FoamXAny* maxValue();
    void maxValue( const ::FoamXServer::FoamXAny& value );
    char* lookupDict();
    void lookupDict( const char* value );
    ::FoamXServer::FoamXAnyList* valueList();
    void valueList( const ::FoamXServer::FoamXAnyList& value );
    char* dictionaryPath();
    void dictionaryPath( const char* value );
    CORBA::Long numElements();
    void numElements( CORBA::Long value );
    ::FoamXServer::StringList* elementLabels();
    void elementLabels( const ::FoamXServer::StringList& value );
    ::FoamXServer::TypeDescriptorList* subTypes();
    ::FoamXServer::ITypeDescriptor_ptr elementType();

    CORBA::Boolean hasDefaultValue();
    void getDefaultValue( IDictionaryEntry_out defaultValue );
    void addSubType( FoamXType type, ITypeDescriptor_out subEntry );
    void removeSubType( ITypeDescriptor_ptr subEntry );
    void validate();

  private:
    void operator=( const ITypeDescriptor_stub& );
};

#ifndef MICO_CONF_NO_POA

class ITypeDescriptor_stub_clp :
  virtual public ITypeDescriptor_stub,
  virtual public PortableServer::StubBase
{
  public:
    ITypeDescriptor_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~ITypeDescriptor_stub_clp ();
    ::FoamXServer::FoamXType type();
    void type( ::FoamXServer::FoamXType value );
    CORBA::Boolean isPrimitiveType();
    CORBA::Boolean isCompoundType();
    char* path();
    void path( const char* value );
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );
    char* comment();
    void comment( const char* value );
    char* category();
    void category( const char* value );
    char* helpURL();
    void helpURL( const char* value );
    char* iconURL();
    void iconURL( const char* value );
    CORBA::Boolean optional();
    void optional( CORBA::Boolean value );
    CORBA::Boolean visible();
    void visible( CORBA::Boolean value );
    CORBA::Boolean editable();
    void editable( CORBA::Boolean value );
    ::FoamXServer::FoamXAny* minValue();
    void minValue( const ::FoamXServer::FoamXAny& value );
    ::FoamXServer::FoamXAny* maxValue();
    void maxValue( const ::FoamXServer::FoamXAny& value );
    char* lookupDict();
    void lookupDict( const char* value );
    ::FoamXServer::FoamXAnyList* valueList();
    void valueList( const ::FoamXServer::FoamXAnyList& value );
    char* dictionaryPath();
    void dictionaryPath( const char* value );
    CORBA::Long numElements();
    void numElements( CORBA::Long value );
    ::FoamXServer::StringList* elementLabels();
    void elementLabels( const ::FoamXServer::StringList& value );
    ::FoamXServer::TypeDescriptorList* subTypes();
    ::FoamXServer::ITypeDescriptor_ptr elementType();

    CORBA::Boolean hasDefaultValue();
    void getDefaultValue( IDictionaryEntry_out defaultValue );
    void addSubType( FoamXType type, ITypeDescriptor_out subEntry );
    void removeSubType( ITypeDescriptor_ptr subEntry );
    void validate();

  protected:
    ITypeDescriptor_stub_clp ();
  private:
    void operator=( const ITypeDescriptor_stub_clp & );
};

#endif // MICO_CONF_NO_POA

typedef IfaceSequenceTmpl< IDictionaryEntry_var,IDictionaryEntry_ptr> DictionaryEntryList;
typedef TSeqVar< IfaceSequenceTmpl< IDictionaryEntry_var,IDictionaryEntry_ptr> > DictionaryEntryList_var;
typedef TSeqOut< IfaceSequenceTmpl< IDictionaryEntry_var,IDictionaryEntry_ptr> > DictionaryEntryList_out;

extern CORBA::TypeCodeConst _tc_DictionaryEntryList;


/*
 * Base class and common definitions for interface IDictionaryEntry
 */

class IDictionaryEntry : 
  virtual public CORBA::Object
{
  public:
    virtual ~IDictionaryEntry();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IDictionaryEntry_ptr _ptr_type;
    typedef IDictionaryEntry_var _var_type;
    #endif

    static IDictionaryEntry_ptr _narrow( CORBA::Object_ptr obj );
    static IDictionaryEntry_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IDictionaryEntry_ptr _duplicate( IDictionaryEntry_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IDictionaryEntry_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::ITypeDescriptor_ptr typeDescriptor() = 0;
    virtual ::FoamXServer::FoamXAny* value() = 0;
    virtual void value( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual ::FoamXServer::DictionaryEntryList* subElements() = 0;
    virtual CORBA::Long selection() = 0;
    virtual void selection( CORBA::Long value ) = 0;

    virtual void setValue( const FoamXAny& value ) = 0;
    virtual CORBA::Long nSubElements() = 0;
    virtual CORBA::Boolean packedList() = 0;
    virtual void addElement( IDictionaryEntry_out subEntry ) = 0;
    virtual void removeElement( IDictionaryEntry_ptr subEntry ) = 0;
    virtual void validate() = 0;
    virtual CORBA::Boolean modified() = 0;
    virtual void save() = 0;

  protected:
    IDictionaryEntry() {};
  private:
    IDictionaryEntry( const IDictionaryEntry& );
    void operator=( const IDictionaryEntry& );
};

extern CORBA::TypeCodeConst _tc_IDictionaryEntry;

// Stub for interface IDictionaryEntry
class IDictionaryEntry_stub:
  virtual public IDictionaryEntry
{
  public:
    virtual ~IDictionaryEntry_stub();
    ::FoamXServer::ITypeDescriptor_ptr typeDescriptor();
    ::FoamXServer::FoamXAny* value();
    void value( const ::FoamXServer::FoamXAny& value );
    ::FoamXServer::DictionaryEntryList* subElements();
    CORBA::Long selection();
    void selection( CORBA::Long value );

    void setValue( const FoamXAny& value );
    CORBA::Long nSubElements();
    CORBA::Boolean packedList();
    void addElement( IDictionaryEntry_out subEntry );
    void removeElement( IDictionaryEntry_ptr subEntry );
    void validate();
    CORBA::Boolean modified();
    void save();

  private:
    void operator=( const IDictionaryEntry_stub& );
};

#ifndef MICO_CONF_NO_POA

class IDictionaryEntry_stub_clp :
  virtual public IDictionaryEntry_stub,
  virtual public PortableServer::StubBase
{
  public:
    IDictionaryEntry_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IDictionaryEntry_stub_clp ();
    ::FoamXServer::ITypeDescriptor_ptr typeDescriptor();
    ::FoamXServer::FoamXAny* value();
    void value( const ::FoamXServer::FoamXAny& value );
    ::FoamXServer::DictionaryEntryList* subElements();
    CORBA::Long selection();
    void selection( CORBA::Long value );

    void setValue( const FoamXAny& value );
    CORBA::Long nSubElements();
    CORBA::Boolean packedList();
    void addElement( IDictionaryEntry_out subEntry );
    void removeElement( IDictionaryEntry_ptr subEntry );
    void validate();
    CORBA::Boolean modified();
    void save();

  protected:
    IDictionaryEntry_stub_clp ();
  private:
    void operator=( const IDictionaryEntry_stub_clp & );
};

#endif // MICO_CONF_NO_POA

enum JobStatus {
  JOB_UNDEFINED = 0,
  JOB_LAUNCHING,
  JOB_RUNNING,
  JOB_STOPPING,
  JOB_SUSPENDED,
  JOB_FINISHED,
  JOB_ABORTED
};

typedef JobStatus& JobStatus_out;

extern CORBA::TypeCodeConst _tc_JobStatus;

enum MessageType {
  M_DIAGNOSTIC = 0,
  M_WARNING,
  M_ERROR
};

typedef MessageType& MessageType_out;

extern CORBA::TypeCodeConst _tc_MessageType;

typedef SequenceTmpl< CORBA::Double,MICO_TID_DEF> DoubleList;
typedef TSeqVar< SequenceTmpl< CORBA::Double,MICO_TID_DEF> > DoubleList_var;
typedef TSeqOut< SequenceTmpl< CORBA::Double,MICO_TID_DEF> > DoubleList_out;

extern CORBA::TypeCodeConst _tc_DoubleList;

typedef SequenceTmpl< CORBA::Float,MICO_TID_DEF> FloatList;
typedef TSeqVar< SequenceTmpl< CORBA::Float,MICO_TID_DEF> > FloatList_var;
typedef TSeqOut< SequenceTmpl< CORBA::Float,MICO_TID_DEF> > FloatList_out;

extern CORBA::TypeCodeConst _tc_FloatList;

typedef SequenceTmpl< CORBA::Long,MICO_TID_DEF> LongList;
typedef TSeqVar< SequenceTmpl< CORBA::Long,MICO_TID_DEF> > LongList_var;
typedef TSeqOut< SequenceTmpl< CORBA::Long,MICO_TID_DEF> > LongList_out;

extern CORBA::TypeCodeConst _tc_LongList;

typedef CORBA::Float Point3[ 3 ];
typedef CORBA::Float Point3_slice;
typedef ArrayFixVar< CORBA::Float,Point3_slice,Point3,3 > Point3_var;
typedef Point3 Point3_out;
enum _dummy_Point3 { _dummy_Point3_0 };
typedef ArrayFixForAny < CORBA::Float,Point3_slice,Point3,3,_dummy_Point3> Point3_forany;
#undef MICO_ARRAY_ARG
#define MICO_ARRAY_ARG CORBA::Float
DECLARE_ARRAY_ALLOC(Point3,MICO_ARRAY_ARG,Point3_slice,3)
DECLARE_ARRAY_DUP(Point3,MICO_ARRAY_ARG,Point3_slice,3)
DECLARE_ARRAY_FREE(Point3,MICO_ARRAY_ARG,Point3_slice,3)
DECLARE_ARRAY_COPY(Point3,MICO_ARRAY_ARG,Point3_slice,3)
extern CORBA::TypeCodeConst _tc_Point3;

struct StringPair;
typedef TVarVar< StringPair > StringPair_var;
typedef TVarOut< StringPair > StringPair_out;


struct StringPair {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef StringPair_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  StringPair();
  ~StringPair();
  StringPair( const StringPair& s );
  StringPair& operator=( const StringPair& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::String_var name;
  CORBA::String_var value;
};

extern CORBA::TypeCodeConst _tc_StringPair;

typedef SequenceTmpl< StringPair,MICO_TID_DEF> StringPairList;
typedef TSeqVar< SequenceTmpl< StringPair,MICO_TID_DEF> > StringPairList_var;
typedef TSeqOut< SequenceTmpl< StringPair,MICO_TID_DEF> > StringPairList_out;

extern CORBA::TypeCodeConst _tc_StringPairList;

struct HostDescriptor;
typedef TVarVar< HostDescriptor > HostDescriptor_var;
typedef TVarOut< HostDescriptor > HostDescriptor_out;


struct HostDescriptor {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef HostDescriptor_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  HostDescriptor();
  ~HostDescriptor();
  HostDescriptor( const HostDescriptor& s );
  HostDescriptor& operator=( const HostDescriptor& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::String_var name;
  CORBA::Boolean alive;
};

extern CORBA::TypeCodeConst _tc_HostDescriptor;

typedef SequenceTmpl< HostDescriptor,MICO_TID_DEF> HostDescriptorList;
typedef TSeqVar< SequenceTmpl< HostDescriptor,MICO_TID_DEF> > HostDescriptorList_var;
typedef TSeqOut< SequenceTmpl< HostDescriptor,MICO_TID_DEF> > HostDescriptorList_out;

extern CORBA::TypeCodeConst _tc_HostDescriptorList;

struct ApplicationDescriptor;
typedef TVarVar< ApplicationDescriptor > ApplicationDescriptor_var;
typedef TVarOut< ApplicationDescriptor > ApplicationDescriptor_out;


struct ApplicationDescriptor {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef ApplicationDescriptor_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  ApplicationDescriptor();
  ~ApplicationDescriptor();
  ApplicationDescriptor( const ApplicationDescriptor& s );
  ApplicationDescriptor& operator=( const ApplicationDescriptor& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::String_var name;
  CORBA::String_var category;
  CORBA::String_var path;
  CORBA::Boolean systemClass;
};

extern CORBA::TypeCodeConst _tc_ApplicationDescriptor;

typedef SequenceTmpl< ApplicationDescriptor,MICO_TID_DEF> ApplicationDescriptorList;
typedef TSeqVar< SequenceTmpl< ApplicationDescriptor,MICO_TID_DEF> > ApplicationDescriptorList_var;
typedef TSeqOut< SequenceTmpl< ApplicationDescriptor,MICO_TID_DEF> > ApplicationDescriptorList_out;

extern CORBA::TypeCodeConst _tc_ApplicationDescriptorList;

struct CaseDescriptor;
typedef TVarVar< CaseDescriptor > CaseDescriptor_var;
typedef TVarOut< CaseDescriptor > CaseDescriptor_out;


struct CaseDescriptor {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef CaseDescriptor_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  CaseDescriptor();
  ~CaseDescriptor();
  CaseDescriptor( const CaseDescriptor& s );
  CaseDescriptor& operator=( const CaseDescriptor& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::String_var rootDir;
  CORBA::String_var rawRootDir;
  CORBA::String_var caseName;
  CORBA::String_var app;
  CORBA::Long nProcs;
  CORBA::Boolean managed;
  CORBA::Boolean locked;
  CORBA::Boolean error;
};

extern CORBA::TypeCodeConst _tc_CaseDescriptor;

typedef SequenceTmpl< CaseDescriptor,MICO_TID_DEF> CaseDescriptorList;
typedef TSeqVar< SequenceTmpl< CaseDescriptor,MICO_TID_DEF> > CaseDescriptorList_var;
typedef TSeqOut< SequenceTmpl< CaseDescriptor,MICO_TID_DEF> > CaseDescriptorList_out;

extern CORBA::TypeCodeConst _tc_CaseDescriptorList;

struct JobID;
typedef TVarVar< JobID > JobID_var;
typedef TVarOut< JobID > JobID_out;


struct JobID {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef JobID_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  JobID();
  ~JobID();
  JobID( const JobID& s );
  JobID& operator=( const JobID& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::String_var hostName;
  CORBA::Long processID;
};

extern CORBA::TypeCodeConst _tc_JobID;

typedef SequenceTmpl< JobID,MICO_TID_DEF> JobIDList;
typedef TSeqVar< SequenceTmpl< JobID,MICO_TID_DEF> > JobIDList_var;
typedef TSeqOut< SequenceTmpl< JobID,MICO_TID_DEF> > JobIDList_out;

extern CORBA::TypeCodeConst _tc_JobIDList;

struct JobDescriptor;
typedef TVarVar< JobDescriptor > JobDescriptor_var;
typedef TVarOut< JobDescriptor > JobDescriptor_out;


struct JobDescriptor {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef JobDescriptor_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  JobDescriptor();
  ~JobDescriptor();
  JobDescriptor( const JobDescriptor& s );
  JobDescriptor& operator=( const JobDescriptor& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  JobID jobID;
  CORBA::Long ppid;
  CORBA::Long pgid;
  CORBA::String_var startDate;
  CORBA::String_var startTime;
  CORBA::String_var userName;
  CORBA::String_var foamVersion;
  CORBA::String_var code;
  CORBA::String_var argList;
  CORBA::String_var currentDir;
  CORBA::String_var rootDir;
  CORBA::String_var caseName;
  CORBA::Long nProcs;
  JobIDList slaves;
  CORBA::Long nCountedProcs;
  CORBA::Double cpuTime;
  CORBA::String_var endDate;
  CORBA::String_var endTime;
  JobStatus status;
};

extern CORBA::TypeCodeConst _tc_JobDescriptor;

typedef SequenceTmpl< JobDescriptor,MICO_TID_DEF> JobDescriptorList;
typedef TSeqVar< SequenceTmpl< JobDescriptor,MICO_TID_DEF> > JobDescriptorList_var;
typedef TSeqOut< SequenceTmpl< JobDescriptor,MICO_TID_DEF> > JobDescriptorList_out;

extern CORBA::TypeCodeConst _tc_JobDescriptorList;

struct DimensionSet;
typedef TFixVar< DimensionSet > DimensionSet_var;
typedef DimensionSet& DimensionSet_out;


struct DimensionSet {
  #ifdef HAVE_TYPEDEF_OVERLOAD
  typedef DimensionSet_var _var_type;
  #endif
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  DimensionSet();
  ~DimensionSet();
  DimensionSet( const DimensionSet& s );
  DimensionSet& operator=( const DimensionSet& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  CORBA::Double mass;
  CORBA::Double length;
  CORBA::Double time;
  CORBA::Double temperature;
  CORBA::Double moles;
  CORBA::Double current;
  CORBA::Double luminousIntensity;
};

extern CORBA::TypeCodeConst _tc_DimensionSet;

typedef SequenceTmpl< DimensionSet,MICO_TID_DEF> DimensionSetList;
typedef TSeqVar< SequenceTmpl< DimensionSet,MICO_TID_DEF> > DimensionSetList_var;
typedef TSeqOut< SequenceTmpl< DimensionSet,MICO_TID_DEF> > DimensionSetList_out;

extern CORBA::TypeCodeConst _tc_DimensionSetList;

struct FoamXSYSError : public CORBA::UserException {
  #ifdef HAVE_EXPLICIT_STRUCT_OPS
  FoamXSYSError();
  ~FoamXSYSError();
  FoamXSYSError( const FoamXSYSError& s );
  FoamXSYSError& operator=( const FoamXSYSError& s );
  #endif //HAVE_EXPLICIT_STRUCT_OPS

  #ifndef HAVE_EXPLICIT_STRUCT_OPS
  FoamXSYSError();
  #endif //HAVE_EXPLICIT_STRUCT_OPS
  FoamXSYSError( ::FoamXServer::ErrorCode _m0, const char* _m1, const char* _m2, const char* _m3, const char* _m4, CORBA::Long _m5 );

  #ifdef HAVE_STD_EH
  FoamXSYSError *operator->() { return this; }
  FoamXSYSError& operator*() { return *this; }
  operator FoamXSYSError*() { return this; }
  #endif // HAVE_STD_EH

  void _throwit() const;
  const char *_repoid() const;
  void _encode( CORBA::DataEncoder &en ) const;
  void _encode_any( CORBA::Any &a ) const;
  CORBA::Exception *_clone() const;
  static FoamXSYSError *_downcast( CORBA::Exception *ex );
  static const FoamXSYSError *_downcast( const CORBA::Exception *ex );
  ErrorCode errorCode;
  CORBA::String_var errorMessage;
  CORBA::String_var hostName;
  CORBA::String_var methodName;
  CORBA::String_var fileName;
  CORBA::Long lineNo;
};

#ifdef HAVE_STD_EH
typedef FoamXSYSError FoamXSYSError_catch;
#else
typedef ExceptVar< FoamXSYSError > FoamXSYSError_var;
typedef TVarOut< FoamXSYSError > FoamXSYSError_out;
typedef FoamXSYSError_var FoamXSYSError_catch;
#endif // HAVE_STD_EH

extern CORBA::TypeCodeConst _tc_FoamXSYSError;



namespace CaseServer
{


/*
 * Base class and common definitions for interface ICaseServer
 */

class ICaseServer : 
  virtual public CORBA::Object
{
  public:
    virtual ~ICaseServer();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef ICaseServer_ptr _ptr_type;
    typedef ICaseServer_var _var_type;
    #endif

    static ICaseServer_ptr _narrow( CORBA::Object_ptr obj );
    static ICaseServer_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static ICaseServer_ptr _duplicate( ICaseServer_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static ICaseServer_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual CORBA::Boolean managed() = 0;
    virtual void managed( CORBA::Boolean value ) = 0;
    virtual char* caseRoot() = 0;
    virtual char* caseName() = 0;
    virtual ::FoamXServer::CaseServer::IApplication_ptr application() = 0;
    virtual ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties() = 0;
    virtual ::FoamXServer::StringList* availableTimeSteps() = 0;
    virtual CORBA::Boolean meshDefined() = 0;
    virtual ::FoamXServer::StringList* patchNames() = 0;

    virtual char* getTime() = 0;
    virtual void setTime( const char* timeName, CORBA::Long timeIndex ) = 0;
    virtual void readMesh() = 0;
    virtual void importMesh( const char* hostName, const char* rootDir, const char* caseName ) = 0;
    virtual void getFieldValues( const char* fieldName, IGeometricField_out fieldValues ) = 0;
    virtual void addPatch( const char* patchName, const char* patchPhysicalType ) = 0;
    virtual void deletePatch( const char* patchName ) = 0;
    virtual void deleteAllPatches() = 0;
    virtual void setPatchPhysicalType( const char* patchName, const char* patchPhysicalType ) = 0;
    virtual void getPatchPhysicalType( const char* patchName, CORBA::String_out patchPhysicalType ) = 0;
    virtual void getDictionary( const char* dictionaryName, CORBA::Boolean forceRead, IDictionaryEntry_out dictRoot ) = 0;
    virtual void readFile( const char* name, CORBA::String_out contents ) = 0;
    virtual void writeFile( const char* name, const char* contents ) = 0;
    virtual CORBA::Long fileModificationDate( const char* fileName ) = 0;
    virtual CORBA::Long runCase( const char* arguments ) = 0;
    virtual void killCase() = 0;
    virtual void validate() = 0;
    virtual CORBA::Boolean modified() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    ICaseServer() {};
  private:
    ICaseServer( const ICaseServer& );
    void operator=( const ICaseServer& );
};

extern CORBA::TypeCodeConst _tc_ICaseServer;

// Stub for interface ICaseServer
class ICaseServer_stub:
  virtual public ICaseServer
{
  public:
    virtual ~ICaseServer_stub();
    CORBA::Boolean managed();
    void managed( CORBA::Boolean value );
    char* caseRoot();
    char* caseName();
    ::FoamXServer::CaseServer::IApplication_ptr application();
    ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties();
    ::FoamXServer::StringList* availableTimeSteps();
    CORBA::Boolean meshDefined();
    ::FoamXServer::StringList* patchNames();

    char* getTime();
    void setTime( const char* timeName, CORBA::Long timeIndex );
    void readMesh();
    void importMesh( const char* hostName, const char* rootDir, const char* caseName );
    void getFieldValues( const char* fieldName, IGeometricField_out fieldValues );
    void addPatch( const char* patchName, const char* patchPhysicalType );
    void deletePatch( const char* patchName );
    void deleteAllPatches();
    void setPatchPhysicalType( const char* patchName, const char* patchPhysicalType );
    void getPatchPhysicalType( const char* patchName, CORBA::String_out patchPhysicalType );
    void getDictionary( const char* dictionaryName, CORBA::Boolean forceRead, IDictionaryEntry_out dictRoot );
    void readFile( const char* name, CORBA::String_out contents );
    void writeFile( const char* name, const char* contents );
    CORBA::Long fileModificationDate( const char* fileName );
    CORBA::Long runCase( const char* arguments );
    void killCase();
    void validate();
    CORBA::Boolean modified();
    void save();
    void close();

  private:
    void operator=( const ICaseServer_stub& );
};

#ifndef MICO_CONF_NO_POA

class ICaseServer_stub_clp :
  virtual public ICaseServer_stub,
  virtual public PortableServer::StubBase
{
  public:
    ICaseServer_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~ICaseServer_stub_clp ();
    CORBA::Boolean managed();
    void managed( CORBA::Boolean value );
    char* caseRoot();
    char* caseName();
    ::FoamXServer::CaseServer::IApplication_ptr application();
    ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties();
    ::FoamXServer::StringList* availableTimeSteps();
    CORBA::Boolean meshDefined();
    ::FoamXServer::StringList* patchNames();

    char* getTime();
    void setTime( const char* timeName, CORBA::Long timeIndex );
    void readMesh();
    void importMesh( const char* hostName, const char* rootDir, const char* caseName );
    void getFieldValues( const char* fieldName, IGeometricField_out fieldValues );
    void addPatch( const char* patchName, const char* patchPhysicalType );
    void deletePatch( const char* patchName );
    void deleteAllPatches();
    void setPatchPhysicalType( const char* patchName, const char* patchPhysicalType );
    void getPatchPhysicalType( const char* patchName, CORBA::String_out patchPhysicalType );
    void getDictionary( const char* dictionaryName, CORBA::Boolean forceRead, IDictionaryEntry_out dictRoot );
    void readFile( const char* name, CORBA::String_out contents );
    void writeFile( const char* name, const char* contents );
    CORBA::Long fileModificationDate( const char* fileName );
    CORBA::Long runCase( const char* arguments );
    void killCase();
    void validate();
    CORBA::Boolean modified();
    void save();
    void close();

  protected:
    ICaseServer_stub_clp ();
  private:
    void operator=( const ICaseServer_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IFoamProperties
 */

class IFoamProperties : 
  virtual public CORBA::Object
{
  public:
    virtual ~IFoamProperties();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IFoamProperties_ptr _ptr_type;
    typedef IFoamProperties_var _var_type;
    #endif

    static IFoamProperties_ptr _narrow( CORBA::Object_ptr obj );
    static IFoamProperties_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IFoamProperties_ptr _duplicate( IFoamProperties_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IFoamProperties_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::StringList* availableModules() = 0;
    virtual ::FoamXServer::StringList* rootDirectories() = 0;
    virtual ::FoamXServer::StringList* rawRootDirectories() = 0;
    virtual ::FoamXServer::StringList* foamTypes() = 0;
    virtual ::FoamXServer::StringList* geometryTypes() = 0;
    virtual ::FoamXServer::StringList* patchTypes() = 0;
    virtual ::FoamXServer::StringList* patchFieldTypes() = 0;
    virtual ::FoamXServer::ApplicationDescriptorList* applicationes() = 0;
    virtual ::FoamXServer::ApplicationDescriptorList* utilities() = 0;

    virtual void addRootDirectory( const char* rawRootDir ) = 0;
    virtual void deleteRootDirectory( const char* rawRootDir ) = 0;
    virtual void getFoamType( const char* foamTypeName, ITypeDescriptor_out typeDesc ) = 0;
    virtual void getGeometryType( const char* geometryTypeName, IGeometryDescriptor_out geometryDesc ) = 0;
    virtual void getPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc ) = 0;
    virtual void findPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc ) = 0;
    virtual void getPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc ) = 0;
    virtual void findPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc ) = 0;
    virtual void getFoamControlDict( IDictionaryEntry_out controlDict ) = 0;
    virtual void getApplication( const char* appName, IApplication_out app ) = 0;
    virtual void addApplication( const char* appName, IApplication_out app ) = 0;
    virtual void deleteApplication( const char* appName ) = 0;
    virtual void cloneApplication( const char* appNameSrc, const char* appNameDest, const char* appDestPath, IApplication_out app ) = 0;
    virtual void getUtilityControlDict( const char* utilityName, const char* rootDir, const char* caseName, IDictionaryEntry_out controlDict ) = 0;
    virtual void validate() = 0;
    virtual void saveSystemProperties() = 0;
    virtual void saveUserProperties() = 0;

  protected:
    IFoamProperties() {};
  private:
    IFoamProperties( const IFoamProperties& );
    void operator=( const IFoamProperties& );
};

extern CORBA::TypeCodeConst _tc_IFoamProperties;

// Stub for interface IFoamProperties
class IFoamProperties_stub:
  virtual public IFoamProperties
{
  public:
    virtual ~IFoamProperties_stub();
    ::FoamXServer::StringList* availableModules();
    ::FoamXServer::StringList* rootDirectories();
    ::FoamXServer::StringList* rawRootDirectories();
    ::FoamXServer::StringList* foamTypes();
    ::FoamXServer::StringList* geometryTypes();
    ::FoamXServer::StringList* patchTypes();
    ::FoamXServer::StringList* patchFieldTypes();
    ::FoamXServer::ApplicationDescriptorList* applicationes();
    ::FoamXServer::ApplicationDescriptorList* utilities();

    void addRootDirectory( const char* rawRootDir );
    void deleteRootDirectory( const char* rawRootDir );
    void getFoamType( const char* foamTypeName, ITypeDescriptor_out typeDesc );
    void getGeometryType( const char* geometryTypeName, IGeometryDescriptor_out geometryDesc );
    void getPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc );
    void findPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc );
    void getPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc );
    void findPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc );
    void getFoamControlDict( IDictionaryEntry_out controlDict );
    void getApplication( const char* appName, IApplication_out app );
    void addApplication( const char* appName, IApplication_out app );
    void deleteApplication( const char* appName );
    void cloneApplication( const char* appNameSrc, const char* appNameDest, const char* appDestPath, IApplication_out app );
    void getUtilityControlDict( const char* utilityName, const char* rootDir, const char* caseName, IDictionaryEntry_out controlDict );
    void validate();
    void saveSystemProperties();
    void saveUserProperties();

  private:
    void operator=( const IFoamProperties_stub& );
};

#ifndef MICO_CONF_NO_POA

class IFoamProperties_stub_clp :
  virtual public IFoamProperties_stub,
  virtual public PortableServer::StubBase
{
  public:
    IFoamProperties_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IFoamProperties_stub_clp ();
    ::FoamXServer::StringList* availableModules();
    ::FoamXServer::StringList* rootDirectories();
    ::FoamXServer::StringList* rawRootDirectories();
    ::FoamXServer::StringList* foamTypes();
    ::FoamXServer::StringList* geometryTypes();
    ::FoamXServer::StringList* patchTypes();
    ::FoamXServer::StringList* patchFieldTypes();
    ::FoamXServer::ApplicationDescriptorList* applicationes();
    ::FoamXServer::ApplicationDescriptorList* utilities();

    void addRootDirectory( const char* rawRootDir );
    void deleteRootDirectory( const char* rawRootDir );
    void getFoamType( const char* foamTypeName, ITypeDescriptor_out typeDesc );
    void getGeometryType( const char* geometryTypeName, IGeometryDescriptor_out geometryDesc );
    void getPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc );
    void findPatchType( const char* patchTypeName, IPatchDescriptor_out patchDesc );
    void getPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc );
    void findPatchFieldType( const char* patchFieldTypeName, ITypeDescriptor_out patchFieldDesc );
    void getFoamControlDict( IDictionaryEntry_out controlDict );
    void getApplication( const char* appName, IApplication_out app );
    void addApplication( const char* appName, IApplication_out app );
    void deleteApplication( const char* appName );
    void cloneApplication( const char* appNameSrc, const char* appNameDest, const char* appDestPath, IApplication_out app );
    void getUtilityControlDict( const char* utilityName, const char* rootDir, const char* caseName, IDictionaryEntry_out controlDict );
    void validate();
    void saveSystemProperties();
    void saveUserProperties();

  protected:
    IFoamProperties_stub_clp ();
  private:
    void operator=( const IFoamProperties_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IApplication
 */

class IApplication : 
  virtual public CORBA::Object
{
  public:
    virtual ~IApplication();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IApplication_ptr _ptr_type;
    typedef IApplication_var _var_type;
    #endif

    static IApplication_ptr _narrow( CORBA::Object_ptr obj );
    static IApplication_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IApplication_ptr _duplicate( IApplication_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IApplication_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* category() = 0;
    virtual void category( const char* value ) = 0;
    virtual ::FoamXServer::StringList* modules() = 0;
    virtual void modules( const ::FoamXServer::StringList& value ) = 0;
    virtual CORBA::Boolean systemClass() = 0;
    virtual ::FoamXServer::StringList* fields() = 0;
    virtual ::FoamXServer::StringList* patchPhysicalTypes() = 0;
    virtual ::FoamXServer::StringList* dictionaries() = 0;

    virtual void getField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void findField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void addField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void deleteField( const char* fieldName ) = 0;
    virtual void getPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void findPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void addPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void deletePatchPhysicalType( const char* patchPhysicalTypeName ) = 0;
    virtual void getDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor ) = 0;
    virtual void addDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor ) = 0;
    virtual void deleteDictionary( const char* dictName ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;

  protected:
    IApplication() {};
  private:
    IApplication( const IApplication& );
    void operator=( const IApplication& );
};

extern CORBA::TypeCodeConst _tc_IApplication;

// Stub for interface IApplication
class IApplication_stub:
  virtual public IApplication
{
  public:
    virtual ~IApplication_stub();
    char* name();
    void name( const char* value );
    char* description();
    void description( const char* value );
    char* category();
    void category( const char* value );
    ::FoamXServer::StringList* modules();
    void modules( const ::FoamXServer::StringList& value );
    CORBA::Boolean systemClass();
    ::FoamXServer::StringList* fields();
    ::FoamXServer::StringList* patchPhysicalTypes();
    ::FoamXServer::StringList* dictionaries();

    void getField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void findField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void addField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void deleteField( const char* fieldName );
    void getPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void findPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void addPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void deletePatchPhysicalType( const char* patchPhysicalTypeName );
    void getDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor );
    void addDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor );
    void deleteDictionary( const char* dictName );
    void validate();
    void save();

  private:
    void operator=( const IApplication_stub& );
};

#ifndef MICO_CONF_NO_POA

class IApplication_stub_clp :
  virtual public IApplication_stub,
  virtual public PortableServer::StubBase
{
  public:
    IApplication_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IApplication_stub_clp ();
    char* name();
    void name( const char* value );
    char* description();
    void description( const char* value );
    char* category();
    void category( const char* value );
    ::FoamXServer::StringList* modules();
    void modules( const ::FoamXServer::StringList& value );
    CORBA::Boolean systemClass();
    ::FoamXServer::StringList* fields();
    ::FoamXServer::StringList* patchPhysicalTypes();
    ::FoamXServer::StringList* dictionaries();

    void getField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void findField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void addField( const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor );
    void deleteField( const char* fieldName );
    void getPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void findPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void addPatchPhysicalType( const char* patchPhysicalTypeName, IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor );
    void deletePatchPhysicalType( const char* patchPhysicalTypeName );
    void getDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor );
    void addDictionary( const char* dictName, ITypeDescriptor_out dictTypeDescriptor );
    void deleteDictionary( const char* dictName );
    void validate();
    void save();

  protected:
    IApplication_stub_clp ();
  private:
    void operator=( const IApplication_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IGeometricFieldDescriptor
 */

class IGeometricFieldDescriptor : 
  virtual public CORBA::Object
{
  public:
    virtual ~IGeometricFieldDescriptor();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IGeometricFieldDescriptor_ptr _ptr_type;
    typedef IGeometricFieldDescriptor_var _var_type;
    #endif

    static IGeometricFieldDescriptor_ptr _narrow( CORBA::Object_ptr obj );
    static IGeometricFieldDescriptor_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IGeometricFieldDescriptor_ptr _duplicate( IGeometricFieldDescriptor_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IGeometricFieldDescriptor_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::ITypeDescriptor_ptr typeDescriptor() = 0;
    virtual ::FoamXServer::ITypeDescriptor_ptr fieldTypeDescriptor() = 0;
    virtual ::FoamXServer::CaseServer::IGeometryDescriptor_ptr geometryDescriptor() = 0;
    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* fieldTypeName() = 0;
    virtual void fieldTypeName( const char* value ) = 0;
    virtual char* geometryTypeName() = 0;
    virtual void geometryTypeName( const char* value ) = 0;
    virtual ::FoamXServer::DimensionSet dimensions() = 0;
    virtual void dimensions( const ::FoamXServer::DimensionSet& value ) = 0;

  protected:
    IGeometricFieldDescriptor() {};
  private:
    IGeometricFieldDescriptor( const IGeometricFieldDescriptor& );
    void operator=( const IGeometricFieldDescriptor& );
};

extern CORBA::TypeCodeConst _tc_IGeometricFieldDescriptor;

// Stub for interface IGeometricFieldDescriptor
class IGeometricFieldDescriptor_stub:
  virtual public IGeometricFieldDescriptor
{
  public:
    virtual ~IGeometricFieldDescriptor_stub();
    ::FoamXServer::ITypeDescriptor_ptr typeDescriptor();
    ::FoamXServer::ITypeDescriptor_ptr fieldTypeDescriptor();
    ::FoamXServer::CaseServer::IGeometryDescriptor_ptr geometryDescriptor();
    char* name();
    void name( const char* value );
    char* description();
    void description( const char* value );
    char* fieldTypeName();
    void fieldTypeName( const char* value );
    char* geometryTypeName();
    void geometryTypeName( const char* value );
    ::FoamXServer::DimensionSet dimensions();
    void dimensions( const ::FoamXServer::DimensionSet& value );

  private:
    void operator=( const IGeometricFieldDescriptor_stub& );
};

#ifndef MICO_CONF_NO_POA

class IGeometricFieldDescriptor_stub_clp :
  virtual public IGeometricFieldDescriptor_stub,
  virtual public PortableServer::StubBase
{
  public:
    IGeometricFieldDescriptor_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IGeometricFieldDescriptor_stub_clp ();
    ::FoamXServer::ITypeDescriptor_ptr typeDescriptor();
    ::FoamXServer::ITypeDescriptor_ptr fieldTypeDescriptor();
    ::FoamXServer::CaseServer::IGeometryDescriptor_ptr geometryDescriptor();
    char* name();
    void name( const char* value );
    char* description();
    void description( const char* value );
    char* fieldTypeName();
    void fieldTypeName( const char* value );
    char* geometryTypeName();
    void geometryTypeName( const char* value );
    ::FoamXServer::DimensionSet dimensions();
    void dimensions( const ::FoamXServer::DimensionSet& value );

  protected:
    IGeometricFieldDescriptor_stub_clp ();
  private:
    void operator=( const IGeometricFieldDescriptor_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IPatchDescriptor
 */

class IPatchDescriptor : 
  virtual public CORBA::Object
{
  public:
    virtual ~IPatchDescriptor();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IPatchDescriptor_ptr _ptr_type;
    typedef IPatchDescriptor_var _var_type;
    #endif

    static IPatchDescriptor_ptr _narrow( CORBA::Object_ptr obj );
    static IPatchDescriptor_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IPatchDescriptor_ptr _duplicate( IPatchDescriptor_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IPatchDescriptor_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;

  protected:
    IPatchDescriptor() {};
  private:
    IPatchDescriptor( const IPatchDescriptor& );
    void operator=( const IPatchDescriptor& );
};

extern CORBA::TypeCodeConst _tc_IPatchDescriptor;

// Stub for interface IPatchDescriptor
class IPatchDescriptor_stub:
  virtual public IPatchDescriptor
{
  public:
    virtual ~IPatchDescriptor_stub();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );

  private:
    void operator=( const IPatchDescriptor_stub& );
};

#ifndef MICO_CONF_NO_POA

class IPatchDescriptor_stub_clp :
  virtual public IPatchDescriptor_stub,
  virtual public PortableServer::StubBase
{
  public:
    IPatchDescriptor_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IPatchDescriptor_stub_clp ();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );

  protected:
    IPatchDescriptor_stub_clp ();
  private:
    void operator=( const IPatchDescriptor_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IGeometryDescriptor
 */

class IGeometryDescriptor : 
  virtual public CORBA::Object
{
  public:
    virtual ~IGeometryDescriptor();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IGeometryDescriptor_ptr _ptr_type;
    typedef IGeometryDescriptor_var _var_type;
    #endif

    static IGeometryDescriptor_ptr _narrow( CORBA::Object_ptr obj );
    static IGeometryDescriptor_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IGeometryDescriptor_ptr _duplicate( IGeometryDescriptor_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IGeometryDescriptor_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;

  protected:
    IGeometryDescriptor() {};
  private:
    IGeometryDescriptor( const IGeometryDescriptor& );
    void operator=( const IGeometryDescriptor& );
};

extern CORBA::TypeCodeConst _tc_IGeometryDescriptor;

// Stub for interface IGeometryDescriptor
class IGeometryDescriptor_stub:
  virtual public IGeometryDescriptor
{
  public:
    virtual ~IGeometryDescriptor_stub();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );

  private:
    void operator=( const IGeometryDescriptor_stub& );
};

#ifndef MICO_CONF_NO_POA

class IGeometryDescriptor_stub_clp :
  virtual public IGeometryDescriptor_stub,
  virtual public PortableServer::StubBase
{
  public:
    IGeometryDescriptor_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IGeometryDescriptor_stub_clp ();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );

  protected:
    IGeometryDescriptor_stub_clp ();
  private:
    void operator=( const IGeometryDescriptor_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IPatchPhysicalTypeDescriptor
 */

class IPatchPhysicalTypeDescriptor : 
  virtual public CORBA::Object
{
  public:
    virtual ~IPatchPhysicalTypeDescriptor();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IPatchPhysicalTypeDescriptor_ptr _ptr_type;
    typedef IPatchPhysicalTypeDescriptor_var _var_type;
    #endif

    static IPatchPhysicalTypeDescriptor_ptr _narrow( CORBA::Object_ptr obj );
    static IPatchPhysicalTypeDescriptor_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IPatchPhysicalTypeDescriptor_ptr _duplicate( IPatchPhysicalTypeDescriptor_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IPatchPhysicalTypeDescriptor_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* patchType() = 0;
    virtual void patchType( const char* value ) = 0;
    virtual char* parentType() = 0;
    virtual void parentType( const char* value ) = 0;
    virtual ::FoamXServer::StringPairList* patchFieldTypes() = 0;
    virtual void patchFieldTypes( const ::FoamXServer::StringPairList& value ) = 0;

  protected:
    IPatchPhysicalTypeDescriptor() {};
  private:
    IPatchPhysicalTypeDescriptor( const IPatchPhysicalTypeDescriptor& );
    void operator=( const IPatchPhysicalTypeDescriptor& );
};

extern CORBA::TypeCodeConst _tc_IPatchPhysicalTypeDescriptor;

// Stub for interface IPatchPhysicalTypeDescriptor
class IPatchPhysicalTypeDescriptor_stub:
  virtual public IPatchPhysicalTypeDescriptor
{
  public:
    virtual ~IPatchPhysicalTypeDescriptor_stub();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );
    char* patchType();
    void patchType( const char* value );
    char* parentType();
    void parentType( const char* value );
    ::FoamXServer::StringPairList* patchFieldTypes();
    void patchFieldTypes( const ::FoamXServer::StringPairList& value );

  private:
    void operator=( const IPatchPhysicalTypeDescriptor_stub& );
};

#ifndef MICO_CONF_NO_POA

class IPatchPhysicalTypeDescriptor_stub_clp :
  virtual public IPatchPhysicalTypeDescriptor_stub,
  virtual public PortableServer::StubBase
{
  public:
    IPatchPhysicalTypeDescriptor_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IPatchPhysicalTypeDescriptor_stub_clp ();
    char* name();
    void name( const char* value );
    char* displayName();
    void displayName( const char* value );
    char* description();
    void description( const char* value );
    char* patchType();
    void patchType( const char* value );
    char* parentType();
    void parentType( const char* value );
    ::FoamXServer::StringPairList* patchFieldTypes();
    void patchFieldTypes( const ::FoamXServer::StringPairList& value );

  protected:
    IPatchPhysicalTypeDescriptor_stub_clp ();
  private:
    void operator=( const IPatchPhysicalTypeDescriptor_stub_clp & );
};

#endif // MICO_CONF_NO_POA


/*
 * Base class and common definitions for interface IGeometricField
 */

class IGeometricField : 
  virtual public CORBA::Object
{
  public:
    virtual ~IGeometricField();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IGeometricField_ptr _ptr_type;
    typedef IGeometricField_var _var_type;
    #endif

    static IGeometricField_ptr _narrow( CORBA::Object_ptr obj );
    static IGeometricField_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IGeometricField_ptr _duplicate( IGeometricField_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IGeometricField_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* name() = 0;

    virtual void getInternalFieldValue( IDictionaryEntry_out internalFieldValue ) = 0;
    virtual void getPatchFieldParameters( const char* patchName, IDictionaryEntry_out patchFieldValue ) = 0;
    virtual CORBA::Boolean modified() = 0;

  protected:
    IGeometricField() {};
  private:
    IGeometricField( const IGeometricField& );
    void operator=( const IGeometricField& );
};

extern CORBA::TypeCodeConst _tc_IGeometricField;

// Stub for interface IGeometricField
class IGeometricField_stub:
  virtual public IGeometricField
{
  public:
    virtual ~IGeometricField_stub();
    char* name();

    void getInternalFieldValue( IDictionaryEntry_out internalFieldValue );
    void getPatchFieldParameters( const char* patchName, IDictionaryEntry_out patchFieldValue );
    CORBA::Boolean modified();

  private:
    void operator=( const IGeometricField_stub& );
};

#ifndef MICO_CONF_NO_POA

class IGeometricField_stub_clp :
  virtual public IGeometricField_stub,
  virtual public PortableServer::StubBase
{
  public:
    IGeometricField_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IGeometricField_stub_clp ();
    char* name();

    void getInternalFieldValue( IDictionaryEntry_out internalFieldValue );
    void getPatchFieldParameters( const char* patchName, IDictionaryEntry_out patchFieldValue );
    CORBA::Boolean modified();

  protected:
    IGeometricField_stub_clp ();
  private:
    void operator=( const IGeometricField_stub_clp & );
};

#endif // MICO_CONF_NO_POA

}


namespace CasePostServer
{


/*
 * Base class and common definitions for interface ICasePostServer
 */

class ICasePostServer : 
  virtual public CORBA::Object
{
  public:
    virtual ~ICasePostServer();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef ICasePostServer_ptr _ptr_type;
    typedef ICasePostServer_var _var_type;
    #endif

    static ICasePostServer_ptr _narrow( CORBA::Object_ptr obj );
    static ICasePostServer_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static ICasePostServer_ptr _duplicate( ICasePostServer_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static ICasePostServer_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual char* caseRoot() = 0;
    virtual char* caseName() = 0;
    virtual CORBA::Long nProcs() = 0;
    virtual ::FoamXServer::StringList* availableTimeSteps() = 0;

    virtual void setTime( const char* timeName, CORBA::Long timeIndex ) = 0;
    virtual void setTimeSlave() = 0;
    virtual StringList* getPatchNames() = 0;
    virtual void getPatchNamesSlave() = 0;
    virtual StringList* getFieldNames( const char* type ) = 0;
    virtual void getMeshBb( Point3_out min, Point3_out max ) = 0;
    virtual void getMeshBbSlave() = 0;
    virtual void getPatchMesh( const char* patchName, CORBA::Double creaseAngle, FloatList_out points, LongList_out edges ) = 0;
    virtual void getPatchMeshSlave() = 0;
    virtual void getCutMesh( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges ) = 0;
    virtual void getCutMeshSlave() = 0;
    virtual void getCutMeshOutline( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges ) = 0;
    virtual void getCutMeshOutlineSlave() = 0;
    virtual void getTriPatch( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values ) = 0;
    virtual void getTriPatchSlave() = 0;
    virtual void getTriPatchVec( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values ) = 0;
    virtual void getTriPatchVecSlave() = 0;
    virtual void cutPlane( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values ) = 0;
    virtual void cutPlaneSlave() = 0;
    virtual void cutPlaneVec( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values ) = 0;
    virtual void cutPlaneVecSlave() = 0;
    virtual void close() = 0;

  protected:
    ICasePostServer() {};
  private:
    ICasePostServer( const ICasePostServer& );
    void operator=( const ICasePostServer& );
};

extern CORBA::TypeCodeConst _tc_ICasePostServer;

// Stub for interface ICasePostServer
class ICasePostServer_stub:
  virtual public ICasePostServer
{
  public:
    virtual ~ICasePostServer_stub();
    char* caseRoot();
    char* caseName();
    CORBA::Long nProcs();
    ::FoamXServer::StringList* availableTimeSteps();

    void setTime( const char* timeName, CORBA::Long timeIndex );
    void setTimeSlave();
    StringList* getPatchNames();
    void getPatchNamesSlave();
    StringList* getFieldNames( const char* type );
    void getMeshBb( Point3_out min, Point3_out max );
    void getMeshBbSlave();
    void getPatchMesh( const char* patchName, CORBA::Double creaseAngle, FloatList_out points, LongList_out edges );
    void getPatchMeshSlave();
    void getCutMesh( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges );
    void getCutMeshSlave();
    void getCutMeshOutline( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges );
    void getCutMeshOutlineSlave();
    void getTriPatch( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void getTriPatchSlave();
    void getTriPatchVec( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void getTriPatchVecSlave();
    void cutPlane( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void cutPlaneSlave();
    void cutPlaneVec( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void cutPlaneVecSlave();
    void close();

  private:
    void operator=( const ICasePostServer_stub& );
};

#ifndef MICO_CONF_NO_POA

class ICasePostServer_stub_clp :
  virtual public ICasePostServer_stub,
  virtual public PortableServer::StubBase
{
  public:
    ICasePostServer_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~ICasePostServer_stub_clp ();
    char* caseRoot();
    char* caseName();
    CORBA::Long nProcs();
    ::FoamXServer::StringList* availableTimeSteps();

    void setTime( const char* timeName, CORBA::Long timeIndex );
    void setTimeSlave();
    StringList* getPatchNames();
    void getPatchNamesSlave();
    StringList* getFieldNames( const char* type );
    void getMeshBb( Point3_out min, Point3_out max );
    void getMeshBbSlave();
    void getPatchMesh( const char* patchName, CORBA::Double creaseAngle, FloatList_out points, LongList_out edges );
    void getPatchMeshSlave();
    void getCutMesh( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges );
    void getCutMeshSlave();
    void getCutMeshOutline( const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out edges );
    void getCutMeshOutlineSlave();
    void getTriPatch( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void getTriPatchSlave();
    void getTriPatchVec( const char* fieldName, const char* patchName, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void getTriPatchVecSlave();
    void cutPlane( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void cutPlaneSlave();
    void cutPlaneVec( const char* fieldName, const Point3 basePoint, const Point3 normal, FloatList_out points, LongList_out triFaces, FloatList_out values );
    void cutPlaneVecSlave();
    void close();

  protected:
    ICasePostServer_stub_clp ();
  private:
    void operator=( const ICasePostServer_stub_clp & );
};

#endif // MICO_CONF_NO_POA

}


namespace CaseBrowser
{


/*
 * Base class and common definitions for interface ICaseBrowser
 */

class ICaseBrowser : 
  virtual public CORBA::Object
{
  public:
    virtual ~ICaseBrowser();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef ICaseBrowser_ptr _ptr_type;
    typedef ICaseBrowser_var _var_type;
    #endif

    static ICaseBrowser_ptr _narrow( CORBA::Object_ptr obj );
    static ICaseBrowser_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static ICaseBrowser_ptr _duplicate( ICaseBrowser_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static ICaseBrowser_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties() = 0;
    virtual ::FoamXServer::CaseDescriptorList* cases() = 0;
    virtual ::FoamXServer::JobDescriptorList* runningJobs() = 0;
    virtual ::FoamXServer::JobDescriptorList* finishedJobs() = 0;

    virtual void getEnv( const char* envName, CORBA::String_out hostName ) = 0;
    virtual void getHostName( CORBA::String_out hostName ) = 0;
    virtual void getUserName( CORBA::String_out userName ) = 0;
    virtual CORBA::Long fileModificationDate( const char* fileName ) = 0;
    virtual void readFile( const char* fileName, CORBA::String_out contents ) = 0;
    virtual void writeFile( const char* fileName, const char* contents ) = 0;
    virtual CORBA::Long invokeUtility( const char* hostName, const char* utilityName, const StringList& arguments, const char* logName, CORBA::Boolean backGround ) = 0;
    virtual void refreshCaseList() = 0;
    virtual void addToCaseList( const char* rootDir ) = 0;
    virtual void openCase( const CaseDescriptor& caseDesc ) = 0;
    virtual void newCase( const char* rootDir, const char* caseName, const char* app ) = 0;
    virtual void importCase( const char* rootDir, const char* caseName, const char* app ) = 0;
    virtual void deleteCase( const CaseDescriptor& caseDesc ) = 0;
    virtual void cloneCase( const CaseDescriptor& caseDesc, const char* newCaseRootDir, const char* newCaseName, const char* newAppClassName, const char* timeSel ) = 0;
    virtual CORBA::Boolean getCaseServerReference( const char* rootDir, const char* caseName, CaseServer::ICaseServer_out caseObj ) = 0;
    virtual void openCasePost( const CaseDescriptor& caseDesc, CORBA::Long nProcs ) = 0;
    virtual CORBA::Boolean getCasePostServerReference( const char* rootDir, const char* caseName, CORBA::Long nProcs, CasePostServer::ICasePostServer_out caseObj ) = 0;
    virtual CORBA::Boolean caseLocked( const CaseDescriptor& caseDesc ) = 0;
    virtual void unlockCase( const char* rootDir, const char* caseName ) = 0;
    virtual void unlockCaseDescriptor( const CaseDescriptor& caseDesc ) = 0;
    virtual void addCase( const char* rootDir, const char* rawRootDir, const char* caseName, const char* app ) = 0;
    virtual void caseOpen( const char* rootDir, const char* caseName ) = 0;
    virtual CORBA::Boolean isCaseInError( const CaseDescriptor& caseDesc ) = 0;
    virtual void caseIsInError( const CaseDescriptor& caseDesc ) = 0;
    virtual void refreshJobsLists() = 0;
    virtual void checkRunningJobs() = 0;
    virtual void purgeRunningJobs() = 0;
    virtual void purgeFinishedJob( const JobID& jobID ) = 0;
    virtual void purgeFinishedJobs( CORBA::Long nDays ) = 0;
    virtual void kill( const JobID& jobID ) = 0;
    virtual void suspend( const JobID& jobID ) = 0;
    virtual void cont( const JobID& jobID ) = 0;
    virtual void end( const JobID& jobID, const char* rootDir, const char* caseName, CORBA::Boolean now ) = 0;
    virtual void setStatus( const JobID& jobID, JobStatus jobStatus ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    ICaseBrowser() {};
  private:
    ICaseBrowser( const ICaseBrowser& );
    void operator=( const ICaseBrowser& );
};

extern CORBA::TypeCodeConst _tc_ICaseBrowser;

// Stub for interface ICaseBrowser
class ICaseBrowser_stub:
  virtual public ICaseBrowser
{
  public:
    virtual ~ICaseBrowser_stub();
    ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties();
    ::FoamXServer::CaseDescriptorList* cases();
    ::FoamXServer::JobDescriptorList* runningJobs();
    ::FoamXServer::JobDescriptorList* finishedJobs();

    void getEnv( const char* envName, CORBA::String_out hostName );
    void getHostName( CORBA::String_out hostName );
    void getUserName( CORBA::String_out userName );
    CORBA::Long fileModificationDate( const char* fileName );
    void readFile( const char* fileName, CORBA::String_out contents );
    void writeFile( const char* fileName, const char* contents );
    CORBA::Long invokeUtility( const char* hostName, const char* utilityName, const StringList& arguments, const char* logName, CORBA::Boolean backGround );
    void refreshCaseList();
    void addToCaseList( const char* rootDir );
    void openCase( const CaseDescriptor& caseDesc );
    void newCase( const char* rootDir, const char* caseName, const char* app );
    void importCase( const char* rootDir, const char* caseName, const char* app );
    void deleteCase( const CaseDescriptor& caseDesc );
    void cloneCase( const CaseDescriptor& caseDesc, const char* newCaseRootDir, const char* newCaseName, const char* newAppClassName, const char* timeSel );
    CORBA::Boolean getCaseServerReference( const char* rootDir, const char* caseName, CaseServer::ICaseServer_out caseObj );
    void openCasePost( const CaseDescriptor& caseDesc, CORBA::Long nProcs );
    CORBA::Boolean getCasePostServerReference( const char* rootDir, const char* caseName, CORBA::Long nProcs, CasePostServer::ICasePostServer_out caseObj );
    CORBA::Boolean caseLocked( const CaseDescriptor& caseDesc );
    void unlockCase( const char* rootDir, const char* caseName );
    void unlockCaseDescriptor( const CaseDescriptor& caseDesc );
    void addCase( const char* rootDir, const char* rawRootDir, const char* caseName, const char* app );
    void caseOpen( const char* rootDir, const char* caseName );
    CORBA::Boolean isCaseInError( const CaseDescriptor& caseDesc );
    void caseIsInError( const CaseDescriptor& caseDesc );
    void refreshJobsLists();
    void checkRunningJobs();
    void purgeRunningJobs();
    void purgeFinishedJob( const JobID& jobID );
    void purgeFinishedJobs( CORBA::Long nDays );
    void kill( const JobID& jobID );
    void suspend( const JobID& jobID );
    void cont( const JobID& jobID );
    void end( const JobID& jobID, const char* rootDir, const char* caseName, CORBA::Boolean now );
    void setStatus( const JobID& jobID, JobStatus jobStatus );
    void validate();
    void save();
    void close();

  private:
    void operator=( const ICaseBrowser_stub& );
};

#ifndef MICO_CONF_NO_POA

class ICaseBrowser_stub_clp :
  virtual public ICaseBrowser_stub,
  virtual public PortableServer::StubBase
{
  public:
    ICaseBrowser_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~ICaseBrowser_stub_clp ();
    ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties();
    ::FoamXServer::CaseDescriptorList* cases();
    ::FoamXServer::JobDescriptorList* runningJobs();
    ::FoamXServer::JobDescriptorList* finishedJobs();

    void getEnv( const char* envName, CORBA::String_out hostName );
    void getHostName( CORBA::String_out hostName );
    void getUserName( CORBA::String_out userName );
    CORBA::Long fileModificationDate( const char* fileName );
    void readFile( const char* fileName, CORBA::String_out contents );
    void writeFile( const char* fileName, const char* contents );
    CORBA::Long invokeUtility( const char* hostName, const char* utilityName, const StringList& arguments, const char* logName, CORBA::Boolean backGround );
    void refreshCaseList();
    void addToCaseList( const char* rootDir );
    void openCase( const CaseDescriptor& caseDesc );
    void newCase( const char* rootDir, const char* caseName, const char* app );
    void importCase( const char* rootDir, const char* caseName, const char* app );
    void deleteCase( const CaseDescriptor& caseDesc );
    void cloneCase( const CaseDescriptor& caseDesc, const char* newCaseRootDir, const char* newCaseName, const char* newAppClassName, const char* timeSel );
    CORBA::Boolean getCaseServerReference( const char* rootDir, const char* caseName, CaseServer::ICaseServer_out caseObj );
    void openCasePost( const CaseDescriptor& caseDesc, CORBA::Long nProcs );
    CORBA::Boolean getCasePostServerReference( const char* rootDir, const char* caseName, CORBA::Long nProcs, CasePostServer::ICasePostServer_out caseObj );
    CORBA::Boolean caseLocked( const CaseDescriptor& caseDesc );
    void unlockCase( const char* rootDir, const char* caseName );
    void unlockCaseDescriptor( const CaseDescriptor& caseDesc );
    void addCase( const char* rootDir, const char* rawRootDir, const char* caseName, const char* app );
    void caseOpen( const char* rootDir, const char* caseName );
    CORBA::Boolean isCaseInError( const CaseDescriptor& caseDesc );
    void caseIsInError( const CaseDescriptor& caseDesc );
    void refreshJobsLists();
    void checkRunningJobs();
    void purgeRunningJobs();
    void purgeFinishedJob( const JobID& jobID );
    void purgeFinishedJobs( CORBA::Long nDays );
    void kill( const JobID& jobID );
    void suspend( const JobID& jobID );
    void cont( const JobID& jobID );
    void end( const JobID& jobID, const char* rootDir, const char* caseName, CORBA::Boolean now );
    void setStatus( const JobID& jobID, JobStatus jobStatus );
    void validate();
    void save();
    void close();

  protected:
    ICaseBrowser_stub_clp ();
  private:
    void operator=( const ICaseBrowser_stub_clp & );
};

#endif // MICO_CONF_NO_POA

}


namespace HostBrowser
{


/*
 * Base class and common definitions for interface IHostBrowser
 */

class IHostBrowser : 
  virtual public CORBA::Object
{
  public:
    virtual ~IHostBrowser();

    #ifdef HAVE_TYPEDEF_OVERLOAD
    typedef IHostBrowser_ptr _ptr_type;
    typedef IHostBrowser_var _var_type;
    #endif

    static IHostBrowser_ptr _narrow( CORBA::Object_ptr obj );
    static IHostBrowser_ptr _narrow( CORBA::AbstractBase_ptr obj );
    static IHostBrowser_ptr _duplicate( IHostBrowser_ptr _obj )
    {
      CORBA::Object::_duplicate (_obj);
      return _obj;
    }

    static IHostBrowser_ptr _nil()
    {
      return 0;
    }

    virtual void *_narrow_helper( const char *repoid );

    virtual ::FoamXServer::HostDescriptorList* hosts() = 0;

    virtual void refreshHostList() = 0;
    virtual CORBA::Boolean isHostAlive( const char* hostName ) = 0;
    virtual void hostIsAlive( const char* hostName ) = 0;
    virtual void hostIsDead( const char* hostName ) = 0;
    virtual void openCaseBrowser( const char* hostName ) = 0;
    virtual CORBA::Boolean getCaseBrowserReference( const char* hostName, CaseBrowser::ICaseBrowser_out browserObj ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    IHostBrowser() {};
  private:
    IHostBrowser( const IHostBrowser& );
    void operator=( const IHostBrowser& );
};

extern CORBA::TypeCodeConst _tc_IHostBrowser;

// Stub for interface IHostBrowser
class IHostBrowser_stub:
  virtual public IHostBrowser
{
  public:
    virtual ~IHostBrowser_stub();
    ::FoamXServer::HostDescriptorList* hosts();

    void refreshHostList();
    CORBA::Boolean isHostAlive( const char* hostName );
    void hostIsAlive( const char* hostName );
    void hostIsDead( const char* hostName );
    void openCaseBrowser( const char* hostName );
    CORBA::Boolean getCaseBrowserReference( const char* hostName, CaseBrowser::ICaseBrowser_out browserObj );
    void validate();
    void save();
    void close();

  private:
    void operator=( const IHostBrowser_stub& );
};

#ifndef MICO_CONF_NO_POA

class IHostBrowser_stub_clp :
  virtual public IHostBrowser_stub,
  virtual public PortableServer::StubBase
{
  public:
    IHostBrowser_stub_clp (PortableServer::POA_ptr, CORBA::Object_ptr);
    virtual ~IHostBrowser_stub_clp ();
    ::FoamXServer::HostDescriptorList* hosts();

    void refreshHostList();
    CORBA::Boolean isHostAlive( const char* hostName );
    void hostIsAlive( const char* hostName );
    void hostIsDead( const char* hostName );
    void openCaseBrowser( const char* hostName );
    CORBA::Boolean getCaseBrowserReference( const char* hostName, CaseBrowser::ICaseBrowser_out browserObj );
    void validate();
    void save();
    void close();

  protected:
    IHostBrowser_stub_clp ();
  private:
    void operator=( const IHostBrowser_stub_clp & );
};

#endif // MICO_CONF_NO_POA

}
}


#ifndef MICO_CONF_NO_POA



namespace POA_FoamXServer
{

class ITypeDescriptor : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~ITypeDescriptor ();
    FoamXServer::ITypeDescriptor_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static ITypeDescriptor * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::FoamXType type() = 0;
    virtual void type( ::FoamXServer::FoamXType value ) = 0;
    virtual CORBA::Boolean isPrimitiveType() = 0;
    virtual CORBA::Boolean isCompoundType() = 0;
    virtual char* path() = 0;
    virtual void path( const char* value ) = 0;
    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* comment() = 0;
    virtual void comment( const char* value ) = 0;
    virtual char* category() = 0;
    virtual void category( const char* value ) = 0;
    virtual char* helpURL() = 0;
    virtual void helpURL( const char* value ) = 0;
    virtual char* iconURL() = 0;
    virtual void iconURL( const char* value ) = 0;
    virtual CORBA::Boolean optional() = 0;
    virtual void optional( CORBA::Boolean value ) = 0;
    virtual CORBA::Boolean visible() = 0;
    virtual void visible( CORBA::Boolean value ) = 0;
    virtual CORBA::Boolean editable() = 0;
    virtual void editable( CORBA::Boolean value ) = 0;
    virtual ::FoamXServer::FoamXAny* minValue() = 0;
    virtual void minValue( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual ::FoamXServer::FoamXAny* maxValue() = 0;
    virtual void maxValue( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual char* lookupDict() = 0;
    virtual void lookupDict( const char* value ) = 0;
    virtual ::FoamXServer::FoamXAnyList* valueList() = 0;
    virtual void valueList( const ::FoamXServer::FoamXAnyList& value ) = 0;
    virtual char* dictionaryPath() = 0;
    virtual void dictionaryPath( const char* value ) = 0;
    virtual CORBA::Long numElements() = 0;
    virtual void numElements( CORBA::Long value ) = 0;
    virtual ::FoamXServer::StringList* elementLabels() = 0;
    virtual void elementLabels( const ::FoamXServer::StringList& value ) = 0;
    virtual ::FoamXServer::TypeDescriptorList* subTypes() = 0;
    virtual ::FoamXServer::ITypeDescriptor_ptr elementType() = 0;

    virtual CORBA::Boolean hasDefaultValue() = 0;
    virtual void getDefaultValue( ::FoamXServer::IDictionaryEntry_out defaultValue ) = 0;
    virtual void addSubType( ::FoamXServer::FoamXType type, ::FoamXServer::ITypeDescriptor_out subEntry ) = 0;
    virtual void removeSubType( ::FoamXServer::ITypeDescriptor_ptr subEntry ) = 0;
    virtual void validate() = 0;

  protected:
    ITypeDescriptor () {};

  private:
    ITypeDescriptor (const ITypeDescriptor &);
    void operator= (const ITypeDescriptor &);
};

class IDictionaryEntry : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IDictionaryEntry ();
    FoamXServer::IDictionaryEntry_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IDictionaryEntry * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::ITypeDescriptor_ptr typeDescriptor() = 0;
    virtual ::FoamXServer::FoamXAny* value() = 0;
    virtual void value( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual ::FoamXServer::DictionaryEntryList* subElements() = 0;
    virtual CORBA::Long selection() = 0;
    virtual void selection( CORBA::Long value ) = 0;

    virtual void setValue( const ::FoamXServer::FoamXAny& value ) = 0;
    virtual CORBA::Long nSubElements() = 0;
    virtual CORBA::Boolean packedList() = 0;
    virtual void addElement( ::FoamXServer::IDictionaryEntry_out subEntry ) = 0;
    virtual void removeElement( ::FoamXServer::IDictionaryEntry_ptr subEntry ) = 0;
    virtual void validate() = 0;
    virtual CORBA::Boolean modified() = 0;
    virtual void save() = 0;

  protected:
    IDictionaryEntry () {};

  private:
    IDictionaryEntry (const IDictionaryEntry &);
    void operator= (const IDictionaryEntry &);
};



namespace CaseServer
{

class ICaseServer : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~ICaseServer ();
    FoamXServer::CaseServer::ICaseServer_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static ICaseServer * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual CORBA::Boolean managed() = 0;
    virtual void managed( CORBA::Boolean value ) = 0;
    virtual char* caseRoot() = 0;
    virtual char* caseName() = 0;
    virtual ::FoamXServer::CaseServer::IApplication_ptr application() = 0;
    virtual ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties() = 0;
    virtual ::FoamXServer::StringList* availableTimeSteps() = 0;
    virtual CORBA::Boolean meshDefined() = 0;
    virtual ::FoamXServer::StringList* patchNames() = 0;

    virtual char* getTime() = 0;
    virtual void setTime( const char* timeName, CORBA::Long timeIndex ) = 0;
    virtual void readMesh() = 0;
    virtual void importMesh( const char* hostName, const char* rootDir, const char* caseName ) = 0;
    virtual void getFieldValues( const char* fieldName, ::FoamXServer::CaseServer::IGeometricField_out fieldValues ) = 0;
    virtual void addPatch( const char* patchName, const char* patchPhysicalType ) = 0;
    virtual void deletePatch( const char* patchName ) = 0;
    virtual void deleteAllPatches() = 0;
    virtual void setPatchPhysicalType( const char* patchName, const char* patchPhysicalType ) = 0;
    virtual void getPatchPhysicalType( const char* patchName, CORBA::String_out patchPhysicalType ) = 0;
    virtual void getDictionary( const char* dictionaryName, CORBA::Boolean forceRead, ::FoamXServer::IDictionaryEntry_out dictRoot ) = 0;
    virtual void readFile( const char* name, CORBA::String_out contents ) = 0;
    virtual void writeFile( const char* name, const char* contents ) = 0;
    virtual CORBA::Long fileModificationDate( const char* fileName ) = 0;
    virtual CORBA::Long runCase( const char* arguments ) = 0;
    virtual void killCase() = 0;
    virtual void validate() = 0;
    virtual CORBA::Boolean modified() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    ICaseServer () {};

  private:
    ICaseServer (const ICaseServer &);
    void operator= (const ICaseServer &);
};

class IFoamProperties : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IFoamProperties ();
    FoamXServer::CaseServer::IFoamProperties_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IFoamProperties * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::StringList* availableModules() = 0;
    virtual ::FoamXServer::StringList* rootDirectories() = 0;
    virtual ::FoamXServer::StringList* rawRootDirectories() = 0;
    virtual ::FoamXServer::StringList* foamTypes() = 0;
    virtual ::FoamXServer::StringList* geometryTypes() = 0;
    virtual ::FoamXServer::StringList* patchTypes() = 0;
    virtual ::FoamXServer::StringList* patchFieldTypes() = 0;
    virtual ::FoamXServer::ApplicationDescriptorList* applicationes() = 0;
    virtual ::FoamXServer::ApplicationDescriptorList* utilities() = 0;

    virtual void addRootDirectory( const char* rawRootDir ) = 0;
    virtual void deleteRootDirectory( const char* rawRootDir ) = 0;
    virtual void getFoamType( const char* foamTypeName, ::FoamXServer::ITypeDescriptor_out typeDesc ) = 0;
    virtual void getGeometryType( const char* geometryTypeName, ::FoamXServer::CaseServer::IGeometryDescriptor_out geometryDesc ) = 0;
    virtual void getPatchType( const char* patchTypeName, ::FoamXServer::CaseServer::IPatchDescriptor_out patchDesc ) = 0;
    virtual void findPatchType( const char* patchTypeName, ::FoamXServer::CaseServer::IPatchDescriptor_out patchDesc ) = 0;
    virtual void getPatchFieldType( const char* patchFieldTypeName, ::FoamXServer::ITypeDescriptor_out patchFieldDesc ) = 0;
    virtual void findPatchFieldType( const char* patchFieldTypeName, ::FoamXServer::ITypeDescriptor_out patchFieldDesc ) = 0;
    virtual void getFoamControlDict( ::FoamXServer::IDictionaryEntry_out controlDict ) = 0;
    virtual void getApplication( const char* appName, ::FoamXServer::CaseServer::IApplication_out app ) = 0;
    virtual void addApplication( const char* appName, ::FoamXServer::CaseServer::IApplication_out app ) = 0;
    virtual void deleteApplication( const char* appName ) = 0;
    virtual void cloneApplication( const char* appNameSrc, const char* appNameDest, const char* appDestPath, ::FoamXServer::CaseServer::IApplication_out app ) = 0;
    virtual void getUtilityControlDict( const char* utilityName, const char* rootDir, const char* caseName, ::FoamXServer::IDictionaryEntry_out controlDict ) = 0;
    virtual void validate() = 0;
    virtual void saveSystemProperties() = 0;
    virtual void saveUserProperties() = 0;

  protected:
    IFoamProperties () {};

  private:
    IFoamProperties (const IFoamProperties &);
    void operator= (const IFoamProperties &);
};

class IApplication : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IApplication ();
    FoamXServer::CaseServer::IApplication_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IApplication * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* category() = 0;
    virtual void category( const char* value ) = 0;
    virtual ::FoamXServer::StringList* modules() = 0;
    virtual void modules( const ::FoamXServer::StringList& value ) = 0;
    virtual CORBA::Boolean systemClass() = 0;
    virtual ::FoamXServer::StringList* fields() = 0;
    virtual ::FoamXServer::StringList* patchPhysicalTypes() = 0;
    virtual ::FoamXServer::StringList* dictionaries() = 0;

    virtual void getField( const char* fieldName, ::FoamXServer::CaseServer::IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void findField( const char* fieldName, ::FoamXServer::CaseServer::IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void addField( const char* fieldName, ::FoamXServer::CaseServer::IGeometricFieldDescriptor_out fieldDescriptor ) = 0;
    virtual void deleteField( const char* fieldName ) = 0;
    virtual void getPatchPhysicalType( const char* patchPhysicalTypeName, ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void findPatchPhysicalType( const char* patchPhysicalTypeName, ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void addPatchPhysicalType( const char* patchPhysicalTypeName, ::FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out patchPhysicalTypeDescriptor ) = 0;
    virtual void deletePatchPhysicalType( const char* patchPhysicalTypeName ) = 0;
    virtual void getDictionary( const char* dictName, ::FoamXServer::ITypeDescriptor_out dictTypeDescriptor ) = 0;
    virtual void addDictionary( const char* dictName, ::FoamXServer::ITypeDescriptor_out dictTypeDescriptor ) = 0;
    virtual void deleteDictionary( const char* dictName ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;

  protected:
    IApplication () {};

  private:
    IApplication (const IApplication &);
    void operator= (const IApplication &);
};

class IGeometricFieldDescriptor : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IGeometricFieldDescriptor ();
    FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IGeometricFieldDescriptor * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::ITypeDescriptor_ptr typeDescriptor() = 0;
    virtual ::FoamXServer::ITypeDescriptor_ptr fieldTypeDescriptor() = 0;
    virtual ::FoamXServer::CaseServer::IGeometryDescriptor_ptr geometryDescriptor() = 0;
    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* fieldTypeName() = 0;
    virtual void fieldTypeName( const char* value ) = 0;
    virtual char* geometryTypeName() = 0;
    virtual void geometryTypeName( const char* value ) = 0;
    virtual ::FoamXServer::DimensionSet dimensions() = 0;
    virtual void dimensions( const ::FoamXServer::DimensionSet& value ) = 0;

  protected:
    IGeometricFieldDescriptor () {};

  private:
    IGeometricFieldDescriptor (const IGeometricFieldDescriptor &);
    void operator= (const IGeometricFieldDescriptor &);
};

class IPatchDescriptor : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IPatchDescriptor ();
    FoamXServer::CaseServer::IPatchDescriptor_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IPatchDescriptor * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;

  protected:
    IPatchDescriptor () {};

  private:
    IPatchDescriptor (const IPatchDescriptor &);
    void operator= (const IPatchDescriptor &);
};

class IGeometryDescriptor : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IGeometryDescriptor ();
    FoamXServer::CaseServer::IGeometryDescriptor_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IGeometryDescriptor * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;

  protected:
    IGeometryDescriptor () {};

  private:
    IGeometryDescriptor (const IGeometryDescriptor &);
    void operator= (const IGeometryDescriptor &);
};

class IPatchPhysicalTypeDescriptor : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IPatchPhysicalTypeDescriptor ();
    FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IPatchPhysicalTypeDescriptor * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* name() = 0;
    virtual void name( const char* value ) = 0;
    virtual char* displayName() = 0;
    virtual void displayName( const char* value ) = 0;
    virtual char* description() = 0;
    virtual void description( const char* value ) = 0;
    virtual char* patchType() = 0;
    virtual void patchType( const char* value ) = 0;
    virtual char* parentType() = 0;
    virtual void parentType( const char* value ) = 0;
    virtual ::FoamXServer::StringPairList* patchFieldTypes() = 0;
    virtual void patchFieldTypes( const ::FoamXServer::StringPairList& value ) = 0;

  protected:
    IPatchPhysicalTypeDescriptor () {};

  private:
    IPatchPhysicalTypeDescriptor (const IPatchPhysicalTypeDescriptor &);
    void operator= (const IPatchPhysicalTypeDescriptor &);
};

class IGeometricField : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IGeometricField ();
    FoamXServer::CaseServer::IGeometricField_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IGeometricField * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* name() = 0;

    virtual void getInternalFieldValue( ::FoamXServer::IDictionaryEntry_out internalFieldValue ) = 0;
    virtual void getPatchFieldParameters( const char* patchName, ::FoamXServer::IDictionaryEntry_out patchFieldValue ) = 0;
    virtual CORBA::Boolean modified() = 0;

  protected:
    IGeometricField () {};

  private:
    IGeometricField (const IGeometricField &);
    void operator= (const IGeometricField &);
};

}


namespace CasePostServer
{

class ICasePostServer : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~ICasePostServer ();
    FoamXServer::CasePostServer::ICasePostServer_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static ICasePostServer * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual char* caseRoot() = 0;
    virtual char* caseName() = 0;
    virtual CORBA::Long nProcs() = 0;
    virtual ::FoamXServer::StringList* availableTimeSteps() = 0;

    virtual void setTime( const char* timeName, CORBA::Long timeIndex ) = 0;
    virtual void setTimeSlave() = 0;
    virtual ::FoamXServer::StringList* getPatchNames() = 0;
    virtual void getPatchNamesSlave() = 0;
    virtual ::FoamXServer::StringList* getFieldNames( const char* type ) = 0;
    virtual void getMeshBb( ::FoamXServer::Point3_out min, ::FoamXServer::Point3_out max ) = 0;
    virtual void getMeshBbSlave() = 0;
    virtual void getPatchMesh( const char* patchName, CORBA::Double creaseAngle, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out edges ) = 0;
    virtual void getPatchMeshSlave() = 0;
    virtual void getCutMesh( const ::FoamXServer::Point3 basePoint, const ::FoamXServer::Point3 normal, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out edges ) = 0;
    virtual void getCutMeshSlave() = 0;
    virtual void getCutMeshOutline( const ::FoamXServer::Point3 basePoint, const ::FoamXServer::Point3 normal, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out edges ) = 0;
    virtual void getCutMeshOutlineSlave() = 0;
    virtual void getTriPatch( const char* fieldName, const char* patchName, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out triFaces, ::FoamXServer::FloatList_out values ) = 0;
    virtual void getTriPatchSlave() = 0;
    virtual void getTriPatchVec( const char* fieldName, const char* patchName, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out triFaces, ::FoamXServer::FloatList_out values ) = 0;
    virtual void getTriPatchVecSlave() = 0;
    virtual void cutPlane( const char* fieldName, const ::FoamXServer::Point3 basePoint, const ::FoamXServer::Point3 normal, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out triFaces, ::FoamXServer::FloatList_out values ) = 0;
    virtual void cutPlaneSlave() = 0;
    virtual void cutPlaneVec( const char* fieldName, const ::FoamXServer::Point3 basePoint, const ::FoamXServer::Point3 normal, ::FoamXServer::FloatList_out points, ::FoamXServer::LongList_out triFaces, ::FoamXServer::FloatList_out values ) = 0;
    virtual void cutPlaneVecSlave() = 0;
    virtual void close() = 0;

  protected:
    ICasePostServer () {};

  private:
    ICasePostServer (const ICasePostServer &);
    void operator= (const ICasePostServer &);
};

}


namespace CaseBrowser
{

class ICaseBrowser : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~ICaseBrowser ();
    FoamXServer::CaseBrowser::ICaseBrowser_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static ICaseBrowser * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::CaseServer::IFoamProperties_ptr foamProperties() = 0;
    virtual ::FoamXServer::CaseDescriptorList* cases() = 0;
    virtual ::FoamXServer::JobDescriptorList* runningJobs() = 0;
    virtual ::FoamXServer::JobDescriptorList* finishedJobs() = 0;

    virtual void getEnv( const char* envName, CORBA::String_out hostName ) = 0;
    virtual void getHostName( CORBA::String_out hostName ) = 0;
    virtual void getUserName( CORBA::String_out userName ) = 0;
    virtual CORBA::Long fileModificationDate( const char* fileName ) = 0;
    virtual void readFile( const char* fileName, CORBA::String_out contents ) = 0;
    virtual void writeFile( const char* fileName, const char* contents ) = 0;
    virtual CORBA::Long invokeUtility( const char* hostName, const char* utilityName, const ::FoamXServer::StringList& arguments, const char* logName, CORBA::Boolean backGround ) = 0;
    virtual void refreshCaseList() = 0;
    virtual void addToCaseList( const char* rootDir ) = 0;
    virtual void openCase( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void newCase( const char* rootDir, const char* caseName, const char* app ) = 0;
    virtual void importCase( const char* rootDir, const char* caseName, const char* app ) = 0;
    virtual void deleteCase( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void cloneCase( const ::FoamXServer::CaseDescriptor& caseDesc, const char* newCaseRootDir, const char* newCaseName, const char* newAppClassName, const char* timeSel ) = 0;
    virtual CORBA::Boolean getCaseServerReference( const char* rootDir, const char* caseName, ::FoamXServer::CaseServer::ICaseServer_out caseObj ) = 0;
    virtual void openCasePost( const ::FoamXServer::CaseDescriptor& caseDesc, CORBA::Long nProcs ) = 0;
    virtual CORBA::Boolean getCasePostServerReference( const char* rootDir, const char* caseName, CORBA::Long nProcs, ::FoamXServer::CasePostServer::ICasePostServer_out caseObj ) = 0;
    virtual CORBA::Boolean caseLocked( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void unlockCase( const char* rootDir, const char* caseName ) = 0;
    virtual void unlockCaseDescriptor( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void addCase( const char* rootDir, const char* rawRootDir, const char* caseName, const char* app ) = 0;
    virtual void caseOpen( const char* rootDir, const char* caseName ) = 0;
    virtual CORBA::Boolean isCaseInError( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void caseIsInError( const ::FoamXServer::CaseDescriptor& caseDesc ) = 0;
    virtual void refreshJobsLists() = 0;
    virtual void checkRunningJobs() = 0;
    virtual void purgeRunningJobs() = 0;
    virtual void purgeFinishedJob( const ::FoamXServer::JobID& jobID ) = 0;
    virtual void purgeFinishedJobs( CORBA::Long nDays ) = 0;
    virtual void kill( const ::FoamXServer::JobID& jobID ) = 0;
    virtual void suspend( const ::FoamXServer::JobID& jobID ) = 0;
    virtual void cont( const ::FoamXServer::JobID& jobID ) = 0;
    virtual void end( const ::FoamXServer::JobID& jobID, const char* rootDir, const char* caseName, CORBA::Boolean now ) = 0;
    virtual void setStatus( const ::FoamXServer::JobID& jobID, ::FoamXServer::JobStatus jobStatus ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    ICaseBrowser () {};

  private:
    ICaseBrowser (const ICaseBrowser &);
    void operator= (const ICaseBrowser &);
};

}


namespace HostBrowser
{

class IHostBrowser : virtual public PortableServer::StaticImplementation
{
  public:
    virtual ~IHostBrowser ();
    FoamXServer::HostBrowser::IHostBrowser_ptr _this ();
    bool dispatch (CORBA::StaticServerRequest_ptr);
    virtual void invoke (CORBA::StaticServerRequest_ptr);
    virtual CORBA::Boolean _is_a (const char *);
    virtual CORBA::InterfaceDef_ptr _get_interface ();
    virtual CORBA::RepositoryId _primary_interface (const PortableServer::ObjectId &, PortableServer::POA_ptr);

    virtual void * _narrow_helper (const char *);
    static IHostBrowser * _narrow (PortableServer::Servant);
    virtual CORBA::Object_ptr _make_stub (PortableServer::POA_ptr, CORBA::Object_ptr);

    virtual ::FoamXServer::HostDescriptorList* hosts() = 0;

    virtual void refreshHostList() = 0;
    virtual CORBA::Boolean isHostAlive( const char* hostName ) = 0;
    virtual void hostIsAlive( const char* hostName ) = 0;
    virtual void hostIsDead( const char* hostName ) = 0;
    virtual void openCaseBrowser( const char* hostName ) = 0;
    virtual CORBA::Boolean getCaseBrowserReference( const char* hostName, ::FoamXServer::CaseBrowser::ICaseBrowser_out browserObj ) = 0;
    virtual void validate() = 0;
    virtual void save() = 0;
    virtual void close() = 0;

  protected:
    IHostBrowser () {};

  private:
    IHostBrowser (const IHostBrowser &);
    void operator= (const IHostBrowser &);
};

}
}


#endif // MICO_CONF_NO_POA

void operator<<=( CORBA::Any &a, const ::FoamXServer::FoamXType &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::FoamXType &e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXType;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::FoamXAny &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::FoamXAny *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::FoamXAny &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::FoamXAny *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXAny;

void operator<<=( CORBA::Any &a, const ::FoamXServer::ErrorCode &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::ErrorCode &e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_ErrorCode;

void operator<<=( CORBA::Any &a, const ::FoamXServer::FoamXError &e );
void operator<<=( CORBA::Any &a, ::FoamXServer::FoamXError *e );
CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::FoamXError &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, const ::FoamXServer::FoamXError *&e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXError;

void operator<<=( CORBA::Any &a, const ::FoamXServer::FoamXIOError &e );
void operator<<=( CORBA::Any &a, ::FoamXServer::FoamXIOError *e );
CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::FoamXIOError &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, const ::FoamXServer::FoamXIOError *&e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXIOError;

void operator<<=( CORBA::Any &a, const ::FoamXServer::ValidationError &e );
void operator<<=( CORBA::Any &a, ::FoamXServer::ValidationError *e );
CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::ValidationError &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, const ::FoamXServer::ValidationError *&e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_ValidationError;

void operator<<=( CORBA::Any &a, const FoamXServer::ITypeDescriptor_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::ITypeDescriptor_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::ITypeDescriptor_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_ITypeDescriptor;

void operator<<=( CORBA::Any &a, const FoamXServer::IDictionaryEntry_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::IDictionaryEntry_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::IDictionaryEntry_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_IDictionaryEntry;

void operator<<=( CORBA::Any &a, const ::FoamXServer::JobStatus &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::JobStatus &e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobStatus;

void operator<<=( CORBA::Any &a, const ::FoamXServer::MessageType &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::MessageType &e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_MessageType;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::Point3_forany &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::Point3_forany &_s );
void operator<<=( CORBA::Any &_a, const ::FoamXServer::StringPair &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::StringPair *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::StringPair &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::StringPair *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_StringPair;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::HostDescriptor &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::HostDescriptor *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::HostDescriptor &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::HostDescriptor *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_HostDescriptor;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::ApplicationDescriptor &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::ApplicationDescriptor *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::ApplicationDescriptor &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::ApplicationDescriptor *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_ApplicationDescriptor;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::CaseDescriptor &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::CaseDescriptor *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::CaseDescriptor &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::CaseDescriptor *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseDescriptor;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::JobID &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::JobID *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::JobID &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::JobID *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobID;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::JobDescriptor &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::JobDescriptor *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::JobDescriptor &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::JobDescriptor *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_JobDescriptor;

void operator<<=( CORBA::Any &_a, const ::FoamXServer::DimensionSet &_s );
void operator<<=( CORBA::Any &_a, ::FoamXServer::DimensionSet *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, ::FoamXServer::DimensionSet &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const ::FoamXServer::DimensionSet *&_s );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_DimensionSet;

void operator<<=( CORBA::Any &a, const ::FoamXServer::FoamXSYSError &e );
void operator<<=( CORBA::Any &a, ::FoamXServer::FoamXSYSError *e );
CORBA::Boolean operator>>=( const CORBA::Any &a, ::FoamXServer::FoamXSYSError &e );

CORBA::Boolean operator>>=( const CORBA::Any &a, const ::FoamXServer::FoamXSYSError *&e );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_FoamXSYSError;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::ICaseServer_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::ICaseServer_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::ICaseServer_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_ICaseServer;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IFoamProperties_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IFoamProperties_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IFoamProperties_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IFoamProperties;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IApplication_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IApplication_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IApplication_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IApplication;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometricFieldDescriptor;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IPatchDescriptor_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IPatchDescriptor_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IPatchDescriptor_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IPatchDescriptor;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IGeometryDescriptor_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IGeometryDescriptor_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IGeometryDescriptor_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometryDescriptor;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IPatchPhysicalTypeDescriptor;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseServer::IGeometricField_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseServer::IGeometricField_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseServer::IGeometricField_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseServer_IGeometricField;

void operator<<=( CORBA::Any &a, const FoamXServer::CasePostServer::ICasePostServer_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CasePostServer::ICasePostServer_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CasePostServer::ICasePostServer_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CasePostServer_ICasePostServer;

void operator<<=( CORBA::Any &a, const FoamXServer::CaseBrowser::ICaseBrowser_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::CaseBrowser::ICaseBrowser_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::CaseBrowser::ICaseBrowser_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_CaseBrowser_ICaseBrowser;

void operator<<=( CORBA::Any &a, const FoamXServer::HostBrowser::IHostBrowser_ptr obj );
void operator<<=( CORBA::Any &a, FoamXServer::HostBrowser::IHostBrowser_ptr* obj_ptr );
CORBA::Boolean operator>>=( const CORBA::Any &a, FoamXServer::HostBrowser::IHostBrowser_ptr &obj );

extern CORBA::StaticTypeInfo *_marshaller_FoamXServer_HostBrowser_IHostBrowser;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::FoamXAny,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_FoamXAny;

void operator<<=( CORBA::Any &_a, const StringSequenceTmpl<CORBA::String_var> &_s );
void operator<<=( CORBA::Any &_a, StringSequenceTmpl<CORBA::String_var> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, StringSequenceTmpl<CORBA::String_var> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const StringSequenceTmpl<CORBA::String_var> *&_s );

void operator<<=( CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> &_s );
void operator<<=( CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::ITypeDescriptor_var,FoamXServer::ITypeDescriptor_ptr> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_ITypeDescriptor;

void operator<<=( CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> &_s );
void operator<<=( CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const IfaceSequenceTmpl< FoamXServer::IDictionaryEntry_var,FoamXServer::IDictionaryEntry_ptr> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_IDictionaryEntry;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< CORBA::Double,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< CORBA::Double,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< CORBA::Double,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< CORBA::Double,MICO_TID_DEF> *&_s );

void operator<<=( CORBA::Any &_a, const SequenceTmpl< CORBA::Float,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< CORBA::Float,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< CORBA::Float,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< CORBA::Float,MICO_TID_DEF> *&_s );

void operator<<=( CORBA::Any &_a, const SequenceTmpl< CORBA::Long,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< CORBA::Long,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< CORBA::Long,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< CORBA::Long,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__a3_float;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::StringPair,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_StringPair;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::HostDescriptor,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_HostDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::ApplicationDescriptor,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_ApplicationDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::CaseDescriptor,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_CaseDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobID,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_JobID;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::JobDescriptor,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_JobDescriptor;

void operator<<=( CORBA::Any &_a, const SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> &_s );
void operator<<=( CORBA::Any &_a, SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> *_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> &_s );
CORBA::Boolean operator>>=( const CORBA::Any &_a, const SequenceTmpl< FoamXServer::DimensionSet,MICO_TID_DEF> *&_s );

extern CORBA::StaticTypeInfo *_marshaller__seq_FoamXServer_DimensionSet;

#endif
