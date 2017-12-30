/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::token

Description
    A token holds an item read from Istream.

SourceFiles
    tokenI.H
    token.C
    tokenIO.C

\*---------------------------------------------------------------------------*/

#ifndef token_H
#define token_H

#include "label.H"
#include "uLabel.H"
#include "scalar.H"
#include "word.H"
#include "InfoProxy.H"
#include "refCount.H"
#include "typeInfo.H"

#define NoHashTableC
#include "runTimeSelectionTables.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class token;

Istream& operator>>(Istream& is, token& t);
Ostream& operator<<(Ostream& os, const token& t);


/*---------------------------------------------------------------------------*\
                            Class token Declaration
\*---------------------------------------------------------------------------*/

class token
{
public:

    //- Enumeration defining the types of token.
    //  Since these values are also used to tag content in Pstream,
    //  the maximum number of types is limited to 30.
    enum tokenType
    {
        UNDEFINED,        //!< An undefined token-type

        // Fundamental types
        FLAG,             //!< stream flag (1-byte bitmask)
        PUNCTUATION,      //!< single character punctuation
        LABEL,            //!< label (integer) type
        FLOAT_SCALAR,     //!< float (single-precision) type
        DOUBLE_SCALAR,    //!< double (double-precision) type

        // Pointer types
        WORD,             //!< Contents represent a Foam::word
        STRING,           //!< Contents represent a Foam::string
        VARIABLE,         //!< Contents are a Foam::string representing a
                          //!< dictionary \c $variable expansion
        VERBATIMSTRING,   //!< Contents are a Foam::string representing verbatim
                          //!< content
        COMPOUND,         //!< Compound type such as List\<label\> etc.

        ERROR             //!< A token error encountered
    };


    //- Stream or output control flags (1-byte width)
    enum flagType
    {
        NO_FLAG = 0,   //!< No flags
        ASCII   = 1,   //!< ASCII-mode stream
        BINARY  = 2,   //!< BINARY-mode stream
    };


    //- Standard punctuation tokens (a character)
    enum punctuationToken
    {
        NULL_TOKEN     = '\0',
        SPACE          = ' ',
        TAB            = '\t',
        NL             = '\n',

        END_STATEMENT  = ';', //!< End entry [#isseparator]
        BEGIN_LIST     = '(', //!< Begin list [#isseparator]
        END_LIST       = ')', //!< End list [#isseparator]
        BEGIN_SQR      = '[', //!< Begin dimensions [#isseparator]
        END_SQR        = ']', //!< End dimensions [#isseparator]
        BEGIN_BLOCK    = '{', //!< Begin block [#isseparator]
        END_BLOCK      = '}', //!< End block [#isseparator]
        COLON          = ':', //!< Colon [#isseparator]
        COMMA          = ',', //!< Comma [#isseparator]
        HASH           = '#',
        ATSYM          = '@',

        BEGIN_STRING   = '"',
        END_STRING     = BEGIN_STRING,

        ASSIGN         = '=', //!< Assigment/equals [#isseparator]
        ADD            = '+', //!< Addition [#isseparator]
        SUBTRACT       = '-', //!< Substract or start of negative number
        MULTIPLY       = '*', //!< Multiply [#isseparator]
        DIVIDE         = '/'  //!< Divide [#isseparator]
    };


    //- Abstract base class for complex tokens
    class compound
    :
        public refCount
    {
        // Private data

            bool empty_;

    public:

        //- Runtime type information
        TypeName("compound");

        //- Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            compound,
            Istream,
            (Istream& is),
            (is)
        );


        // Constructors

            //- Construct null
            compound()
            :
                empty_(false)
            {}

            //- No default copy construct
            compound(const compound&) = delete;


        // Selectors

            //- Select null constructed
            static autoPtr<compound> New(const word& type, Istream& is);


        //- Destructor
        virtual ~compound();


        // Member Functions

            //- Return true if name is a known compound type
            static bool isCompound(const word& name);

            bool empty() const
            {
                return empty_;
            }

            bool& empty()
            {
                return empty_;
            }

            virtual label size() const = 0;

            virtual void write(Ostream& os) const = 0;


        // Operators

            //- No default assign operator
            compound& operator=(const compound&) = delete;

            //- Output operator
            friend Ostream& operator<<(Ostream& os, const compound& ct);
    };


    //- A templated class for holding compound tokens
    template<class T>
    class Compound
    :
        public token::compound,
        public T
    {
    public:

        //- Runtime type information
        TypeName("Compound<T>");

        Compound(Istream& is)
        :
            T(is)
        {}

        label size() const
        {
            return T::size();
        }

        void write(Ostream& os) const
        {
            operator<<(os, static_cast<const T&>(*this));
        }
    };


    //- Static undefined token
    static const token undefinedToken;


private:

    //- A %union of token types
    union content
    {
        // Fundamental values. Largest first for any {} initialization.
        int64_t int64Val;
        int32_t int32Val;

        int flagVal;   // bitmask - stored as int, not enum
        punctuationToken punctuationVal;
        label labelVal;
        floatScalar floatVal;
        doubleScalar doubleVal;

        // Pointers
        word* wordPtr;
        string* stringPtr;
        mutable compound* compoundPtr;
    };


    // Private data

        //- The data content (as a union).
        //  For memory alignment this should appear as the first member.
        content data_;

        //- The token type
        tokenType type_;

        //- Line number in the file the token was read from
        label lineNumber_;


    // Private Member Functions

        //- Set as UNDEFINED and zero the union content without any checking
        inline void setUndefined();

        //- Clear any allocated storage (word or string) and set to UNDEFINED
        inline void clear();

        // Parse error, expected 'expected', found ...
        void parseError(const char* expected) const;


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        inline token();

        //- Copy construct
        inline token(const token& t);

        //- Move construct. The original token is left as UNDEFINED.
        inline token(token&& t);

        //- Construct punctuation character token
        inline explicit token(punctuationToken p);

        //- Construct label token
        inline explicit token(const label val);

        //- Construct floatScalar token
        inline explicit token(const floatScalar val);

        //- Construct doubleScalar token
        inline explicit token(const doubleScalar val);

        //- Construct word token by copying word contents
        inline explicit token(const word& w);

        //- Construct string token by copying string contents
        inline explicit token(const string& str);


        //- Construct punctuation character token
        inline token(punctuationToken p, const label lineNumber);

        //- Construct label token
        inline token(const label val, const label lineNumber);

        //- Construct floatScalar token
        inline token(const floatScalar val, const label lineNumber);

        //- Construct doubleScalar token
        inline token(const doubleScalar val, const label lineNumber);

        //- Construct word token by copying word contents
        inline token(const word& w, const label lineNumber);

        //- Construct string token by copying string contents
        inline token(const string& str, const label lineNumber);

        //- Construct from Istream
        token(Istream& is);


    //- Destructor
    inline ~token();


    // Static Member Functions

        //- Create a token with stream flags, no sanity check
        //
        //  \param bitmask the flags to set
        inline static token flag(int bitmask);

        //- True if the character is a punctuation separator (eg, in ISstream).
        //  Since it could also start a number, SUBTRACT is not included as
        //  a separator.
        //
        //  \param c the character to test, passed as int for consistency with
        //      isdigit, isspace etc.
        inline static bool isseparator(int c);


    // Member functions

    // Status

        //- Return the name of the token type
        word name() const;

        //- Return the token type
        inline tokenType type() const;

        //- Change the token type, for similar types.
        //  This can be used to change between string-like variants
        //  (eg, STRING, VARIABLE, VERBATIMSTRING)
        //  To change types entirely (eg, STRING to DOUBLE_SCALAR),
        //  use the corresponding assignment operator.
        //
        //  \return true if the change was successful or no change was required
        inline bool setType(const tokenType variant);

        //- The line number for the token
        inline label lineNumber() const;

        //- The line number for the token
        inline label& lineNumber();

        //- True if token is not UNDEFINED or ERROR
        inline bool good() const;

        //- True if token is UNDEFINED
        inline bool undefined() const;

        //- True if token is ERROR
        inline bool error() const;

        //- True if token is FLAG
        inline bool isFlag() const;

        //- True if token is PUNCTUATION
        inline bool isPunctuation() const;

        //- True if token is PUNCTUATION and isseparator
        inline bool isSeparator() const;

        //- True if token is LABEL
        inline bool isLabel() const;

        //- True if token is FLOAT_SCALAR
        inline bool isFloatScalar() const;

        //- True if token is DOUBLE_SCALAR
        inline bool isDoubleScalar() const;

        //- True if token is FLOAT_SCALAR or DOUBLE_SCALAR
        inline bool isScalar() const;

        //- True if token is LABEL, FLOAT_SCALAR or DOUBLE_SCALAR
        inline bool isNumber() const;

        //- True if token is WORD
        inline bool isWord() const;

        //- True if token is STRING, VARIABLE or VERBATIMSTRING
        inline bool isString() const;

        //- True if token is VARIABLE
        inline bool isVariable() const;

        //- True if token is COMPOUND
        inline bool isCompound() const;


    // Access

        //- Return flag bitmask
        //  Report FatalIOError and return NO_FLAG if token is not FLAG
        inline int flagToken() const;

        //- Return punctuation character.
        //  Report FatalIOError and return \b \\0 if token is not PUNCTUATION
        inline punctuationToken pToken() const;

        //- Return label value.
        //  Report FatalIOError and return \b 0 if token is not LABEL
        inline label labelToken() const;

        //- Return float value.
        //  Report FatalIOError and return \b 0.0 if token is not FLOAT_SCALAR
        inline floatScalar floatScalarToken() const;

        //- Return double value.
        //  Report FatalIOError and return \b 0.0 if token is not DOUBLE_SCALAR
        inline doubleScalar doubleScalarToken() const;

        //- Return float or double value.
        //  Report FatalIOError and return \b 0.0 if token is not a
        //  FLOAT_SCALAR or DOUBLE_SCALAR
        inline scalar scalarToken() const;

        //- Return label, float or double value.
        //  Report FatalIOError and return \b 0.0 if token is not a
        //  LABEL, FLOAT_SCALAR or DOUBLE_SCALAR
        inline scalar number() const;

        //- Return const reference to the word contents.
        //  Report FatalIOError and return \b "" if token is not a WORD
        inline const word& wordToken() const;

        //- Return const reference to the string contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  STRING, VARIABLE or VERBATIMSTRING
        inline const string& stringToken() const;

        //- Read access for compound token
        inline const compound& compoundToken() const;

        //- Return reference to compound token and decrease its internal
        //- refCound accordingly.
        //  The Istream is used for reference error messages only.
        compound& transferCompoundToken(const Istream& is);


    // Edit

        //- Clear token and set to be in an error state.
        inline void setBad();

        //- Swap token contents: type, data, line-number
        inline void swap(token& tok);


    // Info

        //- Return info proxy for printing token information to a stream
        InfoProxy<token> info() const
        {
            return *this;
        }


    // Member operators

    // Assignment

        //- Copy assign
        inline void operator=(const token& tok);

        //- Move assign
        inline void operator=(token&& tok);

        //- Copy assign from punctuation
        inline void operator=(const punctuationToken p);

        //- Copy assign from label
        inline void operator=(const label val);

        //- Copy assign from float
        inline void operator=(const floatScalar val);

        //- Copy assign from double
        inline void operator=(const doubleScalar val);

        //- Copy assign from word
        inline void operator=(const word& w);

        //- Copy assign from string
        inline void operator=(const string& str);

        //- Move assign from word
        inline void operator=(word&& w);

        //- Move assign from string
        inline void operator=(string&& str);

        //- Transfer word pointer to the token
        //  \deprecated in favour of using move assign from word
        inline void operator=(word* wordPtr);

        //- Transfer string pointer to the token
        //  \deprecated in favour of using move assign from string
        inline void operator=(string* stringPtr);

        //- Assign compound with reference counting to token
        inline void operator=(compound* compoundPtr);


    // Equality

        inline bool operator==(const token& t) const;

        inline bool operator==(const punctuationToken p) const;
        inline bool operator==(const label val) const;
        inline bool operator==(const floatScalar val) const;
        inline bool operator==(const doubleScalar val) const;
        inline bool operator==(const word& w) const;
        inline bool operator==(const string& str) const;


    // Inequality

        inline bool operator!=(const token& t) const;
        inline bool operator!=(const punctuationToken p) const;
        inline bool operator!=(const label val) const;
        inline bool operator!=(const floatScalar val) const;
        inline bool operator!=(const doubleScalar val) const;
        inline bool operator!=(const word& w) const;
        inline bool operator!=(const string& str) const;


    // IOstream operators

        friend Istream& operator>>(Istream& is, token& t);
        friend Ostream& operator<<(Ostream& os, const token& t);

        friend Ostream& operator<<(Ostream&, const punctuationToken&);
        friend ostream& operator<<(ostream&, const punctuationToken&);

        friend ostream& operator<<(ostream&, const InfoProxy<token>&);
};


Ostream& operator<<(Ostream&, const token::punctuationToken&);
ostream& operator<<(ostream&, const token::punctuationToken&);
Ostream& operator<<(Ostream&, const token::compound&);

ostream& operator<<(ostream&, const InfoProxy<token>&);

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<token>& ip);

#define defineCompoundTypeName(Type, Name)                                     \
    defineTemplateTypeNameAndDebugWithName(token::Compound<Type>, #Type, 0);

#define addCompoundToRunTimeSelectionTable(Type, Name)                         \
    token::compound::addIstreamConstructorToTable<token::Compound<Type>>       \
        add##Name##IstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tokenI.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
