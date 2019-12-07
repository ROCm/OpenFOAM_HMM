/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

// Forward Declarations
class token;
Ostream& operator<<(Ostream& os, const token& tok);


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
        UNDEFINED = 0,    //!< An undefined token-type

        // Fundamental types
        FLAG,             //!< stream flag (1-byte bitmask)
        PUNCTUATION,      //!< single character punctuation
        BOOL,             //!< boolean type
        LABEL,            //!< label (integer) type
        FLOAT_SCALAR,     //!< float (single-precision) type
        DOUBLE_SCALAR,    //!< double (double-precision) type

        // Pointer types
        WORD,             //!< A Foam::word
        STRING,           //!< A string
        VARIABLE,         //!< A dictionary \c \$variable (string variant)
        VERBATIMSTRING,   //!< Verbatim string content
        COMPOUND,         //!< Compound type such as \c List\<label\> etc.

        ERROR             //!< A token error encountered
    };


    //- Stream or output control flags (1-byte width)
    enum flagType
    {
        NO_FLAG = 0,   //!< No flags
        ASCII   = 1,   //!< ASCII-mode stream
        BINARY  = 2    //!< BINARY-mode stream
    };


    //- Standard punctuation tokens (a character)
    enum punctuationToken : char
    {
        NULL_TOKEN     = '\0',  //!< Nul character
        SPACE          = ' ',   //!< Space [isspace]
        TAB            = '\t',  //!< Tab [isspace]
        NL             = '\n',  //!< Newline [isspace]

        END_STATEMENT  = ';',   //!< End entry [#isseparator]
        BEGIN_LIST     = '(',   //!< Begin list [#isseparator]
        END_LIST       = ')',   //!< End list [#isseparator]
        BEGIN_SQR      = '[',   //!< Begin dimensions [#isseparator]
        END_SQR        = ']',   //!< End dimensions [#isseparator]
        BEGIN_BLOCK    = '{',   //!< Begin block [#isseparator]
        END_BLOCK      = '}',   //!< End block [#isseparator]
        COLON          = ':',   //!< Colon [#isseparator]
        COMMA          = ',',   //!< Comma [#isseparator]
        HASH           = '#',   //!< Hash - directive or verbatim string
        DOLLAR         = '$',   //!< Dollar - start variable
        ATSYM          = '@',   //!< At
        SQUOTE         = '\'',  //!< Single quote
        DQUOTE         = '"',   //!< Double quote

        ASSIGN         = '=',   //!< Assignment/equals [#isseparator]
        ADD            = '+',   //!< Addition [#isseparator]
        SUBTRACT       = '-',   //!< Subtract or start of negative number
        MULTIPLY       = '*',   //!< Multiply [#isseparator]
        DIVIDE         = '/',   //!< Divide [#isseparator]

        BEGIN_STRING   = DQUOTE, //!< Begin string with double quote
        END_STRING     = DQUOTE  //!< End string with double quote
    };


    //- Abstract base class for complex tokens
    class compound
    :
        public refCount
    {
        bool empty_;

        //- No copy construct
        compound(const compound&) = delete;

        //- No copy assignment
        compound& operator=(const compound&) = delete;

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


        // Selectors

            //- Select null constructed
            static autoPtr<compound> New(const word& type, Istream& is);

            //- Return true if name is a known (registered) compound type
            static bool isCompound(const word& name);


        //- Destructor
        virtual ~compound() = default;


        // Member Functions

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


    //- An undefined token
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


    // Private Data

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

        // Parse error, expected 'expected', found ...
        void parseError(const char* expected) const;


public:

    // Static Data Members

        //- The type name is "token"
        static constexpr const char* const typeName = "token";


    // Constructors

        //- Construct null, initialized to an UNDEFINED token.
        inline constexpr token() noexcept;

        //- Copy construct
        inline token(const token& t);

        //- Move construct. The original token is left as UNDEFINED.
        inline token(token&& t);

        //- Construct punctuation character token
        inline explicit token(punctuationToken p, label lineNumber=0);

        //- Construct label token
        inline explicit token(const label val, label lineNumber=0);

        //- Construct float token
        inline explicit token(const floatScalar val, label lineNumber=0);

        //- Construct double token
        inline explicit token(const doubleScalar val, label lineNumber=0);

        //- Copy construct word token
        inline explicit token(const word& w, label lineNumber=0);

        //- Copy construct string token
        inline explicit token(const string& str, label lineNumber=0);

        //- Move construct word token
        inline explicit token(word&& w, label lineNumber=0);

        //- Move construct string token
        inline explicit token(string&& str, label lineNumber=0);

        //- Construct from Istream
        explicit token(Istream& is);


    //- Destructor
    inline ~token();


    // Static Member Functions

        //- Create a bool token.
        inline static token boolean(bool on);

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


    // Member Functions

    // Status

        //- Return the name of the token type
        word name() const;

        //- Return the token type
        inline tokenType type() const;

        //- Change the token type, for similar types.
        //  This can be used to change between string-like variants
        //  (eg, STRING, VARIABLE, etc)
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

        //- Token is UNDEFINED
        inline bool undefined() const;

        //- Token is ERROR
        inline bool error() const;

        //- Token is BOOL
        inline bool isBool() const;

        //- Token is FLAG
        inline bool isFlag() const;

        //- Token is PUNCTUATION
        inline bool isPunctuation() const;

        //- Token is PUNCTUATION and isseparator
        inline bool isSeparator() const;

        //- Token is LABEL
        inline bool isLabel() const;

        //- Token is FLOAT_SCALAR
        inline bool isFloatScalar() const;

        //- Token is DOUBLE_SCALAR
        inline bool isDoubleScalar() const;

        //- Token is FLOAT_SCALAR or DOUBLE_SCALAR
        inline bool isScalar() const;

        //- Token is LABEL, FLOAT_SCALAR or DOUBLE_SCALAR
        inline bool isNumber() const;

        //- Token is WORD
        inline bool isWord() const;

        //- Token is STRING, VARIABLE or VERBATIM string
        inline bool isString() const;

        //- Token is VARIABLE
        inline bool isVariable() const;

        //- Token is VERBATIM string
        inline bool isVerbatim() const;

        //- Token is WORD, STRING, VARIABLE or VERBATIM
        inline bool isStringType() const;

        //- Token is COMPOUND
        inline bool isCompound() const;


    // Access

        //- Return boolean token value.
        //  Report FatalIOError and return false if token is not BOOL or LABEL
        inline bool boolToken() const;

        //- Return flag bitmask value.
        //  Report FatalIOError and return NO_FLAG if token is not FLAG
        inline int flagToken() const;

        //- Return punctuation character.
        //  Report FatalIOError and return \b \\0 if token is not PUNCTUATION
        inline punctuationToken pToken() const;

        //- Return label value.
        //  Report FatalIOError and return \b 0 if token is not LABEL
        inline label labelToken() const;

        //- Return float value.
        //  Report FatalIOError and return \b 0 if token is not FLOAT_SCALAR
        inline floatScalar floatScalarToken() const;

        //- Return double value.
        //  Report FatalIOError and return \b 0 if token is not DOUBLE_SCALAR
        inline doubleScalar doubleScalarToken() const;

        //- Return float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  FLOAT_SCALAR or DOUBLE_SCALAR
        inline scalar scalarToken() const;

        //- Return label, float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  LABEL, FLOAT_SCALAR or DOUBLE_SCALAR
        inline scalar number() const;

        //- Return const reference to the word contents.
        //  Report FatalIOError and return \b "" if token is not a WORD
        inline const word& wordToken() const;

        //- Return const reference to the string contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  STRING, VARIABLE, VERBATIM or an upcast WORD
        inline const string& stringToken() const;

        //- Read access for compound token
        inline const compound& compoundToken() const;

        //- Return reference to compound token and decrease its internal
        //- refCound accordingly.
        //  The Istream is used for reference error messages only.
        compound& transferCompoundToken(const Istream& is);


    // Edit

        //- Reset token to UNDEFINED and clear any allocated storage
        inline void reset();

        //- Clear token and set to be ERROR.
        inline void setBad();

        //- Swap token contents: type, data, line-number
        inline void swap(token& tok);


    // Info

        //- Return info proxy for printing token information to a stream
        InfoProxy<token> info() const
        {
            return *this;
        }


    // Member Operators

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

        //- Assign compound with reference counting to token
        inline void operator=(compound* compoundPtr);


    // Equality

        inline bool operator==(const token& tok) const;
        inline bool operator==(const punctuationToken p) const;
        inline bool operator==(const label val) const;
        inline bool operator==(const floatScalar val) const;
        inline bool operator==(const doubleScalar val) const;
        inline bool operator==(const std::string& s) const;


    // Inequality

        inline bool operator!=(const token& tok) const;
        inline bool operator!=(const punctuationToken p) const;
        inline bool operator!=(const label val) const;
        inline bool operator!=(const floatScalar val) const;
        inline bool operator!=(const doubleScalar val) const;
        inline bool operator!=(const std::string& s) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const token& tok);

        friend Ostream& operator<<(Ostream& os, const punctuationToken& pt);
        friend ostream& operator<<(ostream& os, const punctuationToken& pt);

        friend ostream& operator<<(ostream& os, const InfoProxy<token>& ct);


    // Housekeeping

        //- Deprecated(2017-11) transfer word pointer to the token
        //  \deprecated(2017-11) - use move assign from word
        void operator=(word*) = delete;

        //- Deprecated(2017-11) transfer string pointer to the token
        //  \deprecated(2017-11) - use move assign from string
        void operator=(string*) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IOstream Operators

Istream& operator>>(Istream& is, token& tok);
Ostream& operator<<(Ostream& os, const token::punctuationToken& pt);
ostream& operator<<(ostream& os, const token::punctuationToken& pt);
Ostream& operator<<(Ostream& os, const token::compound& ct);

ostream& operator<<(ostream& os, const InfoProxy<token>& ip);

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<token>& ip);


// Handling of compound types

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
