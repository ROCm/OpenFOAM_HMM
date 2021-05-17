/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    //  Since the enumeration is used to tag content in Pstream, it is of
    //  type \c char and shall have values that do not overlap with regular
    //  punctuation characters.
    enum tokenType : char
    {
        UNDEFINED = '\0',     //!< An undefined token-type
        ERROR     = '\x80',   //!< Token error encountered

        // Fundamental types
        FLAG,             //!< stream flag (1-byte bitmask)
        PUNCTUATION,      //!< single character punctuation
        BOOL,             //!< boolean type
        LABEL,            //!< label (integer) type
        FLOAT,            //!< float (single-precision) type
        DOUBLE,           //!< double (double-precision) type

        // Pointer types
        WORD,             //!< Foam::word
        STRING,           //!< Foam::string (usually double-quoted)
        COMPOUND,         //!< Compound type such as \c List\<label\> etc.

        DIRECTIVE,        //!< Word-variant: dictionary \c \#directive
                          //!< stored with sigil
        EXPRESSION,       //!< String-variant: math expression for evaluation
                          //!< stored with delimiters
        VARIABLE,         //!< String-variant: dictionary \c \$variable
                          //!< stored with sigil
        VERBATIM,         //!< String-variant: verbatim string content
                          //!< stored without delimiters

        // Aliases
        FLOAT_SCALAR = FLOAT,
        DOUBLE_SCALAR = DOUBLE,
        VERBATIMSTRING = VERBATIM
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
        TAB            = '\t',  //!< Tab [isspace]
        NL             = '\n',  //!< Newline [isspace]
        SPACE          = ' ',   //!< Space [isspace]

        COLON          = ':',   //!< Colon [#isseparator]
        SEMICOLON      = ';',   //!< Semicolon [#isseparator]
        COMMA          = ',',   //!< Comma [#isseparator]
        HASH           = '#',   //!< Hash - directive or start verbatim string
        DOLLAR         = '$',   //!< Dollar - start variable or expression
        QUESTION       = '?',   //!< Question mark (eg, ternary)
        ATSYM          = '@',   //!< The 'at' symbol
        SQUOTE         = '\'',  //!< Single quote
        DQUOTE         = '"',   //!< Double quote

        ASSIGN         = '=',   //!< Assignment/equals [#isseparator]
        PLUS           = '+',   //!< Addition [#isseparator]
        MINUS          = '-',   //!< Subtract or start of negative number
        MULTIPLY       = '*',   //!< Multiply [#isseparator]
        DIVIDE         = '/',   //!< Divide [#isseparator]

        LPAREN         = '(',   //!< Left parenthesis [#isseparator]
        RPAREN         = ')',   //!< Right parenthesis [#isseparator]
        LSQUARE        = '[',   //!< Left square bracket [#isseparator]
        RSQUARE        = ']',   //!< Right square bracket [#isseparator]
        LBRACE         = '{',   //!< Left brace [#isseparator]
        RBRACE         = '}',   //!< Right brace [#isseparator]

        // With semantically meaning

        ADD            = PLUS,      //!< Addition [#isseparator]
        SUBTRACT       = MINUS,     //!< Subtract or start of negative number
        END_STATEMENT  = SEMICOLON, //!< End entry [#isseparator]
        BEGIN_LIST     = LPAREN,    //!< Begin list [#isseparator]
        END_LIST       = RPAREN,    //!< End list [#isseparator]
        BEGIN_SQR      = LSQUARE,   //!< Begin dimensions [#isseparator]
        END_SQR        = RSQUARE,   //!< End dimensions [#isseparator]
        BEGIN_BLOCK    = LBRACE,    //!< Begin block [#isseparator]
        END_BLOCK      = RBRACE,    //!< End block [#isseparator]

        BEGIN_STRING   = DQUOTE,    //!< Begin string with double quote
        END_STRING     = DQUOTE     //!< End string with double quote
    };


    //- Abstract base class for complex tokens
    class compound
    :
        public refCount
    {
        //- Has compound token already been transferred
        bool moved_;

        //- No copy construct
        compound(const compound&) = delete;

        //- No copy assignment
        compound& operator=(const compound&) = delete;

    public:

        //- Declare type-name, virtual type (with debug switch)
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

            //- Default construct
            constexpr compound() noexcept
            :
                moved_(false)
            {}

            //- Construct compound from Istream
            static autoPtr<compound> New(const word& type, Istream& is);


        //- Destructor
        virtual ~compound() noexcept = default;


        // Member Functions

            //- Test if name is a known (registered) compound type
            static bool isCompound(const word& name);

            //- Get compound transferred status
            bool moved() const noexcept
            {
                return moved_;
            }

            //- Set compound transferred status
            void moved(bool b) noexcept
            {
                moved_ = b;
            }

            //- The size of the underlying content
            virtual label size() const = 0;

            //- Redirect write to underlying content
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

        //- Declare type-name, virtual type (with debug switch)
        TypeName("Compound<T>");

        // Constructors

            //- Copy construct
            explicit Compound(const T& val)
            :
                T(val)
            {}

            //- Move construct
            explicit Compound(T&& val)
            :
                T(std::move(val))
            {}

            //- Read construct from Istream
            explicit Compound(Istream& is)
            :
                T(is)
            {}


        // Member Functions

            //- The size of the underlying content
            virtual label size() const
            {
                return T::size();
            }

            //- Redirect write to underlying content
            virtual void write(Ostream& os) const
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
        label line_;


    // Private Member Functions

        //- Set as UNDEFINED and zero the union content without any checking
        inline void setUndefined() noexcept;

        // Parse error, expected 'expected', found ...
        void parseError(const char* expected) const;


public:

    // Static Data Members

        //- The type name is "token"
        static constexpr const char* const typeName = "token";


    // Constructors

        //- Default construct, initialized to an UNDEFINED token.
        inline constexpr token() noexcept;

        //- Copy construct
        inline token(const token& t);

        //- Move construct. The original token is left as UNDEFINED.
        inline token(token&& t) noexcept;

        //- Construct punctuation character token
        inline explicit token(punctuationToken p, label lineNum=0) noexcept;

        //- Construct label token
        inline explicit token(const label val, label lineNum=0) noexcept;

        //- Construct float token
        inline explicit token(const floatScalar val, label lineNum=0) noexcept;

        //- Construct double token
        inline explicit token(const doubleScalar val, label lineNum=0) noexcept;

        //- Copy construct word token
        inline explicit token(const word& w, label lineNum=0);

        //- Copy construct string token
        inline explicit token(const string& str, label lineNum=0);

        //- Move construct word token
        inline explicit token(word&& w, label lineNum=0);

        //- Move construct string token
        inline explicit token(string&& str, label lineNum=0);

        //- Construct from a compound pointer, taking ownership
        inline explicit token(token::compound* ptr, label lineNum=0);

        //- Construct from Istream
        explicit token(Istream& is);


    //- Destructor
    inline ~token();


    // Static Functions

        //- Create a bool token.
        inline static token boolean(bool on) noexcept;

        //- Create a token with stream flags, no sanity check
        //
        //  \param bitmask the flags to set
        inline static token flag(int bitmask) noexcept;

        //- True if the character is a punctuation separator (eg, in ISstream).
        //  Since it could also start a number, SUBTRACT is not included as
        //  a separator.
        //
        //  \param c the character to test, passed as int for consistency with
        //      isdigit, isspace etc.
        inline static bool isseparator(int c) noexcept;


    // Member Functions

    // Status

        //- Return the name of the token type
        word name() const;

        //- Return the token type
        inline tokenType type() const noexcept;

        //- Change the token type, for similar types.
        //  This can be used to change between string-like variants
        //  (eg, STRING, VARIABLE, etc)
        //  To change types entirely (eg, STRING to DOUBLE),
        //  use the corresponding assignment operator.
        //
        //  \return true if the change was successful or no change was required
        inline bool setType(const tokenType tokType) noexcept;

        //- The line number for the token
        inline label lineNumber() const noexcept;

        //- Change token line number, return old value
        inline label lineNumber(const label lineNum) noexcept;

        //- True if token is not UNDEFINED or ERROR
        inline bool good() const noexcept;

        //- Token is UNDEFINED
        inline bool undefined() const noexcept;

        //- Token is ERROR
        inline bool error() const noexcept;

        //- Token is BOOL
        inline bool isBool() const noexcept;

        //- Token is FLAG
        inline bool isFlag() const noexcept;

        //- Token is PUNCTUATION
        inline bool isPunctuation() const noexcept;

        //- True if token is PUNCTUATION and equal to parameter
        inline bool isPunctuation(const punctuationToken p) const noexcept;

        //- Token is PUNCTUATION and isseparator
        inline bool isSeparator() const noexcept;

        //- Token is LABEL
        inline bool isLabel() const noexcept;

        //- True if token is LABEL and equal to parameter
        inline bool isLabel(const label val) const noexcept;

        //- Token is FLOAT
        inline bool isFloat() const noexcept;

        //- Token is DOUBLE
        inline bool isDouble() const noexcept;

        //- Token is FLOAT or DOUBLE
        inline bool isScalar() const noexcept;

        //- Token is LABEL, FLOAT or DOUBLE
        inline bool isNumber() const noexcept;

        //- Token is word-variant (WORD, DIRECTIVE)
        inline bool isWord() const noexcept;

        //- Token is word-variant and equal to parameter
        inline bool isWord(const std::string& s) const;

        //- Token is DIRECTIVE (word variant)
        inline bool isDirective() const noexcept;

        //- Token is (quoted) STRING (string variant)
        inline bool isQuotedString() const noexcept;

        //- Token is string-variant (STRING, EXPRESSION, VARIABLE, VERBATIM)
        inline bool isString() const noexcept;

        //- Token is EXPRESSION (string variant)
        inline bool isExpression() const noexcept;

        //- Token is VARIABLE (string variant)
        inline bool isVariable() const noexcept;

        //- Token is VERBATIM string (string variant)
        inline bool isVerbatim() const noexcept;

        //- Token is word-variant or string-variant
        //- (WORD, DIRECTIVE, STRING, EXPRESSION, VARIABLE, VERBATIM)
        inline bool isStringType() const noexcept;

        //- Token is COMPOUND
        inline bool isCompound() const noexcept;


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
        //  Report FatalIOError and return \b 0 if token is not FLOAT
        inline floatScalar floatToken() const;

        //- Return double value.
        //  Report FatalIOError and return \b 0 if token is not DOUBLE
        inline doubleScalar doubleToken() const;

        //- Return float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  FLOAT or DOUBLE
        inline scalar scalarToken() const;

        //- Return label, float or double value.
        //  Report FatalIOError and return \b 0 if token is not a
        //  LABEL, FLOAT or DOUBLE
        inline scalar number() const;

        //- Return const reference to the word contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  WORD or DIRECTIVE
        inline const word& wordToken() const;

        //- Return const reference to the string contents.
        //  Report FatalIOError and return \b "" if token is not a
        //  STRING, EXPRESSION, VARIABLE, VERBATIM
        //  or an upcast WORD or DIRECTIVE
        inline const string& stringToken() const;

        //- Read access for compound token
        inline const compound& compoundToken() const;

        //- Return reference to compound and mark internally as \em released.
        compound& transferCompoundToken();

        //- Return reference to compound and mark internally as \em released.
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

        //- Copy assign from word content
        inline void operator=(const word& w);

        //- Copy assign from string content
        inline void operator=(const string& str);

        //- Move assign from word content
        inline void operator=(word&& w);

        //- Move assign from string content
        inline void operator=(string&& str);

        //- Assign compound with reference counting to token
        inline void operator=(token::compound* ptr);

        //- Move assign from compound pointer
        inline void operator=(autoPtr<token::compound>&& ptr);


    // Equality

        inline bool operator==(const token& tok) const;
        inline bool operator==(const punctuationToken p) const noexcept;
        inline bool operator==(const label val) const noexcept;
        inline bool operator==(const floatScalar val) const noexcept;
        inline bool operator==(const doubleScalar val) const noexcept;
        inline bool operator==(const std::string& s) const;


    // Inequality

        inline bool operator!=(const token& tok) const;
        inline bool operator!=(const punctuationToken p) const noexcept;
        inline bool operator!=(const label val) const noexcept;
        inline bool operator!=(const floatScalar val) const noexcept;
        inline bool operator!=(const doubleScalar val) const noexcept;
        inline bool operator!=(const std::string& s) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const token& tok);

        friend Ostream& operator<<(Ostream& os, const punctuationToken& pt);
        friend ostream& operator<<(ostream& os, const punctuationToken& pt);

        friend ostream& operator<<(ostream& os, const InfoProxy<token>& ct);


    // Housekeeping

        //- Write access for the token line number
        //  \deprecated(2021-03) - use lineNumber(label)
        label& lineNumber() noexcept { return line_; }

        //- Token is FLOAT
        //  \deprecated(2020-01) - isFloat()
        bool isFloatScalar() const { return isFloat(); };

        //- Token is DOUBLE
        //  \deprecated(2020-01) - isDouble()
        bool isDoubleScalar() const { return isDouble(); }

        //- Return float value.
        //  \deprecated(2020-01) - floatToken()
        floatScalar floatScalarToken() const { return floatToken(); }

        //- Return double value.
        //  \deprecated(2020-01) - doubleToken()
        doubleScalar doubleScalarToken() const { return doubleToken(); }

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

//- Define compound using \a Type for its name
#define defineCompoundTypeName(Type, UnusedTag)                                \
    defineTemplateTypeNameAndDebugWithName(token::Compound<Type>, #Type, 0);

//- Define compound using \a Name for its name
#define defineNamedCompoundTypeName(Type, Name)                                \
    defineTemplateTypeNameAndDebugWithName(token::Compound<Type>, #Name, 0);

//- Add compound to selection table, lookup using typeName
#define addCompoundToRunTimeSelectionTable(Type, Tag)                          \
    token::compound::addIstreamConstructorToTable<token::Compound<Type>>       \
        add##Tag##IstreamConstructorToTable_;

//- Add compound to selection table, lookup as \a Name
#define addNamedCompoundToRunTimeSelectionTable(Type, Tag, Name)               \
    token::compound::addIstreamConstructorToTable<token::Compound<Type>>       \
        add##Tag##IstreamConstructorToTable_(#Name);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tokenI.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
