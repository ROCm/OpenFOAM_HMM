/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    A token holds items read from Istream.

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

    //- Enumeration defining the types of token
    enum tokenType
    {
        UNDEFINED,

        PUNCTUATION,
        WORD,
        VARIABLE,
        STRING,
        VERBATIMSTRING,
        LABEL,
        FLOAT_SCALAR,
        DOUBLE_SCALAR,
        COMPOUND,

        ERROR
    };


    //- Standard punctuation tokens
    enum punctuationToken
    {
        NULL_TOKEN     = '\0',
        SPACE          = ' ',
        TAB            = '\t',
        NL             = '\n',

        END_STATEMENT  = ';',
        BEGIN_LIST     = '(',
        END_LIST       = ')',
        BEGIN_SQR      = '[',
        END_SQR        = ']',
        BEGIN_BLOCK    = '{',
        END_BLOCK      = '}',
        COLON          = ':',
        COMMA          = ',',
        HASH           = '#',

        BEGIN_STRING   = '"',
        END_STRING     = BEGIN_STRING,

        ASSIGN         = '=',
        ADD            = '+',
        SUBTRACT       = '-',
        MULTIPLY       = '*',
        DIVIDE         = '/'
    };


    //- Abstract base class for complex tokens
    class compound
    :
        public refCount
    {
        // Private data

            bool empty_;


        // Private Member Functions

            //- Disallow default bitwise copy construct
            compound(const compound&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const compound&) = delete;


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
            static autoPtr<compound> New(const word& type, Istream&);


        //- Destructor
        virtual ~compound();


        // Member Functions

            // Access

                //- Return true if name is a compound type
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

            // Write

                virtual void write(Ostream& os) const = 0;


        // IOstream Operators

            friend Ostream& operator<<(Ostream&, const compound&);
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
    static token undefinedToken;


private:

    // Private data

        //- The token type
        tokenType type_;

        //- Anonymous Union of token types
        union
        {
            punctuationToken punctuationToken_;
            word* wordTokenPtr_;
            string* stringTokenPtr_;
            label labelToken_;
            floatScalar floatScalarToken_;
            doubleScalar doubleScalarToken_;
            mutable compound* compoundTokenPtr_;
        };

        //- Line number in the file this token was read from
        label lineNumber_;


    // Private Member Functions

        //- Clear any allocated storage (word or string)
        inline void clear();

        // Parse error, expected 'expected', found ...
        void parseError(const char* expected) const;


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        inline token();

        //- Construct as copy
        inline token(const token& t);

        //- Construct punctuation character token
        inline explicit token(punctuationToken p);

        //- Construct word token
        inline explicit token(const word& w);

        //- Construct string token
        inline explicit token(const string& str);

        //- Construct label token
        inline explicit token(const label val);

        //- Construct floatScalar token
        inline explicit token(const floatScalar val);

        //- Construct doubleScalar token
        inline explicit token(const doubleScalar val);


        //- Construct punctuation character token
        inline token(punctuationToken p, const label lineNumber);

        //- Construct word token
        inline token(const word& w, const label lineNumber);

        //- Construct string token
        inline token(const string& str, const label lineNumber);

        //- Construct label token
        inline token(const label val, const label lineNumber);

        //- Construct floatScalar token
        inline token(const floatScalar val, const label lineNumber);

        //- Construct doubleScalar token
        inline token(const doubleScalar val, const label lineNumber);

        //- Construct from Istream
        token(Istream& is);


    //- Destructor
    inline ~token();


    // Member functions

        // Access

            inline tokenType type() const;
            inline tokenType& type();

            inline bool good() const;
            inline bool undefined() const;
            inline bool error() const;

            inline bool isPunctuation() const;
            inline punctuationToken pToken() const;

            inline bool isWord() const;
            inline const word& wordToken() const;

            inline bool isVariable() const;

            inline bool isString() const;
            inline const string& stringToken() const;

            inline bool isLabel() const;
            inline label labelToken() const;

            inline bool isFloatScalar() const;
            inline floatScalar floatScalarToken() const;

            inline bool isDoubleScalar() const;
            inline doubleScalar doubleScalarToken() const;

            inline bool isScalar() const;
            inline scalar scalarToken() const;

            inline bool isNumber() const;
            inline scalar number() const;

            inline bool isCompound() const;
            inline const compound& compoundToken() const;
            compound& transferCompoundToken(const Istream& is);

            inline label lineNumber() const;
            inline label& lineNumber();


        // Edit

            //- Set bad
            inline void setBad();


        // Info

            //- Return info proxy.
            //  Used to print token information to a stream
            InfoProxy<token> info() const
            {
                return *this;
            }


    // Member operators

        // Assignment

            inline void operator=(const token& t);

            inline void operator=(const punctuationToken p);

            inline void operator=(word* wPtr);
            inline void operator=(const word& w);

            inline void operator=(string* strPtr);
            inline void operator=(const string& str);

            inline void operator=(const label val);
            inline void operator=(const floatScalar val);
            inline void operator=(const doubleScalar val);

            inline void operator=(compound* compPtr);


        // Equality

            inline bool operator==(const token& t) const;
            inline bool operator==(const punctuationToken p) const;
            inline bool operator==(const word& w) const;
            inline bool operator==(const string& str) const;
            inline bool operator==(const label val) const;
            inline bool operator==(const floatScalar val) const;
            inline bool operator==(const doubleScalar val) const;


        // Inequality

            inline bool operator!=(const token& t) const;
            inline bool operator!=(const punctuationToken p) const;
            inline bool operator!=(const word& w) const;
            inline bool operator!=(const string& str) const;
            inline bool operator!=(const label val) const;
            inline bool operator!=(const floatScalar val) const;
            inline bool operator!=(const doubleScalar val) const;


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
