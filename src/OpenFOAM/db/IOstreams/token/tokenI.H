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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::token::clear()
{
    if (type_ == tokenType::WORD)
    {
        delete wordTokenPtr_;
    }
    else if
    (
        type_ == tokenType::STRING
     || type_ == tokenType::VARIABLE
     || type_ == tokenType::VERBATIMSTRING
    )
    {
        delete stringTokenPtr_;
    }
    else if (type_ == tokenType::COMPOUND)
    {
        if (compoundTokenPtr_->unique())
        {
            delete compoundTokenPtr_;
        }
        else
        {
            compoundTokenPtr_->refCount::operator--();
        }
    }

    type_ = tokenType::UNDEFINED;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::token::token()
:
    type_(tokenType::UNDEFINED),
    lineNumber_(0)
{}


inline Foam::token::token(const token& t)
:
    type_(t.type_),
    lineNumber_(t.lineNumber_)
{
    switch (type_)
    {
        case tokenType::UNDEFINED:
        break;

        case tokenType::PUNCTUATION:
            punctuationToken_ = t.punctuationToken_;
        break;

        case tokenType::WORD:
            wordTokenPtr_ = new word(*t.wordTokenPtr_);
        break;

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
            stringTokenPtr_ = new string(*t.stringTokenPtr_);
        break;

        case tokenType::LABEL:
            labelToken_ = t.labelToken_;
        break;

        case tokenType::FLOAT_SCALAR:
            floatScalarToken_ = t.floatScalarToken_;
        break;

        case tokenType::DOUBLE_SCALAR:
            doubleScalarToken_ = t.doubleScalarToken_;
        break;

        case tokenType::COMPOUND:
            compoundTokenPtr_ = t.compoundTokenPtr_;
            compoundTokenPtr_->refCount::operator++();
        break;

        case tokenType::ERROR:
        break;
    }
}


inline Foam::token::token(punctuationToken p)
:
    token(p, 0)
{}


inline Foam::token::token(const word& w)
:
    token(w, 0)
{}


inline Foam::token::token(const string& str)
:
    token(str, 0)
{}


inline Foam::token::token(const label val)
:
    token(val, 0)
{}


inline Foam::token::token(const floatScalar val)
:
    token(val, 0)
{}


inline Foam::token::token(const doubleScalar val)
:
    token(val, 0)
{}


inline Foam::token::token(punctuationToken p, label lineNumber)
:
    type_(tokenType::PUNCTUATION),
    punctuationToken_(p),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const word& w, label lineNumber)
:
    type_(tokenType::WORD),
    wordTokenPtr_(new word(w)),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const string& str, label lineNumber)
:
    type_(tokenType::STRING),
    stringTokenPtr_(new string(str)),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const label val, label lineNumber)
:
    type_(tokenType::LABEL),
    labelToken_(val),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const floatScalar val, label lineNumber)
:
    type_(tokenType::FLOAT_SCALAR),
    floatScalarToken_(val),
    lineNumber_(lineNumber)
{}


inline Foam::token::token(const doubleScalar val, label lineNumber)
:
    type_(tokenType::DOUBLE_SCALAR),
    doubleScalarToken_(val),
    lineNumber_(lineNumber)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::token::~token()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::token::tokenType  Foam::token::type() const
{
    return type_;
}

inline Foam::token::tokenType&  Foam::token::type()
{
    return type_;
}

inline bool Foam::token::good() const
{
    return (type_ != tokenType::ERROR && type_ != tokenType::UNDEFINED);
}

inline bool Foam::token::undefined() const
{
    return (type_ == tokenType::UNDEFINED);
}

inline bool Foam::token::error() const
{
    return (type_ == tokenType::ERROR);
}

inline bool Foam::token::isPunctuation() const
{
    return (type_ == tokenType::PUNCTUATION);
}

inline Foam::token::punctuationToken  Foam::token::pToken() const
{
    if (type_ == tokenType::PUNCTUATION)
    {
        return punctuationToken_;
    }
    else
    {
        parseError("punctuation character");
        return NULL_TOKEN;
    }
}

inline bool Foam::token::isWord() const
{
    return (type_ == tokenType::WORD);
}

inline const Foam::word& Foam::token::wordToken() const
{
    if (type_ == tokenType::WORD)
    {
        return *wordTokenPtr_;
    }
    else
    {
        parseError(word::typeName);
        return word::null;
    }
}

inline bool Foam::token::isVariable() const
{
    return (type_ == tokenType::VARIABLE);
}

inline bool Foam::token::isString() const
{
    return
    (
        type_ == tokenType::STRING
     || type_ == tokenType::VARIABLE
     || type_ == tokenType::VERBATIMSTRING
    );
}

inline const Foam::string& Foam::token::stringToken() const
{
    if
    (
        type_ == tokenType::STRING
     || type_ == tokenType::VARIABLE
     || type_ == tokenType::VERBATIMSTRING
    )
    {
        return *stringTokenPtr_;
    }
    else
    {
        parseError(string::typeName);
        return string::null;
    }
}

inline bool Foam::token::isLabel() const
{
    return (type_ == tokenType::LABEL);
}

inline Foam::label Foam::token::labelToken() const
{
    if (type_ == tokenType::LABEL)
    {
        return labelToken_;
    }
    else
    {
        parseError(pTraits<label>::typeName);
        return 0;
    }
}

inline bool Foam::token::isFloatScalar() const
{
    return (type_ == tokenType::FLOAT_SCALAR);
}

inline Foam::floatScalar Foam::token::floatScalarToken() const
{
    if (type_ == tokenType::FLOAT_SCALAR)
    {
        return floatScalarToken_;
    }
    else
    {
        parseError("floatScalar");
        return 0.0;
    }
}


inline bool Foam::token::isDoubleScalar() const
{
    return (type_ == tokenType::DOUBLE_SCALAR);
}

inline Foam::doubleScalar Foam::token::doubleScalarToken() const
{
    if (type_ == tokenType::DOUBLE_SCALAR)
    {
        return doubleScalarToken_;
    }
    else
    {
        parseError("doubleScalar");
        return 0.0;
    }
}


inline bool Foam::token::isScalar() const
{
    return
    (
        type_ == tokenType::FLOAT_SCALAR
     || type_ == tokenType::DOUBLE_SCALAR
    );
}

inline Foam::scalar Foam::token::scalarToken() const
{
    if (type_ == tokenType::FLOAT_SCALAR)
    {
        return floatScalarToken_;
    }
    else if (type_ == tokenType::DOUBLE_SCALAR)
    {
        return doubleScalarToken_;
    }
    else
    {
        parseError(pTraits<scalar>::typeName);
        return 0.0;
    }
}

inline bool Foam::token::isNumber() const
{
    return (type_ == tokenType::LABEL || isScalar());
}

inline Foam::scalar Foam::token::number() const
{
    if (type_ == tokenType::LABEL)
    {
        return labelToken_;
    }
    else if (isScalar())
    {
        return scalarToken();
    }
    else
    {
        parseError("number (label or scalar)");
        return 0.0;
    }
}

inline bool Foam::token::isCompound() const
{
    return (type_ == tokenType::COMPOUND);
}

inline const Foam::token::compound& Foam::token::compoundToken() const
{
    if (type_ == tokenType::COMPOUND)
    {
        return *compoundTokenPtr_;
    }
    else
    {
        parseError("compound");
        return *compoundTokenPtr_;
    }
}


inline Foam::label Foam::token::lineNumber() const
{
    return lineNumber_;
}

inline Foam::label& Foam::token::lineNumber()
{
    return lineNumber_;
}


inline void Foam::token::setBad()
{
    clear();
    type_ = tokenType::ERROR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::token::operator=(const token& t)
{
    clear();
    type_ = t.type_;

    switch (type_)
    {
        case tokenType::UNDEFINED:
        break;

        case tokenType::PUNCTUATION:
            punctuationToken_ = t.punctuationToken_;
        break;

        case tokenType::WORD:
            wordTokenPtr_ = new word(*t.wordTokenPtr_);
        break;

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
            stringTokenPtr_ = new string(*t.stringTokenPtr_);
        break;

        case tokenType::LABEL:
            labelToken_ = t.labelToken_;
        break;

        case tokenType::FLOAT_SCALAR:
            floatScalarToken_ = t.floatScalarToken_;
        break;

        case tokenType::DOUBLE_SCALAR:
            doubleScalarToken_ = t.doubleScalarToken_;
        break;

        case tokenType::COMPOUND:
            compoundTokenPtr_ = t.compoundTokenPtr_;
            compoundTokenPtr_->refCount::operator++();
        break;

        case tokenType::ERROR:
        break;
    }

    lineNumber_ = t.lineNumber_;
}

inline void Foam::token::operator=(const punctuationToken p)
{
    clear();
    type_ = tokenType::PUNCTUATION;
    punctuationToken_ = p;
}

inline void Foam::token::operator=(word* wPtr)
{
    clear();
    type_ = tokenType::WORD;
    wordTokenPtr_ = wPtr;
}

inline void Foam::token::operator=(const word& w)
{
    operator=(new word(w));
}

inline void Foam::token::operator=(string* strPtr)
{
    clear();
    type_ = tokenType::STRING;
    stringTokenPtr_ = strPtr;
}

inline void Foam::token::operator=(const string& str)
{
    operator=(new string(str));
}

inline void Foam::token::operator=(const label val)
{
    clear();
    type_ = tokenType::LABEL;
    labelToken_ = val;
}

inline void Foam::token::operator=(const floatScalar val)
{
    clear();
    type_ = tokenType::FLOAT_SCALAR;
    floatScalarToken_ = val;
}

inline void Foam::token::operator=(const doubleScalar val)
{
    clear();
    type_ = tokenType::DOUBLE_SCALAR;
    doubleScalarToken_ = val;
}

inline void Foam::token::operator=(Foam::token::compound* compPtr)
{
    clear();
    type_ = tokenType::COMPOUND;
    compoundTokenPtr_ = compPtr;
}


inline bool Foam::token::operator==(const token& t) const
{
    if (type_ != t.type_)
    {
        return false;
    }

    switch (type_)
    {
        case tokenType::UNDEFINED:
            return true;

        case tokenType::PUNCTUATION:
            return punctuationToken_ == t.punctuationToken_;

        case tokenType::WORD:
            return *wordTokenPtr_ == *t.wordTokenPtr_;

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
            return *stringTokenPtr_ == *t.stringTokenPtr_;

        case tokenType::LABEL:
            return labelToken_ == t.labelToken_;

        case tokenType::FLOAT_SCALAR:
            return equal(floatScalarToken_, t.floatScalarToken_);

        case tokenType::DOUBLE_SCALAR:
            return equal(doubleScalarToken_, t.doubleScalarToken_);

        case tokenType::COMPOUND:
            return compoundTokenPtr_ == t.compoundTokenPtr_;

        case tokenType::ERROR:
            return true;
    }

    return false;
}

inline bool Foam::token::operator==(const punctuationToken p) const
{
    return (type_ == tokenType::PUNCTUATION && punctuationToken_ == p);
}

inline bool Foam::token::operator==(const word& w) const
{
    return (type_ == tokenType::WORD && wordToken() == w);
}

inline bool Foam::token::operator==(const string& str) const
{
    return
    (
        (
            type_ == tokenType::STRING
         || type_ == tokenType::VARIABLE
         || type_ == tokenType::VERBATIMSTRING
        )
     && stringToken() == str
    );
}

inline bool Foam::token::operator==(const label val) const
{
    return
    (
        type_ == tokenType::LABEL
     && labelToken_ == val
    );
}

inline bool Foam::token::operator==(const floatScalar val) const
{
    return
    (
        type_ == tokenType::FLOAT_SCALAR
     && equal(floatScalarToken_, val)
    );
}

inline bool Foam::token::operator==(const doubleScalar val) const
{
    return
    (
        type_ == tokenType::DOUBLE_SCALAR
     && equal(doubleScalarToken_, val)
    );
}

inline bool Foam::token::operator!=(const token& t) const
{
    return !operator==(t);
}

inline bool Foam::token::operator!=(const punctuationToken p) const
{
    return !operator==(p);
}

inline bool Foam::token::operator!=(const word& w) const
{
    return !operator==(w);
}

inline bool Foam::token::operator!=(const string& str) const
{
    return !operator==(str);
}

inline bool Foam::token::operator!=(const label val) const
{
    return !operator==(val);
}

inline bool Foam::token::operator!=(const floatScalar val) const
{
    return !operator==(val);
}

inline bool Foam::token::operator!=(const doubleScalar val) const
{
    return !operator==(val);
}


// ************************************************************************* //
