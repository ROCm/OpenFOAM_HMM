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

\*---------------------------------------------------------------------------*/

#include <algorithm>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::token Foam::token::flag(int bitmask)
{
    token tok;
    tok.type_ = tokenType::FLAG;
    tok.data_.flagVal = bitmask;

    return tok;
}


inline bool Foam::token::isseparator(int c)
{
    switch (c)
    {
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::ADD :
        // Excluded token::SUBTRACT since it could start a number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            return true;
        }

        default:
            break;
    }

    return false;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::token::setUndefined()
{
    type_ = tokenType::UNDEFINED;
    data_.int64Val = 0; // bit-wise zero for union content
    // leave lineNumber untouched - may still be needed
}


inline void Foam::token::clear()
{
    switch (type_)
    {
        case tokenType::WORD:
        {
            delete data_.wordPtr;
            break;
        }

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
        {
            delete data_.stringPtr;
            break;
        }

        case tokenType::COMPOUND:
        {
            if (data_.compoundPtr->unique())
            {
                delete data_.compoundPtr;
            }
            else
            {
                data_.compoundPtr->refCount::operator--();
            }
            break;
        }

        default:
            break;
    }

    setUndefined();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::token::token()
:
    data_(), // bit-wise zero for union content
    type_(tokenType::UNDEFINED),
    lineNumber_(0)
{}


inline Foam::token::token(const token& tok)
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    lineNumber_(tok.lineNumber_)
{
    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        case tokenType::WORD:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
            break;
        }

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
        {
            data_.stringPtr = new string(*tok.data_.stringPtr);
            break;
        }

        case tokenType::COMPOUND:
        {
            // Identical pointers, but increase the refCount
            data_.compoundPtr = tok.data_.compoundPtr;
            data_.compoundPtr->refCount::operator++();
            break;
        }

        default:
            break;
    }
}


inline Foam::token::token(token&& tok)
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    lineNumber_(tok.lineNumber_)
{
    tok.setUndefined(); // zero the union content without any checking
    tok.lineNumber_ = 0;
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
    data_(),
    type_(tokenType::PUNCTUATION),
    lineNumber_(lineNumber)
{
    data_.punctuationVal = p;
}


inline Foam::token::token(const label val, label lineNumber)
:
    data_(),
    type_(tokenType::LABEL),
    lineNumber_(lineNumber)
{
    data_.labelVal = val;
}


inline Foam::token::token(const floatScalar val, label lineNumber)
:
    data_(),
    type_(tokenType::FLOAT_SCALAR),
    lineNumber_(lineNumber)
{
    data_.floatVal = val;
}


inline Foam::token::token(const doubleScalar val, label lineNumber)
:
    data_(),
    type_(tokenType::DOUBLE_SCALAR),
    lineNumber_(lineNumber)
{
    data_.doubleVal = val;
}


inline Foam::token::token(const word& w, label lineNumber)
:
    data_(),
    type_(tokenType::WORD),
    lineNumber_(lineNumber)
{
    data_.wordPtr = new word(w);
}


inline Foam::token::token(const string& str, label lineNumber)
:
    data_(),
    type_(tokenType::STRING),
    lineNumber_(lineNumber)
{
    data_.stringPtr = new string(str);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::token::~token()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::token::swap(token& tok)
{
    std::swap(data_, tok.data_);
    std::swap(type_, tok.type_);
    std::swap(lineNumber_, tok.lineNumber_);
}


inline Foam::token::tokenType Foam::token::type() const
{
    return type_;
}


inline bool Foam::token::setType(token::tokenType variant)
{
    if (type_ == variant)
    {
        // No change required
        return true;
    }

    switch (variant)
    {
        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
        {
            switch (type_)
            {
                // could also go from WORD to STRING etc - to be decided
                case tokenType::STRING:
                case tokenType::VARIABLE:
                case tokenType::VERBATIMSTRING:
                    type_ = variant;
                    return true;
                    break;

                default:
                    return false;
                    break;
            }
        }

        default:
            break;
    }

    return false;
}


inline Foam::label Foam::token::lineNumber() const
{
    return lineNumber_;
}


inline Foam::label& Foam::token::lineNumber()
{
    return lineNumber_;
}


inline bool Foam::token::good() const
{
    return (type_ != tokenType::UNDEFINED && type_ != tokenType::ERROR);
}


inline bool Foam::token::undefined() const
{
    return (type_ == tokenType::UNDEFINED);
}


inline bool Foam::token::error() const
{
    return (type_ == tokenType::ERROR);
}


inline bool Foam::token::isFlag() const
{
    return (type_ == tokenType::FLAG);
}


inline int Foam::token::flagToken() const
{
    if (type_ == tokenType::FLAG)
    {
        return data_.flagVal;
    }

    parseError("flag bitmask");
    return NO_FLAG;
}


inline bool Foam::token::isPunctuation() const
{
    return (type_ == tokenType::PUNCTUATION);
}


inline Foam::token::punctuationToken Foam::token::pToken() const
{
    if (type_ == tokenType::PUNCTUATION)
    {
        return data_.punctuationVal;
    }

    parseError("punctuation character");
    return punctuationToken::NULL_TOKEN;
}


inline bool Foam::token::isSeparator() const
{
    return
    (
        type_ == tokenType::PUNCTUATION
     && isseparator(data_.punctuationVal)
    );
}


inline bool Foam::token::isLabel() const
{
    return (type_ == tokenType::LABEL);
}


inline Foam::label Foam::token::labelToken() const
{
    if (type_ == tokenType::LABEL)
    {
        return data_.labelVal;
    }

    parseError(pTraits<label>::typeName);
    return 0;
}


inline bool Foam::token::isFloatScalar() const
{
    return (type_ == tokenType::FLOAT_SCALAR);
}


inline Foam::floatScalar Foam::token::floatScalarToken() const
{
    if (type_ == tokenType::FLOAT_SCALAR)
    {
        return data_.floatVal;
    }

    parseError("floatScalar");
    return 0.0;
}


inline bool Foam::token::isDoubleScalar() const
{
    return (type_ == tokenType::DOUBLE_SCALAR);
}


inline Foam::doubleScalar Foam::token::doubleScalarToken() const
{
    if (type_ == tokenType::DOUBLE_SCALAR)
    {
        return data_.doubleVal;
    }

    parseError("doubleScalar");
    return 0.0;
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
        return data_.floatVal;
    }
    else if (type_ == tokenType::DOUBLE_SCALAR)
    {
        return data_.doubleVal;
    }

    parseError(pTraits<scalar>::typeName);
    return 0.0;
}


inline bool Foam::token::isNumber() const
{
    return (type_ == tokenType::LABEL || isScalar());
}


inline Foam::scalar Foam::token::number() const
{
    if (isLabel())
    {
        return labelToken();
    }
    if (isScalar())
    {
        return scalarToken();
    }

    parseError("number (label or scalar)");
    return 0.0;
}


inline bool Foam::token::isWord() const
{
    return (type_ == tokenType::WORD);
}


inline const Foam::word& Foam::token::wordToken() const
{
    if (type_ == tokenType::WORD)
    {
        return *data_.wordPtr;
    }

    parseError(word::typeName);
    return word::null;
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
        return *data_.stringPtr;
    }

    parseError(string::typeName);
    return string::null;
}


inline bool Foam::token::isCompound() const
{
    return (type_ == tokenType::COMPOUND);
}


inline const Foam::token::compound& Foam::token::compoundToken() const
{
    if (type_ == tokenType::COMPOUND)
    {
        return *data_.compoundPtr;
    }

    parseError("compound");
    return *data_.compoundPtr;  // This is questionable.
}


inline void Foam::token::setBad()
{
    clear();
    type_ = tokenType::ERROR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::token::operator=(const token& tok)
{
    clear();

    type_ = tok.type_;
    data_ = tok.data_;  // bit-wise copy of union content
    lineNumber_ = tok.lineNumber_;

    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        case tokenType::WORD:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
        }
        break;

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
        {
            data_.stringPtr = new string(*tok.data_.stringPtr);
        }
        break;

        case tokenType::COMPOUND:
        {
            // Identical pointers, but increase the refCount
            data_.compoundPtr = tok.data_.compoundPtr;
            data_.compoundPtr->refCount::operator++();
        }
        break;

        default:
            break;
    }
}


inline void Foam::token::operator=(token&& tok)
{
    clear();
    lineNumber_ = 0;
    swap(tok);
}


inline void Foam::token::operator=(const punctuationToken p)
{
    clear();
    type_ = tokenType::PUNCTUATION;
    data_.punctuationVal = p;
}


inline void Foam::token::operator=(const label val)
{
    clear();
    type_ = tokenType::LABEL;
    data_.labelVal = val;
}


inline void Foam::token::operator=(const floatScalar val)
{
    clear();
    type_ = tokenType::FLOAT_SCALAR;
    data_.floatVal = val;
}


inline void Foam::token::operator=(const doubleScalar val)
{
    clear();
    type_ = tokenType::DOUBLE_SCALAR;
    data_.doubleVal = val;
}


inline void Foam::token::operator=(word* wordPtr)
{
    clear();
    type_ = tokenType::WORD;
    data_.wordPtr = wordPtr;
}


inline void Foam::token::operator=(string* stringPtr)
{
    clear();
    type_ = tokenType::STRING;
    data_.stringPtr = stringPtr;
}


inline void Foam::token::operator=(const word& w)
{
    clear();
    type_ = tokenType::WORD;
    data_.wordPtr = new word(w);
}


inline void Foam::token::operator=(const string& str)
{
    clear();
    type_ = tokenType::STRING;
    data_.stringPtr = new string(str);
}


inline void Foam::token::operator=(word&& w)
{
    clear();
    type_ = tokenType::WORD;
    data_.wordPtr = new word;
    data_.wordPtr->swap(w);
}


inline void Foam::token::operator=(string&& s)
{
    clear();
    type_ = tokenType::STRING;
    data_.stringPtr = new string;
    data_.stringPtr->swap(s);
}


inline void Foam::token::operator=(Foam::token::compound* compoundPtr)
{
    clear();
    type_ = tokenType::COMPOUND;
    data_.compoundPtr = compoundPtr;
}


inline bool Foam::token::operator==(const token& tok) const
{
    if (type_ != tok.type_)
    {
        return false;
    }

    switch (type_)
    {
        case tokenType::UNDEFINED:
            return true;

        case tokenType::FLAG:
            return data_.flagVal == tok.data_.flagVal;

        case tokenType::PUNCTUATION:
            return data_.punctuationVal == tok.data_.punctuationVal;

        case tokenType::LABEL:
            return data_.labelVal == tok.data_.labelVal;

        case tokenType::FLOAT_SCALAR:
            return equal(data_.floatVal, tok.data_.floatVal);

        case tokenType::DOUBLE_SCALAR:
            return equal(data_.doubleVal, tok.data_.doubleVal);

        case tokenType::WORD:
            return *data_.wordPtr == *tok.data_.wordPtr;

        case tokenType::STRING:
        case tokenType::VARIABLE:
        case tokenType::VERBATIMSTRING:
            return *data_.stringPtr == *tok.data_.stringPtr;

        case tokenType::COMPOUND:
            return data_.compoundPtr == tok.data_.compoundPtr;

        case tokenType::ERROR:
            return true;
    }

    return false;
}


inline bool Foam::token::operator==(const punctuationToken p) const
{
    return (type_ == tokenType::PUNCTUATION && data_.punctuationVal == p);
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
     && data_.labelVal == val
    );
}


inline bool Foam::token::operator==(const floatScalar val) const
{
    return
    (
        type_ == tokenType::FLOAT_SCALAR
     && equal(data_.floatVal, val)
    );
}


inline bool Foam::token::operator==(const doubleScalar val) const
{
    return
    (
        type_ == tokenType::DOUBLE_SCALAR
     && equal(data_.doubleVal, val)
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


inline bool Foam::token::operator!=(const word& w) const
{
    return !operator==(w);
}


inline bool Foam::token::operator!=(const string& str) const
{
    return !operator==(str);
}


// ************************************************************************* //
