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

\*---------------------------------------------------------------------------*/

#include <algorithm>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::token Foam::token::boolean(bool on) noexcept
{
    token tok;
    tok.type_ = tokenType::BOOL;
    tok.data_.labelVal = on;

    return tok;
}


inline Foam::token Foam::token::flag(int bitmask) noexcept
{
    token tok;
    tok.type_ = tokenType::FLAG;
    tok.data_.flagVal = bitmask;

    return tok;
}


inline bool Foam::token::isseparator(int c) noexcept
{
    // NOTE: keep synchronized with ISstream::read(token&)

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
        case token::PLUS :
        // Excluded token::MINUS since it could start a number
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

inline void Foam::token::setUndefined() noexcept
{
    type_ = tokenType::UNDEFINED;
    data_.int64Val = 0; // bit-wise zero for union content
    // leave lineNumber untouched - may still be needed
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::token::token() noexcept
:
    data_(), // bit-wise zero for union content
    type_(tokenType::UNDEFINED),
    line_(0)
{}


inline Foam::token::token(const token& tok)
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    line_(tok.line_)
{
    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
            break;
        }

        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
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


inline Foam::token::token(token&& tok) noexcept
:
    data_(tok.data_), // bit-wise copy of union content
    type_(tok.type_),
    line_(tok.line_)
{
    tok.setUndefined(); // zero the union content without any checking
    tok.line_ = 0;
}


inline Foam::token::token(punctuationToken p, label lineNum) noexcept
:
    data_(),
    type_(tokenType::PUNCTUATION),
    line_(lineNum)
{
    data_.punctuationVal = p;
}


inline Foam::token::token(const label val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::LABEL),
    line_(lineNum)
{
    data_.labelVal = val;
}


inline Foam::token::token(const floatScalar val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::FLOAT),
    line_(lineNum)
{
    data_.floatVal = val;
}


inline Foam::token::token(const doubleScalar val, label lineNum) noexcept
:
    data_(),
    type_(tokenType::DOUBLE),
    line_(lineNum)
{
    data_.doubleVal = val;
}


inline Foam::token::token(const word& w, label lineNum)
:
    data_(),
    type_(tokenType::WORD),
    line_(lineNum)
{
    data_.wordPtr = new word(w);
}


inline Foam::token::token(const string& str, label lineNum)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    data_.stringPtr = new string(str);
}


inline Foam::token::token(word&& w, label lineNum)
:
    data_(),
    type_(tokenType::WORD),
    line_(lineNum)
{
    data_.wordPtr = new word(std::move(w));
}


inline Foam::token::token(string&& str, label lineNum)
:
    data_(),
    type_(tokenType::STRING),
    line_(lineNum)
{
    data_.stringPtr = new string(std::move(str));
}


inline Foam::token::token(token::compound* ptr, label lineNum)
:
    data_(),
    type_(tokenType::COMPOUND),
    line_(lineNum)
{
    data_.compoundPtr = ptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::token::~token()
{
    reset();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::token::reset()
{
    switch (type_)
    {
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            delete data_.wordPtr;
            break;
        }

        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
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


inline void Foam::token::swap(token& tok)
{
    if (this == &tok)
    {
        return;  // Self-swap is a no-op
    }

    std::swap(data_, tok.data_);
    std::swap(type_, tok.type_);
    std::swap(line_, tok.line_);
}


inline Foam::token::tokenType Foam::token::type() const noexcept
{
    return type_;
}


inline bool Foam::token::setType(token::tokenType tokType) noexcept
{
    if (type_ == tokType)
    {
        // No change required
        return true;
    }

    switch (tokType)
    {
        case tokenType::BOOL:
        case tokenType::LABEL:
        {
            switch (type_)
            {
                case tokenType::BOOL:
                case tokenType::LABEL:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            switch (type_)
            {
                case tokenType::WORD:
                case tokenType::DIRECTIVE:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
        {
            switch (type_)
            {
                // could also go from WORD to STRING etc - to be decided
                case tokenType::STRING:
                case tokenType::EXPRESSION:
                case tokenType::VARIABLE:
                case tokenType::VERBATIM:
                    type_ = tokType;
                    return true;
                    break;

                default:
                    break;
            }
        }
        break;

        default:
            break;
    }

    return false;
}


inline Foam::label Foam::token::lineNumber() const noexcept
{
    return line_;
}


inline Foam::label Foam::token::lineNumber(const label lineNum) noexcept
{
    label old(line_);
    line_ = lineNum;
    return old;
}


inline bool Foam::token::good() const noexcept
{
    return (type_ != tokenType::UNDEFINED && type_ != tokenType::ERROR);
}


inline bool Foam::token::undefined() const noexcept
{
    return (type_ == tokenType::UNDEFINED);
}


inline bool Foam::token::error() const noexcept
{
    return (type_ == tokenType::ERROR);
}


inline bool Foam::token::isBool() const noexcept
{
    return (type_ == tokenType::BOOL);
}


inline bool Foam::token::boolToken() const
{
    if (type_ == tokenType::BOOL || type_ == tokenType::LABEL)
    {
        return data_.labelVal;
    }

    parseError("bool");
    return false;
}


inline bool Foam::token::isFlag() const noexcept
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


inline bool Foam::token::isPunctuation() const noexcept
{
    return (type_ == tokenType::PUNCTUATION);
}


inline bool Foam::token::isPunctuation(const punctuationToken p) const noexcept
{
    return
    (
        type_ == tokenType::PUNCTUATION
     && data_.punctuationVal == p
    );
}


inline bool Foam::token::isSeparator() const noexcept
{
    return
    (
        type_ == tokenType::PUNCTUATION
     && isseparator(data_.punctuationVal)
    );
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


inline bool Foam::token::isLabel() const noexcept
{
    return (type_ == tokenType::LABEL);
}


inline bool Foam::token::isLabel(const label val) const noexcept
{
    return
    (
        type_ == tokenType::LABEL
     && data_.labelVal == val
    );
}


inline Foam::label Foam::token::labelToken() const
{
    if (type_ == tokenType::LABEL)
    {
        return data_.labelVal;
    }

    parseError("label");
    return 0;
}


inline bool Foam::token::isFloat() const noexcept
{
    return (type_ == tokenType::FLOAT);
}


inline Foam::floatScalar Foam::token::floatToken() const
{
    if (type_ == tokenType::FLOAT)
    {
        return data_.floatVal;
    }

    parseError("float");
    return 0;
}


inline bool Foam::token::isDouble() const noexcept
{
    return (type_ == tokenType::DOUBLE);
}


inline Foam::doubleScalar Foam::token::doubleToken() const
{
    if (type_ == tokenType::DOUBLE)
    {
        return data_.doubleVal;
    }

    parseError("double");
    return 0;
}


inline bool Foam::token::isScalar() const noexcept
{
    return
    (
        type_ == tokenType::FLOAT
     || type_ == tokenType::DOUBLE
    );
}


inline Foam::scalar Foam::token::scalarToken() const
{
    if (type_ == tokenType::FLOAT)
    {
        return data_.floatVal;
    }
    else if (type_ == tokenType::DOUBLE)
    {
        return data_.doubleVal;
    }

    parseError("scalar");
    return 0;
}


inline bool Foam::token::isNumber() const noexcept
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
    return 0;
}


inline bool Foam::token::isWord() const noexcept
{
    return
    (
        type_ == tokenType::WORD
     || type_ == tokenType::DIRECTIVE
    );
}


inline bool Foam::token::isWord(const std::string& s) const
{
    return (isWord() && s == *data_.wordPtr);
}


inline bool Foam::token::isDirective() const noexcept
{
    return (type_ == tokenType::DIRECTIVE);
}


inline const Foam::word& Foam::token::wordToken() const
{
    if
    (
        type_ == tokenType::WORD
     || type_ == tokenType::DIRECTIVE
    )
    {
        return *data_.wordPtr;
    }

    parseError("word");
    return word::null;
}


inline bool Foam::token::isQuotedString() const noexcept
{
    return (type_ == tokenType::STRING);
}


inline bool Foam::token::isString() const noexcept
{
    return
    (
        type_ == tokenType::STRING
     || type_ == tokenType::EXPRESSION
     || type_ == tokenType::VARIABLE
     || type_ == tokenType::VERBATIM
    );
}


inline bool Foam::token::isExpression() const noexcept
{
    return (type_ == tokenType::EXPRESSION);
}


inline bool Foam::token::isVariable() const noexcept
{
    return (type_ == tokenType::VARIABLE);
}


inline bool Foam::token::isVerbatim() const noexcept
{
    return (type_ == tokenType::VERBATIM);
}


inline bool Foam::token::isStringType() const noexcept
{
    return (isWord() || isString());
}


inline const Foam::string& Foam::token::stringToken() const
{
    if
    (
        type_ == tokenType::STRING
     || type_ == tokenType::EXPRESSION
     || type_ == tokenType::VARIABLE
     || type_ == tokenType::VERBATIM
    )
    {
        return *data_.stringPtr;
    }
    else if
    (
        type_ == tokenType::WORD
     || type_ == tokenType::DIRECTIVE
    )
    {
        // Foam::word derives from Foam::string, no need to cast.
        return *data_.wordPtr;
    }

    parseError("string");
    return string::null;
}


inline bool Foam::token::isCompound() const noexcept
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
    reset();
    type_ = tokenType::ERROR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::token::operator=(const token& tok)
{
    if (this == &tok)
    {
        return;  // Self-assignment is a no-op
    }

    reset();

    type_ = tok.type_;
    data_ = tok.data_;  // bit-wise copy of union content
    line_ = tok.line_;

    // Fundamental: values already handled by bit-wise copy
    // Pointer: duplicate content or increase refCount

    switch (type_)
    {
        case tokenType::WORD:
        case tokenType::DIRECTIVE:
        {
            data_.wordPtr = new word(*tok.data_.wordPtr);
        }
        break;

        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
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
    if (this == &tok)
    {
        return;  // Self-assignment is a no-op
    }

    reset();
    line_ = 0;
    swap(tok);
}


inline void Foam::token::operator=(const punctuationToken p)
{
    reset();
    type_ = tokenType::PUNCTUATION;
    data_.punctuationVal = p;
}


inline void Foam::token::operator=(const label val)
{
    reset();
    type_ = tokenType::LABEL;
    data_.labelVal = val;
}


inline void Foam::token::operator=(const floatScalar val)
{
    reset();
    type_ = tokenType::FLOAT;
    data_.floatVal = val;
}


inline void Foam::token::operator=(const doubleScalar val)
{
    reset();
    type_ = tokenType::DOUBLE;
    data_.doubleVal = val;
}


inline void Foam::token::operator=(const word& w)
{
    reset();
    type_ = tokenType::WORD;
    data_.wordPtr = new word(w);
}


inline void Foam::token::operator=(const string& str)
{
    reset();
    type_ = tokenType::STRING;
    data_.stringPtr = new string(str);
}


inline void Foam::token::operator=(word&& w)
{
    reset();
    type_ = tokenType::WORD;
    data_.wordPtr = new word(std::move(w));
}


inline void Foam::token::operator=(string&& s)
{
    reset();
    type_ = tokenType::STRING;
    data_.stringPtr = new string(std::move(s));
}


inline void Foam::token::operator=(Foam::token::compound* ptr)
{
    reset();
    type_ = tokenType::COMPOUND;
    data_.compoundPtr = ptr;
}


inline void Foam::token::operator=(autoPtr<token::compound>&& ptr)
{
    reset();
    type_ = tokenType::COMPOUND;
    data_.compoundPtr = ptr.release();
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

        case tokenType::BOOL:
            return data_.labelVal == tok.data_.labelVal;

        case tokenType::FLAG:
            return data_.flagVal == tok.data_.flagVal;

        case tokenType::PUNCTUATION:
            return data_.punctuationVal == tok.data_.punctuationVal;

        case tokenType::LABEL:
            return data_.labelVal == tok.data_.labelVal;

        case tokenType::FLOAT:
            return equal(data_.floatVal, tok.data_.floatVal);

        case tokenType::DOUBLE:
            return equal(data_.doubleVal, tok.data_.doubleVal);

        case tokenType::WORD:
        case tokenType::DIRECTIVE:
            return *data_.wordPtr == *tok.data_.wordPtr;

        case tokenType::STRING:
        case tokenType::EXPRESSION:
        case tokenType::VARIABLE:
        case tokenType::VERBATIM:
            return *data_.stringPtr == *tok.data_.stringPtr;

        case tokenType::COMPOUND:
            return data_.compoundPtr == tok.data_.compoundPtr;

        case tokenType::ERROR:
            return true;
    }

    return false;
}


inline bool Foam::token::operator==(const punctuationToken p) const noexcept
{
    return isPunctuation(p);
}


inline bool Foam::token::operator==(const std::string& s) const
{
    return
    (
        isWord()
      ? s == *data_.wordPtr
      : isString() && s == *data_.stringPtr
    );
}


inline bool Foam::token::operator==(const label val) const noexcept
{
    return isLabel(val);
}


inline bool Foam::token::operator==(const floatScalar val) const noexcept
{
    return
    (
        type_ == tokenType::FLOAT
     && equal(data_.floatVal, val)
    );
}


inline bool Foam::token::operator==(const doubleScalar val) const noexcept
{
    return
    (
        type_ == tokenType::DOUBLE
     && equal(data_.doubleVal, val)
    );
}


inline bool Foam::token::operator!=(const token& tok) const
{
    return !operator==(tok);
}


inline bool Foam::token::operator!=(const punctuationToken p) const noexcept
{
    return !isPunctuation(p);
}


inline bool Foam::token::operator!=(const label val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const floatScalar val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const doubleScalar val) const noexcept
{
    return !operator==(val);
}


inline bool Foam::token::operator!=(const std::string& s) const
{
    return !operator==(s);
}


// ************************************************************************* //
