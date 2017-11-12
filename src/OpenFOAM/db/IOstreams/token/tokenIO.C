/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "error.H"
#include "token.H"
#include "scalar.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
template<class OS>
static OS& printTokenInfo(OS& os, const token& tok)
{
    os  << "on line " << tok.lineNumber() << ": ";

    switch (tok.type())
    {
        case token::tokenType::UNDEFINED:
            os  << "undefined token";
        break;

        case token::tokenType::FLAG:
            os  << "flag '" << int(tok.flagToken()) << '\'';
        break;

        case token::tokenType::PUNCTUATION:
            os  << "punctuation '" << tok.pToken() << '\'';
        break;

        case token::tokenType::LABEL:
            os  << "label " << tok.labelToken();
        break;

        case token::tokenType::FLOAT_SCALAR:
            os  << "float " << tok.floatScalarToken();
        break;

        case token::tokenType::DOUBLE_SCALAR:
            os  << "double " << tok.doubleScalarToken();
        break;

        case token::tokenType::WORD:
            os  << "word '" << tok.wordToken() << '\'';
        break;

        case token::tokenType::STRING:
            os  << "string " << tok.stringToken();
        break;

        case token::tokenType::VARIABLE:
            os  << "variable " << tok.stringToken();
        break;

        case token::tokenType::VERBATIMSTRING:
            os  << "verbatim string " << tok.stringToken();
        break;

        case token::tokenType::COMPOUND:
        {
            if (tok.compoundToken().empty())
            {
                os  << "empty ";
            }
            os  << "compound of type "
                << tok.compoundToken().type();
        }
        break;

        case token::tokenType::ERROR:
            os  << "error";
        break;

        default:
            os  << "unknown token type '" << int(tok.type()) << '\'';
            break;
    }

    return os;
}
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::token::token(Istream& is)
:
    token()
{
    is.read(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::token::name() const
{
    switch (type_)
    {
        case token::tokenType::UNDEFINED: return "undefined";
        case token::tokenType::FLAG: return "flag";
        case token::tokenType::PUNCTUATION: return "punctuation";
        case token::tokenType::LABEL: return "label";
        case token::tokenType::FLOAT_SCALAR: return "float";
        case token::tokenType::DOUBLE_SCALAR: return "double";
        case token::tokenType::WORD: return "word";
        case token::tokenType::STRING: return "string";
        case token::tokenType::VERBATIMSTRING: return "verbatim";
        case token::tokenType::VARIABLE: return "variable";
        case token::tokenType::COMPOUND: return "compound";
        case token::tokenType::ERROR: return "error";

        default:
            break;
    }

    return "unknown(" + std::to_string(int(type_)) + ")";
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, token& tok)
{
    tok.clear();
    return is.read(tok);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token& tok)
{
    switch (tok.type_)
    {
        case token::tokenType::UNDEFINED:
            os << "UNDEFINED";
            WarningInFunction
                << "Undefined token" << endl;
        break;

        case token::tokenType::FLAG:
            // Swallow the flag
        break;

        case token::tokenType::PUNCTUATION:
            os << tok.data_.punctuationVal;
        break;

        case token::tokenType::LABEL:
            os << tok.data_.labelVal;
        break;

        case token::tokenType::FLOAT_SCALAR:
            os << tok.data_.floatVal;
        break;

        case token::tokenType::DOUBLE_SCALAR:
            os << tok.data_.doubleVal;
        break;

        case token::tokenType::WORD:
            os << *tok.data_.wordPtr;
        break;

        case token::tokenType::STRING:
        case token::tokenType::VERBATIMSTRING:
            os << *tok.data_.stringPtr;
        break;

        case token::tokenType::VARIABLE:
            // Behaviour differs according to stream type
            os.write(tok);
        break;

        case token::tokenType::COMPOUND:
            os << *tok.data_.compoundPtr;
        break;

        case token::tokenType::ERROR:
            os << "ERROR";
            WarningInFunction
                << "Error token" << endl;
        break;

        default:
            os << "UNKNOWN";
            SeriousErrorInFunction
                << "Unknown token"
                << endl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


ostream& Foam::operator<<(ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::punctuationToken& pt)
{
    return os << char(pt);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token::compound& ct)
{
    os << ct.type() << token::SPACE;
    ct.write(os);

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ostream& Foam::operator<<(ostream& os, const InfoProxy<token>& ip)
{
    return printTokenInfo(os, ip.t_);
}


template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<token>& ip)
{
    return printTokenInfo(os, ip.t_);
}


// ************************************************************************* //
