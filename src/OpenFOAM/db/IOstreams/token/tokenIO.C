/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "error.H"
#include "token.H"

#include "IOstreams.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::token::token(Istream& is)
:
    type_(tokenType::UNDEFINED)
{
    is.read(*this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, token& t)
{
    t.clear();
    return is.read(t);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const token& t)
{
    switch (t.type_)
    {
        case token::tokenType::UNDEFINED:
            os << "UNDEFINED";
            WarningInFunction
                << "Undefined token" << endl;
        break;

        case token::tokenType::PUNCTUATION:
            os << t.punctuationToken_;
        break;

        case token::tokenType::WORD:
            os << *t.wordTokenPtr_;
        break;

        case token::tokenType::STRING:
        case token::tokenType::VERBATIMSTRING:
            os << *t.stringTokenPtr_;
        break;

        case token::tokenType::VARIABLE:
            // Behaviour differs according to stream type
            os.write(t);
        break;

        case token::tokenType::LABEL:
            os << t.labelToken_;
        break;

        case token::tokenType::FLOAT_SCALAR:
            os << t.floatScalarToken_;
        break;

        case token::tokenType::DOUBLE_SCALAR:
            os << t.doubleScalarToken_;
        break;

        case token::tokenType::COMPOUND:
            os << *t.compoundTokenPtr_;
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
    const token& t = ip.t_;

    os  << "on line " << t.lineNumber();

    switch (t.type())
    {
        case token::tokenType::UNDEFINED:
            os  << " an undefined token";
        break;

        case token::tokenType::PUNCTUATION:
            os  << " the punctuation token " << '\'' << t.pToken() << '\'';
        break;

        case token::tokenType::WORD:
            os  << " the word " << '\'' << t.wordToken() << '\'';
        break;

        case token::tokenType::STRING:
            os  << " the string " << t.stringToken();
        break;

        case token::tokenType::VARIABLE:
            os  << " the variable " << t.stringToken();
        break;

        case token::tokenType::VERBATIMSTRING:
            os  << " the verbatim string " << t.stringToken();
        break;

        case token::tokenType::LABEL:
            os  << " the label " << t.labelToken();
        break;

        case token::tokenType::FLOAT_SCALAR:
            os  << " the floatScalar " << t.floatScalarToken();
        break;

        case token::tokenType::DOUBLE_SCALAR:
            os  << " the doubleScalar " << t.doubleScalarToken();
        break;

        case token::tokenType::COMPOUND:
        {
            if (t.compoundToken().empty())
            {
                os  << " the empty compound of type "
                    << t.compoundToken().type();
            }
            else
            {
                os  << " the compound of type "
                    << t.compoundToken().type();
            }
        }
        break;

        case token::tokenType::ERROR:
            os  << " an error";
        break;

        default:
            os  << " an unknown token type " << '\'' << int(t.type()) << '\'';
    }

    return os;
}


template<>
Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<token>& ip)
{
    const token& t = ip.t_;

    os  << "on line " << t.lineNumber();

    switch (t.type())
    {
        case token::tokenType::UNDEFINED:
            os  << " an undefined token";
        break;

        case token::tokenType::PUNCTUATION:
            os  << " the punctuation token " << '\'' << t.pToken() << '\'';
        break;

        case token::tokenType::WORD:
            os  << " the word " << '\'' << t.wordToken() << '\'';
        break;

        case token::tokenType::STRING:
            os  << " the string " << t.stringToken();
        break;

        case token::tokenType::VARIABLE:
            os  << " the variable " << t.stringToken();
        break;

        case token::tokenType::VERBATIMSTRING:
            os  << " the verbatim string " << t.stringToken();
        break;

        case token::tokenType::LABEL:
            os  << " the label " << t.labelToken();
        break;

        case token::tokenType::FLOAT_SCALAR:
            os  << " the floatScalar " << t.floatScalarToken();
        break;

        case token::tokenType::DOUBLE_SCALAR:
            os  << " the doubleScalar " << t.doubleScalarToken();
        break;

        case token::tokenType::COMPOUND:
        {
            if (t.compoundToken().empty())
            {
                os  << " the empty compound of type "
                    << t.compoundToken().type();
            }
            else
            {
                os  << " the compound of type "
                    << t.compoundToken().type();
            }
        }
        break;

        case token::tokenType::ERROR:
            os  << " an error";
        break;

        default:
            os  << " an unknown token type "  << '\'' << int(t.type()) << '\'';
    }

    return os;
}


// ************************************************************************* //
