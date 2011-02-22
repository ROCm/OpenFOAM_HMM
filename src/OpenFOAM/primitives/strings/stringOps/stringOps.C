/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "stringOps.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::string Foam::stringOps::expand
(
    const string& original,
    const HashTable<string, word, string::hash>& mapping,
    const char sigil
)
{
    string str(original);
    return inplaceExpand(str, mapping);
}


Foam::string& Foam::stringOps::inplaceExpand
(
    string& s,
    const HashTable<string, word, string::hash>& mapping,
    const char sigil
)
{
    string::size_type begVar = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find(sigil, begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            // Find end of first occurrence
            string::size_type endVar = begVar;
            string::size_type delim = 0;

            if (s[begVar+1] == '{')
            {
                endVar = s.find('}', begVar);
                delim = 1;
            }
            else
            {
                string::iterator iter = s.begin() + begVar + 1;

                while
                (
                    iter != s.end()
                 && (isalnum(*iter) || *iter == '_')
                )
                {
                    ++iter;
                    ++endVar;
                }
            }

            if (endVar != string::npos && endVar != begVar)
            {
                string varName = s.substr
                (
                    begVar + 1 + delim,
                    endVar - begVar - 2*delim
                );

                HashTable<string, word, string::hash>::const_iterator fnd =
                    mapping.find(varName);

                if (fnd != HashTable<string, word, string::hash>::end())
                {
                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        *fnd
                    );
                    begVar += (*fnd).size();
                }
                else
                {
                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        ""
                    );
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            ++begVar;
        }
    }

    return s;
}


Foam::string Foam::stringOps::expandEnv
(
    const string& original,
    const bool recurse,
    const bool allowEmptyVar
)
{
    string str(original);
    return inplaceExpandEnv(str, recurse, allowEmptyVar);
}


// Expand all occurences of environment variables and initial tilde sequences
Foam::string& Foam::stringOps::inplaceExpandEnv
(
    string& s,
    const bool recurse,
    const bool allowEmptyVar
)
{
    string::size_type begVar = 0;

    // Expand $VARS
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find('$', begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            // Find end of first occurrence
            string::size_type endVar = begVar;
            string::size_type delim = 0;

            if (s[begVar+1] == '{')
            {
                endVar = s.find('}', begVar);
                delim = 1;
            }
            else
            {
                string::iterator iter = s.begin() + begVar + 1;

                while
                (
                    iter != s.end()
                 && (isalnum(*iter) || *iter == '_')
                )
                {
                    ++iter;
                    ++endVar;
                }
            }

            if (endVar != string::npos && endVar != begVar)
            {
                string varName = s.substr
                (
                    begVar + 1 + delim,
                    endVar - begVar - 2*delim
                );

                string varValue = getEnv(varName);

                if (varValue.size())
                {
                    if (recurse)
                    {
                        varValue.expand(recurse, allowEmptyVar);
                    }
                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        varValue
                    );
                    begVar += varValue.size();
                }
                else if (allowEmptyVar)
                {
                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        ""
                    );
                }
                else
                {
                    FatalErrorIn("string::expand(const bool, const bool)")
                        << "Unknown variable name " << varName << '.'
                        << exit(FatalError);
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            ++begVar;
        }
    }

    if (!s.empty())
    {
        if (s[0] == '~')
        {
            // Expand initial ~
            //   ~/        => home directory
            //   ~OpenFOAM => site/user OpenFOAM configuration directory
            //   ~user     => home directory for specified user

            word user;
            fileName file;

            if ((begVar = s.find('/')) != string::npos)
            {
                user = s.substr(1, begVar - 1);
                file = s.substr(begVar + 1);
            }
            else
            {
                user = s.substr(1);
            }

            // NB: be a bit lazy and expand ~unknownUser as an
            // empty string rather than leaving it untouched.
            // otherwise add extra test
            if (user == "OpenFOAM")
            {
                s = findEtcFile(file);
            }
            else
            {
                s = home(user)/file;
            }
        }
        else if (s[0] == '.')
        {
            // Expand a lone '.' and an initial './' into cwd
            if (s.size() == 1)
            {
                s = cwd();
            }
            else if (s[1] == '/')
            {
                s.std::string::replace(0, 1, cwd());
            }
        }
    }

    return s;
}


Foam::string Foam::stringOps::trimLeft(const string& s)
{
    if (!s.empty())
    {
        string::size_type beg = 0;
        while (isspace(s[beg]))
        {
            ++beg;
        }

        if (beg)
        {
            return s.substr(beg);
        }
    }

    return s;
}


Foam::string& Foam::stringOps::inplaceTrimLeft(string& s)
{
    if (!s.empty())
    {
        string::size_type beg = 0;
        while (isspace(s[beg]))
        {
            ++beg;
        }

        if (beg)
        {
            s.erase(0, beg);
        }
    }

    return s;
}


Foam::string Foam::stringOps::trimRight(const string& s)
{
    notImplemented("string stringOps::trimRight(const string&)");
    return s;
}


Foam::string& Foam::stringOps::inplaceTrimRight(string& s)
{
    notImplemented("string& stringOps::inplaceTrimRight(string&)");
    return s;
}


Foam::string Foam::stringOps::trim(const string& original)
{
    notImplemented("string stringOps::trim(const string&)");
    return original;
}


Foam::string& Foam::stringOps::inplaceTrim(string& s)
{
    notImplemented("string& stringOps::inplaceTrim(string&)");
    return s;
}

// ************************************************************************* //
