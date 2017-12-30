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

#include "stringOps.H"
#include "typeInfo.H"
#include "OSspecific.H"
#include "etcFiles.H"
#include "StringStream.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
// Standard handling of "~/", "./" etc.
static void standardExpansions(std::string& s)
{
    if (s.empty())
    {
        return;
    }

    if (s[0] == '.')
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
    else if (s[0] == '~')
    {
        // Expand initial ~
        //   ~/        => home directory
        //   ~OpenFOAM => site/user OpenFOAM configuration directory
        //   ~user     => home directory for specified user

        string user;
        fileName file;

        const auto slash = s.find('/');
        if (slash == std::string::npos)
        {
            user = s.substr(1);
        }
        else
        {
            user = s.substr(1, slash - 1);
            file = s.substr(slash + 1);
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
}
}


//! \cond fileScope
//  Find the type/position of the ":-" or ":+" alternative values
//  Returns 0, '-', '+' corresponding to not-found or ':-' or ':+'
static inline int findParameterAlternative
(
    const std::string& s,
    std::string::size_type& pos,
    std::string::size_type endPos
)
{
    while (pos != std::string::npos)
    {
        pos = s.find(':', pos);
        if (pos != std::string::npos)
        {
            if (pos < endPos)
            {
                // in-range: check for '+' or '-' following the ':'
                const int altType = s[pos+1];
                if (altType == '+' || altType == '-')
                {
                    return altType;
                }

                ++pos;    // unknown/unsupported - continue at next position
            }
            else
            {
                // out-of-range: abort
                pos = std::string::npos;
            }
        }
    }

    return 0;
}
//! \endcond


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::string::size_type Foam::stringOps::count
(
    const std::string& str,
    const char c
)
{
    std::string::size_type n = 0;

    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        if (*iter == c)
        {
            ++n;
        }
    }

    return n;
}


std::string::size_type Foam::stringOps::count(const char* str, const char c)
{
    if (!str)
    {
        return 0;
    }

    std::string::size_type n = 0;

    for (const char *iter = str; *iter; ++iter)
    {
        if (*iter == c)
        {
            ++n;
        }
    }

    return n;
}


Foam::string Foam::stringOps::expand
(
    const string& original,
    const HashTable<string, word, string::hash>& mapping,
    const char sigil
)
{
    string s(original);
    inplaceExpand(s, mapping);
    return s;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const HashTable<string, word, string::hash>& mapping,
    const char sigil
)
{
    std::string::size_type varBeg = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            // Find end of first occurrence
            std::string::size_type varEnd = varBeg;
            std::string::size_type delim = 0;

            // The type/position of the ":-" or ":+" alternative values
            int altType = 0;
            auto altPos = std::string::npos;

            if (s[varBeg+1] == '{')
            {
                varEnd = s.find('}', varBeg);
                delim = 1;

                // check for ${parameter:-word} or ${parameter:+word}
                if (varEnd != std::string::npos)
                {
                    altPos = varBeg;
                    altType = findParameterAlternative(s, altPos, varEnd);
                }
            }
            else
            {
                string::const_iterator iter = s.cbegin() + varBeg + 1;

                // more generous in accepting keywords than for env variables
                while
                (
                    iter != s.cend()
                 &&
                    (
                        std::isalnum(*iter)
                     || *iter == '.'
                     || *iter == ':'
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++varEnd;
                }
            }

            if (varEnd == std::string::npos)
            {
                // likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (varEnd == varBeg)
            {
                // parsed '${}' or $badChar  - skip over
                varBeg = varEnd + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        varBeg + 1 + delim,
                        (
                            (altPos == std::string::npos ? varEnd : altPos)
                          - varBeg - 2*delim
                        )
                    ),
                    false
                );

                std::string altValue;
                if (altPos != std::string::npos)
                {
                    // had ":-" or ":+" alternative value
                    altValue = s.substr
                    (
                        altPos + 2,
                        varEnd - altPos - 2*delim
                    );
                }


                auto fnd = mapping.cfind(varName);

                if (fnd.found())
                {
                    if (altPos != std::string::npos && altType == '+')
                    {
                        // was found, use ":+" alternative
                        s.std::string::replace
                        (
                            varBeg,
                            varEnd - varBeg + 1,
                            altValue
                        );
                        varBeg += altValue.size();
                    }
                    else
                    {
                        // was found, use value
                        s.std::string::replace
                        (
                            varBeg,
                            varEnd - varBeg + 1,
                            *fnd
                        );
                        varBeg += (*fnd).size();
                    }
                }
                else if (altPos != std::string::npos && altType == '-')
                {
                    // was not found, use ":-" alternative
                    s.std::string::replace
                    (
                        varBeg,
                        varEnd - varBeg + 1,
                        altValue
                    );
                    varBeg += altValue.size();
                }
                else
                {
                    // substitute with nothing, also for ":+" alternative
                    s.std::string::erase(varBeg, varEnd - varBeg + 1);
                }
            }
        }
        else
        {
            ++varBeg;
        }
    }
}


Foam::string Foam::stringOps::expand
(
    const string& original,
    const dictionary& dict,
    const char sigil
)
{
    string s(original);
    inplaceExpand(s, dict, sigil);
    return s;
}


Foam::string Foam::stringOps::getVariable
(
    const word& name,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty
)
{
    string value;

    const entry* eptr = dict.lookupScopedEntryPtr(name, true, false);

    if (eptr)
    {
        OStringStream buf;
        // Force floating point numbers to be printed with at least
        // some decimal digits.
        buf << fixed;
        buf.precision(IOstream::defaultPrecision());

        // Fails for non-primitiveEntry
        dynamicCast<const primitiveEntry>(*eptr).write(buf, true);

        value = buf.str();
    }
    else if (allowEnvVars)
    {
        value = getEnv(name);

        if (value.empty() && !name.empty())
        {
            // The type/position of the ":-" or ":+" alternative values
            std::string::size_type altPos = 0;

            // check for parameter:-word or parameter:+word
            const int altType =
                findParameterAlternative(name, altPos, name.size()-1);

            if (altType)
            {
                value = getEnv
                (
                    // var-name
                    word(name.substr(0, altPos), false)
                );

                // ":-" or ":+" alternative value
                if (value.empty() ? (altType == '-') : (altType == '+'))
                {
                    // alternative
                    value = name.substr(altPos + 2);
                }
            }
        }
    }

    if (!allowEmpty && value.empty())
    {
        if (allowEnvVars)
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Cannot find dictionary or environment variable "
                << name << exit(FatalIOError);
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Cannot find dictionary variable "
                << name << exit(FatalIOError);
        }
    }

    return value;
}


Foam::string Foam::stringOps::expand
(
    const string& s,
    std::string::size_type& index,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty
)
{
    string out;

    while (index < s.size())
    {
        if (s[index] == '$' && s[index+1] == '{')
        {
            // Recurse to parse variable name
            index += 2;
            string val = expand(s, index, dict, allowEnvVars, allowEmpty);
            out.append(val);
        }
        else if (s[index] == '}')
        {
            return getVariable(out, dict, allowEnvVars, allowEmpty);
        }
        else
        {
            out.append(string(s[index]));
        }
        ++index;
    }

    return out;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty,
    const char sigil
)
{
    std::string::size_type varBeg = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            if (s[varBeg+1] == '{')
            {
                // Recursive variable expansion mode
                auto stringStart = varBeg;
                varBeg += 2;
                string varValue
                (
                    expand
                    (
                        s,
                        varBeg,
                        dict,
                        allowEnvVars,
                        allowEmpty
                    )
                );

                s.std::string::replace
                (
                    stringStart,
                    varBeg - stringStart + 1,
                    varValue
                );

                varBeg = stringStart+varValue.size();
            }
            else
            {
                std::string::const_iterator iter = s.cbegin() + varBeg + 1;
                std::string::size_type varEnd = varBeg;

                // more generous in accepting keywords than for env variables
                while
                (
                    iter != s.cend()
                 &&
                    (
                        std::isalnum(*iter)
                     || *iter == '.'
                     || *iter == ':'
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++varEnd;
                }

                const word varName
                (
                    s.substr
                    (
                        varBeg + 1,
                        varEnd - varBeg
                    ),
                    false
                );

                string varValue
                (
                    getVariable
                    (
                        varName,
                        dict,
                        allowEnvVars,
                        allowEmpty
                    )
                );

                s.std::string::replace
                (
                    varBeg,
                    varName.size()+1,
                    varValue
                );
                varBeg += varValue.size();
            }
        }
        else
        {
            ++varBeg;
        }
    }

    // Standard handling of "~/", "./" etc.
    standardExpansions(s);
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const dictionary& dict,
    const char sigil
)
{
    std::string::size_type varBeg = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            // Find end of first occurrence
            std::string::size_type varEnd = varBeg;
            std::string::size_type delim = 0;

            if (s[varBeg+1] == '{')
            {
                varEnd = s.find('}', varBeg);
                delim = 1;
            }
            else
            {
                string::const_iterator iter = s.cbegin() + varBeg + 1;

                // more generous in accepting keywords than for env variables
                while
                (
                    iter != s.cend()
                 &&
                    (
                        std::isalnum(*iter)
                     || *iter == '.'
                     || *iter == ':'
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++varEnd;
                }
            }

            if (varEnd == std::string::npos)
            {
                // likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (varEnd == varBeg)
            {
                // parsed '${}' or $badChar  - skip over
                varBeg = varEnd + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        varBeg + 1 + delim,
                        varEnd - varBeg - 2*delim
                    ),
                    false
                );


                // Lookup in the dictionary without wildcards.
                // See note in primitiveEntry
                const entry* eptr = dict.lookupScopedEntryPtr
                (
                    varName,
                    true,
                    false
                );

                // if defined - copy its entries
                if (eptr)
                {
                    OStringStream buf;
                    // Force floating point numbers to be printed with at least
                    // some decimal digits.
                    buf << fixed;
                    buf.precision(IOstream::defaultPrecision());
                    if (eptr->isDict())
                    {
                        eptr->dict().write(buf, false);
                    }
                    else
                    {
                        // Fail for non-primitiveEntry
                        dynamicCast<const primitiveEntry>
                        (
                            *eptr
                        ).write(buf, true);
                    }

                    s.std::string::replace
                    (
                        varBeg,
                        varEnd - varBeg + 1,
                        buf.str()
                    );
                    varBeg += buf.str().size();
                }
                else
                {
                    // not defined - leave original string untouched
                    varBeg = varEnd + 1;
                }
            }
        }
        else
        {
            ++varBeg;
        }
    }
}


Foam::string Foam::stringOps::expand
(
    const string& original,
    const bool allowEmpty
)
{
    string s(original);
    inplaceExpand(s, allowEmpty);
    return s;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const bool allowEmpty
)
{
    std::string::size_type varBeg = 0;

    // Expand $VARS
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find('$', varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            // Find end of first occurrence
            std::string::size_type varEnd = varBeg;
            std::string::size_type delim = 0;

            // The type/position of the ":-" or ":+" alternative values
            int altType = 0;
            std::string::size_type altPos = std::string::npos;

            if (s[varBeg+1] == '{')
            {
                varEnd = s.find('}', varBeg);
                delim = 1;

                // check for ${parameter:-word} or ${parameter:+word}
                if (varEnd != std::string::npos)
                {
                    altPos = varBeg;
                    altType = findParameterAlternative(s, altPos, varEnd);
                }
            }
            else
            {
                string::const_iterator iter = s.cbegin() + varBeg + 1;

                while
                (
                    iter != s.cend()
                 && (std::isalnum(*iter) || *iter == '_')
                )
                {
                    ++iter;
                    ++varEnd;
                }
            }


            if (varEnd == std::string::npos)
            {
                // likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (varEnd == varBeg)
            {
                // parsed '${}' or $badChar  - skip over
                varBeg = varEnd + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        varBeg + 1 + delim,
                        (
                            (altPos == std::string::npos ? varEnd : altPos)
                          - varBeg - 2*delim
                        )
                    ),
                    false
                );

                std::string altValue;
                if (altPos != std::string::npos)
                {
                    // had ":-" or ":+" alternative value
                    altValue = s.substr
                    (
                        altPos + 2,
                        varEnd - altPos - 2*delim
                    );
                }

                const string varValue = getEnv(varName);
                if (varValue.size())
                {
                    if (altPos != std::string::npos && altType == '+')
                    {
                        // was found, use ":+" alternative
                        s.std::string::replace
                        (
                            varBeg,
                            varEnd - varBeg + 1,
                            altValue
                        );
                        varBeg += altValue.size();
                    }
                    else
                    {
                        // was found, use value
                        s.std::string::replace
                        (
                            varBeg,
                            varEnd - varBeg + 1,
                            varValue
                        );
                        varBeg += varValue.size();
                    }
                }
                else if (altPos != std::string::npos)
                {
                    // use ":-" or ":+" alternative values
                    if (altType == '-')
                    {
                        // was not found, use ":-" alternative
                        s.std::string::replace
                        (
                            varBeg,
                            varEnd - varBeg + 1,
                            altValue
                        );
                        varBeg += altValue.size();
                    }
                    else
                    {
                        // was not found, ":+" alternative implies
                        // substitute with nothing
                        s.std::string::erase(varBeg, varEnd - varBeg + 1);
                    }
                }
                else if (allowEmpty)
                {
                    s.std::string::erase(varBeg, varEnd - varBeg + 1);
                }
                else
                {
                    FatalErrorInFunction
                        << "Unknown variable name '" << varName << "'"
                        << exit(FatalError);
                }
            }
        }
        else
        {
            ++varBeg;
        }
    }

    // Standard handling of "~/", "./" etc.
    standardExpansions(s);
}


bool Foam::stringOps::inplaceReplaceVar(std::string& s, const word& varName)
{
    if (s.empty() || varName.empty())
    {
        return false;
    }

    const string content(getEnv(varName));
    if (content.empty())
    {
        return false;
    }

    const std::string::size_type i = s.find(content);
    if (i == std::string::npos)
    {
        return false;
    }

    s.replace(i, content.size(), string("${" + varName + "}"));
    return true;
}


Foam::string Foam::stringOps::trimLeft(const string& s)
{
    if (!s.empty())
    {
        std::string::size_type beg = 0;
        while (beg < s.size() && std::isspace(s[beg]))
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


void Foam::stringOps::inplaceTrimLeft(std::string& s)
{
    if (!s.empty())
    {
        std::string::size_type beg = 0;
        while (beg < s.size() && std::isspace(s[beg]))
        {
            ++beg;
        }

        if (beg)
        {
            s.erase(0, beg);
        }
    }
}


Foam::string Foam::stringOps::trimRight(const string& s)
{
    if (!s.empty())
    {
        auto n = s.size();
        while (n && std::isspace(s[n-1]))
        {
            --n;
        }

        if (n < s.size())
        {
            return s.substr(0, n);
        }
    }

    return s;
}


void Foam::stringOps::inplaceTrimRight(std::string& s)
{
    if (!s.empty())
    {
        auto n = s.size();
        while (n && std::isspace(s[n-1]))
        {
            --n;
        }

        s.resize(n);
    }
}


Foam::string Foam::stringOps::trim(const string& original)
{
    return trimLeft(trimRight(original));
}


void Foam::stringOps::inplaceTrim(std::string& s)
{
    inplaceTrimRight(s);
    inplaceTrimLeft(s);
}


Foam::string Foam::stringOps::lower(const string& original)
{
    string s(original);
    inplaceLower(s);
    return s;
}


void Foam::stringOps::inplaceLower(std::string& s)
{
    for (auto iter = s.begin(); iter != s.end(); ++iter)
    {
        *iter = static_cast<std::string::value_type>
        (
            std::tolower(static_cast<unsigned char>(*iter))
        );
    }
}


Foam::string Foam::stringOps::upper(const string& original)
{
    string s(original);
    inplaceUpper(s);
    return s;
}


void Foam::stringOps::inplaceUpper(std::string& s)
{
    for (auto iter = s.begin(); iter != s.end(); ++iter)
    {
        *iter = static_cast<std::string::value_type>
        (
            std::toupper(static_cast<unsigned char>(*iter))
        );
    }
}


// ************************************************************************* //
