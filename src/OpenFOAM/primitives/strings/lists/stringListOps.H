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

InNamspace
    Foam

Description
    Operations on lists of strings.

SourceFiles
    stringListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef stringListOps_H
#define stringListOps_H

#include "labelList.H"
#include "stringList.H"
#include "regExp.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Single-string match for one of the regular expressions
    //  \deprecated use stringOps::match instead (deprecated NOV-2017)
    inline bool findStrings(const wordRes& matcher, const std::string& text)
    {
        return matcher(text);
    }


    //- Extract list indices
    //  The unary match predicate has the following signature:
    //  \code
    //  bool operator()(const std::string& text);
    //  \endcode
    //
    //  \return List indices for matching strings
    template<class UnaryMatchPredicate, class StringType>
    labelList findMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        const UList<StringType>& lst,
        const bool invert=false
    );


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const regExp& matcher,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, lst, invert);
    }


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const char* re,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        return findMatchingStrings(matcher, lst, invert);
    }


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const std::string& re,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        return findMatchingStrings(matcher, lst, invert);
    }


    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordRe& matcher,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, lst, invert);
    }


    //- Return list indices for strings matching one of the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordRes& matcher,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, lst, invert);
    }


  // Subsetting multi-string matches (similar to ListOp):

    //- Extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //  \code
    //    subsetMatchingStrings<regExp, stringList>(myRegExp, lst);
    //  \endcode
    template<class UnaryMatchPredicate, class StringListType>
    StringListType subsetMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        const StringListType& lst,
        const bool invert=false
    );


    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const regExp& matcher,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, lst, invert);
    }


    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const char* re,
        const StringListType& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        return subsetMatchingStrings(matcher, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const std::string& re,
        const StringListType& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        return subsetMatchingStrings(matcher, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordRe& matcher,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordRes& matcher,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, lst, invert);
    }


    //- Inplace extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //    inplaceSubsetMatchingStrings<regExp, stringList>(myRegExp, lst);
    template<class UnaryMatchPredicate, class StringListType>
    void inplaceSubsetMatchingStrings
    (
        const UnaryMatchPredicate& matcher,
        StringListType& lst,
        const bool invert=false
    );

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const regExp& matcher,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const char* re,
        StringListType& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const std::string& re,
        StringListType& lst,
        const bool invert=false
    )
    {
        const regExp matcher(re);
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordRe& matcher,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordRes& matcher,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stringListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
