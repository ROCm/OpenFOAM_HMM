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

#include "NamedEnum.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Enum, int nEnum>
template<class StringType>
Foam::List<StringType> Foam::NamedEnum<Enum, nEnum>::getNamesList()
{
    List<StringType> lst(nEnum);

    label count = 0;
    for (int enumi=0; enumi < nEnum; ++enumi)
    {
        if (names[enumi] && names[enumi][0])
        {
            lst[count++] = names[enumi];
        }
    }

    lst.setSize(count);
    return lst;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Enum, int nEnum>
Foam::NamedEnum<Enum, nEnum>::NamedEnum()
:
    table_type(2*nEnum)
{
    for (int enumi=0; enumi < nEnum; ++enumi)
    {
        if (names[enumi] && names[enumi][0])
        {
            insert(names[enumi], enumi);
        }
        else
        {
            // Bad name - generate error message
            stringList goodNames(enumi);

            for (int i = 0; i < enumi; ++i)
            {
                goodNames[i] = names[i];
            }

            FatalErrorInFunction
                << "Illegal enumeration name at position " << enumi << nl
                << "after entries " << goodNames << nl
                << "Possibly your NamedEnum<Enum, nEnum>::names array"
                << " is not of size " << nEnum << endl
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Enum, int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::read(Istream& is) const
{
    const word enumName(is);
    table_type::const_iterator iter = find(enumName);

    if (!iter.found())
    {
        FatalIOErrorInFunction(is)
            << enumName << " is not in enumeration: "
            << sortedToc() << exit(FatalIOError);
    }

    return Enum(iter.object());
}


template<class Enum, int nEnum>
void Foam::NamedEnum<Enum, nEnum>::write(const Enum e, Ostream& os) const
{
    os  << names[int(e)];
}


template<class Enum, int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::lookup
(
    const word& key,
    const dictionary& dict
) const
{
    const word enumName(dict.lookup(key));
    table_type::const_iterator iter = find(enumName);

    if (!iter.found())
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << sortedToc() << exit(FatalIOError);
    }

    return Enum(iter.object());
}


template<class Enum, int nEnum>
Enum Foam::NamedEnum<Enum, nEnum>::lookupOrDefault
(
    const word& key,
    const dictionary& dict,
    const enum_type deflt
) const
{
    if (dict.found(key))
    {
        return lookup(key, dict);
    }
    else
    {
        return deflt;
    }
}


template<class Enum, int nEnum>
Foam::List<Enum> Foam::NamedEnum<Enum, nEnum>::enums()
{
    List<Enum> lst(nEnum);

    label count = 0;
    for (int enumi = 0; enumi < nEnum; ++enumi)
    {
        if (names[enumi] && names[enumi][0])
        {
            lst[count++] = Enum(enumi);
        }
    }

    lst.setSize(count);
    return lst;
}


template<class Enum, int nEnum>
Foam::stringList Foam::NamedEnum<Enum, nEnum>::strings()
{
    return getNamesList<string>();
}


template<class Enum, int nEnum>
Foam::wordList Foam::NamedEnum<Enum, nEnum>::words()
{
    return getNamesList<word>();
}


// ************************************************************************* //
