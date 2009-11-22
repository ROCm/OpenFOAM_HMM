/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Distribution.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::Distribution()
:
    List< Map<label> >(pTraits<Type>::nComponents),
    binWidth_(pTraits<Type>::one)
{}


template<class Type>
Foam::Distribution<Type>::Distribution(const Type& binWidth)
:
    List< Map<label> >(pTraits<Type>::nComponents),
    binWidth_(binWidth)
{}


// template<class Type>
// Foam::Distribution<Type>::Distribution
// (
//     const cmptType& binWidth
// )
// :
//     List< Map<label> >(pTraits<Type>::nComponents),
//     binWidth_(binWidth*pTraits<Type>::one)
// {}


template<class Type>
Foam::Distribution<Type>::Distribution(const Distribution<Type>& d)
:
    List< Map<label> >(static_cast< const List< Map<label> >& >(d)),
    binWidth_(d.binWidth())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::~Distribution()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Distribution<Type>::add(const Type& valueToAdd)
{
    for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<label>& cmptDistribution = (*this)[cmpt];

        Map<label>::iterator iter(cmptDistribution.begin());

        label n =
            label(component(valueToAdd, cmpt)/component(binWidth_, cmpt))
          - label(neg(component(valueToAdd, cmpt)/component(binWidth_, cmpt)));

        iter = cmptDistribution.find(n);

        if (iter == cmptDistribution.end())
        {
            cmptDistribution.insert(n,1);
        }
        else
        {
            cmptDistribution[n]++;
        }

        if (cmptDistribution[n] < 0)
        {
            FatalErrorIn("Distribution::add(const scalar valueToAdd)")
                << "Accumulated Distribution value has become negative: "
                << "bin = " << (0.5 + scalar(n))*component(binWidth_, cmpt)
                << ", value = " << cmptDistribution[n]
                << ". This is most likely to be because too many samples "
                << "have been added to a bin and the label has 'rolled round'"
                << abort(FatalError);
        }
    }
}


template<class Type>
void Foam::Distribution<Type>::insertMissingKeys()
{
    for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<label>& cmptDistribution = (*this)[cmpt];

        Map<label>::iterator iter(cmptDistribution.begin());

        List<label> keys = cmptDistribution.sortedToc();

        if (keys.size())
        {
            for (label k = keys[0]; k < keys[keys.size()-1]; k++)
            {
                iter = cmptDistribution.find(k);

                if (iter == cmptDistribution.end())
                {
                    cmptDistribution.insert(k,0);
                }
            }
        }
    }
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::raw()
{
    List< List < Pair<scalar> > > rawDistributions(pTraits<Type>::nComponents);

    insertMissingKeys();

    for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<label>& cmptDistribution = (*this)[cmpt];

        List<label> keys = cmptDistribution.sortedToc();

        List<Pair<scalar> >& rawDist = rawDistributions[cmpt];

        rawDist.setSize(keys.size());

        forAll(keys,k)
        {
            label key = keys[k];

            rawDist[k].first() = (0.5 + scalar(key))*component(binWidth_, cmpt);

            rawDist[k].second() = scalar(cmptDistribution[key]);
        }
    }

    return rawDistributions;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Distribution<Type>::operator=
(
    const Distribution<Type>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::Distribution<Type>::operator="
            "(const Foam::Distribution<Type>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    List< Map<label> >::operator=(rhs);

    binWidth_ = rhs.binWidth();
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Distribution<Type>& d
)
{
    os  << d.binWidth_
        << static_cast<const List< Map<label> >& >(d);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const Distribution&)"
    );

    return os;
}


// ************************************************************************* //
