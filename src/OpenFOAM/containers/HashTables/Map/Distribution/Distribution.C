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
#include "OFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::Distribution()
:
    List< Map<scalar> >(pTraits<Type>::nComponents),
    binWidth_(pTraits<Type>::one)
{}


template<class Type>
Foam::Distribution<Type>::Distribution(const Type& binWidth)
:
    List< Map<scalar> >(pTraits<Type>::nComponents),
    binWidth_(binWidth)
{}


template<class Type>
Foam::Distribution<Type>::Distribution(const Distribution<Type>& d)
:
    List< Map<scalar> >(static_cast< const List< Map<scalar> >& >(d)),
    binWidth_(d.binWidth())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Distribution<Type>::~Distribution()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Distribution<Type>::totalWeight(direction cmpt) const
{
    const Map<scalar>& cmptDistribution = (*this)[cmpt];

    scalar sumOfWeights = 0.0;

    forAllConstIter(Map<scalar>, cmptDistribution, iter)
    {
        sumOfWeights += iter();
    }

    return sumOfWeights;
}


template<class Type>
inline Type Foam::Distribution<Type>::mean() const
{
    Type meanValue(pTraits<Type>::zero);

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const Map<scalar>& cmptDistribution = (*this)[cmpt];

        scalar totalCmptWeight = totalWeight(cmpt);

        List<label> keys = cmptDistribution.sortedToc();

        forAll(keys,k)
        {
            label key = keys[k];

            setComponent(meanValue, cmpt) +=
                (0.5 + scalar(key))
               *component(binWidth_, cmpt)
               *cmptDistribution[key]
               /totalCmptWeight;
        }
    }

    return meanValue;
}


template<class Type>
inline Type Foam::Distribution<Type>::median()
{
    Type medianValue(pTraits<Type>::zero);

    List< List < Pair<scalar> > > normDistribution = normalised();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        List< Pair<scalar> >& normDist = normDistribution[cmpt];

        if (normDist.size())
        {
            if (normDist.size() == 1)
            {
                setComponent(medianValue, cmpt) = normDist[0].first();
            }
            else if
            (
                normDist.size() > 1
             && normDist[0].second()*component(binWidth_, cmpt) > 0.5
            )
            {
                scalar xk = normDist[1].first();

                scalar xkm1 = normDist[0].first();

                scalar Sk =
                    (normDist[0].second() + normDist[1].second())
                   *component(binWidth_, cmpt);

                scalar Skm1 = normDist[0].second()*component(binWidth_, cmpt);

                setComponent(medianValue, cmpt) =
                    (0.5 - Skm1)*(xk - xkm1)/(Sk - Skm1) + xkm1;
            }
            else
            {
                label lastNonZeroIndex = 0;

                scalar cumulative = 0.0;

                forAll(normDist,nD)
                {
                    if
                    (
                        cumulative
                      + (normDist[nD].second()*component(binWidth_, cmpt))
                      > 0.5
                    )
                    {
                        scalar xk = normDist[nD].first();

                        scalar xkm1 = normDist[lastNonZeroIndex].first();

                        scalar Sk =
                            cumulative
                          + (normDist[nD].second()*component(binWidth_, cmpt));

                        scalar Skm1 = cumulative;

                        setComponent(medianValue, cmpt) =
                            (0.5 - Skm1)*(xk - xkm1)/(Sk - Skm1) + xkm1;

                        break;
                    }
                    else if (mag(normDist[nD].second()) > VSMALL)
                    {
                        cumulative +=
                            normDist[nD].second()*component(binWidth_, cmpt);

                        lastNonZeroIndex = nD;
                    }
                }
            }
        }

    }

    return medianValue;
}


template<class Type>
void Foam::Distribution<Type>::add
(
    const Type& valueToAdd,
    const Type& weight
)
{
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<scalar>& cmptDistribution = (*this)[cmpt];

        Map<scalar>::iterator iter(cmptDistribution.begin());

        label n =
            label(component(valueToAdd, cmpt)/component(binWidth_, cmpt))
          - label(neg(component(valueToAdd, cmpt)/component(binWidth_, cmpt)));

        iter = cmptDistribution.find(n);

        if (iter == cmptDistribution.end())
        {
            cmptDistribution.insert(n, component(weight, cmpt));
        }
        else
        {
            cmptDistribution[n] += component(weight, cmpt);
        }
    }
}


template<class Type>
void Foam::Distribution<Type>::insertMissingKeys()
{
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<scalar>& cmptDistribution = (*this)[cmpt];

        Map<scalar>::iterator iter(cmptDistribution.begin());

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
Distribution<Type>::normalised()
{
    List< List < Pair<scalar> > > normDistribution(pTraits<Type>::nComponents);

    insertMissingKeys();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<scalar>& cmptDistribution = (*this)[cmpt];

        scalar totalCmptWeight = totalWeight(cmpt);

        List<label> keys = cmptDistribution.sortedToc();

        List< Pair<scalar> >& normDist = normDistribution[cmpt];

        normDist.setSize(keys.size());

        forAll(keys,k)
        {
            label key = keys[k];

            normDist[k].first() =
                (0.5 + scalar(key))*component(binWidth_, cmpt);

            normDist[k].second() =
                cmptDistribution[key]
               /totalCmptWeight
               /component(binWidth_, cmpt);
        }
    }

    return normDistribution;
}


template<class Type>
Foam::List< Foam::List< Foam::Pair<Foam::scalar> > >Foam::
Distribution<Type>::raw()
{
    List< List < Pair<scalar> > > rawDistribution(pTraits<Type>::nComponents);

    insertMissingKeys();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        Map<scalar>& cmptDistribution = (*this)[cmpt];

        List<label> keys = cmptDistribution.sortedToc();

        List< Pair<scalar> >& rawDist = rawDistribution[cmpt];

        rawDist.setSize(keys.size());

        forAll(keys,k)
        {
            label key = keys[k];

            rawDist[k].first() = (0.5 + scalar(key))*component(binWidth_, cmpt);

            rawDist[k].second() = cmptDistribution[key];
        }
    }

    return rawDistribution;
}


template<class Type>
void Foam::Distribution<Type>::write
(
    const fileName& filePrefix,
    const List< List< Pair<scalar> > >& pairs
) const
{
    if (pairs.size() != pTraits<Type>::nComponents)
    {
        FatalErrorIn
        (
            "Distribution::write"
            "("
                "const fileName& filePrefix,"
                "const List< List< Pair<scalar> > >& pairs"
            ")"
        )
            << "List of pairs (" << pairs.size()
            << ") is not the same size as the number of components ("
            << pTraits<Type>::nComponents << ")." << nl
            << abort(FatalError);
    }

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const List< Pair<scalar> >& cmptPairs = pairs[cmpt];

        OFstream os(filePrefix + '_' + pTraits<Type>::componentNames[cmpt]);

        forAll(cmptPairs, i)
        {
            os  << cmptPairs[i].first() << ' ' << cmptPairs[i].second() << nl;
        }
    }
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

    List< Map<scalar> >::operator=(rhs);

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
        << static_cast<const List< Map<scalar> >& >(d);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, "
        "const Distribution&)"
    );

    return os;
}


// ************************************************************************* //
