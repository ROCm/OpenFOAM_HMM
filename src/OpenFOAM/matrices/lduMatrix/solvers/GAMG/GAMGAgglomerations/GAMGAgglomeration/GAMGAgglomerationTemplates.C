/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (restrictDistributeMapPtr_.valid())
    {
        const mapDistribute& map = restrictDistributeMapPtr_();
        Field<Type> combinedFf(map.constructSize());
        forAll(ff, i)
        {
            combinedFf[i] = ff[i];
        }
        map.distribute(combinedFf);

        if (combinedFf.size() != fineToCoarse.size())
        {
            FatalErrorIn
            (
                "void GAMGAgglomeration::restrictField"
                "(Field<Type>& cf, const Field<Type>& ff, "
                "const label fineLevelIndex) const"
            )   << "field does not correspond to level " << fineLevelIndex
                << " sizes: field = " << combinedFf.size()
                << " level = " << fineToCoarse.size()
                << abort(FatalError);
        }

        cf = pTraits<Type>::zero;

        forAll(combinedFf, i)
        {
            cf[fineToCoarse[i]] += combinedFf[i];
        }
    }
    else
    {
        if (ff.size() != fineToCoarse.size())
        {
            FatalErrorIn
            (
                "void GAMGAgglomeration::restrictField"
                "(Field<Type>& cf, const Field<Type>& ff, "
                "const label fineLevelIndex) const"
            )   << "field does not correspond to level " << fineLevelIndex
                << " sizes: field = " << ff.size()
                << " level = " << fineToCoarse.size()
                << abort(FatalError);
        }

        cf = pTraits<Type>::zero;

        forAll(ff, i)
        {
            cf[fineToCoarse[i]] += ff[i];
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (faceRestrictDistributeMapPtr_.valid())
    {
        const mapDistribute& map = faceRestrictDistributeMapPtr_();
        Field<Type> combinedFf(map.constructSize());
        forAll(ff, i)
        {
            combinedFf[i] = ff[i];
        }
        map.distribute(combinedFf);

        cf = pTraits<Type>::zero;

        forAll(fineToCoarse, ffacei)
        {
            label cFace = fineToCoarse[ffacei];

            if (cFace >= 0)
            {
                cf[cFace] += combinedFf[ffacei];
            }
        }
    }
    else
    {
        cf = pTraits<Type>::zero;

        forAll(fineToCoarse, ffacei)
        {
            label cFace = fineToCoarse[ffacei];

            if (cFace >= 0)
            {
                cf[cFace] += ff[ffacei];
            }
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    const Field<Type>& cf,
    const label coarseLevelIndex
) const
{
    const labelList& fineToCoarse = restrictAddressing_[coarseLevelIndex];

    if (restrictDistributeMapPtr_.valid())
    {
        const mapDistribute& map = restrictDistributeMapPtr_();
        Field<Type> combinedCf(cf);
        map.reverseDistribute(fineToCoarse.size(), combinedCf);

        forAll(fineToCoarse, i)
        {
            ff[i] = combinedCf[fineToCoarse[i]];
        }
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
}


// ************************************************************************* //
