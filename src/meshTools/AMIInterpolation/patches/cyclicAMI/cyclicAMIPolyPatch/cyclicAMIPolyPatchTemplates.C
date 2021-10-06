/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolateUntransformed
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    if (owner())
    {
        return AMI().interpolateToSource(fld, defaultValues);
    }
    else
    {
        return neighbPatch().AMI().interpolateToTarget(fld, defaultValues);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    if (pTraits<Type>::rank == 0)
    {
        return interpolateUntransformed(fld, defaultValues);
    }
    else
    {
        autoPtr<coordSystem::cylindrical> cs(cylindricalCS());
        if (!cs.valid())
        {
            return interpolateUntransformed(fld, defaultValues);
        }
        else
        {
            const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::interpolate :"
                    << " patch:" << this->name()
                    << " size:" << this->size()
                    << " nbrPatch:" << nbrPp.name()
                    << " size:" << nbrPp.size()
                    << endl;
            }

            if (fld.size() != nbrPp.size())
            {
                FatalErrorInFunction
                    << "Patch:" << this->name()
                    << " size:" << this->size()
                    << " neighbour patch:" << nbrPp.name()
                    << " size:" << nbrPp.size()
                    << " fld size:" << fld.size()
                    << exit(FatalError);
            }


            auto tlocalFld(tmp<Field<Type>>::New(fld.size()));
            Field<Type>& localFld = tlocalFld.ref();

            // Transform to cylindrical coords
            {
                tmp<tensorField> nbrT(cs().R(nbrPp.faceCentres()));
                localFld = Foam::invTransform(nbrT(), fld);
            }

            if (debug&2)
            {
                const vectorField::subField nbrFc(nbrPp.faceCentres());

                Pout<< "On patch:" << this->name()
                    << " size:" << this->size()
                    << " fc:" << gAverage(this->faceCentres())
                    << " getting remote data from:" << nbrPp.name()
                    << " size:" << nbrPp.size()
                    << " fc:" << gAverage(nbrFc)
                    << endl;

                forAll(fld, i)
                {
                    Pout<< "At:" << nbrFc[i] << nl
                        << "    cart:" << fld[i] << nl
                        << "    cyli:" << localFld[i] << nl
                        << endl;
                }
            }


            const tmp<tensorField> T(cs().R(this->faceCentres()));

            List<Type> localDeflt(defaultValues.size());
            if (defaultValues.size() == size())
            {
                // Transform default values into cylindrical coords (using
                // *this faceCentres)
                // We get in UList (why? Copied from cyclicAMI). Convert to
                // Field so we can use transformField routines.
                const SubField<Type> defaultSubFld(defaultValues);
                const Field<Type>& defaultFld(defaultSubFld);
                localDeflt = Foam::invTransform(T(), defaultFld);
            }

            // Do the actual interpolation and interpolate back to cartesian
            // coords
            return Foam::transform
            (
                T,
                interpolateUntransformed(localFld, localDeflt)
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::cyclicAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    //- Commented out for now since called with non-primitives (e.g. wallPoint
    //  from FaceCellWave) - these are missing the pTraits<Type>::rank and
    //  Foam::transform
    /*
    autoPtr<coordSystem::cylindrical> cs(cylindricalCS());

    if (cs.valid() && pTraits<Type>::rank > 0)
    {
        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

        tmp<tensorField> nbrT(cs().R(nbrPp.faceCentres()));

        result = Foam::invTransform(nbrT, result);
        List<Type> localDeflt(defaultValues.size());
        if (defaultValues.size() == nbrT().size())
        {
            // We get in UList (why? Copied from cyclicAMI). Convert to
            // Field so we can use transformField routines.
            const SubField<Type> defaultSubFld(defaultValues);
            const Field<Type>& defaultFld(defaultSubFld);
            localDeflt = Foam::invTransform(nbrT, defaultFld);
        }

        // Do actual AMI interpolation
        if (owner())
        {
            AMI().interpolateToSource
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }
        else
        {
            neighbPatch().AMI().interpolateToTarget
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }

        // Transform back. Result is now at *this
        const vectorField::subField fc(this->faceCentres());
        result = Foam::transform(cs().R(fc), result);
    }
    else
    */
    {
        if (owner())
        {
            AMI().interpolateToSource
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
        else
        {
            neighbPatch().AMI().interpolateToTarget
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
    }
}


// ************************************************************************* //
