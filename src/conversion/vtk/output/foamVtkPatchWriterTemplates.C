/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkPatchWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::foamVtkOutput::patchWriter::write
(
    const UPtrList<const GeometricField<Type, PatchField, volMesh>>& flds
)
{
    forAll(flds, fieldi)
    {
        const auto& fld = flds[fieldi];

        os_ << fld.name() << ' '
            << int(pTraits<Type>::nComponents) << ' '
            << nFaces_ << " float" << nl;

        forAll(patchIDs_, i)
        {
            const auto& pfld = fld.boundaryField()[patchIDs_[i]];

            if (nearCellValue_)
            {
                foamVtkOutput::writeList(format(), pfld.patchInternalField()());
            }
            else
            {
                foamVtkOutput::writeList(format(), pfld);
            }
        }

        format().flush();
    }
}


template<class Type, template<class> class PatchField>
void Foam::foamVtkOutput::patchWriter::write
(
    const UPtrList<const GeometricField<Type, PatchField, pointMesh>>& flds
)
{
    forAll(flds, fieldi)
    {
        const auto& fld = flds[fieldi];

        os_ << fld.name() << ' '
            << int(pTraits<Type>::nComponents) << ' '
            << nPoints_ << " float" << nl;

        forAll(patchIDs_, i)
        {
            const auto& pfld = fld.boundaryField()[patchIDs_[i]];

            foamVtkOutput::writeList(format(), pfld.patchInternalField()());
        }
        format().flush();
    }
}


template<class Type>
void Foam::foamVtkOutput::patchWriter::write
(
    const PrimitivePatchInterpolation<primitivePatch>& pInter,
    const UPtrList<const GeometricField<Type, fvPatchField, volMesh>>& flds
)
{
    forAll(flds, fieldi)
    {
        const auto& fld = flds[fieldi];

        os_ << fld.name() << ' '
            << int(pTraits<Type>::nComponents) << ' '
            << nPoints_ << " float" << nl;

        DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nPoints_);

        forAll(patchIDs_, i)
        {
            const auto& pfld = fld.boundaryField()[patchIDs_[i]];

            if (nearCellValue_)
            {
                auto tfield =
                    pInter.faceToPointInterpolate(pfld.patchInternalField()());

                foamVtkOutput::writeList(format(), tfield());
            }
            else
            {
                auto tfield = pInter.faceToPointInterpolate(pfld);

                foamVtkOutput::writeList(format(), tfield());
            }
        }
        format().flush();
    }
}


// ************************************************************************* //
