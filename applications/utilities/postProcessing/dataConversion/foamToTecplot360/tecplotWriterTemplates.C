/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "tecplotWriter.H"
#include "TECIO.h"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::tecplotWriter::writeField(const Field<Type>& fld) const
{
    typedef typename pTraits<Type>::cmptType cmptType;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        tmp<Field<cmptType>> tcmptFld = fld.component(cmpt);
        const Field<cmptType>& cmptFld = tcmptFld();

        const int32_t size = int32_t(cmptFld.size());

        // Convert to float
        Field<float> floats(size);
        forAll(cmptFld, i)
        {
            floats[i] = float(cmptFld[i]);
        }

        //Pout<< "Writing component:" << cmpt << " of size:" << size
        //    << " floats." << endl;

        if
        (
            tecdat142
            (
                &size,
                floats.cdata(),
                &tecConst_False     //< VIsDouble (0: single, 1: double)
            )
        )
        {
            FatalErrorInFunction
                << "Error in tecdat142."
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::tecplotWriter::writeField(const tmp<Field<Type>>& tfld) const
{
    writeField(tfld());
}


template<class GeoField>
void Foam::tecplotWriter::writeFields
(
    const PtrList<const GeoField>& flds
) const
{
    forAll(flds, i)
    {
        writeField(flds[i]);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::tecplotWriter::getPatchField
(
    const bool nearCellValue,
    const GeometricField<Type, fvPatchField, volMesh>& vfld,
    const label patchi
)
{
    if (nearCellValue)
    {
        return vfld.boundaryField()[patchi].patchInternalField();
    }
    else
    {
        return vfld.boundaryField()[patchi];
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::tecplotWriter::getFaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld,
    const labelUList& faceLabels
)
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    tmp<Field<Type>> tfld(new Field<Type>(faceLabels.size()));
    Field<Type>& fld = tfld.ref();

    forAll(faceLabels, i)
    {
        const label facei  = faceLabels[i];
        const label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            fld[i] = sfld[facei];
        }
        else
        {
            label localFacei = facei - patches[patchi].start();
            fld[i] = sfld.boundaryField()[patchi][localFacei];
        }
    }

    return tfld;
}


template<class GeoField>
Foam::wordList
Foam::tecplotWriter::getNames
(
    const PtrList<const GeoField>& flds
)
{
    wordList names(flds.size());
    forAll(flds, i)
    {
        names[i] = flds[i].name();
    }
    return names;
}


template<class Type>
void Foam::tecplotWriter::getTecplotNames
(
    const wordList& names,
    const int32_t loc,
    string& varNames,
    DynamicList<int32_t>& varLocation
)
{
    const direction nCmpts = pTraits<Type>::nComponents;

    forAll(names, i)
    {
        if (!varNames.empty())
        {
            varNames += " ";
        }

        if (nCmpts == 1)
        {
            varNames += names[i];
            varLocation.append(loc);
        }
        else
        {
            for (direction cmpt = 0; cmpt < nCmpts; ++cmpt)
            {
                varNames +=
                (
                    (cmpt ? " " : string::null)
                  + names[i] + "_" + pTraits<Type>::componentNames[cmpt]
                );
                varLocation.append(loc);
            }
        }
    }
}


template<class GeoField>
void Foam::tecplotWriter::getTecplotNames
(
    const PtrList<GeoField>& flds,
    const int32_t loc,
    string& varNames,
    DynamicList<int32_t>& varLocation
)
{
    getTecplotNames<typename GeoField::value_type>
    (
        getNames(flds),
        loc,
        varNames,
        varLocation
    );
}


// ************************************************************************* //
