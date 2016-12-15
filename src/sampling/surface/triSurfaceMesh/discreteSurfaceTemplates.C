/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "discreteSurface.H"
#include "surfFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::discreteSurface::sampleType
(
    const objectRegistry& obr,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;
    typedef IOField<Type> TmpFieldType;

    if (!mesh().foundObject<VolFieldType>(fieldName))
    {
        return false;
    }

    const VolFieldType& fld = mesh().lookupObject<VolFieldType>(fieldName);
    tmp<Field<Type>> tfield = sampleField(fld);

    // The rest could be moved into a separate helper
    if (isA<surfMesh>(obr))
    {
        const surfMesh& surf = dynamicCast<const surfMesh>(obr);

        SurfFieldType* ptr = surf.lookupObjectRefPtr<SurfFieldType>(fieldName);
        if (!ptr)
        {
            // Doesn't exist or the wrong type
            ptr = new SurfFieldType
            (
                IOobject
                (
                    fieldName,
                    surf.time().timeName(),
                    surf,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                surf,
                dimensioned<Type>("0", fld.dimensions(), Zero)
            );
            ptr->writeOpt() = IOobject::NO_WRITE;

            surf.store(ptr);
        }

        ptr->field() = tfield;
    }
    else
    {
        TmpFieldType* ptr = obr.lookupObjectRefPtr<TmpFieldType>(fieldName);
        if (!ptr)
        {
            // Doesn't exist or the wrong type
            ptr = new TmpFieldType
            (
                IOobject
                (
                    fieldName,
                    obr.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                tfield().size()
            );
            ptr->writeOpt() = IOobject::NO_WRITE;

            obr.store(ptr);
        }

        *ptr = tfield;
    }

    return true;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteSurface::sampleType
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (mesh().foundObject<VolFieldType>(fieldName))
    {
        const VolFieldType& fld = mesh().lookupObject<VolFieldType>(fieldName);
        return sampleField(fld);
    }
    else
    {
        return tmp<Field<Type>>();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteSurface::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));
    Field<Type>& values = tvalues.ref();

    if (!onBoundary())
    {
        // Sample cells

        forAll(sampleElements_, triI)
        {
            values[triI] = vField[sampleElements_[triI]];
        }
    }
    else
    {
        // Sample boundary faces

        const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        const label nBnd = mesh().nFaces()-mesh().nInternalFaces();

        // Create flat boundary field

        Field<Type> bVals(nBnd, Zero);

        forAll(vField.boundaryField(), patchi)
        {
            label bFacei = pbm[patchi].start() - mesh().nInternalFaces();

            SubList<Type>
            (
                bVals,
                vField.boundaryField()[patchi].size(),
                bFacei
            ) = vField.boundaryField()[patchi];
        }

        // Sample in flat boundary field

        forAll(sampleElements_, triI)
        {
            label facei = sampleElements_[triI];
            values[triI] = bVals[facei-mesh().nInternalFaces()];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteSurface::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));
    Field<Type>& values = tvalues.ref();

    if (!onBoundary())
    {
        // Sample cells.

        forAll(sampleElements_, pointi)
        {
            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                sampleElements_[pointi]
            );
        }
    }
    else
    {
        // Sample boundary faces.

        forAll(samplePoints_, pointi)
        {
            label facei = sampleElements_[pointi];

            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                mesh().faceOwner()[facei],
                facei
            );
        }
    }

    return tvalues;
}


// ************************************************************************* //
