/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    const word& fieldName,
    const word& sampleScheme
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;
    typedef IOField<Type> TmpFieldType;

    const auto* volFldPtr = mesh().findObject<VolFieldType>(fieldName);

    if (!volFldPtr)
    {
        return false;
    }

    const auto& volFld = *volFldPtr;
    auto samplerPtr = interpolation<Type>::New(sampleScheme, *volFldPtr);

    tmp<Field<Type>> tfield = sampleOnFaces(*samplerPtr);

    // The rest could be moved into a separate helper
    if (isA<surfMesh>(obr))
    {
        const surfMesh& surf = dynamicCast<const surfMesh>(obr);

        SurfFieldType* ptr = surf.getObjectPtr<SurfFieldType>(fieldName);
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
                dimensioned<Type>(volFld.dimensions(), Zero)
            );
            ptr->writeOpt() = IOobject::NO_WRITE;

            surf.store(ptr);
        }

        ptr->field() = tfield;
    }
    else
    {
        TmpFieldType* ptr = obr.getObjectPtr<TmpFieldType>(fieldName);
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
    const word& fieldName,
    const word& sampleScheme
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const auto* volFldPtr = mesh().findObject<VolFieldType>(fieldName);

    if (volFldPtr)
    {
        auto samplerPtr = interpolation<Type>::New(sampleScheme, *volFldPtr);
        return sampleOnFaces(*samplerPtr);
    }

    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteSurface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const labelList& elements = sampleElements_;

    const auto& vField = sampler.psi();
    const label len = elements.size();

    // One value per face
    auto tvalues = tmp<Field<Type>>::New(len);
    auto& values = tvalues.ref();

    if (!onBoundary())
    {
        // Sample cells

        const faceList& fcs = this->surfFaces();
        const pointField& pts = points();

        for (label i=0; i < len; ++i)
        {
            const label celli = elements[i];
            const point pt = fcs[i].centre(pts);

            values[i] = sampler.interpolate(pt, celli);
        }
    }
    else
    {
        // Sample boundary faces

        const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        const label nBnd = mesh().nBoundaryFaces();

        // Create flat boundary field

        Field<Type> bVals(nBnd, Zero);

        const auto& bField = vField.boundaryField();

        forAll(bField, patchi)
        {
            const label bFacei = pbm[patchi].start() - mesh().nInternalFaces();

            SubList<Type>
            (
                bVals,
                bField[patchi].size(),
                bFacei
            ) = bField[patchi];
        }

        // Sample in flat boundary field

        for (label i=0; i < len; ++i)
        {
            const label bFacei = (elements[i] - mesh().nInternalFaces());
            values[i] = bVals[bFacei];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::discreteSurface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    auto tvalues = tmp<Field<Type>>::New(sampleElements_.size());
    auto& values = tvalues.ref();

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
            const label facei = sampleElements_[pointi];

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
