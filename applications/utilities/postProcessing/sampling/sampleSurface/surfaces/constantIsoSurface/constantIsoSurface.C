/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "constantIsoSurface.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "interpolation.H"
#include "dictionary.H"
#include "meshCutSurface.H"
#include "cellDecompIsoSurfaceCuts.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(constantIsoSurface, 0);

addToRunTimeSelectionTable(surface, constantIsoSurface, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<Type>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    const GeometricField<Type, fvPatchField, volMesh>& vField =
        *cache[fieldName];

    tmp<Field<Type> > tresult(new Field<Type>(faces().size()));
    Field<Type>& result = tresult();

    forAll(result, faceI)
    {
        result[faceI] = vField[cellLabels_[faceI]];
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantIsoSurface::constantIsoSurface
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const word& isoFieldName,
    const scalar isoVal
)
:
    surface(mesh, searchEngine, name),
    isoFieldName_(isoFieldName),
    isoVal_(isoVal),
    points_(0),
    faces_(0),
    cellLabels_(0)
{}


Foam::constantIsoSurface::constantIsoSurface
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    surface(mesh, searchEngine, dict),
    isoFieldName_(dict.lookup("field")),
    isoVal_(readScalar(dict.lookup("value"))),
    points_(0),
    faces_(0),
    cellLabels_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantIsoSurface::~constantIsoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantIsoSurface::correct
(
    const bool meshChanged,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes,
    const fieldsCache<scalar>& scalarCache
)
{
    if (!scalarCache.found(isoFieldName_))
    {
        FatalErrorIn
        (
            "constantIsoSurface::correct(const bool meshChanged,"
            "const volPointInterpolation&, const fieldsCache<scalar>&"
            ", const fieldsCache<vector>&, const fieldsCache<tensor>&)"
        )   << "Field " << isoFieldName_ << " not loaded." << endl
            << "It has to be one of the sampled fields"
            << exit(FatalError);
    }
    const volScalarField& vField = *scalarCache[isoFieldName_];

    const pointScalarField& pField =
        scalarCache.pointField(isoFieldName_, pInterp);

    // Create smooth volField.
    volScalarField smoothVolField(vField);

    const labelListList& cellPoints = vField.mesh().cellPoints();

    forAll(cellPoints, cellI)
    {
        const labelList& cPoints = cellPoints[cellI];

        scalar sum = 0;

        forAll(cPoints, i)
        {
            sum += pField[cPoints[i]];
        }
        smoothVolField[cellI] = sum / cPoints.size();
    }
    

    cellDecompIsoSurfaceCuts isoSurfaceCuts 
    (
        smoothVolField,
        pField,
        isoVal_,
        -0.1
    );

    meshCutSurface isoSurf(isoSurfaceCuts);

    points_ = isoSurf.points();

    // Convert triangles into faces
    faces_.setSize(isoSurf.size());
    cellLabels_.setSize(isoSurf.size());

    forAll(isoSurf, triI)
    {
        face& f = faces_[triI];
        const labelledTri& t = isoSurf[triI];

        f.setSize(3);
        f[0] = t[0];
        f[1] = t[1];
        f[2] = t[2];

        cellLabels_[triI] = t.region();
    }

    Pout<< "Created " << name() << " :"
        << "  isoValue:" << isoVal_
        << "  field:" << isoFieldName_
        << "  faces:" << faces_.size()
        << "  points:" << points_.size() << endl;
}


Foam::tmp<Foam::scalarField> Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<scalar>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<scalar>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::vectorField> Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<vector>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<vector>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::sphericalTensorField> Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<sphericalTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<sphericalTensor>
    (
        fieldName,
        cache,
        pInterp,
        interpolationSchemes
    );
}


Foam::tmp<Foam::symmTensorField> Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<symmTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<symmTensor>
    (
        fieldName,
        cache, 
        pInterp,
        interpolationSchemes
    );
}


Foam::tmp<Foam::tensorField> Foam::constantIsoSurface::interpolate
(
    const word& fieldName,
    const fieldsCache<tensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<tensor>(fieldName, cache, pInterp, interpolationSchemes);
}


// ************************************************************************* //
