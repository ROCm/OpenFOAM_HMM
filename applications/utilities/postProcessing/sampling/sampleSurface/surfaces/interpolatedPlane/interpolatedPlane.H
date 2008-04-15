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

Class
    Foam::interpolatedPlane

Description

SourceFiles
    interpolatedPlane.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedPlane_H
#define interpolatedPlane_H

#include "surface.H"
#include "constantPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class meshSearch;
class dictionary;
class plane;
class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class interpolatedPlane Declaration
\*---------------------------------------------------------------------------*/

class interpolatedPlane
:
    public constantPlane
{

    // Private Member Functions

        //- Interpolate field onto vertices
        template<class Type>
        tmp<Field<Type> > interpolate
        (
            const word& fieldName,
            const fieldsCache<Type>& cache,
            const volPointInterpolation& pInterp,
            const dictionary& interpolationSchemes
        ) const;


public:

    //- Runtime type information
    TypeName("interpolatedPlane");

    // Constructors

        //- Construct from components
        interpolatedPlane
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& name,
            const plane& planeDesc,
            const bool triangulate
        );

        //- Construct from dictionary
        interpolatedPlane
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict          
        );


    // Destructor

        virtual ~interpolatedPlane();


    // Member Functions

        //- correct surface for new mesh done by constantPlane::correct()

        //- interpolate field to surface
        virtual tmp<scalarField> interpolate
        (
            const word&,
            const fieldsCache<scalar>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<vectorField> interpolate
        (
            const word&,
            const fieldsCache<vector>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const word&,
            const fieldsCache<sphericalTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<symmTensorField> interpolate
        (
            const word&,
            const fieldsCache<symmTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<tensorField> interpolate
        (
            const word&,
            const fieldsCache<tensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
