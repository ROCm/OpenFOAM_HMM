/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::sampledTriSurfaceMeshNormal

Description
    Variant of sampledTriSurfaceMesh that samples the surface-normal component
    of a vector field.

    Returns a vector field with the value in the first component and sets
    the other two to zero.

SourceFiles
    sampledTriSurfaceMeshNormal.C
    sampledTriSurfaceMeshNormalTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledTriSurfaceMeshNormal_H
#define sampledTriSurfaceMeshNormal_H

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sampledTriSurfaceMeshNormal Declaration
\*---------------------------------------------------------------------------*/

class sampledTriSurfaceMeshNormal
:
    public sampledTriSurfaceMesh
{
public:

    //- Runtime type information
    TypeName("sampledTriSurfaceMeshNormal");


    // Constructors

        //- Construct from components
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName,
            const samplingSource sampleSource
        );

        //- Construct from dictionary
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        sampledTriSurfaceMeshNormal
        (
            const word& name,
            const polyMesh& mesh,
            const triSurface& surface,
            const word& sampleSourceName
        );


    //- Destructor
    virtual ~sampledTriSurfaceMeshNormal();


    // Member Functions

        //- Sample field on surface
        virtual tmp<scalarField> sample
        (
            const volScalarField&
        ) const
        {
            NotImplemented;
            return tmp<scalarField>();
        }

        //- Sample field on surface
        virtual tmp<vectorField> sample
        (
            const volVectorField&
        ) const;

        //- Sample field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const volSphericalTensorField&
        ) const
        {
            NotImplemented;
            return tmp<sphericalTensorField>();
        }

        //- Sample field on surface
        virtual tmp<symmTensorField> sample
        (
            const volSymmTensorField&
        ) const
        {
            NotImplemented;
            return tmp<symmTensorField>();
        }

        //- Sample field on surface
        virtual tmp<tensorField> sample
        (
            const volTensorField&
        ) const
        {
            NotImplemented;
            return tmp<tensorField>();
        }

        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const
        {
            NotImplemented;
            return tmp<scalarField>();
        }

        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const
        {
            NotImplemented;
            return tmp<sphericalTensorField>();
        }

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const
        {
            NotImplemented;
            return tmp<symmTensorField>();
        }

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const
        {
            NotImplemented;
            return tmp<tensorField>();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
