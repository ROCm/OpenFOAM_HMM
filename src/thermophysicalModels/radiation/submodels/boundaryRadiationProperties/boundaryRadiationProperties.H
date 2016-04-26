/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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
    Foam::radiation::boundaryRadiationProperties

Description
    Boundary radiation properties holder

SourceFiles
    boundaryRadiationProperties.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationProperties_H
#define boundaryRadiationProperties_H

#include "MeshObject.H"
#include "GeometricField.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace radiation
{

/*---------------------------------------------------------------------------*\
                   Class boundaryRadiationProperties Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationProperties
:
    public MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >
{

    // Private data

        //- AutoPtr to volSacalarField properties
        autoPtr<volScalarField> radBoundaryProperties_;


    // Private member functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh, const word) const;


public:

    // Declare name of the class and its debug switch
    TypeName("boundaryRadiationProperties");


    // Constructors

        //- Construct given fvMesh
        boundaryRadiationProperties(const fvMesh&);


    // Member Functions

        //- Access boundary emissivity on patch
        tmp<scalarField> emissivity
        (
            const label patchId,
            const label bandI = 0
        ) const;


        //- Access boundary absorptivity on patch
        tmp<scalarField> absorptivity
        (
            const label patchId,
            const label bandI = 0
        ) const;


        //- Access boundary transmissivity on patch
        tmp<scalarField> transmissivity
        (
            const label patchId,
            const label bandI = 0
        ) const;

        //- Access boundary reflectivity on patch
        tmp<scalarField> reflectivity
        (
            const label patchId,
            const label bandI = 0
        ) const;


        //- Access to radBoundaryProperties
        const volScalarField& radBoundaryProperties() const;


    //- Destructor
    virtual ~boundaryRadiationProperties();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
