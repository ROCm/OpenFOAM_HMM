/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    volSurfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef volSurfaceMapping_H
#define volSurfaceMapping_H

#include "faMesh.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvPatchField;

/*---------------------------------------------------------------------------*\
                      Class volSurfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class volSurfaceMapping
{
    // Private Data

        //- Reference to mesh
        const faMesh& mesh_;


public:

    // Constructors

        //- Construct from mesh
        volSurfaceMapping(const faMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- No copy construct
        volSurfaceMapping(const volSurfaceMapping&) = delete;

        //- No copy assignment
        void operator=(const volSurfaceMapping&) = delete;


    //- Destructor
    ~volSurfaceMapping() = default;


    // Member Functions

        //- Map Boundary field to surface
        template<class Type>
        tmp<Field<Type>> mapToSurface
        (
            const typename
            GeometricField<Type, fvPatchField, volMesh>::Boundary& df
        ) const;

        //- Map vol Field to surface Field
        template<class Type>
        tmp<Field<Type>> mapToSurface(const Field<Type>& f) const;

        //- Map patch internal field to surface
        template<class Type>
        tmp<Field<Type>> mapInternalToSurface
        (
            const typename
            GeometricField<Type, fvPatchField, volMesh>::Boundary& df
        ) const;

        //- Map surface field to volume boundary field
        template<class Type>
        void mapToVolume
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
        ) const;

        //- Map surface tmp field to volume boundary field
        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
            typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
        ) const;

        //- Map surface field to field
        //  Assumes Field faces in the same order as Boundary
        template<class Type>
        void mapToField
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            Field<Type>& f
        ) const;

        //- Map surface field to volume field
        //  Assumes Field faces in the same order as Boundary
        template<class Type>
        void mapToField
        (
            const Field<Type>& af,
            Field<Type>& f
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volSurfaceMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
