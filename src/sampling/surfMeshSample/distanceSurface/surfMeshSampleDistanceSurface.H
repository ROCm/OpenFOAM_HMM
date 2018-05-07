/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::surfMeshSampleDistanceSurface

Description
    Sampling surfFields onto a surfMesh based on a plane.
    The cuttingPlane algorithm 'cuts' the mesh.
    The plane is triangulated by default.

Note
    Does not actually cut until update() called.

SourceFiles
    surfMeshSampleDistanceSurface.C
    surfMeshSampleDistanceSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSampleDistanceSurface_H
#define surfMeshSampleDistanceSurface_H

#include "surfMeshSample.H"
#include "distanceSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfMeshSampleDistanceSurface Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSampleDistanceSurface
:
    public surfMeshSample,
    private distanceSurface
{
    // Private typedefs for convenience
    typedef distanceSurface SurfaceSource;

    // Private data

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Sample field on surface
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;


        //- Sample field on surface.
        template<class Type>
        bool sampleType
        (
            const word& fieldName,
            const word& sampleScheme
        ) const;


public:

    //- Runtime type information
    TypeName("surfMeshSampleDistanceSurface");


    // Constructors

        //- Construct from dictionary
        surfMeshSampleDistanceSurface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshSampleDistanceSurface() = default;


    // Member Functions

        //- The surface is from surfMesh
        using surfMeshSample::surface;

        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Sample the volume field onto surface
        virtual bool sample
        (
            const word& fieldName,
            const word& sampleScheme = "cell"
        ) const;


        //- Write
        virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshSampleDistanceSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
