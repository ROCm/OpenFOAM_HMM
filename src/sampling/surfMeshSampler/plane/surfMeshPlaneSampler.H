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

Class
    Foam::surfMeshPlaneSampler

Description
    Sampling surfFields onto a surfMesh based on a plane.
    The cuttingPlane algorithm 'cuts' the mesh.
    The plane is triangulated by default.

Note
    Does not actually cut until update() called.

SourceFiles
    surfMeshPlaneSampler.C
    surfMeshPlaneSamplerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshPlaneSampler_H
#define surfMeshPlaneSampler_H

#include "surfMeshSampler.H"
#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfMeshPlaneSampler Declaration
\*---------------------------------------------------------------------------*/

class surfMeshPlaneSampler
:
    public surfMeshSampler,
    private cuttingPlane
{
    // Private typedefs for convenience
    typedef cuttingPlane SurfaceSource;

    // Private data

        //- If restricted to zones, name of this zone or a regular expression
        keyType zoneKey_;

        //- Triangulated faces or keep faces as is
        const bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Transfer mesh content from SurfaceSource to surfMesh
        void transferContent();


        //- Sample field on surface.
        template<class Type>
        bool sampleType
        (
            const word& fieldName
        ) const;


        //- Sample field on surface
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


public:

    //- Runtime type information
    TypeName("surfMeshPlaneSampler");


    // Constructors

        //- Construct from components
        surfMeshPlaneSampler
        (
            const word& name,
            const polyMesh& mesh,
            const plane& planeDesc,
            const keyType& zoneKey = word::null,
            const bool triangulate = true
        );

        //- Construct from dictionary
        surfMeshPlaneSampler
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshPlaneSampler();


    // Member Functions

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
        virtual bool sample(const word& fieldName) const;


        //- Write
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshPlaneSamplerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
