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

Class
    Foam::surfMeshSamplePlane

Description
    Sampling surfFields onto a surfMesh based on a plane.
    The cuttingPlane algorithm 'cuts' the mesh.
    The plane is triangulated by default.

    Example of function object partial specification:
    \verbatim
    surfaces
    (
        surface1
        {
            type        plane;
            planeType   pointAndNormal;
            pointAndNormalDict
            {
                ...
            }
        }
    );
    \endverbatim

    Where the sub-entries comprise:
    \table
        Property | Description                             | Required | Default
        type     | plane                                   | yes      |
        planeType | plane description (pointAndNormal etc) | yes      |
        triangulate | triangulate faces                    | no       | true
        bounds   | limit with bounding box                 | no       |
        zone     | limit to cell zone (name or regex)      | no       |
        zones    | limit to cell zones (names, regexs)     | no       |
        coordinateSystem | define plane within given coordinate system | no |
    \endtable

Note
    Does not actually cut until update() called.

SourceFiles
    surfMeshSamplePlane.C
    surfMeshSamplePlaneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSamplePlane_H
#define surfMeshSamplePlane_H

#include "surfMeshSample.H"
#include "cuttingPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfMeshSamplePlane Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSamplePlane
:
    public surfMeshSample,
    private cuttingPlane
{
    // Private typedefs for convenience
    typedef cuttingPlane SurfaceSource;

    // Private data

        //- The zone or zones in which cutting is to occur
        wordRes zoneNames_;

        //- Optional bounding box to trim against
        const boundBox bounds_;

        //- Triangulate faces or not
        const bool triangulate_;

        //- Track if the surface needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Define cell selection from zones and bounding box.
        //  Optionally check and warn if the plane does not intersect
        //  with the bounds of the mesh (or submesh) or if the bounding box
        //  does not overlap with the mesh (or submesh)
        bitSet cellSelection(const bool warn=false) const;


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
    TypeName("surfMeshSamplePlane");


    // Constructors

        //- Construct from components
        surfMeshSamplePlane
        (
            const word& name,
            const polyMesh& mesh,
            const plane& planeDesc,
            const wordRes& zones = wordRes(),
            const bool triangulate = true
        );

        //- Construct from dictionary
        surfMeshSamplePlane
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshSamplePlane() = default;


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
    #include "surfMeshSamplePlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
