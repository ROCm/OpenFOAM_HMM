/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::faceShading

Description

    faceShading uses the transmissivity value in the boundaryRadiationProperties
    in order to evaluate which faces are "hit" by the "direction" vector.
    NOTE: Only transmissivity values of zero are considered for opaque walls.

SourceFiles
    faceShading.C

\*---------------------------------------------------------------------------*/

#ifndef faceShading_H
#define faceShading_H

#include "fvMesh.H"
#include "Time.H"
#include "meshTools.H"
#include "DynamicField.H"
#include "labelIOList.H"
#include "wallPolyPatch.H"
#include "triSurfaceTools.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceShading Declaration
\*---------------------------------------------------------------------------*/

class faceShading
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Direction
        vector direction_;

        //- Faces directly hit by vector direction
        labelList rayStartFaces_;


    // Private members

        //- Calculate ray start faces
        void calculate();


        //- Construct a triSurface from patches and faces on global local index
        triSurface triangulate
        (
            const labelHashSet& includePatches,
            const List<labelHashSet>& includeAllFaces
        );

        //- Write rays
        void writeRays
        (
            const fileName& fName,
            const DynamicField<point>& endCf,
            const pointField& myFc
        );

        //- No copy construct
        faceShading(const faceShading&) = delete;

        //- No copy assignment
        void operator=(const faceShading&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("faceShading");


    // Constructors

        //- Construct from components
        faceShading
        (
            const fvMesh& mesh,
            const vector dir,
            const labelList& hitFaceList
        );

        //- Construct from mesh and vector
        faceShading
        (
            const fvMesh& mesh,
            const vector
        );


    //- Destructor
    ~faceShading();


    // Member Functions

        // Access

            //- const access to direction
            const vector direction() const
            {
                return direction_;
            }

            //- Non-const access to direction
            vector& direction()
            {
                return direction_;
            }

            //- Access to rayStartFaces
            const labelList& rayStartFaces() const
            {
                return rayStartFaces_;
            }


        //- Recalculate rayStartFaces using direction vector
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
