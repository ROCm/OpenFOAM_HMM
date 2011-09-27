/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::AMIInterpolation

Description
    Interpolation class dealing with transfer of data between two
    primitive patches with an arbitrary mesh interface (AMI).

    Based on the algorithm given in:

        Conservative interpolation between volume meshes by local Galerkin
        projection, Farrell PE and Maddison JR, 2011, Comput. Methods Appl.
        Mech Engrg, Volume 200, Issues 1-4, pp 89-100


SourceFiles
    AMIInterpolation.C
    AMIInterpolationName.C

\*---------------------------------------------------------------------------*/

#ifndef AMIInterpolation_H
#define AMIInterpolation_H

#include "className.H"
#include "DynamicList.H"
#include "searchableSurface.H"
#include "boolList.H"
#include "primitivePatch.H"
#include "faceAreaIntersect.H"
#include "treeBoundBox.H"
#include "treeBoundBoxList.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AMIInterpolationName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(AMIInterpolation);


/*---------------------------------------------------------------------------*\
                      Class AMIInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class AMIInterpolation
:
    public AMIInterpolationName
{

    //- Helper class for list
    template<class T>
    class ListPlusEqOp
    {
        public:
        void operator()(List<T>& x, const List<T> y) const
        {
            if (y.size())
            {
                if (x.size())
                {
                    label sz = x.size();
                    x.setSize(sz + y.size());
                    forAll(y, i)
                    {
                        x[sz++] = y[i];
                    }
                }
                else
                {
                    x = y;
                }
            }
        }
    };


    // Private data

        // Source patch

            //- Addresses of target faces per source face
            labelListList srcAddress_;

            //- Weights of target faces per source face
            scalarListList srcWeights_;


        // Target patch

            //- Addresses of source faces per target face
            labelListList tgtAddress_;

            //- Weights of wource faces per target face
            scalarListList tgtWeights_;


        //- Starting face seed index
        label startSeedI_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;

        //- Source map pointer - parallel running only
        autoPtr<mapDistribute> srcMapPtr_;

        //- Target map pointer - parallel running only
        autoPtr<mapDistribute> tgtMapPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        AMIInterpolation(const AMIInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const AMIInterpolation&);


        // Helper functions

            //- Write triangle intersection to OBJ file
            void writeIntersectionOBJ
            (
                const scalar area,
                const face& f1,
                const face& f2,
                const pointField& f1Points,
                const pointField& f2Points
            ) const;

            //- Check that patches are valid
            void checkPatches
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            );


        // Parallel functionality

            //- Return true if faces are spread over multiple domains
            bool distributed
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            );

            label calcOverlappingProcs
            (
                const List<treeBoundBoxList>& procBb,
                const treeBoundBox& bb,
                boolList& overlaps
            );

            void distributePatches
            (
                const mapDistribute& map,
                const primitivePatch& pp,
                const globalIndex& gi,
                List<faceList>& faces,
                List<pointField>& points,
                List<labelList>& tgtFaceIDs
            );

            void distributeAndMergePatches
            (
                const mapDistribute& map,
                const primitivePatch& tgtPatch,
                const globalIndex& gi,
                faceList& tgtFaces,
                pointField& tgtPoints,
                labelList& tgtFaceIDs
            );

            autoPtr<mapDistribute> calcProcMap
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            );


        // Initialisation

            //- Project points to surface
            void projectPointsToSurface
            (
                const searchableSurface& surf,
                pointField& pts
            ) const;


        // Marching front

            //- Find face on target patch that overlaps source face
            label findTargetFace
            (
                const label srcFaceI,
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            ) const;

            //- Add faces neighbouring faceI to the ID list
            void appendNbrFaces
            (
                const label faceI,
                const primitivePatch& patch,
                const DynamicList<label>& visitedFaces,
                DynamicList<label>& faceIDs
            ) const;

            //- Set the source and target seed faces
            void setNextFaces
            (
                label& srcFaceI,
                label& tgtFaceI,
                const primitivePatch& srcPatch0,
                const primitivePatch& tgtPatch0,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces
            );


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFaceI,
                const label tgtFaceI,
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            ) const;

            //- Calculate addressing
            void calcAddressing
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                label srcFaceI = -1,
                label tgtFaceI = -1
            );

            //- Normalise the (area) weights - suppresses numerical error in
            //  weights calculation
            //  NOTE: if area weights are incorrect by 'a significant amount'
            //     normalisation may stabilise the solution, but will introduce
            //     numerical error!
            void normaliseWeights
            (
                const primitivePatch& patch,
                const word& patchName,
                const labelListList& addr,
                scalarListList& wght,
                const bool output
            );


public:

    // Constructors

        //- Construct from components
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode
        );

        //- Construct from components, with projection surface
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const autoPtr<searchableSurface>& surf,
            const faceAreaIntersect::triangulationMode& triMode
        );


    //- Destructor
    ~AMIInterpolation();


    // Member Functions

        // Access

            // Source patch

                //- Return const access to source patch addressing
                inline const labelListList& srcAddress();

                //- Return const access to source patch weights
                inline const scalarListList& srcWeights();


            // Target patch

                //- Return const access to target patch addressing
                inline const labelListList& tgtAddress();

                //- Return const access to target patch weights
                inline const scalarListList& tgtWeights();


        // Manipulation

            //- Update addressing and weights
            void update
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch
            );


        // Evaluation

            //- Interpolate from target to source
            template<class Type>
            tmp<Field<Type> > interpolateToSource(const Field<Type>& fld) const;

            //- Interpolate from target tmp field to source
            template<class Type>
            tmp<Field<Type> > interpolateToSource
            (
                const tmp<Field<Type> >& tFld
            ) const;

            //- Interpolate from source to target
            template<class Type>
            tmp<Field<Type> > interpolateToTarget(const Field<Type>& fld) const;

            //- Interpolate from source tmp field to target
            template<class Type>
            tmp<Field<Type> > interpolateToTarget
            (
                const tmp<Field<Type> >& tFld
            ) const;


        // Checks

            //- Write face connectivity as OBJ file
            void writeFaceConnectivity
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const labelListList& srcAddress
            ) const;

            //- Write weights to VTK file
            void writeWeights
            (
                const scalarListList& weights,
                const primitivePatch& patch,
                const word& folder,
                const word& prefix
            )
            const;

            //- Write patch
            void writePatch
            (
                const primitivePatch& patch,
                const word& folder,
                const word& prefix
            )
            const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIInterpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "AMIInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
