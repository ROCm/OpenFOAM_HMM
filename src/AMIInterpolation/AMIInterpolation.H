/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    primitive patches with an arbitrary mesh interface (AMI)

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
    // Static data

        // Source patch

            //- Addressing
            labelListList srcAddress_;

            //- Weighting factors
            scalarListList srcWeights_;


        // Target patch

            //- Addressing
            labelListList tgtAddress_;

            //- Weighting factors
            scalarListList tgtWeights_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        AMIInterpolation(const AMIInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const AMIInterpolation&);

        //- Write triangle intersection to OBJ file
        void writeIntersectionOBJ
        (
            const scalar area,
            const face& f1,
            const face& f2,
            const pointField& f1Points,
            const pointField& f2Points
        ) const;

        void checkPatches
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        );

        void projectPointsToSurface
        (
            const searchableSurface& surf,
            pointField& pts
        ) const;

        label findTargetFace
        (
            const label srcFaceI,
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        ) const;

        void appendNbrFaces
        (
            const label faceI,
            const primitivePatch& patch,
            const DynamicList<label>& visitedFaces,
            DynamicList<label>& faceIDs
        ) const;

        void setNextFaces
        (
            label& srcFaceI,
            label& tgtFaceI,
            const primitivePatch& srcPatch0,
            const primitivePatch& tgtPatch0,
            const boolList& mapFlag,
            labelList& seedFaces,
            const DynamicList<label>& visitedFaces
        ) const;

        scalar interArea
        (
            const label srcFaceI,
            const label tgtFaceI,
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        ) const;

        void normaliseWeights
        (
            const word& patchType,
            const primitivePatch& patch,
            const pointField& points,
            const List<DynamicList<label> >& addr,
            List<DynamicList<scalar> >& wght
        );


public:

    // Constructors

        //- Construct from components
        AMIInterpolation
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const searchableSurface& surf
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

            //- Calculate addressing
            void calcAddressing
            (
                const primitivePatch& srcPatch,
                const primitivePatch& tgtPatch,
                const searchableSurface& surf
            );

            //- Correct weighting factors for moving mesh
            bool movePoints();


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
