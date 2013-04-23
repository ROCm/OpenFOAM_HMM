/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::faceAreaWeightAMI

Description
    Face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    faceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightAMI_H
#define faceAreaWeightAMI_H

#include "AMIMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class faceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class faceAreaWeightAMI
:
    public AMIMethod<SourcePatch, TargetPatch>
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceAreaWeightAMI(const faceAreaWeightAMI&);

        //- Disallow default bitwise assignment
        void operator=(const faceAreaWeightAMI&);

        // Marching front

            //- Determine overlap contributions for source face srcFaceI
            bool processSourceFace
            (
                const label srcFaceI,
                const label tgtStartFaceI,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label> >& srcAddr,
                List<DynamicList<scalar> >& srcWght,
                List<DynamicList<label> >& tgtAddr,
                List<DynamicList<scalar> >& tgtWght
            );

            //- Attempt to re-evaluate source faces that have not been included
            void restartUncoveredSourceFace
            (
                List<DynamicList<label> >& srcAddr,
                List<DynamicList<scalar> >& srcWght,
                List<DynamicList<label> >& tgtAddr,
                List<DynamicList<scalar> >& tgtWght
            );

            //- Set the source and target seed faces
            void setNextFaces
            (
                label& startSeedI,
                label& srcFaceI,
                label& tgtFaceI,
                const boolList& mapFlag,
                labelList& seedFaces,
                const DynamicList<label>& visitedFaces
            ) const;


        // Evaluation

            //- Area of intersection between source and target faces
            scalar interArea
            (
                const label srcFaceI,
                const label tgtFaceI
            ) const;


public:

    //- Runtime type information
    TypeName("faceAreaWeightAMI");


    // Constructors

        //- Construct from components
        faceAreaWeightAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false
        );


    //- Destructor
    virtual ~faceAreaWeightAMI();


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFaceI = -1,
                label tgtFaceI = -1
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
