/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::partialFaceAreaWeightAMI

Description
    Partial face area weighted Arbitrary Mesh Interface (AMI) method

SourceFiles
    partialFaceAreaWeightAMI.C

\*---------------------------------------------------------------------------*/

#ifndef partialFaceAreaWeightAMI_H
#define partialFaceAreaWeightAMI_H

#include "faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class partialFaceAreaWeightAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class partialFaceAreaWeightAMI
:
    public faceAreaWeightAMI<SourcePatch, TargetPatch>
{

private:

    // Private Member Functions

        //- No copy construct
        partialFaceAreaWeightAMI(const partialFaceAreaWeightAMI&) = delete;

        //- No copy assignment
        void operator=(const partialFaceAreaWeightAMI&) = delete;


protected:

    // Protected Member Functions

        //- Set the source and target seed faces
        virtual void setNextFaces
        (
            label& startSeedi,
            label& srcFacei,
            label& tgtFacei,
            const bitSet& mapFlag,
            labelList& seedFaces,
            const DynamicList<label>& visitedFaces,
            const bool errorOnNotFound = true
        ) const;


public:

    //- Runtime type information
    TypeName("partialFaceAreaWeightAMI");


    // Constructors

        //- Construct from components
        partialFaceAreaWeightAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );


    //- Destructor
    virtual ~partialFaceAreaWeightAMI() = default;


    // Member Functions

        // Access

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;


        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                pointListList& srcCentroids,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            );

            //- Set the face areas for parallel runs
            virtual void setMagSf
            (
                const TargetPatch& tgtPatch,
                const mapDistribute& map,
                scalarList& srcMagSf,
                scalarList& tgtMagSf
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "partialFaceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
