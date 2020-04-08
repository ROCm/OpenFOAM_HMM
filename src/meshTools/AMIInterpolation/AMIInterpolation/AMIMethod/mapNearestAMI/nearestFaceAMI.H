/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::nearestFaceAMI

Description
    Nearest-face Arbitrary Mesh Interface (AMI) method

SourceFiles
    nearestFaceAMI.C

\*---------------------------------------------------------------------------*/

#ifndef nearestFaceAMI_H
#define nearestFaceAMI_H

#include "AMIMethod.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class nearestFaceAMI Declaration
\*---------------------------------------------------------------------------*/

template<class SourcePatch, class TargetPatch>
class nearestFaceAMI
:
    public AMIMethod<SourcePatch, TargetPatch>
{
public:

    typedef Tuple2<pointIndexHit, scalar> nearestAndDist;

    //- Helper class for finding nearest
    class nearestEqOp
    {

    public:

        void operator()(nearestAndDist& x, const nearestAndDist& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second() < x.second())
                {
                    x = y;
                }
            }
        }
    };


private:

    // Private Data

        //- Maximum squared distance
        scalar maxDistance2_;


    // Private Member Functions

        //- No copy construct
        nearestFaceAMI(const nearestFaceAMI&) = delete;

        //- No copy assignment
        void operator=(const nearestFaceAMI&) = delete;

        autoPtr<mapDistribute> calcFaceMap
        (
            const List<nearestAndDist>& localInfo,
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch
        ) const;

        autoPtr<mapDistribute> calcDistributed
        (
            const SourcePatch& src,
            const TargetPatch& tgt,
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght
        ) const;


public:

    //- Runtime type information
    TypeName("nearestFaceAMI");


    // Constructors

        //- Construct from components
        nearestFaceAMI
        (
            const SourcePatch& srcPatch,
            const TargetPatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true
        );


    //- Destructor
    virtual ~nearestFaceAMI() = default;


    // Member Functions

        // Manipulation

            //- Update addressing and weights
            virtual bool calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                pointListList& srcCentroids,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                scalarList& srcMagSf,
                scalarList& tgtMagSf,
                autoPtr<mapDistribute>& srcMapPtr,
                autoPtr<mapDistribute>& tgtMapPtr,
                label srcFacei = -1,
                label tgtFacei = -1
            );

            //- Normalise the weight. Can optionally subset addressing
            //- (e.g. for mapNearest)
            virtual void normaliseWeights
            (
                const bool verbose,
                AMIInterpolation<SourcePatch, TargetPatch>& inter
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nearestFaceAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
