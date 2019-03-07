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

    // Private data

        //- Flag to restart uncovered source faces
        const bool restartUncoveredSourceFace_;

        //- Storage for src-side triangle decomposition
        List<DynamicList<face>> srcTris_;

        //- Storage for tgt-side triangle decomposition
        List<DynamicList<face>> tgtTris_;


protected:

    // Protected Member Functions

        //- No copy construct
        faceAreaWeightAMI(const faceAreaWeightAMI&) = delete;

        //- No copy assignment
        void operator=(const faceAreaWeightAMI&) = delete;

        // Marching front

            //- Calculate addressing, weights and centroids using temporary
            //- storage
            virtual void calcAddressing
            (
                List<DynamicList<label>>& srcAddress,
                List<DynamicList<scalar>>& srcWeights,
                List<DynamicList<point>>& srcCentroids,
                List<DynamicList<label>>& tgtAddress,
                List<DynamicList<scalar>>& tgtWeights,
                label srcFacei,
                label tgtFacei
            );

            //- Determine overlap contributions for source face srcFacei
            virtual bool processSourceFace
            (
                const label srcFacei,
                const label tgtStartFacei,
                DynamicList<label>& nbrFaces,
                DynamicList<label>& visitedFaces,
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

            //- Attempt to re-evaluate source faces that have not been included
            virtual void restartUncoveredSourceFace
            (
                List<DynamicList<label>>& srcAddr,
                List<DynamicList<scalar>>& srcWght,
                List<DynamicList<point>>& srcCtr,
                List<DynamicList<label>>& tgtAddr,
                List<DynamicList<scalar>>& tgtWght
            );

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


        // Evaluation

            //- Area of intersection between source and target faces
            virtual void calcInterArea
            (
                const label srcFacei,
                const label tgtFacei,
                scalar& area,
                vector& centroid
            ) const;

            //- Return true if faces overlap
            virtual bool overlaps
            (
                const label srcFacei,
                const label tgtFacei,
                const scalar threshold
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
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget = false,
            const bool requireMatch = true,
            const bool restartUncoveredSourceFace = true
        );


    //- Destructor
    virtual ~faceAreaWeightAMI() = default;


    // Member Functions

        // Manipulation

            //- Update addressing, weights and centroids
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
    #include "faceAreaWeightAMI.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
