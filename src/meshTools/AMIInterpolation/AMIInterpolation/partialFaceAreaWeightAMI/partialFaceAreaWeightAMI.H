/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

class partialFaceAreaWeightAMI
:
    public faceAreaWeightAMI
{

private:

    // Private Member Functions

        //- No copy assignment
        void operator=(const partialFaceAreaWeightAMI&) = delete;


protected:

    // Protected Member Functions

        //- Set the source and target seed faces
        virtual bool setNextFaces
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

        //- Construct from dictionary
        partialFaceAreaWeightAMI
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        partialFaceAreaWeightAMI
        (
            const bool requireMatch = false,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1,
            const faceAreaIntersect::triangulationMode triMode =
                faceAreaIntersect::tmMesh,
            const bool restartUncoveredSourceFace = true
        );

        //- Construct as copy
        partialFaceAreaWeightAMI(const partialFaceAreaWeightAMI& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return
                autoPtr<AMIInterpolation>(new partialFaceAreaWeightAMI(*this));
        }


    //- Destructor
    virtual ~partialFaceAreaWeightAMI() = default;


    // Member Functions

        //- Flag to indicate that interpolation patches are conformal
        virtual bool conformal() const;

        //- Update addressing, weights and (optional) centroids
        virtual bool calculate
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
