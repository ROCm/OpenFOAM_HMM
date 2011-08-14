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
    Foam::directMappedVariableThicknessWallPolyPatch

Description
    Foam::directMappedVariableThicknessWallPolyPatch

SourceFiles
    directMappedVariableThicknessWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedVariableThicknessWallPolyPatch_H
#define directMappedVariableThicknessWallPolyPatch_H

#include "wallPolyPatch.H"
#include "directMappedWallPolyPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
        Class directMappedVariableThicknessWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class directMappedVariableThicknessWallPolyPatch
:
    public directMappedWallPolyPatch
{

    // Private data

        //- Thickness
        scalarList thickness_;


public:

    //- Runtime type information
    TypeName("directMappedWallVariableThickness");


    // Constructors

        //- Construct from components
        directMappedVariableThicknessWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        directMappedVariableThicknessWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const directMappedPatchBase::sampleMode mode,
            const word& samplePatch,
            const vectorField& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from components. Uniform offset.
        directMappedVariableThicknessWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const directMappedPatchBase::sampleMode mode,
            const word& samplePatch,
            const vector& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        directMappedVariableThicknessWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        directMappedVariableThicknessWallPolyPatch
        (
            const directMappedVariableThicknessWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        directMappedVariableThicknessWallPolyPatch
        (
            const directMappedVariableThicknessWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        directMappedVariableThicknessWallPolyPatch
        (
            const directMappedVariableThicknessWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new directMappedVariableThicknessWallPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new directMappedVariableThicknessWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new directMappedVariableThicknessWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~directMappedVariableThicknessWallPolyPatch();


    // Member functions

        //- Return non const thickness
        scalarList& thickness()
        {
            return thickness_;
        }


        //- Return const thickness
        const scalarList& thickness() const
        {
            return thickness_;
        }


        //- Write the polyPatch data as a dictionary
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
