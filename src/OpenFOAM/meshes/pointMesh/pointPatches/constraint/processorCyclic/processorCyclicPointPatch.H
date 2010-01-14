/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::processorCyclicPointPatch

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides.

    Looking at the creation of the faces on both sides of the processor
    patch they need to be identical on both sides with the normals pointing
    in opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to re-create
    the ordering of patch points on the slave side by reversing all the
    patch faces of the owner.

SourceFiles
    processorCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicPointPatch_H
#define processorCyclicPointPatch_H

#include "processorPointPatch.H"
#include "processorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorCyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicPointPatch
:
    public processorPointPatch
{

        //- Disallow default construct as copy
        processorCyclicPointPatch(const processorCyclicPointPatch&);

        //- Disallow default assignment
        void operator=(const processorCyclicPointPatch&);

public:

    //- Runtime type information
    TypeName(processorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~processorCyclicPointPatch();


    // Member functions

//        //- Return true if running parallel
//        virtual bool coupled() const
//        {
//            if (Pstream::parRun())
//            {
//                return true;
//            }
//            else
//            {
//                return false;
//            }
//        }
//
//        //- Return processor number
//        int myProcNo() const
//        {
//            return procPolyPatch_.myProcNo();
//        }
//
//        //- Return neigbour processor number
//        int neighbProcNo() const
//        {
//            return procPolyPatch_.neighbProcNo();
//        }
//
//        //- Is this a master patch
//        bool isMaster() const
//        {
//            return myProcNo() < neighbProcNo();
//        }
//
//        //- Is this a slave patch
//        bool isSlave() const
//        {
//            return !isMaster();
//        }
//
//        //- Return the underlying processorPolyPatch
//        const processorPolyPatch& procPolyPatch() const
//        {
//            return procPolyPatch_;
//        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
