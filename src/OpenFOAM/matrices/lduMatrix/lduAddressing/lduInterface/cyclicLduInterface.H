/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::cyclicLduInterface

Description
    An abstract base class for cyclic coupled interfaces.

SourceFiles
    cyclicLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLduInterface_H
#define cyclicLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicLduInterface Declaration
\*---------------------------------------------------------------------------*/

class cyclicLduInterface
{
    // Private data

        //- Coupled patch
        const label neighbPatchID_;

//        const cyclicLduInterface* neighbPatch_;

public:

    //- Runtime type information
    TypeName("cyclicLduInterface");

    // Constructors

        //- Construct from components
        cyclicLduInterface(const label neighbPatchID);

    // Destructor

        virtual ~cyclicLduInterface();


    // Member Functions

        // Access

            //- Return processor number
            label neighbPatchID() const
            {
                return neighbPatchID_;
            }

//            //- Return processor number
//            const cyclicLduInterface& neighbPatch() const
//            {
//                if (!neighbPatch_)
//                {
//                    FatalErrorIn("cyclicLduInterface::neighbPatch() const")
//                        << "Not owner." << abort(FatalError);
//                }
//                return neighbPatchID_;
//            }

            //- Return face transformation tensor
            virtual const tensor& forwardT() const = 0;

            //- Return face reverse transformation tensor
            virtual const tensor& reverseT() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
