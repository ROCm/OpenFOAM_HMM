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
    Foam::pimpleLoop

Description
    PIMPLE loop class to formalise the iteration and automate the handling
    of the "finalIteration" mesh data entry.

\*---------------------------------------------------------------------------*/

#ifndef pimpleLoop_H
#define pimpleLoop_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pimpleLoop Declaration
\*---------------------------------------------------------------------------*/

class pimpleLoop
{
    // Private data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Number of PIMPLE correctors
        const int nCorr_;

        //- Current PIMPLE corrector
        int corr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pimpleLoop(const pimpleLoop&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleLoop&);


public:

    // Constructors

        //- Construct from components
        pimpleLoop(fvMesh& mesh, const int nCorr)
        :
            mesh_(mesh),
            nCorr_(nCorr),
            corr_(0)
        {}


    //- Destructor
    ~pimpleLoop()
    {}


    // Member Functions

        bool loop()
        {
            if (finalIter())
            {
                mesh_.data::add("finalIteration", true);
            }

            return corr_ < nCorr_;
        }

        bool finalIter() const
        {
            return corr_ == nCorr_-1;
        }


    // Member Operators

        void operator++(int)
        {
            if (finalIter())
            {
                mesh_.data::remove("finalIteration");
            }

            corr_++;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
