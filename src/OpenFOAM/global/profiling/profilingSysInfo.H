/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::profiling::sysInfo

Description
    General system information

SourceFiles
    profilingSysInfo.C

\*---------------------------------------------------------------------------*/

#ifndef profilingSysInfo_H
#define profilingSysInfo_H

#include "profiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

/*---------------------------------------------------------------------------*\
                     Class profiling::sysInfo Declaration
\*---------------------------------------------------------------------------*/

class profiling::sysInfo
{
    // Private Static Data Members


    // Private Data Members


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sysInfo(const sysInfo&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const sysInfo&) = delete;


protected:

    // Friendship

        friend class profiling;


    // Member Functions

public:


    // Constructors

        //- Construct from components
        sysInfo();


    //- Destructor
    ~sysInfo();


    // Member Functions

    // Access


    // Edit

        //- Update it with a new timing information
        void update();

        //- Write the profiling sys-info, use dictionary format.
        Ostream& write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
