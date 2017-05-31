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
    Foam::ccm::base

Description
    Base functionality common to reader and writer classes

Note
    this class is in development
    - any/all of the class names and members may change

SourceFiles
    ccmBase.C

\*---------------------------------------------------------------------------*/

#ifndef ccmBase_H
#define ccmBase_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class ccmGlobalState;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                          Class ccm::base Declaration
\*---------------------------------------------------------------------------*/

class base
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        base(const base&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const base&) = delete;


protected:

    // Protected Data

        //- Maintain overall global states (error, root-node)
        ccmGlobalState* globalState_;


    // Protected Member Functions

        //- Die with msg if there is an error
        //  Return true if there is no error
        static bool assertNoError(int err, const char *msg);

        //- Die with msg if there is an error
        //  Return true if there is no error
        static bool assertNoError(int err, const std::string& msg);


        //- check global state for errors and die as required
        //  Return true if there is no error
        bool assertNoError(const char* msg) const;

        //- check global state for errors and die as required
        //  Return true if there is no error
        bool assertNoError(const std::string& msg) const;


public:

    // Constructors

        //- Construct null
        base();


    //- Destructor: close file
    ~base();


    // Member Functions

      // Access

        //- Explicity close the file and terminate ccmio access.
        //  Return false if it was already closed.
        bool close();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
