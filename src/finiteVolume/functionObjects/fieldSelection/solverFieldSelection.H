/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::volFieldSelection

Description
    Helper class to manage volume field selections

SourceFiles
    volFieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_solverFieldSelection_H
#define functionObjects_solverFieldSelection_H

#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class solverFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class solverFieldSelection
:
    public volFieldSelection
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        solverFieldSelection(const solverFieldSelection&) = delete;


public:

    //- Construct from object registry
    solverFieldSelection(const objectRegistry& obr);


    //- Destructor
    virtual ~solverFieldSelection();


    // Member Functions

        //- Update the selection using current contents of obr_
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
