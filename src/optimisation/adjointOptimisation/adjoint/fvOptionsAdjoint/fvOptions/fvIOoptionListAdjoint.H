/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::fv::IOoptionListAdjoint

Description
    IOoptionListAdjoint

SourceFiles
    IOoptionListAdjoint.C

\*---------------------------------------------------------------------------*/

#ifndef IOoptionListAdjoint_H
#define IOoptionListAdjoint_H

#include "fvOptionList.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class IOoptionListAdjoint Declaration
\*---------------------------------------------------------------------------*/

class IOoptionListAdjoint
:
    public IOdictionary,
    public optionList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        IOoptionListAdjoint(const IOoptionListAdjoint&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const IOoptionListAdjoint&) = delete;


public:

    // Constructors

        //- Construct from components with list of field names
        IOoptionListAdjoint(const fvMesh& mesh);


        //- Destructor
        virtual ~IOoptionListAdjoint() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
