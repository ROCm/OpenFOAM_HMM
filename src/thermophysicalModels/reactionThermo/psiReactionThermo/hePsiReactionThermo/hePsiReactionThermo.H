/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::hePsiReactionThermo

Description
    Foam::hePsiReactionThermo

SourceFiles
    hePsiReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiReactionThermo_H
#define hePsiReactionThermo_H

#include "heThermo.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hePsiReactionThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class hePsiReactionThermo
:
    public heThermo<psiReactionThermo, MixtureType>
{
    // Private Member Functions

        void calculate();

        //- Construct as copy (not implemented)
        hePsiReactionThermo(const hePsiReactionThermo<MixtureType>&);


public:

    //- Runtime type information
    TypeName("hePsiReactionThermo");


    // Constructors

        //- Construct from mesh
        hePsiReactionThermo(const fvMesh&);


    //- Destructor
    virtual ~hePsiReactionThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "hePsiReactionThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
