/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::RASModels::kEpsilonPhitFBase

Group
    grpRASTurbulence

Description
    Abstract base-class for the \c k-epsilon-phit-f model to provide boundary
    condition access to the \c phit and \c f fields.

See also
    Foam::RASModels::kEpsilonPhitF

SourceFiles
    kEpsilonPhitFBase.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonPhitFBase_H
#define kEpsilonPhitFBase_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kEpsilonPhitFBase Declaration
\*---------------------------------------------------------------------------*/

class kEpsilonPhitFBase
{
public:

    //- Runtime type information
    TypeName("kEpsilonPhitFBase");


    // Constructors

        kEpsilonPhitFBase()
        {}


    //- Destructor
    virtual ~kEpsilonPhitFBase()
    {}


    // Member Functions

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> phit() const = 0;

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> f() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
