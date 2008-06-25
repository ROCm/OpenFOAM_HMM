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
    basic

Description
    Basic sub-grid obstacle drag model.
    Details supplied by J Puttock 2/7/06.

SourceFiles
    basic.C

\*---------------------------------------------------------------------------*/

#ifndef basic_H
#define basic_H

#include "PDRDragModel.H"
#include "XiEqModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRDragModels
{

/*---------------------------------------------------------------------------*\
                          Class basic Declaration
\*---------------------------------------------------------------------------*/

class basic
:
    public PDRDragModel
{
    // Private data

        dimensionedScalar Csu;
        dimensionedScalar Csk;

        volScalarField Aw2_;
        volSymmTensorField CR_;
        volSymmTensorField CT_;


    // Private Member Functions

        //- Disallow copy construct
        basic(const basic&);

        //- Disallow default bitwise assignment
        void operator=(const basic&);


public:

    //- Runtime type information
    TypeName("basic");


    // Constructors

        //- Construct from components
        basic
        (
            const dictionary& PDRProperties,
            const compressible::RASModel& turbulence,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~basic();


    // Member Functions

        //- Return the momentum drag coefficient
        virtual tmp<volSymmTensorField> Dcu() const;

        //- Return the momentum drag turbulence generation rate
        virtual tmp<volScalarField> Gk() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& PDRProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDRDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
