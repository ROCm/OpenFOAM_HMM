/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    XiModel

Description
    Base-class for all Xi models used by the b-Xi combustion model.

SourceFiles
    XiModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiModel_H
#define XiModel_H

#include "IOdictionary.H"
#include "hhuCombustionThermo.H"
#include "RASmodel.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiModel Declaration
\*---------------------------------------------------------------------------*/

class XiModel
{

protected:

    // Protected data

        dictionary XiModelCoeffs_;

        const hhuCombustionThermo& thermo_;
        const compressible::RASmodel& turbulence_;
        const volScalarField& Su_;
        const volScalarField& rho_;
        const volScalarField& b_;
        const surfaceScalarField& phi_;

        //- Flame wrinking field
        volScalarField Xi_;


private:

    // Private Member Functions

        //- Disallow copy construct
        XiModel(const XiModel&);

        //- Disallow default bitwise assignment
        void operator=(const XiModel&);


public:

    //- Runtime type information
    TypeName("XiModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiModel,
            dictionary,
            (
                const dictionary& XiProperties,
                const hhuCombustionThermo& thermo,
                const compressible::RASmodel& turbulence,
                const volScalarField& Su,
                const volScalarField& rho,
                const volScalarField& b,
                const surfaceScalarField& phi
            ),
            (
                XiProperties,
                thermo,
                turbulence,
                Su,
                rho,
                b,
                phi
            )
        );


    // Selectors

        //- Return a reference to the selected Xi model
        static autoPtr<XiModel> New
        (
            const dictionary& XiProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASmodel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        XiModel
        (
            const dictionary& XiProperties,
            const hhuCombustionThermo& thermo,
            const compressible::RASmodel& turbulence,
            const volScalarField& Su,
            const volScalarField& rho,
            const volScalarField& b,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~XiModel();


    // Member Functions

        //- Return the flame-wrinking Xi
        virtual const volScalarField& Xi() const
        {
            return Xi_;
        }

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const
        {
            return turbulence_.muEff();
        }

        //- Add Xi to the multivariateSurfaceInterpolationScheme table
        //  if required
        virtual void addXi
        (
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable&
        )
        {}

        //- Correct the flame-wrinking Xi
        virtual void correct() = 0;

        //- Correct the flame-wrinking Xi using the given convection scheme
        virtual void correct(const fv::convectionScheme<scalar>&)
        {
            correct();
        }

        //- Update properties from given dictionary
        virtual bool read(const dictionary& XiProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
