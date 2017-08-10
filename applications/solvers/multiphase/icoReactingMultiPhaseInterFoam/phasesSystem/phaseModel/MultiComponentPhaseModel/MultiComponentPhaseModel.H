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
    Foam::MultiComponentPhaseModel

Description
    Class which represents a phase with multiple species. Returns the species'
    mass fractions, and their governing equations.

SourceFiles
    MultiComponentPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MultiComponentPhaseModel_H
#define MultiComponentPhaseModel_H

#include "phaseModel.H"
#include "hashedWordList.H"
#include "rhoReactionThermo.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class phaseThermo>
class MultiComponentPhaseModel
:
    public BasePhaseModel
{
protected:


    // Protected data

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Species table
        hashedWordList species_;

        //- Inert species index
        label inertIndex_;

        //- Thermophysical model
        autoPtr<phaseThermo> thermoPtr_;


public:

    // Constructors

        MultiComponentPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName
        );


    //- Destructor
    virtual ~MultiComponentPhaseModel();


    // Member Functions

        // Access

            //- Species table
            const hashedWordList& species() const
            {
                return species_;
            }


        // Thermo

            //- Access to thermo
            virtual const phaseThermo& thermo() const;

            //- Access non-const thermo
            virtual phaseThermo& thermo();

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);

            //- Constant access the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& Y();

            //- Return heat of formation for species
            virtual dimensionedScalar hf(label index) const;

            //- Correct the mass fractions
            //void correctMassFractions();


        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const;

            //- Return the laminar dynamic viscosity on a patch
            virtual tmp<scalarField> mu(const label patchI) const;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Return the laminar kinematic viscosity on a patch
            virtual tmp<scalarField> nu(const label patchI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MultiComponentPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
