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
    Foam::PurePhaseModel

Description
    Class which represents pure phases, i.e. without any species. Returns an
    empty list of mass fractions.

SourceFiles
    PurePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef PurePhaseModel_H
#define PurePhaseModel_H

#include "PtrList.H"
#include "volFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

/*---------------------------------------------------------------------------*\
                           Class PurePhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class phaseThermo>
class PurePhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Empty mass fraction field list
        PtrList<volScalarField> Y_;

        //- Empty thermophysical model Ptr
        autoPtr<phaseThermo> thermoPtr_;


public:

    // Constructors

        PurePhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName
        );


    //- Destructor
    virtual ~PurePhaseModel();


    // Member Functions

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);

        // Thermo

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& Y();

            //- Access to thermo
            virtual const phaseThermo& thermo() const;

             //- Access non-const thermo
            virtual phaseThermo& thermo();

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const;

            //- Access the laminar dynamic viscosity
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return the laminar kinematic viscosity
            virtual tmp<volScalarField> nu() const;

            //- Access the laminar kinematic viscosity
            virtual tmp<scalarField> nu(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PurePhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
