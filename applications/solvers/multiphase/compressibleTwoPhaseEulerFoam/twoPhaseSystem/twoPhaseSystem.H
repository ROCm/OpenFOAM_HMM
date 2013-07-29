/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::twoPhaseSystem

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unsed in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    twoPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSystem_H
#define twoPhaseSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dragModel;
class heatTransferModel;


/*---------------------------------------------------------------------------*\
                      Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public IOdictionary
{
    // Private data

        const fvMesh& mesh_;

        phaseModel phase1_;
        phaseModel phase2_;

        dimensionedScalar Cvm_;
        dimensionedScalar Cl_;

        autoPtr<dragModel> drag1_;
        autoPtr<dragModel> drag2_;

        autoPtr<heatTransferModel> heatTransfer1_;
        autoPtr<heatTransferModel> heatTransfer2_;

        word dispersedPhase_;
        scalar residualPhaseFraction_;
        dimensionedScalar residualSlip_;


public:

    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~twoPhaseSystem()
    {}


    // Member Functions

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const phaseModel& phase1() const
        {
            return phase1_;
        }

        const phaseModel& phase2() const
        {
            return phase2_;
        }

        const phaseModel& otherPhase(const phaseModel& phase) const
        {
            if (&phase == &phase1_)
            {
                return phase2_;
            }
            else
            {
                return phase1_;
            }
        }

        phaseModel& phase1()
        {
            return phase1_;
        }

        phaseModel& phase2()
        {
            return phase2_;
        }

        const dragModel& drag1() const
        {
            return drag1_();
        }

        const dragModel& drag2() const
        {
            return drag2_();
        }

        const dragModel& drag(const phaseModel& phase) const
        {
            if (&phase == &phase1_)
            {
                return drag1_();
            }
            else
            {
                return drag2_();
            }
        }

        scalar residualPhaseFraction() const
        {
            return residualPhaseFraction_;
        }

        const dimensionedScalar& residualSlip() const
        {
            return residualSlip_;
        }

        tmp<volScalarField> dragCoeff() const;
        tmp<volVectorField> liftForce(const volVectorField& U) const;

        const heatTransferModel& heatTransfer1() const
        {
            return heatTransfer1_();
        }

        const heatTransferModel& heatTransfer2() const
        {
            return heatTransfer2_();
        }

        tmp<volScalarField> heatTransferCoeff() const;

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the mixture flux
        tmp<surfaceScalarField> phi() const;

        //- Return the virtual-mass coefficient
        dimensionedScalar Cvm() const
        {
            return Cvm_;
        }

        //- Return the lift coefficient
        dimensionedScalar Cl() const
        {
            return Cl_;
        }

        //- Dummy correct
        void correct()
        {}

        //- Read base phaseProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
