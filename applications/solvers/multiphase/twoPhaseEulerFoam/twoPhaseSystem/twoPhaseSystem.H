/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
#include "liftModel.H"
#include "heatTransferModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dragModel;
class heatTransferModel;
class liftModel;

/*---------------------------------------------------------------------------*\
                      Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase model 1
        phaseModel phase1_;

        //- Phase model 2
        phaseModel phase2_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- 
        volScalarField dgdt_;

        //- Surface tension coefficient
        dimensionedScalar sigma_;

        //- Virtual mass coefficient
        dimensionedScalar Cvm_;

        //- Drag model for phase 1
        autoPtr<dragModel> drag1_;

        //- Drag model for phase 2
        autoPtr<dragModel> drag2_;

        //- Heat transfer model for phase 1 
        autoPtr<heatTransferModel> heatTransfer1_;

        //- Heat transfer model for phase 2
        autoPtr<heatTransferModel> heatTransfer2_;

        //- Lift model for phase 1
        autoPtr<liftModel> lift1_;

        //- Lift model for phase 2
        autoPtr<liftModel> lift2_;

        //- Name of the dispersed phase, or "both"
        word dispersedPhase_;

        //- Residual phase fraction
        scalar residualPhaseFraction_;

        //- Redisual slip
        dimensionedScalar residualSlip_;


    // Private member functions

        //- Return the mixture flux
        tmp<surfaceScalarField> calcPhi() const;


public:

    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~twoPhaseSystem()
    {}


    // Member Functions

        //- Return the drag coefficient
        tmp<volScalarField> dragCoeff() const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> heatTransferCoeff() const;

        //- Return the lift force
        tmp<volVectorField> liftForce(const volVectorField& U) const;

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Solve for the two-phase-fractions
        void solve();

        //- Correct two-phase properties other than turbulence
        void correct();

        //- Correct two-phase turbulence
        void correctTurbulence();

        //- Read base phaseProperties dictionary
        bool read();

        // Access

            //- Return the mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return phase model 1
            const phaseModel& phase1() const
            {
                return phase1_;
            }

            //- Return non-const access to phase model 1
            phaseModel& phase1()
            {
                return phase1_;
            }

            //- Return phase model 2
            const phaseModel& phase2() const
            {
                return phase2_;
            }

            //- Return non-const access to phase model 2
            phaseModel& phase2()
            {
                return phase2_;
            }

            //- Return the phase not given as an argument
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

            //- Return the mixture flux
            const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return non-const access to the the mixture flux
            surfaceScalarField& phi()
            {
                return phi_;
            }

            //- Return
            const volScalarField& dgdt() const
            {
                return dgdt_;
            }

            //- Return non-const access to the
            volScalarField& dgdt()
            {
                return dgdt_;
            }

            //- Return the drag model for phase 1
            const dragModel& drag1() const
            {
                return drag1_();
            }

            //- Return the drag model for phase 2
            const dragModel& drag2() const
            {
                return drag2_();
            }

            //- Return the drag model for the supplied phase
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

            //- Return non-const access to the residual phase fraction
            scalar residualPhaseFraction() const
            {
                return residualPhaseFraction_;
            }

            //- Return the residual slip
            const dimensionedScalar& residualSlip() const
            {
                return residualSlip_;
            }

            //- Return the heat transfer model for phase 1
            const heatTransferModel& heatTransfer1() const
            {
                return heatTransfer1_();
            }

            //- Return the heat transfer model for phase 2
            const heatTransferModel& heatTransfer2() const
            {
                return heatTransfer2_();
            }

            //- Return the lift model for phase 1
            const liftModel& lift1() const
            {
                return lift1_();
            }

            //- Return the lift model for phase 2
            const liftModel& lift2() const
            {
                return lift2_();
            }

            //- Return the surface tension coefficient
            dimensionedScalar sigma() const
            {
                return sigma_;
            }

            //- Return the virtual-mass coefficient
            dimensionedScalar Cvm() const
            {
                return Cvm_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
