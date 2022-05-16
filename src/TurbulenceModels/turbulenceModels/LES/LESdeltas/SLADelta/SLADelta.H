/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::LESModels::DeltaSLA

Description
    Delta formulation that accounts for the orientation of the vorticity vector
    and a flow-sensitised function to detect 2D flow regions. Delta value is
    strongly reduced in these regions to accelerate the transition from RANS
    to LES in hybrid RANS/LES simulations.

    Reference:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2015).
        An enhanced version of DES with rapid transition from RANS to LES in
        separated flows.
        Flow, Turbulence and Combustion, 95, 709–737, 2015.
    \endverbatim

SourceFiles
    SLADelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaSLADelta_H
#define LESModels_DeltaSLADelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                          Class DeltaSLA Declaration
\*---------------------------------------------------------------------------*/

class DeltaSLA
:
    public LESdelta
{
    // Private data

        //- Run-time selectable delta for hmax
        //  Defaults to the maxDeltaxyz model if not supplied
        autoPtr<LESdelta> hmaxPtr_;

        // Model coefficients

            scalar deltaCoeff_;
            bool requireUpdate_;
            bool shielding_;
            scalar FKHmin_;
            scalar FKHmax_;
            scalar a1_;
            scalar a2_;
            scalar epsilon_;
            scalar kappa_;
            scalar Cd1_;
            scalar Cd2_;


    // Private Member Functions

        //- No copy construct
        DeltaSLA(const DeltaSLA&) = delete;

        //- No copy assignment
        void operator=(const DeltaSLA&) = delete;

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("SLADelta");


    // Constructors

        //- Construct from name, turbulenceModel and IOdictionary
        DeltaSLA
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaSLA() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        void read(const dictionary&);

        //- Return the hmax delta field
        const volScalarField& hmax() const
        {
            return hmaxPtr_();
        }

        // Correct values
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
