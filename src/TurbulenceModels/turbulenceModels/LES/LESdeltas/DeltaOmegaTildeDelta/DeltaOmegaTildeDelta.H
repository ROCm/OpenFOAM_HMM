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
    Foam::DeltaOmegaTildeDelta

Description
    Delta formulation that accounts for the orientation of the vorticity
    vector. In "2D-regions" (i.e. xy-plane), delta is of the order
    max(delta_x,delta_y), so that the influence of delta_z is discarded.
    This can help to accelerate the transition from RANS to LES in hybrid
    RANS/LES simulations.

    Reference:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2015).
        An enhanced version of DES with rapid transition from RANS to LES in
        separated flows.
        Flow, Turbulence and Combustion, 95, 709-737, 2015.
    \endverbatim

SourceFiles
    DeltaOmegaTildeDelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaOmegaTildeDelta_H
#define LESModels_DeltaOmegaTildeDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeltaOmegaTildeDelta Declaration
\*---------------------------------------------------------------------------*/

class DeltaOmegaTildeDelta
:
    public LESdelta
{
    // Private data

        //- Run-time selectable delta for hmax
        //  Defaults to the maxDeltaxyz model if not supplied
        autoPtr<LESdelta> hmaxPtr_;

        //- Model coefficient
        scalar deltaCoeff_;

        //- Flag to indicate whether hmax requires updating
        bool requireUpdate_;


    // Private Member Functions

        //- No copy construct
        DeltaOmegaTildeDelta(const DeltaOmegaTildeDelta&) = delete;

        //- No copy assignment
        void operator=(const DeltaOmegaTildeDelta&) = delete;

        //- Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("DeltaOmegaTilde");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaOmegaTildeDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaOmegaTildeDelta() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

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
