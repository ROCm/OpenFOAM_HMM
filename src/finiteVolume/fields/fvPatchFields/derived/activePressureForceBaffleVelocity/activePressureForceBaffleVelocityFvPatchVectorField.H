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
    Foam::activePressureForceBaffleVelocityFvPatchVectorField

Description
    This boundary condition is applied to the flow velocity, to simulate the
    opening of a baffle due to local flow conditions, by merging the behaviours
    of wall and cyclic conditions.

    The baffle joins two mesh regions, where the open fraction determines
    the interpolation weights applied to each cyclic- and neighbour-patch
    contribution.

    Once opened the baffle continues to open at a fixed rate using

        \f[
            x = x_old + dt/DT
        \f]

    where
        \var x      = baffle open fraction [0-1]
        \var x_old  = baffle open fraction on previous evaluation
        \var dt     = simulation time step
        \var DT     = time taken to open the baffle

    The open fraction is then applied to scale the patch areas.

    Example of the boundary condition specification:
    \verbatim
        myPatch
        {
            type            activePressureForceBaffleVelocity;
            p               p;          // name of pressure field (default = p)
            cyclicPatch     cyclic1;    // name of cyclic patch
            orientation     1;          // 1 or -1 used to switch flow direction
            openFraction    0.2;        // current opatch open fraction [0-1]
            openingTime     5.0;        // time taken to open the baffle
            maxOpenFractionDelta 0.1;   // max open fraction change/timestep
            minThresholdValue 0.01;     // minimum open fraction for activation
            forceBased      false;      // force (true) or pressue (false)
        }
    \endverbatim

SourceFiles
    activePressureForceBaffleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef activePressureForceBaffleVelocityFvPatchVectorField_H
#define activePressureForceBaffleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class activePressureForceBaffleVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class activePressureForceBaffleVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the pressure field used to calculate the force
        //  on the active baffle
        word pName_;

        //- Name of the cyclic patch used when the active baffle is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active baffle is open
        label cyclicPatchLabel_;

        //- Orientation (1 or -1) of the active baffle patch.
        //  Used to change the direction of opening without the need for
        //  reordering the patch faces
        label orientation_;

        //- Initial wall patch areas
        vectorField initWallSf_;

        //- Initial cyclic patch areas
        vectorField initCyclicSf_;

        //- Initial neighbour-side cyclic patch areas
        vectorField nbrCyclicSf_;

        //- Current fraction of the active baffle which is open
        scalar openFraction_;

        //- Time taken for the active baffle to open
        scalar openingTime_;

        //- Maximum fractional change to the active baffle openness
        //  per time-step
        scalar maxOpenFractionDelta_;

        label curTimeIndex_;

        //- Minimum value for the active baffle to start opening
        scalar minThresholdValue_;

        //- Force based active baffle
        bool fBased_;

        //- Baffle is activated
        bool baffleActivated_;


public:

    //- Runtime type information
    TypeName("activePressureForceBaffleVelocity");


    // Constructors

        //- Construct from patch and internal field
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new activePressureForceBaffleVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        activePressureForceBaffleVelocityFvPatchVectorField
        (
            const activePressureForceBaffleVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new activePressureForceBaffleVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
