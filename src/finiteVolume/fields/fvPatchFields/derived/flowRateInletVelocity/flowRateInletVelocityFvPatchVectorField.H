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
    Foam::flowRateInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition, derived
    from the flux (volumetric or mass-based), whose direction is assumed
    to be normal to the patch.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, \c phi, and the value of the density field name,
    \c rhoName entry.

    For a mass-based flux:
    - if \c rhoName is a valid density field name, the flow rate is in kg/s
    - if \c rhoName is "none" the flow rate is in m3/s

    For a volumetric-based flux:
    - the flow rate is in m3/s

    /heading Patch usage

    /table
        Property     | Description             | Required    | Default value
        flowRate     | volumetric [m3/s] OR mass flow rate [kg/s] | yes |
    /endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type        flowRateInletVelocity;
        flowRate    0.2;
        rho         rho;
        value       uniform (0 0 0); // placeholder
    }
    \endverbatim

    The \c flowRate entry is a \c DataEntry type, meaning that it can be
    specified as constant, a polynomial fuction of time, and ...

    \note
    - the value is positive into the domain (as an inlet)
    - may not work correctly for transonic inlets
    - strange behaviour with potentialFoam since the U equation is not solved

SeeAlso
    Foam::DataEntry
    Foam::fixedValueFvPatchField

SourceFiles
    flowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateInletVelocityFvPatchVectorField_H
#define flowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class flowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class flowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        autoPtr<DataEntry<scalar> > flowRate_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;


public:

   //- Runtime type information
   TypeName("flowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        flowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        flowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new flowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        flowRateInletVelocityFvPatchVectorField
        (
            const flowRateInletVelocityFvPatchVectorField&,
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
                new flowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
