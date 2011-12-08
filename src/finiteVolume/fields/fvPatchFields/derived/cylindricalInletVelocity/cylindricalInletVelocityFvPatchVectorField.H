/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::cylindricalInletVelocityFvPatchVectorField

Description
    Describes an inlet vector boundary condition in cylindrical coordinates
    given a central axis, central point, rpm, axial and radial velocity.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            cylindricalInletVelocity;
        axis            (0 0 1);
        centre          (0 0 0);
        axialVelocity   constant 30;
        radialVelocity  constant -10;
        rpm             constant 100;
    }
    \endverbatim

    The axialVelocity, radialVelocity and rpm entries are DataEntry types, able
    to describe time varying functions.  The example above gives the usage for
    supplying constant values.

SourceFiles
    cylindricalInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalInletVelocityFvPatchVectorField_H
#define cylindricalInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class cylindricalInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class cylindricalInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Central point
        const vector centre_;

        //- Axis
        const vector axis_;

        //- Axial velocity
        autoPtr<DataEntry<scalar> > axialVelocity_;

        //- Radial velocity
        autoPtr<DataEntry<scalar> > radialVelocity_;

        //- RPM
        autoPtr<DataEntry<scalar> > rpm_;


public:

   //- Runtime type information
   TypeName("cylindricalInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        cylindricalInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cylindricalInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        cylindricalInletVelocityFvPatchVectorField
        (
            const cylindricalInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cylindricalInletVelocityFvPatchVectorField
        (
            const cylindricalInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cylindricalInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cylindricalInletVelocityFvPatchVectorField
        (
            const cylindricalInletVelocityFvPatchVectorField&,
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
                new cylindricalInletVelocityFvPatchVectorField(*this, iF)
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
