/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    ABLInletEpsilonFvPatchScalarField

Description
    Boundary condition specifies a epsilon inlet for the atmospheric boundary
    layer (ABL). This boundaty is to be used in conjunction with
    ABLInletVelocity.

    @verbatim
        epsilon = Ustar^3 / (K(z + z0))

    where:

        Ustar is the frictional velocity
        K is karman's constant
        z is the verical coordinate
        z0 is the surface roughness lenght

    @endverbatim

    Reference:
    D.M. Hargreaves and N.G. Wright
    "On the use of the k-epsilon model in commercial CFD software to model the
     neutral atmospheric boundary layer"
    Journal of Wind Engineering and Industrial Aerodynamics 95(2007) 355-369.

SourceFiles
    ABLInletEpsilonFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ABLInletEpsilonFvPatchScalarField_H
#define ABLInletEpsilonFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
              Class ABLInletEpsilonFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ABLInletEpsilonFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Peak velocity magnitude
        scalar Ustar_;

        //- Direction of the z-coordinate
        vector z_;

         //- Surface roughness lenght
        scalar z0_;

         //- Von Karman constant
        scalar kappa_;

public:

    //- Runtime type information
    TypeName("ABLInletEpsilon");


    // Constructors

        //- Construct from patch and internal field
        ABLInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ABLInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given ABLInletEpsilonFvPatchScalarField
        //  onto a new patch
        ABLInletEpsilonFvPatchScalarField
        (
            const ABLInletEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ABLInletEpsilonFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ABLInletEpsilonFvPatchScalarField
        (
            const ABLInletEpsilonFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new ABLInletEpsilonFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Return max value
        scalar& Ustar()
        {
            return Ustar_;
        }

        //- Return z direction
        vector& z()
        {
            return z_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
