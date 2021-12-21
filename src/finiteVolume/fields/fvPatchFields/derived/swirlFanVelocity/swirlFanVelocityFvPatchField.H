/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::swirlFanVelocityFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition for U across a
    cyclic pressure jump condition and applies a transformation to U.

    The U-jump is specified with a swirl component as follows:
    \verbatim
         Utan = deltaP/rEff/fanEff/(rpm*pi/30.0);

         where

            deltaP : pressure drop across the cyclic.
            rEff   : effective radius
            fanEff : fan efficiency coefficient
            rpm    : RPM of the fan
    \endverbatim

    Alternatively an inner and outer radii can be used instead of rEff. The
    Utan is as follow for r > rInner and r < rOuter
    \verbatim
            Utan = deltaP/r/fanEff/(rpm/pi/30.0);

            where

                r : p - origin, p is the face center
    \endverbatim

    Outside rInner and rOuter, Utan = 0. The input for this mode is:
    \verbatim
            useRealRadius   true;
            rInner          0.005;
            rOuter          0.01;
    \endverbatim

    The radial velocity is zero in the present model.

Usage
    \table
        Property     | Description             | Required    | Default
        patchType    | underlying patch type should be \c cyclic| yes |
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        p            | pressure field name     | no          | p
        origin       | fan centre              | no          | calculated
        rpm          | RPM of the fan          | yes
        fanEff       | Fan efficiency          | no          | 1
        rEff         | Effective radius        | no          | 0
        useRealRadius| Use inner/outer radii   | no          | false
        rInner       | Inner radius            | no          | 0
        rOuter       | Outer radius            | no          | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        cyclicFaces_master
        {
            type            swirlFanVelocity;
            patchType       cyclic;
            jump            uniform (0 0 0);
            value           uniform (0 0 0);
            rpm             1000;
            rEff            0.01;
        }
    }
    \endverbatim

SourceFiles
    swirlFanVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlFanVelocityFvPatchField_H
#define swirlFanVelocityFvPatchField_H

#include "fixedJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class swirlFanVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class swirlFanVelocityFvPatchField
:
     public fixedJumpFvPatchField<vector>
{
    // Private data

        //- Name of the flux field
        const word phiName_;

        //- Name of the pressure field
        const word pName_;

        //- Name of the rho field
        const word rhoName_;

        //- Origin of the rotation
        const vector origin_;

        //- Fan rpm
        autoPtr<Function1<scalar>> rpm_;

        //- Fan efficiency
        scalar fanEff_;

        //- Effective fan radius
        scalar rEff_;

        //- Inner radius
        scalar rInner_;

        //- Outer radius
        scalar rOuter_;

        //- Switch to use effective radius or inner and outer radius
        bool useRealRadius_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("swirlFanVelocity");


    // Constructors

        //- Construct from patch and internal field
        swirlFanVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swirlFanVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given swirlFanVelocityFvPatchField
        //- onto a new patch
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return tmp<fvPatchField<vector>>
            (
                new swirlFanVelocityFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swirlFanVelocityFvPatchField
        (
            const swirlFanVelocityFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new swirlFanVelocityFvPatchField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
