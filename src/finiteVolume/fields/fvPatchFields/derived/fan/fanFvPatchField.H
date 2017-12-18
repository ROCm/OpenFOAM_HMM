/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fanFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base.

    The jump is specified as a \c Function1 type, to enable the use of, e.g.
    contant, polynomial, table values.

    The switch nonDimensional can be used for a non-dimensional table. It needs
    inputs rpm and dm of the fan. It should be used with uniformJump true.

    The nonDimensional flux for the table is calculate as :

        phi = 4.0*mDot/(rho*sqr(PI)*dm^3*omega)
        where:
            dm is the mean diameter.
            omega is rad/sec.

    The nonDimensinal pressure :

        Psi = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop

    The non-dimensional table should be given as Psi = F(phi).

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclic| yes |
        jumpTable    | jump data, e.g. \c csvFile | yes      |
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | none
        uniformJump  | applies a uniform pressure based on the averaged
                       velocity | no | false
        nonDimensional | uses non-dimensional table | no | false
        rpm          | fan rpm for non-dimensional table | no | 0.0
        dm           | mean diameter for non-dimensional table | no | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fan;
        patchType       cyclic;
        jumpTable       csvFile;
        csvFileCoeffs
        {
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            file            "$FOAM_CASE/constant/pressureVsU";
        }
        value           uniform 0;
    }
    \endverbatim

    The above example shows the use of a comma separated (CSV) file to specify
    the jump condition.

Note
    The underlying \c patchType should be set to \c cyclic

See also
    Foam::Function1Types

SourceFiles
    fanFvPatchField.C
    fanFvPatchFields.H
    fanFvPatchFields.C
    fanFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef fanFvPatchField_H
#define fanFvPatchField_H

#include "uniformJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fanFvPatchField
:
    public uniformJumpFvPatchField<Type>
{
    // Private data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Uniform pressure drop
        bool uniformJump_;

        //- Swtich for using non-dimensional curve
        Switch nonDimensional_;

        // Parameters for non-dimensional table

            //- Fan rpm
            scalar rpm_;

            //- Fan mean diameter
            scalar dm_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("fan");


    // Constructors

        //- Construct from patch and internal field
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanFvPatchField onto a new patch
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanFvPatchField
        (
            const fanFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fanFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fanFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
void fanFvPatchField<scalar>::calcFanJump();

template<>
fanFvPatchField<scalar>::fanFvPatchField
(
    const fvPatch&,
    const DimensionedField<scalar, volMesh>&,
    const dictionary&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fanFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
