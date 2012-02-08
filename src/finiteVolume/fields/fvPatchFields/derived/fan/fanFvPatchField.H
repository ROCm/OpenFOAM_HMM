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
    Foam::fanFvPatchField

Description
    Jump boundary condition. Operates on cyclic. Jump specified as a DataEntry.
    (table, fileTable, csv etc.)

    Example of the boundary condition specification:
    \verbatim
        fan_half0
        {
            type            fan;
            patchType       cyclic;
            jump            uniform 0;
            jumpTable       csvFile;
            csvFileCoeffs
            {
                hasHeaderLine   1;
                refColumn       0;
                componentColumns 1(1);
                separator       ",";
                fileName        "$FOAM_CASE/constant/pressureVsU";
            }
            value           uniform 0;
        }
    \endverbatim


    Backwards compatibility: if the 'f' keyword is detected it assumes
    it is a power of the flowrate.

    \verbatim
        fan_half0
        {
            type            fan;
            patchType       cyclic;
            jump            uniform 0;
            f               2(100 -0.1);
            value           uniform 0;
        }
    \endverbatim

SourceFiles
    fanFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fanFvPatchField_H
#define fanFvPatchField_H

#include "fixedJumpFvPatchField.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fanFvPatchField
:
    public fixedJumpFvPatchField<Type>
{
    // Private data

        //- Interpolation table
        autoPtr<DataEntry<Type> > jumpTable_;


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
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
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
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new fanFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
fanFvPatchField<scalar>::fanFvPatchField
(
    const fvPatch&,
    const DimensionedField<scalar, volMesh>&,
    const dictionary&
);
template<>
void fanFvPatchField<scalar>::updateCoeffs();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fanFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
