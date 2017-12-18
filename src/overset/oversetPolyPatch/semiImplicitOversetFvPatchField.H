/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::semiImplicitOversetFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    Boundary condition for use on overset patches. Implements update of
    interpolated cells inside the linear solvers.

SeeAlso

SourceFiles
    semiImplicitOversetFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitOversetFvPatchField_H
#define semiImplicitOversetFvPatchField_H

#include "oversetFvPatch.H"
#include "zeroGradientFvPatchField.H"
#include "LduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class semiImplicitOversetFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class semiImplicitOversetFvPatchField
:
    public LduInterfaceField<Type>,
    public zeroGradientFvPatchField<Type>
{
protected:

    // Protected data

        //- Local reference cast into the overset patch
        const oversetFvPatch& oversetPatch_;

        //- Master patch ID
        mutable label masterPatchID_;


public:

    //- Runtime type information
    TypeName("semiImplicitOverset");


    // Constructors

        //- Construct from patch and internal field
        semiImplicitOversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        semiImplicitOversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given semiImplicitOversetFvPatchField onto a
        //  new patch
        semiImplicitOversetFvPatchField
        (
            const semiImplicitOversetFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        semiImplicitOversetFvPatchField
        (
            const semiImplicitOversetFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new semiImplicitOversetFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        semiImplicitOversetFvPatchField
        (
            const semiImplicitOversetFvPatchField<Type>&,
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
                new semiImplicitOversetFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return local reference cast into the cyclic AMI patch
            const oversetFvPatch& oversetPatch() const
            {
                return oversetPatch_;
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

             //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                scalarField& result,
                const bool add,
                const scalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "semiImplicitOversetFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
