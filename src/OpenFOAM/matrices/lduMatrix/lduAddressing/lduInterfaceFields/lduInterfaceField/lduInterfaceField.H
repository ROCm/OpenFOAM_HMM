/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::lduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    lduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterfaceField_H
#define lduInterfaceField_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"
#include "Pstream.H"
#include "lduAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class lduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class lduInterfaceField
{
    // Private Data

        //- Reference to the coupled patch this field is defined for
        const lduInterface& interface_;

        //- Update index used so that updateInterfaceMatrix is called only once
        //  during the construction of the matrix
        bool updatedMatrix_;


    // Private Member Functions

        //- No copy construct
        lduInterfaceField(const lduInterfaceField&) = delete;

        //- No copy assignment
        void operator=(const lduInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("lduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        explicit lduInterfaceField(const lduInterface& patch)
        :
            interface_(patch),
            updatedMatrix_(false)
        {}


    //- Destructor
    virtual ~lduInterfaceField() = default;


    // Member Functions

        // Access

            //- Return the interface
            const lduInterface& interface() const
            {
                return interface_;
            }

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return type();
            }


        // Coupled interface matrix update

            //- Whether matrix has been updated
            bool updatedMatrix() const
            {
                return updatedMatrix_;
            }

            //- Whether matrix has been updated
            bool& updatedMatrix()
            {
                return updatedMatrix_;
            }

            //- Is all data available
            virtual bool ready() const
            {
                return true;
            }

            //- Initialise neighbour matrix update.
            //- Add/subtract coupled contributions to matrix
            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality.
            //- Add/subtract coupled contributions to matrix
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const = 0;

            //- Add/subtract weighted contributions to internal field
            template<class Type>
            void addToInternalField
            (
                Field<Type>& result,
                const bool add,
                const labelUList& faceCells,
                const scalarField& coeffs,
                const Field<Type>& vals
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lduInterfaceFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
