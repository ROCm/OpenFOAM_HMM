/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::LduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    LduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef LduInterfaceField_H
#define LduInterfaceField_H

#include "lduInterface.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class LduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LduInterfaceField
{
    // Private data

        //- Reference to the coupled patch this field is defined for
        const lduInterface& interface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        LduInterfaceField(const LduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const LduInterfaceField&);


public:

        class Amultiplier
        {
        public:
            
            Amultiplier()
            {}

            virtual void addAmul
            (
                Field<Type>& Apsi,
                const Field<Type>& psi
            ) const = 0;
        };


    //- Runtime type information
    TypeName("LduInterfaceField");


    // Constructors

        //- Construct given interface
        LduInterfaceField(const lduInterface& interface)
        :
            interface_(interface)
        {}


    // Destructor

        virtual ~LduInterfaceField();


    // Member Functions

        // Access

            //- Return the interface
            const lduInterface& interface() const
            {
                return interface_;
            }


        // Coupled interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>& Apsi,
                const Field<Type>& psi,
                const Amultiplier&,
                const Pstream::commsTypes commsType
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& Apsi,
                const Field<Type>& psi,
                const Amultiplier&,
                const Pstream::commsTypes commsType
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
