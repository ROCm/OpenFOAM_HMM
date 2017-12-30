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
    Foam::oversetGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetGAMGInterfaceField_H
#define oversetGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "oversetGAMGInterface.H"
#include "processorLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oversetGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class oversetGAMGInterfaceField
:
    public GAMGInterfaceField
{
protected:

    // Protected data

        //- Local reference cast into the interface
        const oversetGAMGInterface& oversetInterface_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetGAMGInterfaceField(const oversetGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const oversetGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        oversetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );


    //- Destructor
    virtual ~oversetGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return oversetInterface_.size();
            }


        // Interface matrix update

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
