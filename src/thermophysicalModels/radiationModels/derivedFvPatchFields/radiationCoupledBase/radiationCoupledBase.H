/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    radiationCoupledBase

Description
    Common functions for use in temperature coupled boundaries when
    radiation is used. For now only emissivity() : emissivity
    gets supplied from lookup or calculate:

    - 'lookup' : lookup volScalarField with name
    - 'solidThermo' : use basicSolidThermo emissivity()

SourceFiles
    radiationCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCoupledBase_H
#define radiationCoupledBase_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class radiationCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class radiationCoupledBase
{
public:
        //- Type of supplied emissivity
        enum emissivityMethodType
        {
            SOLIDTHERMO,
            LOOKUP
        };

private:

    // Private data

        static const NamedEnum<emissivityMethodType, 2>
            emissivityMethodTypeNames_;

        //- Underlying patch
        const fvPatch& patch_;

        //- How to get emissivity
        const emissivityMethodType method_;

         //- Emissivity
         // Cached locally when is read from dictionary
        scalarField emissivity_;


public:

    // Constructors

        //- Construct from patch and emissivity name
        radiationCoupledBase
        (
            const fvPatch& patch,
            const word& calculationMethod
        );

        //- Construct from patch and dictionary
        radiationCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );


    // Member functions

        //- Method to obtain emissivity
        word emissivityMethod() const
        {
            return emissivityMethodTypeNames_[method_];
        }

/*
        //- Name of thermal conductivity field
        const word& emissivityName() const
        {
            return emissivityName_;
        }
*/
        //- Calculate corresponding emissivity field
        tmp<scalarField> emissivity() const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
