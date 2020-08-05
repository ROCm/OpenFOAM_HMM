/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::reactingEulerHtcModel

Description
    A heat transfer coefficient for reactingEuler solvers

SourceFiles
    reactingEulerHtcModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactingEulerHtcModel_H
#define reactingEulerHtcModel_H

#include "HashSet.H"
#include "volFields.H"
#include "fieldExpression.H"
#include "runTimeSelectionTables.H"
#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class reactingEulerHtcModel Declaration
\*---------------------------------------------------------------------------*/

class reactingEulerHtcModel
:
    public fieldExpression
{
    // Private Data

        //- Heat transfer coefficient model
        autoPtr<heatTransferCoeffModel> htcModelPtr_;


protected:

    // Protected Member Functions

        //- Calculate the heat transfer coefficient field
        //  \return true on success
        virtual bool calc();

        //- Calculate heat flux
        tmp<FieldField<Field, scalar>> q() const;

        //- No copy construct
        reactingEulerHtcModel(const reactingEulerHtcModel&) = delete;

        //- No copy assignment
        void operator=(const reactingEulerHtcModel&) = delete;


public:

    //- Runtime type information
    TypeName("reactingEulerHtcModel");


    // Constructors

        //- Construct from components
        reactingEulerHtcModel
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reactingEulerHtcModel() = default;


    // Member Functions

        //- Read the heatTransferCoeff data
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
