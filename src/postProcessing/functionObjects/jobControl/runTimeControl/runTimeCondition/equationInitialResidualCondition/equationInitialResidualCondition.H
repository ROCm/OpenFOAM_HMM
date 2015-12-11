/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::equationInitialResidualCondition

Description
    Minimum or maximum initial residual run time condition

SourceFiles
    equationInitialResidualCondition.H
    equationInitialResidualCondition.C

\*---------------------------------------------------------------------------*/

#ifndef equationInitialResidualCondition_H
#define equationInitialResidualCondition_H

#include "runTimeCondition.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class equationInitialResidualCondition Declaration
\*---------------------------------------------------------------------------*/

class equationInitialResidualCondition
:
    public runTimeCondition
{
public:

    enum operatingMode
    {
        omMin,
        omMax
    };

    static const NamedEnum<operatingMode, 2> operatingModeNames;


protected:

    // Protected data

        //- Field name
        const wordList fieldNames_;

        //- Value to compare
        const scalar value_;

        //- Start checking from time - always skips first iteration
        scalar timeStart_;

        //- Operating mode
        operatingMode mode_;


public:

    //- Runtime type information
    TypeName("equationInitialResidual");

    //- Constructor
    equationInitialResidualCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        functionObjectState& state
    );

    //- Destructor
    virtual ~equationInitialResidualCondition();


    // Public Member Functions

        //- Apply the condition
        virtual bool apply();

        //- Write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
