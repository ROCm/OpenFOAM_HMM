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
    Foam::equationMaxIterCondition

Description
    Maximum number of equation iterations run time condition

SourceFiles
    equationMaxIterCondition.H
    equationMaxIterCondition.C

\*---------------------------------------------------------------------------*/

#ifndef equationMaxIterCondition_H
#define equationMaxIterCondition_H

#include "runTimeCondition.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class equationMaxIterCondition Declaration
\*---------------------------------------------------------------------------*/

class equationMaxIterCondition
:
    public runTimeCondition
{
protected:

    // Protected data

        //- Field name
        const wordList fieldNames_;

        //- Threshold for maximum number of iterations
        const label threshold_;

        //- Start checking from iteration - always skips first iteration
        label startIter_;


public:

    //- Runtime type information
    TypeName("equationMaxIter");

    //- Constructor
    equationMaxIterCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        functionObjectState& state
    );

    //- Destructor
    virtual ~equationMaxIterCondition();


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
