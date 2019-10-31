/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::functionObjects::runTimeControls::averageCondition

Description
    Average run time condition - satisfied when average does not change by
    more than a given value.

SourceFiles
    averageCondition.H
    averageCondition.C

\*---------------------------------------------------------------------------*/

#ifndef averageCondition_H
#define averageCondition_H

#include "runTimeCondition.H"
#include "Switch.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{

/*---------------------------------------------------------------------------*\
                      Class averageCondition Declaration
\*---------------------------------------------------------------------------*/

class averageCondition
:
    public runTimeCondition
{
public:

    // Public enumerations

        enum class windowType
        {
            NONE,
            APPROXIMATE,
            EXACT
        };

        static const Enum<windowType> windowTypeNames;


protected:

    // Protected data

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Satisfied when difference in mean values is less than this value
        const scalar tolerance_;

        //- Averaging window
        const scalar window_;

        //- Averaging window type
        windowType windowType_;

        //- Average time per field
        List<scalar> totalTime_;

        //- Reset the averaging process on restart flag
        Switch resetOnRestart_;

        //- Number of start-up iterations before allowing satisfied checks
        label nIterStartUp_;

        //- Current iteration count
        label iter_;

    // Protected Member Functions

        //- Templated function to calculate the average
        template<class Type>
        void calc
        (
            const label fieldi,
            bool& satisfied,
            bool& processed
        );


public:

    //- Runtime type information
    TypeName("average");

    //- Constructor
    averageCondition
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        stateFunctionObject& state
    );

    //- Destructor
    virtual ~averageCondition() = default;


    // Public Member Functions

        //- Apply the condition
        virtual bool apply();

        //- Write
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimeControls
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "averageConditionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
