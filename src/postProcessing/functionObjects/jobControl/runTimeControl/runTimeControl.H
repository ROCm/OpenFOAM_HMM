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
    Foam::runTimeControl

Group
    grpJobControlFunctionObjects

Description
    This function object controls when the calculation is terminated based on
    satisfying user-specified conditions.

    Optionally specify a number of write steps before the calculation is
    terminated.  Here, a write is performed each time that all conditons are
    satisfied.

SourceFiles
    runTimeControl.C
    IOrunTimeControl.H

\*---------------------------------------------------------------------------*/

#ifndef runTimeControl_H
#define runTimeControl_H

#include "functionObjectState.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class runTimeCondition;

/*---------------------------------------------------------------------------*\
                       Class runTimeControl Declaration
\*---------------------------------------------------------------------------*/

class runTimeControl
:
    public functionObjectState
{
    // Private data

        //- Reference to the database
        const objectRegistry& obr_;

        //- List of conditions to satisfy
        PtrList<runTimeCondition> conditions_;

        //- Map to define group IDs
        Map<label> groupMap_;

        //- Number of write steps before exiting
        label nWriteStep_;

        //- Current number of steps written
        label writeStepI_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        runTimeControl(const runTimeControl&);

        //- Disallow default bitwise assignment
        void operator=(const runTimeControl&);


public:

    //- Runtime type information
    TypeName("runTimeControl");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        runTimeControl
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~runTimeControl();


    // Member Functions

        //- Return name of the set of runTimeControl
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the runTimeControl data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the runTimeControl and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
