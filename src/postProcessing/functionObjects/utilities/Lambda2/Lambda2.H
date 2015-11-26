/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::Lambda2

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the second largest eigenvalue
    of the sum of the square of the symmetrical and anti-symmetrical parts of
    the velocity gradient tensor.

    Example of function object specification to calculate Lambda2:
    \verbatim
    Lambda2_1
    {
        type        Lambda2;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: Lambda2      | yes         |
        UName        | Name of velocity field  | no          | U
        resultName   | Name of Lambda2 field   | no          | <function name>
        log          | Log to standard output  | no          | yes
    \endtable

SourceFiles
    Lambda2.C
    IOLambda2.H

\*---------------------------------------------------------------------------*/

#ifndef Lambda2_H
#define Lambda2_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class Lambda2 Declaration
\*---------------------------------------------------------------------------*/

class Lambda2
{
    // Private data

        //- Name of this set of Lambda2 objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of velocity field, default is "U"
        word UName_;

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Lambda2(const Lambda2&);

        //- Disallow default bitwise assignment
        void operator=(const Lambda2&);


public:

    //- Runtime type information
    TypeName("Lambda2");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        Lambda2
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~Lambda2();


    // Member Functions

        //- Return name of the set of Lambda2
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the Lambda2 data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the Lambda2 and write
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
