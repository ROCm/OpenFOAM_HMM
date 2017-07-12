/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::abort

Group
    grpUtilitiesFunctionObjects

Description
    Watches for presence of the named file in the $FOAM_CASE directory
    and aborts the calculation if it is present.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

SourceFiles
    abort.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_abort_H
#define functionObjects_abort_H

#include "functionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class abort Declaration
\*---------------------------------------------------------------------------*/

class abort
:
    public functionObject
{
    // Private data

        //- A subset of Time stopAtControls
        static const Enum<Time::stopAtControls> actionNames_;

        //- Reference to the Time
        const Time& time_;

        //- The fully-qualified name of the abort file
        fileName abortFile_;

        //- The type of action
        Time::stopAtControls action_;


    // Private Member Functions

        //- Remove abort file.
        void removeFile() const;

        //- Disallow default bitwise copy construct
        abort(const abort&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const abort&) = delete;


public:

    //- Runtime type information
    TypeName("abort");


    // Constructors

        //- Construct from Time and dictionary
        abort
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~abort();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Execute, check existence of abort file and take action
        virtual bool execute();

        //- Execute, check existence of abort file and take action
        virtual bool write();

        //- Execute at the final time-loop, used for cleanup
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
