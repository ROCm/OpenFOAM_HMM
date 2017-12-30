/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

    The presence of the abort file is only checked on the master process.

    Currently the following action types are supported:
    - noWriteNow
    - writeNow
    - nextWrite

    \heading Function object usage
    \table
        Property  | Description         | Required  | Default value
        type      | Type name: abort    | yes       |
        file      | The abort filename  | no        | $FOAM_CASE/name
        action    | Abort action        | no        | nextWrite
    \endtable

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

        //- Reference to the Time
        const Time& time_;

        //- The fully-qualified name of the abort file
        fileName abortFile_;

        //- The type of action
        Time::stopAtControls action_;

        //- Only trigger action once
        bool triggered_;


    // Private Member Functions

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
            const dictionary& dict
        );


    //- Destructor
    virtual ~abort() = default;


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Check existence of abort file and take action
        virtual bool execute();

        //- No-op
        virtual bool write();

        //- Remove abort file after the final time-loop.
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
