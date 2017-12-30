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
    Foam::functionEntries::includeEntry

Description
    Specify an include file when reading dictionaries, expects a
    single string to follow.

    An example of the \c \#include directive:
    \verbatim
        #include "includeFile"
    \endverbatim

    The usual expansion of environment variables and other constructs
    (eg, the \c ~OpenFOAM/ expansion) is retained.

See also
    fileName, string::expand()

SourceFiles
    includeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeEntry_H
#define includeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeEntry Declaration
\*---------------------------------------------------------------------------*/

class includeEntry
:
    public functionEntry
{
protected:

    // Protected Member Functions

        //- Expand include fileName and return
        static fileName resolveFile
        (
            const fileName& dir,
            const fileName& f,
            const dictionary& dict
        );


public:

    // Static data members

        //- Report to stdout which file is included
        static bool log;


    // Member Functions

        //- Include file in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);

        //- Include file in a primitiveEntry context
        static bool execute
        (
            const dictionary& parentDict,
            primitiveEntry& entry,
            Istream& is
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
