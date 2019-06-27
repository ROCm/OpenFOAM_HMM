/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2016 OpenFOAM Foundation
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
    Foam::functionEntries::includeEtcEntry

Description
    A dictionary directive for including a file found using the
    Foam::findEtcFile() mechanism.

    Specify an etc file to include when reading dictionaries, expects a
    single string to follow.

    An example of the \c \#includeEtc directive:
    \verbatim
        #includeEtc "etcFile"
    \endverbatim

    The usual expansion of environment variables and other constructs is
    retained.

See also
    findEtcFile, fileName, string::expand()

SourceFiles
    includeEtcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeEtcEntry_H
#define includeEtcEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeEtcEntry Declaration
\*---------------------------------------------------------------------------*/

class includeEtcEntry
:
    public functionEntry
{
protected:

    // Protected Member Functions

        //- Expand include fileName and search etc directories for the file
        static fileName resolveEtcFile
        (
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
