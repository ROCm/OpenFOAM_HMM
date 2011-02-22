/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::functionEntries::codeStream

Description
    Dictionary entry that contains C++ OpenFOAM code that is compiled to
    generate the entry itself. So
    - codeStream reads three entries: 'code', 'codeInclude' (optional),
    'codeOptions' (optional)
    and uses those to generate library sources inside constant/codeStream/
    - these get compiled using 'wmake libso'
    - the resulting library is loaded in executed with as arguments
            const dictionary& dict,
            Ostream& os
    where the dictionary is the current dictionary.
    - the code has to write into Ostream which is then used to construct
    the actual dictionary entry.


    E.g. to set the internal field of a field:

    internalField  #codeStream
    {
        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());
            scalarField fld(mesh.nCells(), 12.34);
            fld.writeEntry("", os);
        #};

        //- Optional:
        codeInclude
        #{
            #include "fvCFD.H"
        #};
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude
        #};
    };


    Note the #{ #} syntax which is just a way of inputting strings with embedded
    newlines.

    Limitations:
    - '~' symbol not allowed inside the code sections.
    - probably some other limitations (uses string::expand which expands $, ~)

SourceFiles
    codeStream.C

\*---------------------------------------------------------------------------*/

#ifndef codeStream_H
#define codeStream_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                         Class codeStream Declaration
\*---------------------------------------------------------------------------*/

class codeStream
:
    public functionEntry
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        codeStream(const codeStream&);

        //- Disallow default bitwise assignment
        void operator=(const codeStream&);


public:

    //- Runtime type information
    ClassName("codeStream");


    // Member Functions

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
