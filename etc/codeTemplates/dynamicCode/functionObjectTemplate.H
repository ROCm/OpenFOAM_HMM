/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
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

Description
    Template for use with dynamic code generation of a
    OutputFilter functionObject.

SourceFiles
    functionObjectTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectTemplate_H
#define functionObjectTemplate_H

#include "regionFunctionObject.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class ${typeName}FunctionObject
:
    public functionObjects::regionFunctionObject
{
    // Private data

//{{{ begin codeData
        ${codeData}
//}}} end codeData


    // Private Member Functions

        const fvMesh& mesh() const;

        //- No copy construct
        ${typeName}FunctionObject(const ${typeName}FunctionObject&) = delete;

        //- No copy assignment
        void operator=(const ${typeName}FunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from Time and dictionary
        ${typeName}FunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~${typeName}FunctionObject();


    // Member Functions

        //- Read the system calls
        virtual bool read(const dictionary& dict);

        //- Execute the "executeCalls" at each time-step
        virtual bool execute();

        //- Execute the "endCalls" at the final time-loop
        virtual bool end();

        //- Write, execute the "writeCalls"
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
