/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::continuityError

Group
    grpFieldFunctionObjects

Description
    Calculates the continuity error for a flux field

Usage
    Example of function object specification:
    \verbatim
    continuityError1
    {
        type        continuityError;
        libs        ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        phi         phi;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: continuityError  | yes         |
        writeToFile  | write min/max data to file | no       | yes
        log          | write min/max data to standard output | no | yes
        phi          | name of flux field      | no | phi
    \endtable

    Output data is written to the file \<timeDir\>/continuityError.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    continuityError.C
    continuityErrorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_continuityError_H
#define functionObjects_continuityError_H

#include "Switch.H"
#include "fvMeshFunctionObject.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class continuityError Declaration
\*---------------------------------------------------------------------------*/

class continuityError
:
    public fvMeshFunctionObject,
    public writeFile
{
private:

    // Private Member Functions

        //- No copy construct
        continuityError(const continuityError&) = delete;

        //- No copy assignment
        void operator=(const continuityError&) = delete;


protected:

    // Protected data

        //- Name of the flux field; default = "phi
        word phiName_;

        //- Cumulative error
        scalar cumulative_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("continuityError");


    // Constructors

        //- Construct from Time and dictionary
        continuityError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~continuityError() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the continuityError
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
