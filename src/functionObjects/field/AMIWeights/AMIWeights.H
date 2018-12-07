/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::functionObjects::AMIWeights

Group
    grpFieldFunctionObjects

Description

Usage
    Example of function object specification:
    \verbatim
    AMIWeights1
    {
        type        AMIWeights;
        libs        ("libfieldFunctionObjects.so");
        writeFields yes;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: AMIWeights   | yes         |
        writeFields  | write weights as VTK fields | yes     |
    \endtable

    Output data is written to the file \<timeDir\>/AMIWeights.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    AMIWeights.C
    AMIWeightsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_AMIWeights_H
#define functionObjects_AMIWeights_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class AMIWeights Declaration
\*---------------------------------------------------------------------------*/

class AMIWeights
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Flag to write AMI fields (as VTK files)
        bool writeFields_;

        //- List of AMI patch IDs
        labelList patchIDs_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Helper function to report patch information
        virtual void reportPatch(const cyclicAMIPolyPatch& pp);

        void writeWeightField
        (
            const cyclicAMIPolyPatch& cpp,
            const scalarField& weightSum,
            const word& side
        ) const;

        void writeWeightFields(const cyclicAMIPolyPatch& cpp) const;

        //- No copy construct
        AMIWeights(const AMIWeights&) = delete;

        //- No copy assignment
        void operator=(const AMIWeights&) = delete;


public:

    //- Runtime type information
    TypeName("AMIWeights");


    // Constructors

        //- Construct from Time and dictionary
        AMIWeights
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~AMIWeights() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the AMIWeights
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
