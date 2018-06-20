/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::functionObjects::residuals

Group
    grpUtilitiesFunctionObjects

Description
    Writes out the initial residual for specified fields.

Usage
    Example of function object specification:
    \verbatim
    residuals
    {
        type            residuals;
        libs            ("libutilityFunctionObjects.so");
        ...
        fields          (U p);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description               | Required    | Default value
        type         | Type name: residua  ls    | yes         |
        fields       | List of fields to process | yes         |
        writeFields  | Write the residual fields | no          | no
    \endtable

    Output data is written to the dir postProcessing/residuals/\<timeDir\>/
    For vector/tensor fields, e.g. U, where an equation is solved for each
    component, the largest residual of each component is written.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    residuals.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_residuals_H
#define functionObjects_residuals_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "solverFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class residuals Declaration
\*---------------------------------------------------------------------------*/

class residuals
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Fields to write residuals
        solverFieldSelection fieldSet_;

        //- Flag to write the residual as a vol field
        bool writeFields_;

        //- Initialisation flag
        bool initialised_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(Ostream& os);

        //- Create and store a residual field on the mesh database
        void createField(const word& fieldName);

        //- Write a residual field
        void writeField(const word& fieldName) const;

        //- Output file header information per primitive type value
        template<class Type>
        void writeFileHeader(Ostream& os, const word& fileName) const;

        //- Initialise a residual field
        template<class Type>
        void initialiseField(const word& fieldName);

        //- Calculate the residual
        template<class Type>
        void writeResidual(const word& fieldName);


private:

    // Private member functions

        //- No copy construct
        residuals(const residuals&) = delete;

        //- No copy assignment
        void operator=(const residuals&) = delete;


public:

    //- Runtime type information
    TypeName("residuals");


    // Constructors

        //- Construct from Time and dictionary
        residuals
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~residuals();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the residuals
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "residualsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
