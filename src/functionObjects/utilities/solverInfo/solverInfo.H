/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::functionObjects::solverInfo

Group
    grpUtilitiesFunctionObjects

Description
    Writes solver information for a list of user-specified fields

    Information written to file includes:
    - residual fields
    - solver type
    - initial residual
    - final residual
    - number of solver iterations
    - convergecnce flag

Usage
    Example of function object specification:
    \verbatim
    solverInfo
    {
        type            solverInfo;
        libs            ("libutilityFunctionObjects.so");
        ...
        fields          (U p);
        writeResidualFields yes;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description               | Required    | Default value
        type         | Type name: solverInfo     | yes         |
        fields       | List of fields to process | yes         |
        writeResidualFields | Write the residual fields | no          | no
    \endtable

    Output data is written to the dir postProcessing/solverInfo/\<timeDir\>/

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    solverInfo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_solverInfo_H
#define functionObjects_solverInfo_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "solverFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class solverInfo Declaration
\*---------------------------------------------------------------------------*/

class solverInfo
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Fields to process
        solverFieldSelection fieldSet_;

        //- Flag to write the residual as a vol field
        bool writeResidualFields_;

        //- Names of (result) residual fields
        wordHashSet residualFieldNames_;

        //- Initialisation flag
        bool initialised_;


    // Protected Member Functions

        //- Output file header information
        void writeFileHeader(Ostream& os);

        //- Create and store a residual field on the mesh database
        void createResidualField(const word& fieldName);

        //- Output file header information per primitive type value
        template<class Type>
        void writeFileHeader(Ostream& os, const word& fileName) const;

        //- Initialise a residual field
        template<class Type>
        void initialiseResidualField(const word& fieldName);

        //- Calculate the solver information
        template<class Type>
        void updateSolverInfo(const word& fieldName);


private:

    // Private member functions

        //- No copy construct
        solverInfo(const solverInfo&) = delete;

        //- No copy assignment
        void operator=(const solverInfo&) = delete;


public:

    //- Runtime type information
    TypeName("solverInfo");


    // Constructors

        //- Construct from Time and dictionary
        solverInfo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solverInfo() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the solverInfo
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solverInfoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
