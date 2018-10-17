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
    Foam::functionObjects::dataCloud

Group
    grpLagrangianFunctionObjects

Description
    This functionObject writes a cloud position and in ASCII.

    Example of function object specification:
    \verbatim
    cloudWrite1
    {
        type            dataCloud;
        libs            ("liblagrangianFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        cloud           myCloud;
        field           d;
    }
    \endverbatim

Usage
    \table
        Property     | Description                      | Required | Default
        type         | Type name: dataCloud             | yes   |
        cloud        |                                  | no    | defaultCloud
        clouds       | wordRe list of clouds            | no    |
        field        | Name of the field                | yes   |
        directory    | The output directory name     | no | postProcessing/NAME
        width        | Padding width for file name      | no    | 8
    \endtable

See also
    Foam::functionObjects::vtkCloud

SourceFiles
    dataCloud.C
    dataCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_dataCloud_H
#define functionObjects_dataCloud_H

#include "fvMeshFunctionObject.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class dataCloud Declaration
\*---------------------------------------------------------------------------*/

class dataCloud
:
    public fvMeshFunctionObject
{
    // Private data

        //- The printf format for zero-padding names
        string printf_;

        //- Requested names of clouds to process
        wordRes selectClouds_;

        //- Subset of cloud fields to process
        word fieldName_;

        //- Output directory
        fileName directory_;


    // Private Member Functions

        template<class Type>
        static void writeField
        (
            Ostream& os,
            const vectorField& points,
            const Field<Type>& field
        );

        //- Write to disk
        bool writeCloud(const fileName& outputName, const word& cloudName);

        //- Write from objectRegistry entry
        template<class Type>
        bool writeField
        (
            const fileName& outputName,
            const objectRegistry& obrTmp
        ) const;


        //- No copy construct
        dataCloud(const dataCloud&) = delete;

        //- No copy assignment
        void operator=(const dataCloud&) = delete;


public:

    //- Runtime type information
    TypeName("dataCloud");


    // Constructors

        //- Construct from Time and dictionary
        dataCloud
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dataCloud() = default;


    // Member Functions

        //- Read the dataCloud specification
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dataCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
