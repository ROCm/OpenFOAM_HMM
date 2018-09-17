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
    Foam::functionObjects::vtkCloud

Group
    grpLagrangianFunctionObjects

Description
    This functionObject writes cloud(s) in VTK PolyData format
    (.vtp extension) with the time information saved in a '.series' file.

    Example of function object specification:
    \verbatim
    cloudWrite1
    {
        type            vtkCloud;
        libs            ("liblagrangianFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        format          ascii;

        cloud           myCloud;
        fields          (T U rho);
        width           4;  // file-padding
    }
    \endverbatim

Usage
    \table
        Property     | Description                      | Required    | Default
        type         | Type name: vtkCloud              | yes         |
        writeControl | Output control                   | recommended | timeStep
        cloud        |                                  | no  | defaultCloud
        clouds       | wordRe list of clouds            | no    |
        fields       | wordRe list of fields            | no    |
        prune        | Suppress writing of empty clouds | no    | false
        cellData     | Emit cellData instead of pointData | no  | false
        directory    | The output directory name     | no | postProcessing/NAME
        width        | Padding width for file name      | no    | 8
        format       | Format as ascii or binary        | no    | binary
        precision    | Write precision in ascii         | no | same as IOstream
    \endtable

    The output filename and fields are added to the functionObjectProperties
    information. For the previous example specification:

    \verbatim
    cloudWrite1
    {
        myCloud
        {
            file    "<case>/VTK/myCloud_0001.vtp";
            fields  (T U rho);
        }
    }
    \endverbatim

See also
    Foam::functionObjects::ensightWrite
    Foam::functionObjects::vtkWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    vtkCloud.C
    vtkCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vtkCloud_H
#define functionObjects_vtkCloud_H

#include "fvMeshFunctionObject.H"
#include "foamVtkOutputOptions.H"
#include "foamVtkSeriesWriter.H"
#include "wordRes.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vtkCloud Declaration
\*---------------------------------------------------------------------------*/

class vtkCloud
:
    public fvMeshFunctionObject
{
    // Private data

        //- Writer options
        vtk::outputOptions writeOpts_;

        //- The printf format for zero-padding names
        string printf_;

        //- Write lagrangian as cell data (verts) instead of point data
        bool useVerts_;

        //- Suppress writing of empty clouds
        bool pruneEmpty_;

        //- Requested names of clouds to process
        wordRes selectClouds_;

        //- Subset of cloud fields to process
        wordRes selectFields_;

        //- Output directory
        fileName directory_;

        //- Per cloud output for file series
        HashTable<vtk::seriesWriter, fileName> series_;


    // Private Member Functions

        //- Write a cloud to disk (creates parent directory),
        //- and record on the cloud OutputProperties.
        //  \param file is the output file name, with extension.
        bool writeCloud(const fileName& file, const word& cloudName);

        //- Write vertex (cells) - callable on master only
        void writeVerts
        (
            autoPtr<vtk::formatter>& format,
            const label nTotParcels
        ) const;

        //- Write fields of IOField<Type>
        template<class Type>
        wordList writeFields
        (
            autoPtr<vtk::formatter>& format,
            const objectRegistry& obrTmp,
            const label nTotParcels
        ) const;


        //- No copy construct
        vtkCloud(const vtkCloud&) = delete;

        //- No copy assignment
        void operator=(const vtkCloud&) = delete;


public:

    //- Runtime type information
    TypeName("vtkCloud");


    // Constructors

        //- Construct from Time and dictionary
        vtkCloud
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~vtkCloud() = default;


    // Member Functions

        //- Read the vtkCloud specification
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
    #include "vtkCloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
