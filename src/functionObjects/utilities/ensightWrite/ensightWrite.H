/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::ensightWrite

Group
    grpUtilitiesFunctionObjects

Description
    Writes fields in ensight format.

    Example of function object specification:
    \verbatim
    ensight
    {
        type            ensightWrite;
        libs            ("libutilityFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        format          binary;

        overwrite       true;
        width           12;
        directory       "EnSight";

        fields
        (
            U
            p
        );
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description                  | Required  | Default value
        type         | Type name: ensightWrite      | yes       |
        fields       | Fields to output             | yes       |
        writeControl | Output control               | recommended | timeStep
        directory    | The output directory name    | no        | "ensightWrite"
        overwrite    | Remove existing directory    | no        | false
        format       | ASCII or binary format       | no        | same as simulation
        width        | Mask width for \c data/XXXX  | no        | 8
        noPatches    | Suppress writing patches     | no        | false
        patches      | Select patches to write      | no        |
        faceZones    | Select faceZones to write    | no        |
        consecutive  | Consecutive output numbering | no        | false
    \endtable

    Note that if the \c patches entry is an empty list, this will select all
    patches and suppress writing the internalMesh.
    Consecutive output numbering can be used in conjunction with \c overwrite.

See also
    Foam::functionObjects::vtkWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    ensightWrite.C
    ensightWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ensightWrite_H
#define functionObjects_ensightWrite_H

#include "fvMeshFunctionObject.H"
#include "ensightCase.H"
#include "ensightMesh.H"

#include "wordReList.H"
#include "interpolation.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ensightWrite Declaration
\*---------------------------------------------------------------------------*/

class ensightWrite
:
    public fvMeshFunctionObject
{
    // Private data

        //- Writer options
        ensightMesh::options writeOpts_;
        ensightCase::options caseOpts_;

        //- Name of fields to process
        wordReList selectFields_;

        //- Output directory name
        fileName dirName_;

        //- Consecutive output numbering
        bool consecutive_;

        //- Ensight case handler
        autoPtr<ensightCase> ensCase_;

        //- Ensight mesh handler
        autoPtr<ensightMesh> ensMesh_;


    // Private Member Functions

        //- Ensight case handler
        ensightCase& ensCase()
        {
            return ensCase_();
        }

        //- Ensight mesh handler
        ensightMesh& ensMesh()
        {
            return ensMesh_();
        }


        //- Apply for the volume field type
        template<class Type>
        int writeVolField(const word& inputName, int& state);


        //- Process by trying to apply for various volume field types.
        int process(const word& inputName);


        //- Disallow default bitwise copy construct
        ensightWrite(const ensightWrite&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightWrite&) = delete;

public:

    //- Runtime type information
    TypeName("ensightWrite");


    // Constructors

        //- Construct from runTime and dictionary.
        ensightWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ensightWrite();


    // Member Functions

        //- Read the ensightWrite specification
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();

        //- Execute at the final time-loop, flush case file
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
