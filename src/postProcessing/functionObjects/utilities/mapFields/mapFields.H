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
    Foam::mapFields

Group
    grpUtilitiesFunctionObjects

Description
    Map fields from local mesh to secondary mesh at run-time.

    Example of function object specification to map fields:
    \verbatim
    mapFields1
    {
        type            mapFields;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...

        mapRegion       coarseMesh;
        mapMethod       cellVolumeWeight;
        consistent      yes;

        fields          ("U.*" p);
    }

    \heading Function object usage
    \table
        Property     | Description                | Required  | Default value
        type         | Type name: mapFields       | yes       |
        mapRgion     | Name of region to map to   | yes       |
        mapMethod    | Mapping method             | yes       |
        patchMapMethod | Patch mapping method     | no        | <auto>
        consistent   | Mapping meshes have consistent boundaries | yes |
        fields       | List of field names to map | yes       |
        log          | Log to standard output     | no        | yes
    \endtable


SourceFiles
    mapFields.C
    IOmapFields.H

\*---------------------------------------------------------------------------*/

#ifndef mapFieldsFO_H
#define mapFieldsFO_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class meshToMesh;

/*---------------------------------------------------------------------------*\
                          Class mapFields Declaration
\*---------------------------------------------------------------------------*/

class mapFieldsFO
{
    // Private data

        //- Name of this set of mapFields objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Locally cached map region mesh (map to this mesh)
        autoPtr<fvMesh> mapRegionPtr_;

        //- Mesh-to-mesh interpolation
        autoPtr<meshToMesh> interpPtr_;

        //- List of field names to interpolate
        wordReList fieldNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mapFieldsFO(const mapFieldsFO&);

        //- Disallow default bitwise assignment
        void operator=(const mapFieldsFO&) = delete;

        //- Helper function to create the mesh-to-mesh interpolation
        void createInterpolation(const dictionary& dict);

        //- Helper function to interpolate and write the fied
        template<class Type>
        bool writeFieldType() const;


public:

    //- Runtime type information
    TypeName("mapFieldsFO");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        mapFieldsFO
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~mapFieldsFO();


    // Member Functions

        //- Return name of the set of mapFields
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the mapFields data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the mapFields and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
