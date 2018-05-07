/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::surfMeshSample

Group
    grpUtilitiesFunctionObjects

Description
    An abstract class for surfMesh with sampling.

SourceFiles
    surfMeshSample.C
    surfMeshSampleTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSample_H
#define surfMeshSample_H

#include "surfMesh.H"
#include "surfFields.H"

#include "pointField.H"
#include "word.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "polyMesh.H"
#include "interpolation.H"
#include "error.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class surfMeshSample;

/*---------------------------------------------------------------------------*\
                       Class surfMeshSample Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSample
{
    // Private data

        //- Name for surfMesh lookup
        word name_;

        //- Mesh reference
        const polyMesh& mesh_;


protected:

    // Protected Member Functions

        //- General loop for sampling elements to faces
        template<class Type>
        static tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler,
            const labelUList& elements,
            const faceList& fcs,
            const pointField& pts
        );


        //- Get existing or create new surfMesh
        surfMesh& getOrCreateSurfMesh() const;

        //- Get existing or create new surfField.
        //  Create with same name and dimensions as the 'parent' volField.
        template<class Type>
        DimensionedField<Type, surfGeoMesh>&
        getOrCreateSurfField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;


        // //- Transfer field from object registry to surfField
        // template<class Type>
        // bool transferField
        // (
        //     const objectRegistry& store,
        //     const word& fieldName
        // ) const;


        //- Write the given fields
        template<class Type>
        label writeFields(const wordRes& select) const;

public:

    //- Runtime type information
    TypeName("surfMeshSample");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfMeshSample,
        word,
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        ),
        (name, mesh, dict)
    );


   // iNew helper class

        //- Class used for the PtrLists read-construction
        class iNew
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<surfMeshSample> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return surfMeshSample::New(name, mesh_, dict);
            }
        };


    // Constructors

        //- Construct from name, mesh
        surfMeshSample
        (
            const word& name,
            const polyMesh& mesh
        );

        //- Construct from dictionary
        surfMeshSample
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


        //- Clone
        autoPtr<surfMeshSample> clone() const
        {
            NotImplemented;
            return nullptr;
        }



    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<surfMeshSample> New
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshSample() = default;


    // Member Functions

    // Access

        //- Access to the underlying mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Name of surface
        const word& name() const
        {
            return name_;
        }


        //- Create surfMesh if required.
        void create() const;

        //- Return existing surfMesh.
        virtual const surfMesh& surface() const;


        //- Does the surface need an update?
        virtual bool needsUpdate() const = 0;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire() = 0;

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update() = 0;


    // Edit

        //- Rename
        virtual void rename(const word& newName)
        {
            name_ = newName;
        }


    // Execute

        //- Sample from specified volume field to obtain surface field.
        virtual bool sample
        (
            const word& fieldName,
            const word& sampleScheme = "cell"
        ) const = 0;

        //- Sample from specified volume fields to obtain surface fields.
        virtual label sample
        (
            const UList<word>& fieldNames,
            const word& sampleScheme = "cell"
        ) const;


    // Write

        //- Write specified fields
        virtual label write(const wordRes& fieldSelection) const;

        //- Write
        virtual void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "surfMeshSampleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
