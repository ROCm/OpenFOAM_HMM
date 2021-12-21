/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::probes

Group
    grpUtilitiesFunctionObjects

Description
    Set of locations to sample.

    Call write() to sample and write files.

    Example of function object specification:
    \verbatim
    probes
    {
        type            probes;
        libs            (sampling);

        // Name of the directory for probe data
        name            probes;

        // Write at same frequency as fields
        writeControl    outputTime;
        writeInterval   1;

        // Fields to be probed
        fields          (p U);

        // Optional: do not recalculate cells if mesh moves
        fixedLocations  false;

        // Optional: interpolation scheme to use (default is cell)
        interpolationScheme cellPoint;

        probeLocations
        (
            ( 1e-06 0 0.01 )      // at inlet
            (0.21 -0.20999 0.01)  // at outlet1
            (0.21 0.20999 0.01)   // at outlet2
            (0.21 0 0.01)         // at central block
        );

        // Optional: filter out points that haven't been found. Default
        //           is to include them (with value -VGREAT)
        includeOutOfBounds  true;
    }
    \endverbatim

SourceFiles
    probes.C

\*---------------------------------------------------------------------------*/

#ifndef probes_H
#define probes_H

#include "stateFunctionObject.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "pointField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "wordRes.H"

using namespace Foam::functionObjects;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class probes Declaration
\*---------------------------------------------------------------------------*/

class probes
:
    public stateFunctionObject,
    public pointField
{
protected:

    // Protected classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}
        };


    // Protected data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;


      // Read from dictionary

        //- Names of fields to probe
        wordRes fieldSelection_;

        //- Fixed locations, default = yes
        //  Note: set to false for moving mesh calculations where locations
        //        should move with the mesh
        bool fixedLocations_;

        //- Interpolation scheme name
        //  Note: only possible when fixedLocations_ is true
        word interpolationScheme_;

        //- Include probes that were not found
        bool includeOutOfBounds_;


      // Calculated

        //- Categorized scalar/vector/tensor vol fields
        fieldGroup<scalar> scalarFields_;
        fieldGroup<vector> vectorFields_;
        fieldGroup<sphericalTensor> sphericalTensorFields_;
        fieldGroup<symmTensor> symmTensorFields_;
        fieldGroup<tensor> tensorFields_;

        //- Categorized scalar/vector/tensor surf fields
        fieldGroup<scalar> surfaceScalarFields_;
        fieldGroup<vector> surfaceVectorFields_;
        fieldGroup<sphericalTensor> surfaceSphericalTensorFields_;
        fieldGroup<symmTensor> surfaceSymmTensorFields_;
        fieldGroup<tensor> surfaceTensorFields_;

        //- Cells to be probed (obtained from the locations)
        labelList elementList_;

        //- Faces to be probed
        labelList faceList_;

        //- Processor holding the cell or face (-1 if point not found
        //  on any processor)
        labelList processor_;

        //- Current open files
        HashPtrTable<OFstream> probeFilePtrs_;

        // Additional fields for patchProbes

            //- Patch IDs on which the new probes are located
            labelList patchIDList_;

            //- Original probes location (only used for patchProbes)
            pointField oldPoints_;


    // Protected Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Find cells and faces containing probes
        virtual void findElements(const fvMesh& mesh);

        //- Classify field type and open/close file streams,
        //  returns number of fields to sample
        label prepare();


private:

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

        //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- No copy construct
        probes(const probes&) = delete;

        //- No copy assignment
        void operator=(const probes&) = delete;


public:

    //- Runtime type information
    TypeName("probes");


    // Constructors

        //- Construct from Time and dictionary
        probes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~probes() = default;


    // Member Functions

        //- Return names of fields to probe
        virtual const wordRes& fieldNames() const
        {
            return fieldSelection_;
        }

        //- Return locations to probe
        virtual const pointField& probeLocations() const
        {
            return *this;
        }

        //- Return location for probe i
        virtual const point& probe(const label i) const
        {
            return operator[](i);
        }

        //- Cells to be probed (obtained from the locations)
        const labelList& elements() const
        {
            return elementList_;
        }

        //- Read the probes
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Sample a single vol field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a single scalar field on all sample locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceFields(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "probesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
