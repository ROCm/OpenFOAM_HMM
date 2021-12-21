/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::sampledSets

Description
    Set of sets to sample.
    Call sampledSets.write() to sample and write files.

SourceFiles
    sampledSets.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSets_H
#define sampledSets_H

#include "regionFunctionObject.H"
#include "sampledSet.H"
#include "volFieldsFwd.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "coordSet.H"
#include "writer.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class objectRegistry;
class dictionary;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class sampledSets Declaration
\*---------------------------------------------------------------------------*/

class sampledSets
:
    public functionObjects::regionFunctionObject,
    public PtrList<sampledSet>
{
    // Private Classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:

            //- The set formatter
            autoPtr<writer<Type>> formatter;

            //- Construct null
            fieldGroup() = default;

            //- Reset format and field list
            void clear()
            {
                DynamicList<word>::clear();
                formatter.clear();
            }

            void setFormatter(const word& writeFormat, const dictionary& dict)
            {
                formatter = writer<Type>::New(writeFormat, dict);
            }
        };


        //- Class used for sampling volFields
        template<class Type>
        class volFieldSampler
        :
            public List<Field<Type>>
        {
            //- Name of this collection of values
            const word name_;

        public:

            //- Construct interpolating field to the sampleSets
            volFieldSampler
            (
                const word& interpolationScheme,
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct mapping field to the sampleSets
            volFieldSampler
            (
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct from components
            volFieldSampler
            (
                const List<Field<Type>>& values,
                const word& name
            );

            //- Return the field name
            const word& name() const
            {
                return name_;
            }
        };


    // Static Data Members

        //- Output verbosity
        static bool verbose_;


    // Private Data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Keep the dictionary to recreate sets for moving mesh cases
        dictionary dict_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Output path
        fileName outputPath_;

        //- Mesh search engine
        meshSearch searchEngine_;


    // Read from dictionary

        //- Names of fields to sample
        wordRes fieldSelection_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- Output format to use
        word writeFormat_;

        //- Dictionary containing writer options
        dictionary writeFormatOptions_;


    // Categorized scalar/vector/tensor fields

        fieldGroup<scalar> scalarFields_;
        fieldGroup<vector> vectorFields_;
        fieldGroup<sphericalTensor> sphericalTensorFields_;
        fieldGroup<symmTensor> symmTensorFields_;
        fieldGroup<tensor> tensorFields_;


    // Merging structures

        PtrList<coordSet> masterSampledSets_;
        labelListList indexSets_;


    // Private Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Combine points from all processors. Sort by curveDist and produce
        //- index list. Valid result only on master processor.
        void combineSampledSets
        (
            PtrList<coordSet>& masterSampledSets,
            labelListList& indexSets
        );

        //- Combine values from all processors.
        //  Valid result only on master processor.
        template<class T>
        void combineSampledValues
        (
            const PtrList<volFieldSampler<T>>& sampledFields,
            const labelListList& indexSets,
            PtrList<volFieldSampler<T>>& masterFields
        );

        //- Write set on master, return fileName
        template<class Type>
        fileName writeSampleFile
        (
            const coordSet& masterSampleSet,
            const PtrList<volFieldSampler<Type>>& masterFields,
            const label setI,
            const fileName& timeDir,
            const writer<Type>& formatter
        );

        template<class Type>
        void sampleAndWrite(fieldGroup<Type>& fields);


        //- No copy construct
        sampledSets(const sampledSets&) = delete;

        //- No copy assignment
        void operator=(const sampledSets&) = delete;


public:

    //- Runtime type information
    TypeName("sets");


    // Constructors

        //- Construct from Time and dictionary
        sampledSets
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSets
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~sampledSets() = default;


    // Member Functions

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on);

        //- Read the sampledSets
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Correct for mesh changes
        void correct();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSetsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
