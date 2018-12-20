/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::sampledSurfaces

Description
    Set of surfaces to sample.

    The write() method is used to sample and write files.

    Example of function object specification:

    \verbatim
    surfaces
    {
        type    surfaces;
        libs    ("libsampling.so");

        // Write at same frequency as fields
        writeControl    outputTime;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Scheme to obtain face centre value
        sampleScheme    cell;

        // Scheme to obtain node values
        // (only used if interpolate=true for the surfaces below)
        interpolationScheme cell;

        // Output surface format
        surfaceFormat   vtk;
        formatOptions
        {
            vtk
            {
                precision  10;
            }
        }

        surfaces
        (
            f0surf
            {
                type        sampledTriSurfaceMesh;
                surface     f0surf.obj;
                source      cells;

                // Optional: keep original regions
                keepIds     true;

                // Optional: generate values on points instead of faces
                interpolate true;
            }
        );
    }
    \endverbatim

    Entries:
    \table
        Property | Description                             | Required | Default
        type     | surfaces                                | yes      |
        surfaces | the list of sample surfaces             | recommended |
        fields   | word/regex list of fields to sampled    | yes      |
        sampleScheme | scheme to obtain face centre value  | no       | cell
        interpolationScheme | scheme to obtain node values | yes      |
        surfaceFormat | output surface format              | yes      |
        formatOptions | dictionary of format options       | no       |
    \endtable

Note
    The interpolationScheme is only used if interpolate=true is used by any
    of the surfaces.

See also
    Foam::surfMeshSamplers

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaces_H
#define sampledSurfaces_H

#include "regionFunctionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "mergedSurf.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "wordRes.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;
class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public functionObjects::regionFunctionObject,
    public PtrList<sampledSurface>
{
    // Static data members

        //- Output verbosity
        static bool verbose_;

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;


    // Private data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output path
        fileName outputPath_;


    // Read from dictonary

        //- Names of fields to sample
        wordRes fieldSelection_;

        //- Interpolation/sample scheme to obtain face values
        word sampleFaceScheme_;

        //- Interpolation/sample scheme to obtain node values
        word sampleNodeScheme_;


    // Surfaces

        //- Merged meshed surfaces (parallel only)
        List<mergedSurf> mergedList_;

        //- Track which surfaces have changed
        List<bool> changedGeom_;


    // Calculated

        //- Surface formatter
        autoPtr<surfaceWriter> formatter_;


    // Private Member Functions

        //- Return the surfaces
        const PtrList<sampledSurface>& surfaces() const
        {
            return *this;
        }

        //- Return the surfaces
        PtrList<sampledSurface>& surfaces()
        {
            return *this;
        }

        //- Return number of fields
        label classifyFields();

        //- Write geometry only
        void writeGeometry() const;

        //- Write scalar field with original ids
        void writeOriginalIds();

        //- Write sampled fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            const Field<Type>& values,
            const label surfi,
            const word& fieldName,
            const fileName& outputDir
        );

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
        );

        //- Sample and write all sampled fields
        template<class Type> void sampleAndWrite(const IOobjectList& objects);


        //- No copy construct
        sampledSurfaces(const sampledSurfaces&) = delete;

        //- No copy assignment
        void operator=(const sampledSurfaces&) = delete;


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~sampledSurfaces() = default;


    // Member Functions

        //- Does any of the surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        //- Update the surfaces as required and merge surface points (parallel).
        //  Return false if no surfaces required an update.
        virtual bool update();

        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set tolerance (and return old tolerance)
        static scalar mergeTol(const scalar tol);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
