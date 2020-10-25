/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
        libs    (sampling);

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

        // Optional: registry storage
        store           true

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
        {
            f0surf
            {
                type        meshedSurface;
                surface     f0surf.obj;
                source      cells;

                // Optional: keep original regions
                keepIds     true;

                // Optional: generate values on points instead of faces
                interpolate true;

                // Optional: alternative output type
                surfaceFormat   ensight;

                // Optional: registry storage
                store       true
            }
        }
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type     | Type-name: surfaces                      | yes |
        surfaces | Dictionary or list of sample surfaces    | expected |
        fields   | word/regex list of fields to sample      | yes |
        sampleScheme | scheme to obtain face centre value   | no  | cell
        interpolationScheme | scheme to obtain node values  | no  | cellPoint
        surfaceFormat | output surface format               | yes |
        formatOptions | dictionary of format options        | no  |
        sampleOnExecute | Sample (store) on execution as well | no | false
        store    | Store surface/fields on registry         | no  | false
    \endtable

    Additional per-surface entries:
    \table
        Property | Description                              | Required | Default
        store    | Store surface/fields on registry         | no  |
        surfaceFormat | output surface format               | no  |
        formatOptions | dictionary of format options        | no  |
        surfMeshStore | Store surface/fields as surfMesh (transitional) | no  |
    \endtable

Note
    The interpolationScheme is only used if interpolate=true is used by any
    of the surfaces.

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaces_H
#define sampledSurfaces_H

#include "fvMeshFunctionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polySurface;

/*---------------------------------------------------------------------------*\
                       Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public functionObjects::fvMeshFunctionObject,
    public PtrList<sampledSurface>
{
    // Static Data Members

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_SURF_MESH = 0x4,
            ACTION_ALL = 0xF
        };


    // Private Data

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output verbosity
        bool verbose_;

        //- Perform sample/store actions on execute as well
        bool onExecute_;

        //- Output path
        fileName outputPath_;


    // Read from dictionary

        //- Names of fields to sample
        wordRes fieldSelection_;

        //- Interpolation/sample scheme to obtain face values
        word sampleFaceScheme_;

        //- Interpolation/sample scheme to obtain node values
        word sampleNodeScheme_;


    // Output control

        //- Surface writers (one per surface)
        PtrList<surfaceWriter> writers_;

        //- Per-surface selection of store/write actions
        List<unsigned> actions_;

        //- Cached values of the global number of faces per-surface
        labelList nFaces_;


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

        //- A new surfaceWriter, with per-surface formatOptions
        static autoPtr<surfaceWriter> newWriter
        (
            word writeType,
            const dictionary& formatOptions,
            const dictionary& surfDict
        );


        //- Perform sampling action with store/write
        bool performAction(unsigned request);

        //- Count selected/sampled fields per surface
        void countFields();

        //- Write sampled fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            surfaceWriter& writer,
            const Field<Type>& values,
            const word& fieldName
        );

        //- Sample and store/write a specific volume field
        template<class Type>
        void performAction
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            unsigned request
        );

        //- Sample and store/write a specific surface field
        template<class Type>
        void performAction
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
            unsigned request
        );

        //- Sample and write all applicable sampled fields
        template<class GeoField>
        void performAction
        (
            const IOobjectList& objects,
            unsigned request
        );


        //- Get surface from registry if available.
        //  \return surface or nullptr
        polySurface* getRegistrySurface(const sampledSurface& s) const;

        //- Put surface onto registry, when enabled.
        //  \return surface or nullptr it surface should not be stored
        polySurface* storeRegistrySurface(const sampledSurface& s);

        //- Remove surface from registry.
        //  \return True if surface existed and was removed
        bool removeRegistrySurface(const sampledSurface& s);

        //- Store sampled field onto surface registry if it exists
        template<class Type, class GeoMeshType>
        bool storeRegistryField
        (
            const sampledSurface& s,
            const word& fieldName,
            const dimensionSet& dims,
            Field<Type>&& values
        );

        //- Test surfaces for condition
        template<class Container, class Predicate>
        static bool testAny(const Container& items, const Predicate& pred);

        //- Do any surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  Return false if all surfaces were already marked as expired.
        //  Optionally force expire, even if a surface has been marked as
        //  invariant.
        virtual bool expire(const bool force=false);

        //- Update the surfaces as required.
        //  Return false if no surfaces required an update.
        virtual bool update();

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

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on);

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary& dict);

        //- Sample and store if the sampleOnExecute is enabled.
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh& mesh);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set merge tolerance and return old value
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
