/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::surfMeshSamplers

Description
    Set of surfaces to sample into a surfMesh/surfField.

    The execute() method is used to sample, and the write() method to write.
    It is fairly common to use for sampling only and have the write disabled.

    \verbatim
    surfaces
    {
        type    surfMeshes;
        libs    ("libsampling.so");

        // Sample at every time-step
        executeControl  timeStep;
        executeInterval 1;

        // Disable writing (or write at same frequency as fields)
        writeControl    none;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Optional: pre-defined derived fields to be sampled
        derived         (rhoU pTotal);

        // Optional: create surface immediately on read
        // The default is to create a placeholder without any faces.
        createOnRead    false;

        surfaces
        (
            f0surf
            {
                type        sampledTriSurfaceMesh;
                surface     f0surf.obj;
                source      cells;
            }
        );
    }
    \endverbatim


See also
    Foam::sampledSurfaces

SourceFiles
    surfMeshSamplers.C

\*---------------------------------------------------------------------------*/

#ifndef surfMeshSamplers_H
#define surfMeshSamplers_H

#include "regionFunctionObject.H"
#include "surfMeshSampler.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "wordReList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class surfMeshSamplers Declaration
\*---------------------------------------------------------------------------*/

class surfMeshSamplers
:
    public functionObjects::regionFunctionObject,
    public PtrList<surfMeshSampler>
{
    // Static data members

        //- output verbosity
        static bool verbose_;


    // Private data

        //- Const reference to fvMesh
        const fvMesh& mesh_;


      // Read from dictonary

        //- Names of fields to sample
        wordReList fieldSelection_;

        //- Names of derived fields to create and sample
        wordList derivedNames_;


    // Private Member Functions

        //- Remove items by name from objectRegistry
        static void checkOutNames
        (
            const objectRegistry& registry,
            const UList<word>& names
        );

        //- Access the sampling surfaces
        inline const PtrList<surfMeshSampler>& surfaces() const
        {
            return static_cast<const PtrList<surfMeshSampler>&>(*this);
        }

        //- Access the sampling surfaces
        inline PtrList<surfMeshSampler>& surfaces()
        {
            return static_cast<PtrList<surfMeshSampler>&>(*this);
        }


        //- Filter acceptable fields types
        template<class Type>
        wordList acceptType() const;


        //- Disallow default bitwise copy construct and assignment
        surfMeshSamplers(const surfMeshSamplers&) = delete;
        void operator=(const surfMeshSamplers&) = delete;


public:

    //- Runtime type information
    TypeName("surfMeshes");


    // Constructors

        //- Construct from name, Time and dictionary
        surfMeshSamplers
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        surfMeshSamplers
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfMeshSamplers();


    // Member Functions

        //- Do any of the surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        //- Update the surfaces as required and merge surface points (parallel).
        //  Return false if no surfaces required an update.
        virtual bool update();

        //- set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the surfMeshSamplers dictionary
        virtual bool read(const dictionary&);

        //- Execute, does sampling
        virtual bool execute();

        //- Write sampled values
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfMeshSamplersTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
