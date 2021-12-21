/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fvMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

    NOTE ON USAGE:
    fvMesh contains all the topological and geometric information
    related to the mesh.  It is also responsible for keeping the data
    up-to-date.  This is done by deleting the cell volume, face area,
    cell/face centre, addressing and other derived information as
    required and recalculating it as necessary.  The fvMesh therefore
    reserves the right to delete the derived information upon every
    topological (mesh refinement/morphing) or geometric change (mesh
    motion).  It is therefore unsafe to keep local references to the
    derived data outside of the time loop.

SourceFiles
    fvMesh.C
    fvMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fvMesh_H
#define fvMesh_H

#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "fvBoundaryMesh.H"
#include "surfaceInterpolation.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "data.H"
#include "DimensionedField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "slicedVolFieldsFwd.H"
#include "slicedSurfaceFieldsFwd.H"
#include "className.H"
#include "SolverPerformance.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshLduAddressing;
class volMesh;
template<class Type>
class fvMatrix;

/*---------------------------------------------------------------------------*\
                           Class fvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMesh
:
    public polyMesh,
    public lduMesh,
    public fvSchemes,
    public surfaceInterpolation,    // needs input from fvSchemes
    public fvSolution,
    public data
{
protected:

    // Private data

        //- Boundary mesh
        fvBoundaryMesh boundary_;

    // Demand-driven data

        mutable fvMeshLduAddressing* lduPtr_;

        //- Current time index for cell volumes
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Cell volumes old time level
        mutable void* VPtr_;

        //- Cell volumes old time level
        mutable DimensionedField<scalar, volMesh>* V0Ptr_;

        //- Cell volumes old-old time level
        mutable DimensionedField<scalar, volMesh>* V00Ptr_;

        //- Face area vectors
        mutable slicedSurfaceVectorField* SfPtr_;

        //- Mag face area vectors
        mutable surfaceScalarField* magSfPtr_;

        //- Cell centres
        mutable slicedVolVectorField* CPtr_;

        //- Face centres
        mutable slicedSurfaceVectorField* CfPtr_;

        //- Face motion fluxes
        mutable surfaceScalarField* phiPtr_;


    // Private Member Functions

        // Storage management

            //- Clear geometry but not the old-time cell volumes
            void clearGeomNotOldVol();

            //- Clear geometry like clearGeomNotOldVol but recreate any
            //  geometric demand-driven data that was set
            void updateGeomNotOldVol();

            //- Clear local geometry
            void clearGeom();

            //- Clear local addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Clear local-only storage (geometry, addressing etc)
            void clearOutLocal();

            //- Preserve old volume(s)
            void storeOldVol(const scalarField&);


       // Make geometric data

            void makeSf() const;
            void makeMagSf() const;

            void makeC() const;
            void makeCf() const;


        //- No copy construct
        fvMesh(const fvMesh&) = delete;

        //- No copy assignment
        void operator=(const fvMesh&) = delete;


public:

    // Public typedefs

        typedef fvMesh Mesh;
        typedef fvBoundaryMesh BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("fvMesh");


    // Constructors

        //- Construct from IOobject
        explicit fvMesh(const IOobject& io, const bool doInit=true);

        //- Construct from IOobject or as zero-sized mesh
        //  Boundary is added using addFvPatches() member function
        fvMesh(const IOobject& io, const zero, bool syncPar=true);

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Construct as copy (for dictionaries) and components without
        //  boundary. Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            const fvMesh& baseMesh,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct as copy (for dictionaries) without boundary from cells
        //  rather than owner/neighbour. Boundary is added using addFvPatches()
        //  member function
        fvMesh
        (
            const IOobject& io,
            const fvMesh& baseMesh,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );


    //- Destructor
    virtual ~fvMesh();


    // Member Functions

        // Helpers

            //- Initialise all non-demand-driven data
            virtual bool init(const bool doInit);

            //- Add boundary patches. Constructor helper
            void addFvPatches
            (
                PtrList<polyPatch>& plist,
                const bool validBoundary = true
            );

            //- Add boundary patches. Constructor helper
            void addFvPatches
            (
                const List<polyPatch*>& p,
                const bool validBoundary = true
            );

            //- Update the mesh based on the mesh files saved in time
            //  directories
            virtual readUpdateState readUpdate();


        // Access

            //- Return the top-level database
            const Time& time() const
            {
                return polyMesh::time();
            }

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const
            {
                return true;
            }

            //- Return the object registry - resolve conflict polyMesh/lduMesh
            virtual const objectRegistry& thisDb() const
            {
                return polyMesh::thisDb();
            }

            //- Return reference to name
            //  Note: name() is currently ambiguous due to derivation from
            //  surfaceInterpolation
            const word& name() const
            {
                return polyMesh::name();
            }

            //- Return reference to boundary mesh
            const fvBoundaryMesh& boundary() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const;

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return polyMesh::comm();
            }


            // Overlap

                //- Interpolate interpolationCells only
                virtual void interpolate(volScalarField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volVectorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volSphericalTensorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volSymmTensorField&) const
                {}

                //- Interpolate interpolationCells only
                virtual void interpolate(volTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(scalarField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(vectorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(sphericalTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(symmTensorField&) const
                {}

                //- Interpolate interpolationCells only. No bcs.
                virtual void interpolate(tensorField&) const
                {}

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<scalar> solve
                (
                    fvMatrix<scalar>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<vector> solve
                (
                    fvMatrix<vector>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<sphericalTensor> solve
                (
                    fvMatrix<sphericalTensor>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<symmTensor> solve
                (
                    fvMatrix<symmTensor>&,
                    const dictionary&
                ) const;

                //- Solve returning the solution statistics given convergence
                //- tolerance. Use the given solver controls
                virtual SolverPerformance<tensor> solve
                (
                    fvMatrix<tensor>&,
                    const dictionary&
                ) const;


            //- Internal face owner. Note bypassing virtual mechanism so
            //  e.g. relaxation always gets done using original addressing
            const labelUList& owner() const
            {
                return fvMesh::lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return fvMesh::lduAddr().upperAddr();
            }

            //- Return cell volumes
            const DimensionedField<scalar, volMesh>& V() const;

            //- Return old-time cell volumes
            const DimensionedField<scalar, volMesh>& V0() const;

            //- Return old-old-time cell volumes
            const DimensionedField<scalar, volMesh>& V00() const;

            //- Return sub-cycle cell volumes
            tmp<DimensionedField<scalar, volMesh>> Vsc() const;

            //- Return sub-cycle old-time cell volumes
            tmp<DimensionedField<scalar, volMesh>> Vsc0() const;

            //- Return cell face area vectors
            const surfaceVectorField& Sf() const;

            //- Return cell face area magnitudes
            const surfaceScalarField& magSf() const;

            //- Return cell face motion fluxes
            const surfaceScalarField& phi() const;

            //- Return cell centres as volVectorField
            const volVectorField& C() const;

            //- Return face centres as surfaceVectorField
            const surfaceVectorField& Cf() const;

            //- Return face deltas as surfaceVectorField
            tmp<surfaceVectorField> delta() const;

            //- Return a labelType of valid component indicators
            //  1 : valid (solved)
            // -1 : invalid (not solved)
            template<class Type>
            typename pTraits<Type>::labelType validComponents() const;


        // Edit

            //- Clear all geometry and addressing
            void clearOut();

            //- Update mesh corresponding to the given map
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Avoid masking surfaceInterpolation method
            using surfaceInterpolation::movePoints;

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const pointField&);

            //- Update all geometric data. This gets redirected up from
            //- primitiveMesh level
            virtual void updateGeom();

            //- Map all fields in time using given map.
            virtual void mapFields(const mapPolyMesh& mpm);

            //- Remove boundary patches. Warning: fvPatchFields hold ref to
            //- these fvPatches.
            void removeFvBoundary();

            //- Return cell face motion fluxes
            surfaceScalarField& setPhi();

            //- Return old-time cell volumes
            DimensionedField<scalar, volMesh>& setV0();


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool valid
            ) const;

            //- Write mesh using IO settings from time
            virtual bool write(const bool valid = true) const;


    // Member Operators

        //- Compares addresses
        bool operator!=(const fvMesh& rhs) const;

        //- Compares addresses
        bool operator==(const fvMesh& rhs) const;
};


template<>
typename pTraits<sphericalTensor>::labelType
fvMesh::validComponents<sphericalTensor>() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshTemplates.C"
    #include "fvPatchFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
