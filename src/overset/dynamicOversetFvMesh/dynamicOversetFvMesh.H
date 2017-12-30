/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::dynamicOversetFvMesh

Description
    dynamicFvMesh with support for overset meshes.

SourceFiles
    dynamicOversetFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicOversetFvMesh_H
#define dynamicOversetFvMesh_H

#include "dynamicMotionSolverFvMesh.H"
#include "labelIOList.H"
#include "fvMeshPrimitiveLduAddressing.H"
#include "lduInterfaceFieldPtrsList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;
class lduPrimitiveProcessorInterface;

/*---------------------------------------------------------------------------*\
                   Class dynamicOversetFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicOversetFvMesh
:
    public dynamicMotionSolverFvMesh
{
protected:

    // Protected data

        //- Select base addressing (false) or locally stored extended
        //  lduAddressing (true)
        mutable bool active_;

        //- Extended addressing (extended with local interpolation stencils)
        mutable autoPtr<fvMeshPrimitiveLduAddressing> lduPtr_;

        //- Added (processor)lduInterfaces for remote bits of stencil.
        //PtrList<const lduInterface> remoteStencilInterfaces_;
        mutable PtrList<const lduPrimitiveProcessorInterface>
            remoteStencilInterfaces_;

        //- Interfaces for above mesh. Contains both original and
        //  above added processorLduInterfaces
        mutable lduInterfacePtrsList allInterfaces_;

        //- Corresponding faces (in above lduPtr) to the stencil
        mutable labelListList stencilFaces_;

        //- From old to new face labels
        mutable labelList reverseFaceMap_;


   // Protected Member Functions

        //- Calculate the extended lduAddressing
        virtual bool updateAddressing() const;

        //- Debug: print matrix
        template<class Type>
        void write(Ostream&, const fvMatrix<Type>&, const lduAddressing&) const;

        //- Explicit interpolation of acceptor cells from donor cells
        template<class T>
        void interpolate(Field<T>& psi) const;

        //- Explicit interpolation of acceptor cells from donor cells with
        //  boundary condition handling
        template<class GeoField>
        void interpolate(GeoField& psi) const;

        //- Freeze values at holes
        template<class Type>
        void freezeHoles(fvMatrix<Type>&) const;

        //- Get scalar interfaces of certain type
        template<class GeoField, class PatchType>
        lduInterfaceFieldPtrsList scalarInterfaces(const GeoField& psi) const;

        //- Add interpolation to matrix (coefficients)
        template<class Type>
        void addInterpolation(fvMatrix<Type>&) const;

        //- Solve given dictionary with settings
        template<class Type>
        SolverPerformance<Type> solve(fvMatrix<Type>&, const dictionary&) const;

        //- Debug: correct coupled bc
        template<class GeoField>
        static void correctCoupledBoundaryConditions(GeoField& fld);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicOversetFvMesh(const dynamicOversetFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicOversetFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicOversetFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicOversetFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicOversetFvMesh();


    // Member Functions


        // Extended addressing

            //- Return extended ldu addressing
            const fvMeshPrimitiveLduAddressing& primitiveLduAddr() const;

            //- Return ldu addressing. If active: is (extended)
            //  primitiveLduAddr
            virtual const lduAddressing& lduAddr() const;

            //- Return old to new face addressing
            const labelList& reverseFaceMap() const
            {
                return reverseFaceMap_;
            }

            //- Return true if using extended addressing
            bool active() const
            {
                return active_;
            }

            //- Enable/disable extended addressing
            void active(const bool f) const
            {
                active_ = f;

                if (active_)
                {
                    DebugInfo<< "Switching to extended addressing with nFaces:"
                        << primitiveLduAddr().lowerAddr().size()
                        << endl;
                }
                else
                {
                    DebugInfo<< "Switching to base addressing with nFaces:"
                        << fvMesh::lduAddr().lowerAddr().size()
                        << endl;
                }
            }


        // Overset

            // Explicit interpolation

                virtual void interpolate(scalarField& psi) const
                {
                    interpolate<scalar>(psi);
                }

                virtual void interpolate(vectorField& psi) const
                {
                    interpolate<vector>(psi);
                }

                virtual void interpolate(sphericalTensorField& psi) const
                {
                    interpolate<sphericalTensor>(psi);
                }

                virtual void interpolate(symmTensorField& psi) const
                {
                    interpolate<symmTensor>(psi);
                }

                virtual void interpolate(tensorField& psi) const
                {
                    interpolate<tensor>(psi);
                }

                virtual void interpolate(volScalarField& psi) const
                {
                    interpolate<volScalarField>(psi);
                }

                virtual void interpolate(volVectorField& psi) const
                {
                    interpolate<volVectorField>(psi);
                }

                virtual void interpolate(volSphericalTensorField& psi) const
                {
                    interpolate<volSphericalTensorField>(psi);
                }

                virtual void interpolate(volSymmTensorField& psi) const
                {
                    interpolate<volSymmTensorField>(psi);
                }

                virtual void interpolate(volTensorField& psi) const
                {
                    interpolate<volTensorField>(psi);
                }


            // Implicit interpolation (matrix manipulation)

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<scalar> solve
                (
                    fvMatrix<scalar>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<scalar>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<vector> solve
                (
                    fvMatrix<vector>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<vector>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<symmTensor> solve
                (
                    fvMatrix<symmTensor>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<symmTensor>(m, dict);
                }

                //- Solve returning the solution statistics given convergence
                //  tolerance. Use the given solver controls
                virtual SolverPerformance<tensor> solve
                (
                    fvMatrix<tensor>& m,
                    const dictionary& dict
                ) const
                {
                    return solve<tensor>(m, dict);
                }


        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool valid
        ) const;

        //- Debug: check halo swap is ok
        template<class GeoField>
        static void checkCoupledBC(const GeoField& fld);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dynamicOversetFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
