/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "AtomicAccumulator.H"



#ifdef USE_ROCTX
#include <roctracer/roctx.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::surfaceIntegrate_A");
    #endif

    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field<Type>& issf = ssf;

    #if 0
    forAll(owner, facei)
    #else
    const label loop_len = owner.size();
    #pragma omp target teams distribute parallel for if(target:loop_len>20000)
    for (label facei = 0; facei < loop_len; ++facei)
    #endif
    {
            
        atomicAccumulator(ivf[owner[facei]]) += issf[facei];
        atomicAccumulator(ivf[neighbour[facei]]) -= issf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        #if 0
        forAll(mesh.boundary()[patchi], facei)
        #else
        const label loop_len = mesh.boundary()[patchi].size();
        #pragma omp target teams distribute parallel for if(target:loop_len>20000)
        for (label facei = 0; facei < loop_len; ++facei)
        #endif
        {
            atomicAccumulator(ivf[pFaceCells[facei]]) += pssf[facei];
        }
    }

    ivf /= mesh.Vsc();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::surfaceIntegrate_B");
    #endif

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions()/dimVol, Zero),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf);
    vf.correctBoundaryConditions();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{

    #ifdef USE_ROCTX
    roctxRangePush("fvc::surfaceIntegrate_call");
    #endif

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        fvc::surfaceIntegrate(tssf())
    );
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    #ifdef USE_ROCTX
    roctxRangePush("fvc::surfaceIntegrate_clear");
    #endif

    tssf.clear();
    
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
    
    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceSum
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::surfaceSum");
    #endif

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceSum("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions(), Zero),
            extrapolatedCalculatedFvPatchField<Type>::typeName
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    #if 0
    forAll(owner, facei)
    #else
    const label loop_len = owner.size();
    #pragma omp target teams distribute parallel for if(target:loop_len>20000)
    for (label facei = 0; facei < loop_len; ++facei)
    #endif
    {
        atomicAccumulator(vf[owner[facei]]) += ssf[facei];
        atomicAccumulator(vf[neighbour[facei]]) += ssf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        #if 0
        forAll(mesh.boundary()[patchi], facei)
        #else
        const label loop_len = mesh.boundary()[patchi].size();
        #pragma omp target teams distribute parallel for if(target:loop_len>20000)
	for (label facei = 0; facei < loop_len; ++facei)
        #endif
        {   
            atomicAccumulator(vf[pFaceCells[facei]]) += pssf[facei];
        }
    }

    vf.correctBoundaryConditions();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tvf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfaceSum
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf = surfaceSum(tssf());
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
