/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "cellLimitedGrad.H"
#include "gaussGrad.H"

#ifdef USE_ROCTX
#include <roctx.h>
#endif



  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<scalar>& limiter,
    Field<vector>& gIf
) const
{
    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_A");
    #endif

    gIf *= limiter;

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

}


template<class Type, class Limiter>
void Foam::fv::cellLimitedGrad<Type, Limiter>::limitGradient
(
    const Field<vector>& limiter,
    Field<tensor>& gIf
) const
{
    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_B");
    #endif


    forAll(gIf, celli)
    {
        gIf[celli] = tensor
        (
            cmptMultiply(limiter[celli], gIf[celli].x()),
            cmptMultiply(limiter[celli], gIf[celli].y()),
            cmptMultiply(limiter[celli], gIf[celli].z())
        );
    }
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
}


template<class Type, class Limiter>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::cellLimitedGrad<Type, Limiter>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_C");
    #endif

    const fvMesh& mesh = vsf.mesh();

    tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
    > tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (k_ < SMALL)
    {
        return tGrad;
    }

    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    Field<Type> maxVsf(vsf.primitiveField());
    Field<Type> minVsf(vsf.primitiveField());

    #if 1
    forAll(owner, facei)
    #else
    //LG1 AMD: make sure we do not have any race conditions  
    #pragma omp target teams distribute parallel for
    for (label facei = 0; facei < owner.size(); ++facei)
    #endif
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const Type& vsfOwn = vsf[own];
        const Type& vsfNei = vsf[nei];

        maxVsf[own] = Foam::max(maxVsf[own], vsfNei);
        minVsf[own] = Foam::min(minVsf[own], vsfNei);

        maxVsf[nei] = Foam::max(maxVsf[nei], vsfOwn);
        minVsf[nei] = Foam::min(minVsf[nei], vsfOwn);
    }


    const auto& bsf = vsf.boundaryField();

    //printf(" bsf.size()=%d\n",bsf.size());
    #if 1
    forAll(bsf, patchi)
    #else
    //LG2 AMD problem with creating target region here - need investigation
    #pragma omp target teams distribute 
    for (label patchi=0; patchi < bsf.size(); ++patchi)
    #endif
    {
        const fvPatchField<Type>& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const Field<Type> psfNei(psf.patchNeighbourField());

            //printf("patchi=%d :  pOwner.size()=%d\n",patchi,pOwner.size());
            #if 1
            forAll(pOwner, pFacei)
            #else
            #pragma omp target teams distribute parallel for
            for (label pFacei=0; pFacei < pOwner.size(); ++pFacei)            
            #endif 
            {
                const label own = pOwner[pFacei];
                const Type& vsfNei = psfNei[pFacei];

                //atomic MAX/MIN should solve the race condition issue 
                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            #if 1
            forAll(pOwner, pFacei)
            #else
            #pragma omp target teams distribute parallel for
            for (label pFacei=0; pFacei < pOwner.size(); ++pFacei)            
            #endif             
            {
                const label own = pOwner[pFacei];
                const Type& vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const Field<Type> maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;
    }


    // Create limiter initialized to 1
    // Note: the limiter is not permitted to be > 1
    Field<Type> limiter(vsf.primitiveField().size(), pTraits<Type>::one);

    #if 1
    forAll(owner, facei)
    #else
    //LG2 AMD race conditions in parallel implementation ?
    #pragma omp target teams distribute parallel for
    for (label facei=0; facei < owner.size(); ++facei)    
    #endif
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitFace
        (
            limiter[own],
            maxVsf[own],
            minVsf[own],
            (Cf[facei] - C[own]) & g[own]
        );

        // neighbour side
        limitFace
        (
            limiter[nei],
            maxVsf[nei],
            minVsf[nei],
            (Cf[facei] - C[nei]) & g[nei]
        );
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    limitGradient(limiter, g);
    g.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, g);

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tGrad;
}


// ************************************************************************* //
