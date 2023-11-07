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
#include <roctracer/roctx.h>
#endif

#ifdef USE_OMP
  #include <omp.h>
  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif 
#endif

#include <type_traits>


#define USM_Cell_Limit_Grad

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


    //forAll(gIf, celli)
    #pragma omp target teams distribute parallel for if(target:gIf.size() > 20000) //LG4 tested, OK
    for (label celli=0; celli < gIf.size(); ++celli)
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
/*
template< class T >
void test_type(){
   if (std::is_same_v<T, double>) printf("test_type: T==double\n");
   if (std::is_same_v<T, float>) printf("test_type: T==float\n");
   if (std::is_same_v<T,Foam::Vector<double>>) printf("test_type: T==Foam::Vector<double>\n");
   if (std::is_same_v<T,Foam::Vector<float>>) printf("test_type: T==Foam::Vector<float>\n");
}
*/
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
        #ifdef USE_ROCTX
        roctxRangePop();
        #endif
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

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_min_max");
    #endif

    if constexpr ( std::is_same_v<Type,Foam::Vector<double>> ) {   

      #pragma omp target teams distribute parallel for if(target:owner.size() > 10000) 
      for (label facei = 0; facei < owner.size(); ++facei){
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const Foam::Vector<double>& vsfOwn = vsf[own];
        const Foam::Vector<double>& vsfNei = vsf[nei];
      
	//maxVsf[own] = Foam::max(maxVsf[own], vsfNei);
        for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
 	  double& var = setComponent(maxVsf[own],cmpt);	
          #pragma omp atomic compare            
          if (var < (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);	
	}

	//minVsf[own] = Foam::min(minVsf[own], vsfNei);
	for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
          double& var = setComponent(minVsf[own],cmpt);
          #pragma omp atomic compare
          if (var > (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);
        }
        
	//maxVsf[nei] = Foam::max(maxVsf[nei], vsfOwn);
        for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
          double& var = setComponent(maxVsf[nei],cmpt);
          #pragma omp atomic compare
          if (var < (double) component(vsfOwn,cmpt)) var = (double) component(vsfOwn,cmpt);
        }

	//minVsf[nei] = Foam::min(minVsf[nei], vsfOwn);
        for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
          double& var = setComponent(minVsf[nei],cmpt);
          #pragma omp atomic compare
          if (var > (double) component(vsfOwn,cmpt)) var = (double) component(vsfOwn,cmpt);
        }

      }
    }
    else{
      	    
      for (label facei = 0; facei < owner.size(); ++facei){
        const label own = owner[facei];
        const label nei = neighbour[facei];
        const Type& vsfOwn = vsf[own];
        const Type& vsfNei = vsf[nei];
        
	maxVsf[own] = Foam::max(maxVsf[own], vsfNei);
        minVsf[own] = Foam::min(minVsf[own], vsfNei);
        maxVsf[nei] = Foam::max(maxVsf[nei], vsfOwn);
        minVsf[nei] = Foam::min(minVsf[nei], vsfOwn);
      }
    }
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_boundary_min_max");
    #endif

    const auto& bsf = vsf.boundaryField();


    for (label patchi=0; patchi < bsf.size(); ++patchi)
    {
        const fvPatchField<Type>& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

       // printf("patchi=%d :  pOwner.size()=%d\n",patchi,pOwner.size());

        if (psf.coupled())
        {
            const Field<Type> psfNei(psf.patchNeighbourField());

            if constexpr ( std::is_same_v<Type,Foam::Vector<double>> ) {
              #pragma omp target teams distribute parallel for if(target:owner.size() > 10000) //LG4 testing , possibly OK
               for (label pFacei = 0; pFacei < pOwner.size(); ++pFacei){
		 const label own = pOwner[pFacei];
                 const Type& vsfNei = psfNei[pFacei];

                 //maxVsf[own] = max(maxVsf[own], vsfNei);
		 for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
                    double& var = setComponent(maxVsf[own],cmpt);
                    #pragma omp atomic compare
                    if (var < (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);
                 }
                 //minVsf[own] = min(minVsf[own], vsfNei);
		 for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
                    double& var = setComponent(minVsf[own],cmpt);
                    #pragma omp atomic compare
                    if (var > (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);
                 }
	       }
	    }
	    else {

              forAll(pOwner, pFacei)
              {
                const label own = pOwner[pFacei];
                const Type& vsfNei = psfNei[pFacei];

                //atomic MAX/MIN should solve the race condition issue 
                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
              }
	    }
        }
        else
        {
            if constexpr ( std::is_same_v<Type,Foam::Vector<double>> ) {
              #pragma omp target teams distribute parallel for if(target:owner.size() > 10000) //LG4 testing , possibly OK
              for (label pFacei = 0; pFacei < pOwner.size(); ++pFacei){
                const label own = pOwner[pFacei];
                const Type& vsfNei = psf[pFacei];
                
		//maxVsf[own] = max(maxVsf[own], vsfNei);
                for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
                    double& var = setComponent(maxVsf[own],cmpt);
                    #pragma omp atomic compare
                     if (var < (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);
                 }
		 
                //minVsf[own] = min(minVsf[own], vsfNei);
                for (direction cmpt=0; cmpt<pTraits<Foam::Vector<double>>::nComponents; ++cmpt){
                    double& var = setComponent(minVsf[own],cmpt);
                    #pragma omp atomic compare
                     if (var > (double) component(vsfNei,cmpt)) var = (double) component(vsfNei,cmpt);
                 }
	      }
	    }
	    else
	    {

              forAll(pOwner, pFacei)             
              {
                const label own = pOwner[pFacei];
                const Type& vsfNei = psf[pFacei];
                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
              }
	    }
        }
     }

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_C:update");
    #endif

    maxVsf -= vsf;
    minVsf -= vsf;

    if (k_ < 1.0)
    {
        const Field<Type> maxMinVsf((1.0/k_ - 1.0)*(maxVsf - minVsf));
        maxVsf += maxMinVsf;
        minVsf -= maxMinVsf;
    }
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    // Create limiter initialized to 1
    // Note: the limiter is not permitted to be > 1
    Field<Type> limiter(vsf.primitiveField().size(), pTraits<Type>::one);

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_C:limitFace");
    #endif
   
    #if 0
    forAll(owner, facei)
    #else
    #pragma omp target teams distribute parallel for if(target:owner.size() > 20000)
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

        #if 0
        forAll(pOwner, pFacei)
        #else
        #pragma omp target teams distribute parallel for if(target:owner.size() > 20000)
        for (label pFacei = 0; pFacei < pOwner.size(); ++pFacei)
        #endif
        {
            const label own = pOwner[pFacei]; //it is possible to have same own for different pFacei

            limitFace
            (
                limiter[own],
                maxVsf[own],
                minVsf[own],
                ((pCf[pFacei] - C[own]) & g[own])
            );
        }
    }
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_C:limitGradient");
    #endif

    limitGradient(limiter, g);
    
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    #ifdef USE_ROCTX
    roctxRangePush("fv::cellLimitedGrad_C:correctBoundaryConditions");
    #endif
    g.correctBoundaryConditions();
    gaussGrad<Type>::correctBoundaryConditions(vsf, g);
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tGrad;
}


// ************************************************************************* //
