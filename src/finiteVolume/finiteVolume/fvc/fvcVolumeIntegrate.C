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

#include "fvcVolumeIntegrate.H"
#include "fvMesh.H"
#include "Field.H"

#ifdef USE_ROCTX
#include <roctx.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    printf("in volumeIntegrate line = %d\n", __LINE__);
    return vf.mesh().V()*vf.primitiveField();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::volumeIntegrate_B");
    #endif

    tmp<Field<Type>> tvivf = tvf().mesh().V()*tvf().primitiveField();
    tvf.clear();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tvivf;
}


template<class Type>
tmp<Field<Type>> volumeIntegrate(const DimensionedField<Type, volMesh>& df)
{
    printf("in volumeIntegrate line = %d\n", __LINE__);
    return df.mesh().V()*df.field();
}


template<class Type>
tmp<Field<Type>>
volumeIntegrate(const tmp<DimensionedField<Type, volMesh>>& tdf)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::volumeIntegrate_D");
    #endif

    tmp<Field<Type>> tdidf = tdf().mesh().V()*tdf().field();
    tdf.clear();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return tdidf;
}


template<class Type>
dimensioned<Type>
domainIntegrate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    printf("in domainIntegrate line = %d\n", __LINE__);
    return dimensioned<Type>
    (
        "domainIntegrate(" + vf.name() + ')',
        dimVol*vf.dimensions(),
        gSum(fvc::volumeIntegrate(vf))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::domainIntegrate_B");
    #endif

    dimensioned<Type> integral = domainIntegrate(tvf());
    tvf.clear();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return integral;
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const DimensionedField<Type, volMesh>& df
)
{
    printf("in domainIntegrate line = %d\n", __LINE__);
    return dimensioned<Type>
    (
        "domainIntegrate(" + df.name() + ')',
        dimVol*df.dimensions(),
        gSum(fvc::volumeIntegrate(df))
    );
}


template<class Type>
dimensioned<Type> domainIntegrate
(
    const tmp<DimensionedField<Type, volMesh>>& tdf
)
{
    #ifdef USE_ROCTX
    roctxRangePush("fvc::domainIntegrate_D");
    #endif

    dimensioned<Type> integral = domainIntegrate(tdf());
    tdf.clear();

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return integral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
