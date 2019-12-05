/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#ifndef DimensionedFieldReuseFunctions_H
#define DimensionedFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR, class Type1, class GeoMesh>
struct reuseTmpDimensionedField
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
struct reuseTmpDimensionedField<TypeR, TypeR, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tdf1.isTmp())
        {
            auto& df1 = tdf1.constCast();

            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }

        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


//- Global function forwards to reuseTmpDimensionedField::New
template<class TypeR, class GeoMesh>
tmp<DimensionedField<TypeR, GeoMesh>> New
(
    const tmp<DimensionedField<TypeR, GeoMesh>>& tdf1,
    const word& name,
    const dimensionSet& dimensions
)
{
    return reuseTmpDimensionedField<TypeR, TypeR, GeoMesh>::New
    (
        tdf1,
        name,
        dimensions
    );
}


template<class TypeR, class Type1, class Type12, class Type2, class GeoMesh>
struct reuseTmpTmpDimensionedField
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class Type1, class Type12, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, Type1, Type12, TypeR, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tdf2.isTmp())
        {
            auto& df2 = tdf2.constCast();

            df2.rename(name);
            df2.dimensions().reset(dimensions);
            return tdf2;
        }

        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class Type2, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, Type2, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tdf1.isTmp())
        {
            auto& df1 = tdf1.constCast();

            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }

        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
struct reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, TypeR, GeoMesh>
{
    static tmp<DimensionedField<TypeR, GeoMesh>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh>>& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (tdf1.isTmp())
        {
            auto& df1 = tdf1.constCast();

            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else if (tdf2.isTmp())
        {
            auto& df2 = tdf2.constCast();

            df2.rename(name);
            df2.dimensions().reset(dimensions);
            return tdf2;
        }

        const auto& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh>>::New
        (
            IOobject
            (
                name,
                df1.instance(),
                df1.db()
            ),
            df1.mesh(),
            dimensions
        );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
