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

#ifndef GeometricFieldReuseFunctions_H
#define GeometricFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Test if a tmp GeometricField is reusable, which means that it has been
//- allocated (not a const-reference).
//
//  When GeometricField debug is on, also checks that the patches have
//  reusable types.
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
bool reusable(const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf)
{
    if (tgf.isTmp())
    {
        if (GeometricField<Type, PatchField, GeoMesh>::debug)
        {
            for (const auto& p : tgf().boundaryField())
            {
                if
                (
                    !polyPatch::constraintType(p.patch().type())
                 && !isA<typename PatchField<Type>::Calculated>(p)
                )
                {
                    WarningInFunction
                        << "Attempt to reuse temporary with non-reusable BC "
                        << p.type() << endl;

                    return false;
                }
            }
        }

        return true;
    }

    return false;
}


template
<
    class TypeR,
    class Type1,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpGeometricField
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
struct reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>
{
    //- Since the input and output types are identical, also allow
    //- optional copy assignment of the initial content.
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
        const word& name,
        const dimensionSet& dimensions,
        const bool initCopy = false
    )
    {
        if (reusable(tgf1))
        {
            auto& gf1 = tgf1.constCast();

            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }

        const auto& gf1 = tgf1();

        auto rtgf = tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );

        if (initCopy)
        {
            rtgf.ref() == gf1;
        }

        return rtgf;
    }
};


//- This global function forwards to reuseTmpGeometricField::New
template<class TypeR, template<class> class PatchField, class GeoMesh>
tmp
<
    GeometricField<TypeR, PatchField, GeoMesh>
> New
(
    const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initCopy = false
)
{
    return reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>::New
    (
        tgf1,
        name,
        dimensions,
        initCopy
    );
}


template
<
    class TypeR,
    class Type1,
    class Type12,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );
    }
};


template
<
    class TypeR,
    class Type1,
    class Type12,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
<
    TypeR, Type1, Type12, TypeR, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (reusable(tgf2))
        {
            auto& gf2 = tgf2.constCast();

            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            return tgf2;
        }

        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );
    }
};


template
<
    class TypeR,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
struct reuseTmpTmpGeometricField
<
    TypeR, TypeR, TypeR, Type2, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (reusable(tgf1))
        {
            auto& gf1 = tgf1.constCast();

            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }

        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
struct reuseTmpTmpGeometricField
<
    TypeR, TypeR, TypeR, TypeR, PatchField, GeoMesh
>
{
    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        if (reusable(tgf1))
        {
            auto& gf1 = tgf1.constCast();

            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }
        else if (reusable(tgf2))
        {
            auto& gf2 = tgf2.constCast();

            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            return tgf2;
        }

        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>::New
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        );
    }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
