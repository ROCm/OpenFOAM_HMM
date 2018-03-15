/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
            const auto& gf = tgf();
            const auto& gbf = gf.boundaryField();

            forAll(gbf, patchi)
            {
                if
                (
                    !polyPatch::constraintType(gbf[patchi].patch().type())
                 && !isA<typename PatchField<Type>::Calculated>(gbf[patchi])
                )
                {
                    WarningInFunction
                        << "Attempt to reuse temporary with non-reusable BC "
                        << gbf[patchi].type() << endl;

                    return false;
                }
            }
        }

        return true;
    }

    return false;
}


template<class TypeR, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
(
    const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initRet = false
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

    tmp<GeometricField<TypeR, PatchField, GeoMesh>> rtgf
    (
        new GeometricField<TypeR, PatchField, GeoMesh>
        (
            IOobject
            (
                name,
                gf1.instance(),
                gf1.db()
            ),
            gf1.mesh(),
            dimensions
        )
    );

    if (initRet)
    {
        rtgf.ref() == gf1;
    }

    return rtgf;
}


template
<
    class TypeR,
    class Type1,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
class reuseTmpGeometricField<TypeR, TypeR, PatchField, GeoMesh>
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<TypeR, PatchField, GeoMesh>>& tgf1,
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

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


template
<
    class TypeR,
    class Type1,
    class Type12,
    class Type2,
    template<class> class PatchField,
    class GeoMesh
>
class reuseTmpTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, PatchField, GeoMesh>> New
    (
        const tmp<GeometricField<Type1, PatchField, GeoMesh>>& tgf1,
        const tmp<GeometricField<Type2, PatchField, GeoMesh>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const auto& gf1 = tgf1();

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
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
class reuseTmpTmpGeometricField
    <TypeR, Type1, Type12, TypeR, PatchField, GeoMesh>
{
public:

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

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
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
class reuseTmpTmpGeometricField<TypeR, TypeR, TypeR, Type2, PatchField, GeoMesh>
{
public:

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

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


template<class TypeR, template<class> class PatchField, class GeoMesh>
class reuseTmpTmpGeometricField<TypeR, TypeR, TypeR, TypeR, PatchField, GeoMesh>
{
public:

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

        return tmp<GeometricField<TypeR, PatchField, GeoMesh>>
        (
            new GeometricField<TypeR, PatchField, GeoMesh>
            (
                IOobject
                (
                    name,
                    gf1.instance(),
                    gf1.db()
                ),
                gf1.mesh(),
                dimensions
            )
        );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
