/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "DimensionedFieldReuseFunctions.H"

#define TEMPLATE template<class Type, class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh, direction r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh>>
pow
(
    const DimensionedField<Type, GeoMesh>& df,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    auto tres =
        tmp<DimensionedField<powProductType, GeoMesh>>::New
        (
            IOobject
            (
                "pow(" + df.name() + ',' + name(r) + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            pow(df.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tres.ref().field(), df.field());

    return tres;
}


template<class Type, class GeoMesh, direction r>
tmp<DimensionedField<typename powProduct<Type, r>::type, GeoMesh>>
pow
(
    const tmp<DimensionedField<Type, GeoMesh>>& tdf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    auto tres =
        reuseTmpDimensionedField<powProductType, Type, GeoMesh>::New
        (
            tdf,
            "pow(" + df.name() + ',' + name(r) + ')',
            pow(df.dimensions(), r)
        );

    pow<Type, r, GeoMesh>(tres.ref().field(), df.field());

    tdf.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    auto tres =
        tmp<DimensionedField<outerProductType, GeoMesh>>::New
        (
            IOobject
            (
                "sqr(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            sqr(df.dimensions())
        );

    sqr(tres.ref().field(), df.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename outerProduct<Type, Type>::type, GeoMesh>>
sqr(const tmp<DimensionedField<Type, GeoMesh>>& tdf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    auto tres =
        reuseTmpDimensionedField<outerProductType, Type, GeoMesh>::New
        (
            tdf,
            "sqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    sqr(tres.ref().field(), df.field());

    tdf.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres =
        tmp<DimensionedField<magType, GeoMesh>>::New
        (
            IOobject
            (
                "magSqr(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            sqr(df.dimensions())
        );

    magSqr(tres.ref().field(), df.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
magSqr(const tmp<DimensionedField<Type, GeoMesh>>& tdf)
{
    typedef typename typeOfMag<Type>::type magType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    auto tres =
        reuseTmpDimensionedField<magType, Type, GeoMesh>::New
        (
            tdf,
            "magSqr(" + df.name() + ')',
            sqr(df.dimensions())
        );

    magSqr(tres.ref().field(), df.field());

    tdf.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres =
        tmp<DimensionedField<magType, GeoMesh>>::New
        (
            IOobject
            (
                "mag(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            df.dimensions()
        );

    mag(tres.ref().field(), df.field());

    return tres;
}

template<class Type, class GeoMesh>
tmp<DimensionedField<typename typeOfMag<Type>::type, GeoMesh>>
mag(const tmp<DimensionedField<Type, GeoMesh>>& tdf)
{
    typedef typename typeOfMag<Type>::type magType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    auto tres =
        reuseTmpDimensionedField<magType, Type, GeoMesh>::New
        (
            tdf,
            "mag(" + df.name() + ')',
            df.dimensions()
        );

    mag(tres.ref().field(), df.field());

    tdf.clear();
    return tres;
}


template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >
>
cmptAv(const DimensionedField<Type, GeoMesh>& df)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType cmptType;

    auto tres =
        tmp<DimensionedField<cmptType, GeoMesh>>::New
        (
            IOobject
            (
                "cmptAv(" + df.name() + ')',
                df.instance(),
                df.db()
            ),
            df.mesh(),
            df.dimensions()
        );

    cmptAv(tres.ref().field(), df.field());

    return tres;
}


template<class Type, class GeoMesh>
tmp
<
    DimensionedField
    <
        typename DimensionedField<Type, GeoMesh>::cmptType, GeoMesh
    >
>
cmptAv(const tmp<DimensionedField<Type, GeoMesh>>& tdf)
{
    typedef typename DimensionedField<Type, GeoMesh>::cmptType cmptType;

    const DimensionedField<Type, GeoMesh>& df = tdf();

    auto tres =
        reuseTmpDimensionedField<cmptType, Type, GeoMesh>::New
        (
            tdf,
            "cmptAv(" + df.name() + ')',
            df.dimensions()
        );

    cmptAv(tres.ref().field(), df.field());

    tdf.clear();
    return tres;
}


#define UNARY_REDUCTION_FUNCTION(returnType, func, dfunc)                      \
                                                                               \
template<class Type, class GeoMesh>                                            \
dimensioned<returnType> func                                                   \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& df                                  \
)                                                                              \
{                                                                              \
    return dimensioned<returnType>                                             \
    (                                                                          \
        #func "(" + df.name() + ')',                                           \
        df.dimensions(),                                                       \
        dfunc(df.field())                                                      \
    );                                                                         \
}                                                                              \
                                                                               \
template<class Type, class GeoMesh>                                            \
dimensioned<returnType> func                                                   \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tdf1                           \
)                                                                              \
{                                                                              \
    dimensioned<returnType> res = func(tdf1());                                \
    tdf1.clear();                                                              \
    return res;                                                                \
}

UNARY_REDUCTION_FUNCTION(Type, max, gMax)
UNARY_REDUCTION_FUNCTION(Type, min, gMin)
UNARY_REDUCTION_FUNCTION(Type, sum, gSum)
UNARY_REDUCTION_FUNCTION(Type, average, gAverage)
UNARY_REDUCTION_FUNCTION(MinMax<Type>, minMax, gMinMax)
UNARY_REDUCTION_FUNCTION(scalarMinMax, minMaxMag, gMinMaxMag)

UNARY_REDUCTION_FUNCTION(typename typeOfMag<Type>::type, sumMag, gSumMag)

#undef UNARY_REDUCTION_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clip)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate, transform)

BINARY_OPERATOR(Type, Type, scalar, *, '*', multiply)
BINARY_OPERATOR(Type, scalar, Type, *, '*', multiply)
BINARY_OPERATOR(Type, Type, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, '*', multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, '*', multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, '|', divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, op, opFunc)                                  \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    auto tres =                                                                \
        tmp<DimensionedField<productType, GeoMesh>>::New                       \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                '(' + df1.name() + #op + df2.name() + ')',                     \
                df1.instance(),                                                \
                df1.db()                                                       \
            ),                                                                 \
            df1.mesh(),                                                        \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), df2.field());                \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type1, GeoMesh>& df1,                               \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<productType, Type2, GeoMesh>::New             \
        (                                                                      \
            tdf2,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf2.clear();                                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const DimensionedField<Type2, GeoMesh>& df2                                \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<productType, Type1, GeoMesh>::New             \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf1.clear();                                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Type1, class Type2, class GeoMesh>                              \
tmp<DimensionedField<typename product<Type1, Type2>::type, GeoMesh>>           \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type1, GeoMesh>>& tdf1,                         \
    const tmp<DimensionedField<Type2, GeoMesh>>& tdf2                          \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
                                                                               \
    const DimensionedField<Type1, GeoMesh>& df1 = tdf1();                      \
    const DimensionedField<Type2, GeoMesh>& df2 = tdf2();                      \
                                                                               \
    auto tres =                                                                \
        reuseTmpTmpDimensionedField                                            \
        <productType, Type1, Type1, Type2, GeoMesh>::New                       \
        (                                                                      \
            tdf1,                                                              \
            tdf2,                                                              \
            '(' + df1.name() + #op + df2.name() + ')',                         \
            df1.dimensions() op df2.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), df2.field());                \
                                                                               \
    tdf1.clear();                                                              \
    tdf2.clear();                                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& df1,                                \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    auto tres =                                                                \
        tmp<DimensionedField<productType, GeoMesh>>::New                       \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                '(' + df1.name() + #op + dvs.name() + ')',                     \
                df1.instance(),                                                \
                df1.db()                                                       \
            ),                                                                 \
            df1.mesh(),                                                        \
            df1.dimensions() op dvs.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), dvs.value());                \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const DimensionedField<Type, GeoMesh>& df1,                                \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return df1 op dimensioned<Form>(static_cast<const Form&>(vs));             \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Type, Form>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tdf1,                          \
    const dimensioned<Form>& dvs                                               \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
                                                                               \
    const DimensionedField<Type, GeoMesh>& df1 = tdf1();                       \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<productType, Type, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            '(' + df1.name() + #op + dvs.name() + ')',                         \
            df1.dimensions() op dvs.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), df1.field(), dvs.value());                \
                                                                               \
    tdf1.clear();                                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tdf1,                          \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    return tdf1 op dimensioned<Form>(static_cast<const Form&>(vs));            \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const DimensionedField<Type, GeoMesh>& df1                                 \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
                                                                               \
    auto tres =                                                                \
        tmp<DimensionedField<productType, GeoMesh>>::New                       \
        (                                                                      \
            IOobject                                                           \
            (                                                                  \
                '(' + dvs.name() + #op + df1.name() + ')',                     \
                df1.instance(),                                                \
                df1.db()                                                       \
            ),                                                                 \
            df1.mesh(),                                                        \
            dvs.dimensions() op df1.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), dvs.value(), df1.field());                \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const DimensionedField<Type, GeoMesh>& df1                                 \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op df1;             \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Type, class GeoMesh>                                \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const dimensioned<Form>& dvs,                                              \
    const tmp<DimensionedField<Type, GeoMesh>>& tdf1                           \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
                                                                               \
    const DimensionedField<Type, GeoMesh>& df1 = tdf1();                       \
                                                                               \
    auto tres =                                                                \
        reuseTmpDimensionedField<productType, Type, GeoMesh>::New              \
        (                                                                      \
            tdf1,                                                              \
            '(' + dvs.name() + #op + df1.name() + ')',                         \
            dvs.dimensions() op df1.dimensions()                               \
        );                                                                     \
                                                                               \
    Foam::opFunc(tres.ref().field(), dvs.value(), df1.field());                \
                                                                               \
    tdf1.clear();                                                              \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type, class GeoMesh>   \
tmp<DimensionedField<typename product<Form, Type>::type, GeoMesh>>             \
operator op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const tmp<DimensionedField<Type, GeoMesh>>& tdf1                           \
)                                                                              \
{                                                                              \
    return dimensioned<Form>(static_cast<const Form&>(vs)) op tdf1;            \
}


PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
