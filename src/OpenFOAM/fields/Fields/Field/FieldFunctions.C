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

#include "PstreamReduceOps.H"
#include "FieldReuseFunctions.H"

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

template<class Type>
void component
(
    Field<typename Field<Type>::cmptType>& res,
    const UList<Type>& f,
    const direction d
)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_F_FUNC_S
    (
        cmptType, res, =, Type, f, .component, const direction, d
    )
}


template<class Type>
void T(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_F_FUNC(Type, res, =, Type, f, T)
}


template<class Type, direction r>
void pow
(
    Field<typename powProduct<Type, r>::type>& res,
    const UList<Type>& vf
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    TFOR_ALL_F_OP_FUNC_F_S
    (
        powProductType, res, =, pow, Type, vf, powProductType,
        pTraits<powProductType>::zero
    )
}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const UList<Type>& f,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    auto tres
    (
        tmp<Field<powProductType>>::New(f.size())
    );
    pow<Type, r>(tres.ref(), f);
    return tres;
}

template<class Type, direction r>
tmp<Field<typename powProduct<Type, r>::type>>
pow
(
    const tmp<Field<Type>>& tf,
    typename powProduct<Type, r>::type
)
{
    typedef typename powProduct<Type, r>::type powProductType;
    auto tres = reuseTmp<powProductType, Type>::New(tf);
    pow<Type, r>(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void sqr
(
    Field<typename outerProduct<Type, Type>::type>& res,
    const UList<Type>& vf
)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    TFOR_ALL_F_OP_FUNC_F(outerProductType, res, =, sqr, Type, vf)
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const UList<Type>& f)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    auto tres
    (
        tmp<Field<outerProductType>>::New(f.size())
    );
    sqr(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type>>
sqr(const tmp<Field<Type>>& tf)
{
    typedef typename outerProduct<Type, Type>::type outerProductType;
    auto tres = reuseTmp<outerProductType, Type>::New(tf);
    sqr(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void magSqr
(
    Field<typename typeOfMag<Type>::type>& res,
    const UList<Type>& f
)
{
    typedef typename typeOfMag<Type>::type magType;

    TFOR_ALL_F_OP_FUNC_F(magType, res, =, magSqr, Type, f)
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const UList<Type>& f)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres = tmp<Field<magType>>::New(f.size());
    magSqr(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
magSqr(const tmp<Field<Type>>& tf)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres = reuseTmp<magType, Type>::New(tf);
    magSqr(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void mag
(
    Field<typename typeOfMag<Type>::type>& res,
    const UList<Type>& f
)
{
    typedef typename typeOfMag<Type>::type magType;

    TFOR_ALL_F_OP_FUNC_F(magType, res, =, mag, Type, f)
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const UList<Type>& f)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres = tmp<Field<magType>>::New(f.size());
    mag(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename typeOfMag<Type>::type>>
mag(const tmp<Field<Type>>& tf)
{
    typedef typename typeOfMag<Type>::type magType;

    auto tres = reuseTmp<magType, Type>::New(tf);
    mag(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMax, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMax(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = tmp<Field<cmptType>>::New(f.size());
    cmptMax(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMax(const tmp<Field<Type>>& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = reuseTmp<cmptType, Type>::New(tf);
    cmptMax(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptMin, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMin(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = tmp<Field<cmptType>>::New(f.size());
    cmptMin(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptMin(const tmp<Field<Type>>& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = reuseTmp<cmptType, Type>::New(tf);
    cmptMin(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    TFOR_ALL_F_OP_FUNC_F(cmptType, res, =, cmptAv, Type, f)
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const UList<Type>& f)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = tmp<Field<cmptType>>::New(f.size());
    cmptAv(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<typename Field<Type>::cmptType>> cmptAv(const tmp<Field<Type>>& tf)
{
    typedef typename Field<Type>::cmptType cmptType;
    auto tres = reuseTmp<cmptType, Type>::New(tf);
    cmptAv(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(Type, res, =, cmptMag, Type, f)
}

template<class Type>
tmp<Field<Type>> cmptMag(const UList<Type>& f)
{
    auto tres = tmp<Field<Type>>::New(f.size());
    cmptMag(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<Type>> cmptMag(const tmp<Field<Type>>& tf)
{
    auto tres = New(tf);
    cmptMag(tres.ref(), tf());
    tf.clear();
    return tres;
}


template<class Type>
void cmptMagSqr(Field<Type>& res, const UList<Type>& f)
{
    TFOR_ALL_F_OP_FUNC_F(Type, res, =, cmptMagSqr, Type, f)
}

template<class Type>
tmp<Field<Type>> cmptMagSqr(const UList<Type>& f)
{
    auto tres = tmp<Field<Type>>::New(f.size());
    cmptMagSqr(tres.ref(), f);
    return tres;
}

template<class Type>
tmp<Field<Type>> cmptMagSqr(const tmp<Field<Type>>& tf)
{
    auto tres = New(tf);
    cmptMagSqr(tres.ref(), tf());
    tf.clear();
    return tres;
}


#define TMP_UNARY_FUNCTION(ReturnType, Func)                                   \
                                                                               \
template<class Type>                                                           \
ReturnType Func(const tmp<Field<Type>>& tf1)                                   \
{                                                                              \
    ReturnType res = Func(tf1());                                              \
    tf1.clear();                                                               \
    return res;                                                                \
}

template<class Type>
Type max(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Max, =, max, Type, f, Type, Max)
        return Max;
    }

    return pTraits<Type>::min;
}

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S(Type, Min, =, min, Type, f, Type, Min)
        return Min;
    }

    return pTraits<Type>::max;
}

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f)
{
    typedef typename Foam::typeOfSolve<Type>::type solveType;

    solveType Sum = Zero;

    if (f.size())
    {
        TFOR_ALL_S_OP_FUNC_F(solveType, Sum, +=, solveType, Type, f)
    }

    return Type(Sum);
}

TMP_UNARY_FUNCTION(Type, sum)


// From MinMaxOps.H:
//   - Foam::minMax(const UList<Type>&)
//   - Foam::minMaxMag(const UList<Type>&)

TMP_UNARY_FUNCTION(MinMax<Type>, minMax)
TMP_UNARY_FUNCTION(scalarMinMax, minMaxMag)


template<class Type>
Type maxMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Max(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Max,
            =,
            maxMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Max
        )
        return Max;
    }

    return Zero;
}

TMP_UNARY_FUNCTION(Type, maxMagSqr)

template<class Type>
Type minMagSqr(const UList<Type>& f)
{
    if (f.size())
    {
        Type Min(f[0]);
        TFOR_ALL_S_OP_FUNC_F_S
        (
            Type,
            Min,
            =,
            minMagSqrOp<Type>(),
            Type,
            f,
            Type,
            Min
        )
        return Min;
    }

    return pTraits<Type>::rootMax;
}

TMP_UNARY_FUNCTION(Type, minMagSqr)

template<class Type>
typename scalarProduct<Type, Type>::type
sumProd(const UList<Type>& f1, const UList<Type>& f2)
{
    typedef typename scalarProduct<Type, Type>::type prodType;

    prodType result = Zero;
    if (f1.size() && (f1.size() == f2.size()))
    {
        TFOR_ALL_S_OP_F_OP_F(prodType, result, +=, Type, f1, &&, Type, f2)
    }
    return result;
}


template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2)
{
    Type SumProd = Zero;
    if (f1.size() && (f1.size() == f2.size()))
    {
        TFOR_ALL_S_OP_FUNC_F_F
        (
            Type,
            SumProd,
            +=,
            cmptMultiply,
            Type,
            f1,
            Type,
            f2
        )
    }
    return SumProd;
}


template<class Type>
typename outerProduct1<Type>::type
sumSqr(const UList<Type>& f)
{
    typedef typename outerProduct1<Type>::type prodType;
    prodType result = Zero;
    if (f.size())
    {
        TFOR_ALL_S_OP_FUNC_F(prodType, result, +=, sqr, Type, f)
    }
    return result;
}

template<class Type>
typename outerProduct1<Type>::type
sumSqr(const tmp<Field<Type>>& tf)
{
    typedef typename outerProduct1<Type>::type prodType;
    prodType result = sumSqr(tf());
    tf.clear();
    return result;
}


template<class Type>
typename typeOfMag<Type>::type
sumMag(const UList<Type>& f)
{
    typedef typename typeOfMag<Type>::type magType;
    magType result = Zero;
    if (f.size())
    {
        TFOR_ALL_S_OP_FUNC_F(magType, result, +=, mag, Type, f)
    }
    return result;
}

TMP_UNARY_FUNCTION(typename typeOfMag<Type>::type, sumMag)


template<class Type>
Type sumCmptMag(const UList<Type>& f)
{
    Type result = Zero;
    if (f.size())
    {
        TFOR_ALL_S_OP_FUNC_F(Type, result, +=, cmptMag, Type, f)
    }
    return result;
}

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const UList<Type>& f)
{
    if (f.size())
    {
        Type avrg = sum(f)/f.size();

        return avrg;
    }

    WarningInFunction
        << "empty field, returning zero" << endl;

    return Zero;
}

TMP_UNARY_FUNCTION(Type, average)


// With reduction on ReturnType
#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                       \
                                                                               \
template<class Type>                                                           \
ReturnType gFunc(const UList<Type>& f, const label comm)                       \
{                                                                              \
    ReturnType res = Func(f);                                                  \
    reduce(res, rFunc##Op<ReturnType>(), Pstream::msgType(), comm);            \
    return res;                                                                \
}                                                                              \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(Type, gMaxMagSqr, maxMagSqr, maxMagSqr)
G_UNARY_FUNCTION(Type, gMinMagSqr, minMagSqr, minMagSqr)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

G_UNARY_FUNCTION(MinMax<Type>, gMinMax, minMax, sum)
G_UNARY_FUNCTION(scalarMinMax, gMinMaxMag, minMaxMag, sum)

G_UNARY_FUNCTION(typename outerProduct1<Type>::type, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(typename typeOfMag<Type>::type, gSumMag, sumMag, sum)

#undef G_UNARY_FUNCTION


template<class Type>
typename scalarProduct<Type, Type>::type gSumProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    typedef typename scalarProduct<Type, Type>::type prodType;

    prodType result = sumProd(f1, f2);
    reduce(result, sumOp<prodType>(), Pstream::msgType(), comm);
    return result;
}

template<class Type>
Type gSumCmptProd
(
    const UList<Type>& f1,
    const UList<Type>& f2,
    const label comm
)
{
    Type SumProd = sumCmptProd(f1, f2);
    reduce(SumProd, sumOp<Type>(), Pstream::msgType(), comm);
    return SumProd;
}

template<class Type>
Type gAverage
(
    const UList<Type>& f,
    const label comm
)
{
    label n = f.size();
    Type s = sum(f);
    sumReduce(s, n, Pstream::msgType(), comm);

    if (n > 0)
    {
        Type avrg = s/n;

        return avrg;
    }

    WarningInFunction
        << "empty field, returning zero." << endl;

    return Zero;
}

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION_FS(Type, Type, MinMax<Type>, clip)


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                  \
                                                                               \
template<class Type1, class Type2>                                             \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type1, Type2>::type>& res,                          \
    const UList<Type1>& f1,                                                    \
    const UList<Type2>& f2                                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    TFOR_ALL_F_OP_F_OP_F(productType, res, =, Type1, f1, Op, Type2, f2)        \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2)                    \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    auto tres = tmp<Field<productType>>::New(f1.size());                       \
    OpFunc(tres.ref(), f1, f2);                                                \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2>>& tf2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    auto tres = reuseTmp<productType, Type2>::New(tf2);                        \
    OpFunc(tres.ref(), f1, tf2());                                             \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const UList<Type2>& f2)              \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    auto tres = reuseTmp<productType, Type1>::New(tf1);                        \
    OpFunc(tres.ref(), tf1(), f2);                                             \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type1, class Type2>                                             \
tmp<Field<typename product<Type1, Type2>::type>>                               \
operator Op(const tmp<Field<Type1>>& tf1, const tmp<Field<Type2>>& tf2)        \
{                                                                              \
    typedef typename product<Type1, Type2>::type productType;                  \
    auto tres = reuseTmpTmp<productType, Type1, Type1, Type2>::New(tf1, tf2);  \
    OpFunc(tres.ref(), tf1(), tf2());                                          \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Type, Form>::type>& res,                            \
    const UList<Type>& f1,                                                     \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    TFOR_ALL_F_OP_F_OP_S                                                       \
        (productType, res, =,Type, f1, Op, Form, static_cast<const Form&>(vs)) \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs)     \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    auto tres = tmp<Field<productType>>::New(f1.size());                       \
    OpFunc(tres.ref(), f1, static_cast<const Form&>(vs));                      \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Type, class Form, class Cmpt, direction nCmpt>                  \
tmp<Field<typename product<Type, Form>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const tmp<Field<Type>>& tf1,                                               \
    const VectorSpace<Form,Cmpt,nCmpt>& vs                                     \
)                                                                              \
{                                                                              \
    typedef typename product<Type, Form>::type productType;                    \
    auto tres = reuseTmp<productType, Type>::New(tf1);                         \
    OpFunc(tres.ref(), tf1(), static_cast<const Form&>(vs));                   \
    tf1.clear();                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
void OpFunc                                                                    \
(                                                                              \
    Field<typename product<Form, Type>::type>& res,                            \
    const VectorSpace<Form,Cmpt,nCmpt>& vs,                                    \
    const UList<Type>& f1                                                      \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    TFOR_ALL_F_OP_S_OP_F                                                       \
        (productType, res, =,Form,static_cast<const Form&>(vs), Op, Type, f1)  \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1)     \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    auto tres = tmp<Field<productType>>::New(f1.size());                       \
    OpFunc(tres.ref(), static_cast<const Form&>(vs), f1);                      \
    return tres;                                                               \
}                                                                              \
                                                                               \
template<class Form, class Cmpt, direction nCmpt, class Type>                  \
tmp<Field<typename product<Form, Type>::type>>                                 \
operator Op                                                                    \
(                                                                              \
    const VectorSpace<Form,Cmpt,nCmpt>& vs, const tmp<Field<Type>>& tf1        \
)                                                                              \
{                                                                              \
    typedef typename product<Form, Type>::type productType;                    \
    auto tres = reuseTmp<productType, Type>::New(tf1);                         \
    OpFunc(tres.ref(), static_cast<const Form&>(vs), tf1());                   \
    tf1.clear();                                                               \
    return tres;                                                               \
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
