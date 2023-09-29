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

Description
    High performance macro functions for Field\<Type\> algebra.  These expand
    using either array element access (for vector machines) or pointer
    dereferencing for scalar machines as appropriate.

\*---------------------------------------------------------------------------*/
#ifdef USE_OMP
  #include <omp.h>
  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif 
#endif

#ifndef TARGET_CUT_OFF
#define TARGET_CUT_OFF 10000
#endif

#ifndef FieldM_H
#define FieldM_H

#include "error.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef FULLDEBUG

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const char* op
)
{
    if (f1.size() != f2.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << " and Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << endl
            << " for operation " << op
            << abort(FatalError);
    }
}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>& f1,
    const UList<Type2>& f2,
    const UList<Type3>& f3,
    const char* op
)
{
    if (f1.size() != f2.size() || f1.size() != f3.size())
    {
        FatalErrorInFunction
            << " Field<"<<pTraits<Type1>::typeName<<"> f1("<<f1.size()<<')'
            << ", Field<"<<pTraits<Type2>::typeName<<"> f2("<<f2.size()<<')'
            << " and Field<"<<pTraits<Type3>::typeName<<"> f3("<<f3.size()<<')'
            << endl
            << "    for operation " << op
            << abort(FatalError);
    }
}

#else

template<class Type1, class Type2>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const char*
)
{}

template<class Type1, class Type2, class Type3>
void checkFields
(
    const UList<Type1>&,
    const UList<Type2>&,
    const UList<Type3>&,
    const char*
)
{}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Member function : f1 OP Func f2

#define TFOR_ALL_F_OP_FUNC_F(typeF1, f1, OP, FUNC, typeF2, f2)                 \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2)");                           \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    label loop_len = (f1).size();                                                 \
    /* Loop: f1 OP FUNC(f2) */                                                 \
    /* _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)")  */ \
    for (label i = 0; i < loop_len; ++i)                                          \
    {                                                                          \
        (f1P[i]) OP FUNC(f2P[i]);                                              \
    }


#define TPARALLELFOR_ALL_F_OP_FUNC_F_ARITHMETIC(typeF1, f1, OP, FUNC, typeF2, f2)                 \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2)");                           \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    label loop_len = (f1).size();                                                 \
    /* Loop: f1 OP FUNC(f2) */                                                 \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000) ")  \
    for (label i = 0; i < loop_len; ++i)                                          \
    {                                                                          \
        (f1P[i]) OP FUNC(f2P[i]);                                              \
    }


#define TFOR_ALL_F_OP_F_FUNC(typeF1, f1, OP, typeF2, f2, FUNC)                 \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2" #FUNC);                                \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP f2.FUNC() */                                                \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP (f2P[i]).FUNC();                                           \
    }


#define TPARALLELFOR_ALL_F_OP_F_FUNC(typeF1, f1, OP, typeF2, f2, FUNC)                 \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2" #FUNC);                                \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                              \
    /* Loop: f1 OP f2.FUNC() */                                                \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000) ")  \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]).FUNC();                                           \
    }




// Member function : this field f1 OP FUNC(f2, f3)

#define TFOR_ALL_F_OP_FUNC_F_F(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3)   \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3)");                   \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3) */                                             \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]));                                  \
    }


#define TPARALLELFOR_ALL_F_OP_FUNC_F_F(typeF1, f1, OP, FUNC, typeF2, f2, typeF3, f3)   \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP " " #FUNC "(f2, f3)");                   \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
    label loop_len = (f1).size();                                              \
                                                                               \
    /* Loop: f1 OP FUNC(f2, f3) */                                             \
    /* List_FOR_ALL(f1, i) */                                                  \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (f3P[i]));                                  \
    }



// Member function : s OP FUNC(f1, f2)

#define TFOR_ALL_S_OP_FUNC_F_F(typeS, s, OP, FUNC, typeF1, f1, typeF2, f2)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "s " #OP " " #FUNC "(f1, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: s OP FUNC(f1, f2) */                                              \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (s) OP FUNC((f1P[i]), (f2P[i]));                                       \
    }


// Member function : this f1 OP FUNC(f2, s)

#define TFOR_ALL_F_OP_FUNC_F_S(typeF1, f1, OP, FUNC, typeF2, f2, typeS, s)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2, s)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP FUNC(f2, s) */                                              \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (s));                                       \
    }


#define TPARALLELFOR_ALL_F_OP_FUNC_F_S(typeF1, f1, OP, FUNC, typeF2, f2, typeS, s)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(f2, s)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                               \
    /* Loop: f1 OP FUNC(f2, s) */                                              \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((f2P[i]), (s));                                       \
    }


// Member function : s1 OP FUNC(f, s2)

#define TFOR_ALL_S_OP_FUNC_F_S(typeS1, s1, OP, FUNC, typeF, f, typeS2, s2)     \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* Loop: s1 OP FUNC(f, s2) */                                              \
    List_FOR_ALL(f, i)                                                         \
    {                                                                          \
        (s1) OP FUNC((fP[i]), (s2));                                           \
    }


#define TPARALLELFOR_ALL_S_OP_MIN_REDUCTION_F_S_ARITHMETIC(typeS1, s, OP, FUNC, typeF, f)     \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
    label loop_len = (f).size();                                               \
    typeS1 min_val = (s);                                                      \
    /* Loop: s OP FUNC(f, s) */                                                \
    _Pragma("omp target teams distribute parallel for reduction(min:min_val) map(tofrom:min_val) if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        /*(s1) OP FUNC((fP[i]), (s2));     */                                  \
        min_val = (fP[i]) < min_val ?  (fP[i]) : min_val;                      \
    }                                                                          \
    (s) = min_val;


#define TPARALLELFOR_ALL_S_OP_MAX_REDUCTION_F_S_ARITHMETIC(typeS1, s, OP, FUNC, typeF, f)     \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
    label loop_len = (f).size();                                               \
    typeS1 max_val = (s);                                                      \
    /* Loop: s OP FUNC(f, s) */                                                \
    _Pragma("omp target teams distribute parallel for reduction(max:max_val) map(tofrom:max_val) if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        /*(s1) OP FUNC((fP[i]), (s2));     */                                  \
        max_val = (fP[i]) > max_val ?  (fP[i]) : max_val;                      \
    }                                                                          \
    (s) = max_val;

// Member function : this f1 OP FUNC(s, f2)

#define TFOR_ALL_F_OP_FUNC_S_F(typeF1, f1, OP, FUNC, typeS, s, typeF2, f2)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(s, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP FUNC((s), (f2P[i]));                                       \
    }



#define TPARALLELFOR_ALL_F_OP_FUNC_S_F(typeF1, f1, OP, FUNC, typeS, s, typeF2, f2)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " " #FUNC "(s, f2)");                        \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                               \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
    /* List_FOR_ALL(f1, i)   */                                                \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP FUNC((s), (f2P[i]));                                       \
    }



// Member function : this f1 OP FUNC(s1, s2)

#define TFOR_ALL_F_OP_FUNC_S_S(typeF1, f1, OP, FUNC, typeS1, s1, typeS2, s2)   \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
                                                                               \
    /* Loop: f1 OP FUNC(s1, s2) */                                             \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP FUNC((s1), (s2));                                          \
    }


// Member function : this f1 OP f2 FUNC(s)

#define TFOR_ALL_F_OP_F_FUNC_S(typeF1, f1, OP, typeF2, f2, FUNC, typeS, s)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2 " #FUNC "(s)");                         \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP f2 FUNC(s) */                                               \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP (f2P[i]) FUNC((s));                                        \
    }



#define TPARALLELFOR_ALL_F_OP_F_FUNC_S(typeF1, f1, OP, typeF2, f2, FUNC, typeS, s)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2 " #FUNC "(s)");                         \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                               \
    /* Loop: f1 OP f2 FUNC(s) */                                               \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)")  \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP (f2P[i]) FUNC((s));                                        \
    }




// Member operator : this field f1 OP1 f2 OP2 f3

#define TFOR_ALL_F_OP_F_OP_F(typeF1, f1, OP1, typeF2, f2, OP2, typeF3, f3)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP1 " f2 " #OP2 " f3");                     \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
    label f1_sz = (f1).size();                                                 \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
    for (label i = 0; i < f1_sz; ++i)                                          \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (f3P[i]);                                    \
    }


#define TPARALLELFOR_ALL_F_OP_F_OP_F(typeF1, f1, OP1, typeF2, f2, OP2, typeF3, f3)     \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, f3, "f1 " #OP1 " f2 " #OP2 " f3");                     \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    List_CONST_ACCESS(typeF3, f3, f3P);                                        \
    label loop_len = (f1).size();                                              \
                                                                               \
    /* Loop: f1 OP1 f2 OP2 f3 */                                               \
     _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)")  \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (f3P[i]);                                    \
    }





// Member operator : this field f1 OP1 www OP2 f2

#define TFOR_ALL_F_OP_S_OP_F(typeF1, f1, OP1, typeS, s, OP2, typeF2, f2)       \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " s " #OP2 " f2");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: f1 OP1 s OP2 f2 */                                                \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP1 (s) OP2 (f2P[i]);                                         \
    }



#define TPARALLELFOR_ALL_F_OP_S_OP_F(typeF1, f1, OP1, typeS, s, OP2, typeF2, f2)       \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " s " #OP2 " f2");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                              \
                                                                               \
    /* Loop: f1 OP1 s OP2 f2 */                                                \
    /* List_FOR_ALL(f1, i)   */                                                \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)")  \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (s) OP2 (f2P[i]);                                         \
    }



// Member operator : this field f1 OP1 f2 OP2 s

#define TFOR_ALL_F_OP_F_OP_S(typeF1, f1, OP1, typeF2, f2, OP2, typeS, s)       \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " f2 " #OP2 " s");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop f1 OP1 s OP2 f2 */                                                 \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (s);                                         \
    }


#define TPARALLELFOR_ALL_F_OP_F_OP_S(typeF1, f1, OP1, typeF2, f2, OP2, typeS, s)       \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP1 " f2 " #OP2 " s");                          \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
    label loop_len = (f1).size();                                              \
                                                                               \
    /* Loop f1 OP1 s OP2 f2 */                                                 \
    /*List_FOR_ALL(f1, i) */                                                   \
    _Pragma("omp target teams distribute parallel for if(target:loop_len > 20000)")  \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        (f1P[i]) OP1 (f2P[i]) OP2 (s);                                         \
    }



// Member operator : this field f1 OP f2

#define TFOR_ALL_F_OP_F(typeF1, f1, OP, typeF2, f2)                            \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, "f1 " #OP " f2");                                      \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    label f1_sz = (f1).size();                                                 \
    /* Loop: f1 OP f2 */                                                       \
    _Pragma("omp target teams distribute parallel for if(target:f1_sz > 20000)") \
    for (label i = 0; i < f1_sz; ++i)                                          \
    {                                                                          \
        (f1P[i]) OP (f2P[i]);                                                  \
    }

// Member operator : this field f1 OP1 OP2 f2

#define TFOR_ALL_F_OP_OP_F(typeF1, f1, OP1, OP2, typeF2, f2)                   \
                                                                               \
    /* Check fields have same size */                                          \
    checkFields(f1, f2, #OP1 " " #OP2 " f2");                                  \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF1, f1, f1P);                                              \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    label f1_sz = (f1).size();                                                 \
    /* Loop: f1 OP1 OP2 f2 */                                                  \
    /* _Pragma("omp target teams distribute parallel for if(target:f1_sz > 0)") */  \
    for (label i = 0; i < f1_sz; ++i)                                          \
    {                                                                          \
        (f1P[i]) OP1 OP2 (f2P[i]);                                             \
    }


// Member operator : this field f OP s

#define TFOR_ALL_F_OP_S(typeF, f, OP, typeS, s)                                \
                                                                               \
    /* Field access */                                                         \
    List_ACCESS(typeF, f, fP);                                                 \
                                                                               \
    label f1_sz = (f).size();                                                  \
    /* Loop: f OP s */                                                         \
    _Pragma("omp target teams distribute parallel for if(target:f1_sz > 20000)")  \
    for (label i = 0; i < f1_sz; ++i)                                          \
    {                                                                          \
        (fP[i]) OP (s);                                                        \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Friend operator function : s OP f, allocates storage for s

#define TFOR_ALL_S_OP_F(typeS, s, OP, typeF, f)                                \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
    label loop_len = (f).size();                                                  \
                                                                               \
    /* Loop: s OP f */                                                         \
    /* _Pragma("omp target teams distribute parallel for if(target:loop_len > 2000)") */ \
    for (label i = 0; i < loop_len; ++i)                                          \
    {                                                                          \
        (s) OP (fP[i]);                                                        \
    }


// Friend operator function : s OP1 f1 OP2 f2, allocates storage for s

#define TFOR_ALL_S_OP_F_OP_F(typeS, s, OP1, typeF1, f1, OP2, typeF2, f2)       \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF1, f1, f1P);                                        \
    List_CONST_ACCESS(typeF2, f2, f2P);                                        \
                                                                               \
    /* Loop: s OP f */                                                         \
    List_FOR_ALL(f1, i)                                                        \
    {                                                                          \
        (s) OP1 (f1P[i]) OP2 (f2P[i]);                                         \
    }


// Friend operator function : s OP FUNC(f), allocates storage for s
#define TFOR_ALL_S_OP_FUNC_F(typeS, s, OP, FUNC, typeF, f)                     \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
                                                                               \
    /* Loop: s OP FUNC(f) */                                                   \
    List_FOR_ALL(f, i)                                                         \
    {                                                                          \
        (s) OP FUNC(fP[i]);                                                    \
    }


#define TPARALLELFOR_ALL_S_OP_FUNC_REDUCTION_F_ARITHMETIC(typeS, s, OP, FUNC, typeF, f)     \
                                                                               \
    /* Field access */                                                         \
    List_CONST_ACCESS(typeF, f, fP);                                           \
    label loop_len = (f).size();                                               \
    typeS sum = 0.0;                                                           \
    /* Loop: s OP FUNC(f) */                                                   \
    _Pragma("omp target teams distribute parallel for reduction(+:sum) map(tofrom:sum) if(target:loop_len > 20000)") \
    for (label i = 0; i < loop_len; ++i)                                       \
    {                                                                          \
        sum OP FUNC(fP[i]);                                                    \
    }                                                                          \
    (s) OP sum;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
