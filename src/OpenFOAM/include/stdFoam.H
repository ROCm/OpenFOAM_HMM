/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Namespace
    stdFoam

Description
    Global macros and templates used by OpenFOAM and some standard
    C++ headers.

    Some of the templates defined here correspond to useful
    std templates that are part of future C++ standards, or that
    are in a state of change. Defining them here provides some additional
    control over which definitions are used within the OpenFOAM code-base.

SeeAlso
    - http://en.cppreference.com/w/cpp/iterator/begin
    - http://en.cppreference.com/w/cpp/iterator/end
    - http://en.cppreference.com/w/cpp/iterator/rbegin
    - http://en.cppreference.com/w/cpp/iterator/rend

\*---------------------------------------------------------------------------*/

#ifndef stdFoam_H
#define stdFoam_H

#include <initializer_list>
#include <memory>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compile-time warning for use of deprecated methods (compiler-dependent).
// Use within the class declaration.

#if (__cplusplus >= 201402L)
# define FOAM_DEPRECATED(since) [[deprecated("Since " #since)]]
# define FOAM_DEPRECATED_FOR(since, replacement) [[deprecated("Since " #since "; use " #replacement)]]
#elif defined(__GNUC__)
# define FOAM_DEPRECATED(since) __attribute__((__deprecated__("Since " #since)))
# define FOAM_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))
#else
# define FOAM_DEPRECATED(since)
# define FOAM_DEPRECATED_FOR(since, replacement)
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Namespace for OpenFOAM
namespace Foam
{
    //- Implementation details for various OpenFOAM classes
    namespace Detail {}

    //- Additional OpenFOAM modules
    namespace Module {}

    // Standard items

    //- Allow use of std::unique_ptr directly as Foam::unique_ptr
    using std::unique_ptr;
}


namespace stdFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward iteration

//- Return iterator to the beginning of the container \a c.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c.
//  Definition as per std::cbegin C++17
template<class C>
constexpr auto cbegin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return iterator to the end of the container \a c.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(const C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c.
//  Definition as per std::cend C++17
template<class C>
constexpr auto cend(const C& c) -> decltype(c.end())
{
    return c.end();
}


// Reverse iteration

//- Return reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::rbegin C++17
template<class C>
constexpr auto rbegin(C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return const_reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::rbegin C++17
template<class C>
constexpr auto rbegin(const C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return const_reverse_iterator to the reverse-begin of container \a c.
//  Definition as per std::crbegin C++17
template<class C>
constexpr auto crbegin(const C& c) -> decltype(c.rbegin())
{
    return c.rbegin();
}

//- Return reverse_iterator to reverse-end of container \a c.
//  Definition as per std::rend C++17
template<class C>
constexpr auto rend(C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return const_reverse_iterator to reverse-end of container \a c.
//  Definition as per std::rend C++17
template<class C>
constexpr auto rend(const C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return const_reverse_iterator to reverse-end of container \a c.
//  Definition as per std::crend C++17
template<class C>
constexpr auto crend(const C& c) -> decltype(c.rend())
{
    return c.rend();
}

//- Return the lesser of the parameters.
//  Definition as per std::min C++14
template<class T>
constexpr inline const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}

//- Return the greater of the parameters.
//  Definition as per std::max C++14
template<class T>
constexpr inline const T& max(const T& a, const T& b)
{
    return (a < b) ? b : a;
}

} // End namespace stdFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Iterate across all elements in the \a container object.
//  \par Usage
//  \code
//  forAllIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIters, forAllIter, forAllConstIters
#define forAllIters(container,iter)                                            \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::begin(container);                                 \
        iter != stdFoam::end(container);                                       \
        ++iter                                                                 \
    )


//- Iterate across all elements of the \a container object with const access.
//  \par Usage
//  \code
//  forAllConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIters, forAllIter, forAllConstIter
#define forAllConstIters(container,iter)                                       \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::cbegin(container);                                \
        iter != stdFoam::cend(container);                                      \
        ++iter                                                                 \
    )


//- Reverse iterate across elements in the \a container object of type
//  \a Container.
//  \par Usage
//  \code
//  forAllReverseIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstReverseIters
#define forAllReverseIters(container,iter)                                     \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::rbegin(container);                                \
        iter != stdFoam::rend(container);                                      \
        ++iter                                                                 \
    )


//- Reverse iterate across elements of \a container object with const access.
//  \par Usage
//  \code
//  forAllReverseConstIters(container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllReverseIters
#define forAllConstReverseIters(container,iter)                                \
    for                                                                        \
    (                                                                          \
        auto iter = stdFoam::crbegin(container);                               \
        iter != stdFoam::crend(container);                                     \
        ++iter                                                                 \
    )


//- Loop across all elements in \a list
// \par Usage
// \code
// forAll(anyList, i)
// {
//      statements;
// }
// \endcode
// \sa forAllReverse
#define forAll(list, i) \
    for (Foam::label i=0; i<(list).size(); ++i)


//- Reverse loop across all elements in \a list
//  \par Usage
//  \code
//  forAllReverse(anyList, i)
//  {
//       statements;
//  }
//  \endcode
//  \sa forAll
#define forAllReverse(list, i) \
    for (Foam::label i=(list).size()-1; i>=0; --i)


// Compatibility macros for pre C++11

//- Iterate across all elements in the \a container object
//  of type \a Container.
//  \par Usage
//  \code
//  forAllIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllConstIter
#define forAllIter(Container,container,iter)                                   \
    for                                                                        \
    (                                                                          \
        Container::iterator iter = (container).begin();                        \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )


//- Iterate across all elements in the \a container object
//  of type \a Container with const access.
//  \par Usage
//  \code
//  forAllConstIter(ContainerType, container, iter)
//  {
//      statements;
//  }
//  \endcode
//  \sa forAllIter
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).cbegin();                 \
        iter != (container).cend();                                            \
        ++iter                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
