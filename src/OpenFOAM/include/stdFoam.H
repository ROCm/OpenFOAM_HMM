/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Some global templates and macros used by OpenFOAM and some standard
    C++ headers.

    Some of the templates defined here correspond to useful
    std templates that are part of future C++ standards, or that
    are in a state of change. Defining them here provides some additional
    control over which definitions are used within the OpenFOAM code-base.

SeeAlso
    - http://en.cppreference.com/w/cpp/iterator/end
    - http://en.cppreference.com/w/cpp/iterator/begin

\*---------------------------------------------------------------------------*/

#ifndef StdFoam_H
#define StdFoam_H

#include <initializer_list>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace stdFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return iterator to the beginning of the container \a c or array.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c or array.
//  Definition as per std::begin C++17
template<class C>
constexpr auto begin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return const_iterator to the beginning of the container \a c or array.
//  Definition as per std::cbegin C++17
template<class C>
constexpr auto cbegin(const C& c) -> decltype(c.begin())
{
    return c.begin();
}

//- Return iterator to the end of the container \a c or array.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c or array.
//  Definition as per std::end C++17
template<class C>
constexpr auto end(const C& c) -> decltype(c.end())
{
    return c.end();
}

//- Return const_iterator to the end of the container \a c or array.
//  Definition as per std::cend C++17
template<class C>
constexpr auto cend(const C& c) -> decltype(c.end())
{
    return c.end();
}


} // End namespace stdFoam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Iterate across all elements in the \a container object of type
//  \a Container.
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
