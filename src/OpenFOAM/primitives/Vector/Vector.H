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

Class
    Foam::Vector

Description
    Templated 3D Vector derived from VectorSpace adding construction from
    3 components, element access using x(), y() and z() member functions and
    the inner-product (dot-product) and cross product operators.

    A centre() member function which returns the Vector for which it is called
    is defined so that point which is a typedef to Vector\<scalar\> behaves as
    other shapes in the shape hierarchy.

SourceFiles
    VectorI.H

\*---------------------------------------------------------------------------*/

#ifndef Vector_H
#define Vector_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class List;

/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Vector
:
    public VectorSpace<Vector<Cmpt>, Cmpt, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Vector<label> labelType;


    // Member Constants

        //- Rank of Vector is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y, Z };


    // Generated Methods

        //- Default construct
        Vector() = default;

        //- Copy construct
        Vector(const Vector&) = default;

        //- Copy assignment
        Vector& operator=(const Vector&) = default;


    // Constructors

        //- Construct initialized to zero
        inline Vector(const Foam::zero);

        //- Copy construct from VectorSpace of the same rank
        template<class Cmpt2>
        inline Vector(const VectorSpace<Vector<Cmpt2>, Cmpt2, 3>& vs);

        //- Construct from three components
        inline Vector(const Cmpt& vx, const Cmpt& vy, const Cmpt& vz);

        //- Construct from Istream
        inline explicit Vector(Istream& is);


    // Member Functions

        //- Access to the vector x component
        inline const Cmpt& x() const;

        //- Access to the vector y component
        inline const Cmpt& y() const;

        //- Access to the vector z component
        inline const Cmpt& z() const;

        //- Access to the vector x component
        inline Cmpt& x();

        //- Access to the vector y component
        inline Cmpt& y();

        //- Access to the vector z component
        inline Cmpt& z();


        //- Normalise the vector by its magnitude
        inline Vector<Cmpt>& normalise();


        //- Return *this (used for point which is a typedef to Vector<scalar>.
        inline const Vector<Cmpt>& centre
        (
            const Foam::List<Vector<Cmpt>>&
        ) const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Vector<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Vector<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Vector<Cmpt>> : is_contiguous_scalar<Cmpt> {};


template<class Cmpt>
class typeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class typeOfSolve<Vector<Cmpt>>
{
public:

    typedef Vector<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
