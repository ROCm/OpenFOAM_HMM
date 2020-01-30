/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::Barycentric2D

Description
    Templated 2D Barycentric derived from VectorSpace.
    Has 3 components, one of which is redundant.

SourceFiles
    Barycentric2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Barycentric2D_H
#define Barycentric2D_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Barycentric2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Barycentric2D
:
    public VectorSpace<Barycentric2D<Cmpt>, Cmpt, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Barycentric2D<label> labelType;


    // Member Constants

        //- Rank of Barycentric2D is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { A, B, C };


    // Generated Methods: copy construct/assignment

        //- Default construct
        Barycentric2D() = default;


    // Constructors

        //- Construct initialized to zero
        inline Barycentric2D(const Foam::zero);

        //- Construct from components
        inline Barycentric2D
        (
            const Cmpt& va,
            const Cmpt& vb,
            const Cmpt& vc
        );


    // Member Functions

    // Access

        inline const Cmpt& a() const;
        inline const Cmpt& b() const;
        inline const Cmpt& c() const;


    // Edit

        inline Cmpt& a();
        inline Cmpt& b();
        inline Cmpt& c();


    // Tests

        //- True if any coordinates are negative
        inline bool outside() const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Barycentric2D<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Barycentric2D<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Barycentric2D<Cmpt>> : is_contiguous_scalar<Cmpt>{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Barycentric2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
