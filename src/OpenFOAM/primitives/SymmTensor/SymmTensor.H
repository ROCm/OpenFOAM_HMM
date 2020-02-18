/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::SymmTensor

Description
    A templated (3 x 3) symmetric tensor of objects of \<T\>, effectively
    containing 6 elements, derived from VectorSpace.

See also
    Test-SymmTensor.C

SourceFiles
    SymmTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef SymmTensor_H
#define SymmTensor_H

#include "contiguous.H"
#include "Vector.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SymmTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor
:
    public VectorSpace<SymmTensor<Cmpt>, Cmpt, 6>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef SymmTensor<label> labelType;


    // Member Constants

        //- Rank of SymmTensor is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const SymmTensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YY, YZ, ZZ };


    // Generated Methods

        //- Default construct
        SymmTensor() = default;

        //- Copy construct
        SymmTensor(const SymmTensor&) = default;

        //- Copy assignment
        SymmTensor& operator=(const SymmTensor&) = default;


    // Constructors

        //- Construct initialized to zero
        inline SymmTensor(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline SymmTensor(const VectorSpace<SymmTensor<Cmpt2>, Cmpt2, 6>&);

        //- Construct given SphericalTensor
        inline SymmTensor(const SphericalTensor<Cmpt>&);

        //- Construct given the six components
        inline SymmTensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
                            const Cmpt tyy, const Cmpt tyz,
                                            const Cmpt tzz
        );

        //- Construct from Istream
        inline explicit SymmTensor(Istream& is);


    // Member Functions

        // Component access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& xz() const;
            inline const Cmpt& yx() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& yz() const;
            inline const Cmpt& zx() const;
            inline const Cmpt& zy() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& xz();
            inline Cmpt& yx();
            inline Cmpt& yy();
            inline Cmpt& yz();
            inline Cmpt& zx();
            inline Cmpt& zy();
            inline Cmpt& zz();


        // Diagonal access and manipulation

            //- Extract the diagonal as a vector
            inline Vector<Cmpt> diag() const;

            //- Set values of the diagonal
            inline void diag(const Vector<Cmpt>& v);


    // Tensor Operations

        //- Return non-Hermitian transpose
        inline const SymmTensor<Cmpt>& T() const;


    // Member Operators

        //- Inherit VectorSpace assignment operators
        using SymmTensor::vsType::operator=;

        //- Assign to given SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SymmTensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SymmTensor<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SymmTensor<Cmpt>> : is_contiguous_scalar<Cmpt> {};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 2>
{
public:

    typedef SymmTensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSolve<SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SymmTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
