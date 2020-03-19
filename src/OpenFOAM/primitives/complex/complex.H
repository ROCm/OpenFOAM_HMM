/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::complex

Description
    A complex number, similar to the C++ complex type.

SourceFiles
    complexI.H
    complex.C

\*---------------------------------------------------------------------------*/

#ifndef complex_H
#define complex_H

#include "scalar.H"
#include "word.H"
#include "zero.H"
#include "contiguous.H"
#include <complex>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class complex;

inline scalar magSqr(const complex&);
inline scalar mag(const complex&);
inline complex sqr(const complex&);
inline const complex& min(const complex&, const complex&);
inline const complex& max(const complex&, const complex&);
inline complex limit(const complex&, const complex&);
inline const complex& sum(const complex&);
inline complex operator-(const complex&);
inline complex operator+(const complex&, const complex&);
inline complex operator+(const complex&, const scalar);
inline complex operator+(const scalar, const complex&);
inline complex operator-(const complex&, const complex&);
inline complex operator-(const complex&, const scalar);
inline complex operator-(const scalar, const complex&);
inline complex operator*(const complex&, const complex&);
inline complex operator*(const complex&, const scalar);
inline complex operator*(const scalar, const complex&);
inline complex operator/(const complex&, const complex&);
inline complex operator/(const complex&, const scalar);
inline complex operator/(const scalar, const complex&);


/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/

class complex
{
    // Private Data

        //- Real and imaginary parts
        scalar re, im;


public:

    // Generated Methods

        //- Copy construct
        complex(const complex&) = default;

        //- Copy assignment
        complex& operator=(const complex&) = default;

        //- Move construct
        complex(complex&&) = default;

        //- Move assignment
        complex& operator=(complex&&) = default;


    // Constructors

        //- Default construct, as zero-initialized
        inline constexpr complex() noexcept;

        //- Construct zero-initialized from zero class
        inline constexpr complex(const Foam::zero) noexcept;

        //- Construct from real component
        inline explicit constexpr complex(const scalar r) noexcept;

        //- Construct from real and imaginary parts
        inline constexpr complex(const scalar r, const scalar i) noexcept;

        //- Implicit construct from std::complex
        inline complex(const std::complex<float>& c);

        //- Implicit construct from std::complex
        inline complex(const std::complex<double>& c);

        //- Construct from Istream
        explicit complex(Istream& is);


    // Member Functions

    // STL getter/setter

        //- Real part of complex number - STL naming
        constexpr scalar real() const
        {
            return re;
        }

        //- Imaginary part of complex number - STL naming
        constexpr scalar imag() const
        {
            return im;
        }

        //- Set real part of complex number - STL naming
        inline void real(scalar val);

        //- Set imaginary part of complex number - STL naming
        inline void imag(scalar val);


    // Access

        //- Real part of complex number
        inline scalar Re() const;

        //- Imaginary part of complex number
        inline scalar Im() const;


    // Edit

        //- Real part of complex number
        inline scalar& Re();

        //- Imaginary part of complex number
        inline scalar& Im();


    // Operations

        //- Complex conjugate
        inline complex conjugate() const;


    // Member Operators

        //- Implicit conversion to std::complex
        operator std::complex<scalar>() const
        {
            return std::complex<scalar>(re, im);
        }


        //- Assign zero
        inline void operator=(const Foam::zero);

        //- Assign scalar (imag = zero)
        inline void operator=(const scalar s);

        inline void operator+=(const complex& c);
        inline void operator+=(const scalar s);

        inline void operator-=(const complex& c);
        inline void operator-=(const scalar s);

        inline void operator*=(const complex& c);
        inline void operator*=(const scalar s);

        inline void operator/=(const complex& c);
        inline void operator/=(const scalar s);

        inline bool operator==(const complex& c) const;
        inline bool operator!=(const complex& c) const;


    // Friend Functions

        friend scalar magSqr(const complex& c);
        friend scalar mag(const complex& c);
        friend complex sqr(const complex& c);

        //- sgn() https://en.wikipedia.org/wiki/Sign_function#Complex_signum
        friend complex sign(const complex& c);

        //- csgn() https://en.wikipedia.org/wiki/Sign_function#Complex_signum
        friend scalar csign(const complex& c);

        friend const complex& min(const complex& c1, const complex& c2);
        friend const complex& max(const complex& c1, const complex& c2);
        friend complex limit(const complex& c1, const complex& c2);
        friend const complex& sum(const complex& c);


    // Friend Operators

        friend complex operator-(const complex& c);

        friend complex operator+(const complex& c1, const complex& c2);
        friend complex operator+(const complex& c, const scalar s);
        friend complex operator+(const scalar s, const complex& c);

        friend complex operator-(const complex& c1, const complex& c2);
        friend complex operator-(const complex& c, const scalar s);
        friend complex operator-(const scalar s, const complex& c);

        friend complex operator*(const complex& c1, const complex& c2);
        friend complex operator*(const complex& c, const scalar s);
        friend complex operator*(const scalar s, const complex& c);

        friend complex operator/(const complex& c1, const complex& c2);
        friend complex operator/(const complex& c, const scalar s);
        friend complex operator/(const scalar s, const complex& c);
};


/*---------------------------------------------------------------------------*\
                      Class pTraits<complex> Declaration
\*---------------------------------------------------------------------------*/

// Template specialisation for pTraits<complex>
template<>
class pTraits<complex>
{
    complex p_;

public:

    // Typedefs

        //- Component type
        typedef complex cmptType;

        //- Magnitude type
        typedef scalar magType;

        //- Equivalent type of labels used for valid component indexing
        typedef label labelType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of complex is 0
        static constexpr direction rank = 0;

        //- Number of components in complex is 2
        static constexpr direction nComponents = 2;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];

        static const complex zero;      //!< complex (0,0)
        static const complex one;       //!< complex (1,0)
        static const complex min;       //!< complex (-VGREAT,-VGREAT)
        static const complex max;       //!< complex (VGREAT,VGREAT)
        static const complex rootMin;   //!< complex (-ROOTVGREAT, -ROOTVGREAT)
        static const complex rootMax;   //!< complex (ROOTVGREAT, ROOTVGREAT)


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(const complex& val);

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Access to the value
        operator complex() const
        {
            return p_;
        }

        //- Access to the value
        operator complex&()
        {
            return p_;
        }
};


/*---------------------------------------------------------------------------*\
                               Namespace Detail
\*---------------------------------------------------------------------------*/

namespace Detail
{
    // Helper functions for complex, in Detail namespace to avoid possible
    // name collisions (could change in the future)

    //- The 'conjugate' of non-complex returns itself (pass-through)
    //- it does not return a complex!
    template<class T>
    typename std::enable_if
    <
        !std::is_same<complex, T>::value,
        const T&
    >::type conj(const T& val)
    {
        return val;
    }

    //- The conjugate of a complex number
    template<class T>
    typename std::enable_if
    <
        std::is_same<complex, T>::value,
        complex
    >::type conj(const T& val)
    {
        return val.conjugate();
    }

} // End namespace Detail


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for complex
template<> struct is_contiguous<complex> : std::true_type {};

//- Contiguous scalar data for complex
template<> struct is_contiguous_scalar<complex> : std::true_type {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, complex& c);
Ostream& operator<<(Ostream& os, const complex& c);

//- Complex conjugate
inline complex operator~(const complex& c);


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return string representation of complex
word name(const complex& c);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
