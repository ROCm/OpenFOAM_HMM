/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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
inline complex operator+(const complex&, const complex&);
inline complex operator-(const complex&);
inline complex operator-(const complex&, const complex&);
inline complex operator*(const complex&, const complex&);
inline complex operator/(const complex&, const complex&);
inline complex operator*(const scalar, const complex&);
inline complex operator*(const complex&, const scalar);
inline complex operator/(const complex&, const scalar);
inline complex operator/(const scalar, const complex&);

Istream& operator>>(Istream& is, complex& c);
Ostream& operator<<(Ostream& os, const complex& c);


/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/

class complex
{
    // Private Data

        //- Real and imaginary parts
        scalar re, im;


public:

    //- Component type
    typedef complex cmptType;


    // Static Data Members

        //- The type name is "complex"
        static const char* const typeName;

        //- A complex zero (0,0)
        static const complex zero;

        //- A complex one (1,0)
        static const complex one;


    // Constructors

        //- Construct null as zero-initialized
        inline constexpr complex() noexcept;

        //- Default copy constructor
        complex(const complex&) = default;

        //- Construct zero-initialized from zero class
        inline constexpr complex(const Foam::zero) noexcept;

        //- Construct from real component
        inline explicit constexpr complex(const scalar r) noexcept;

        //- Construct from real and imaginary parts
        inline constexpr complex(const scalar r, const scalar i) noexcept;

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

        //- Copy assignment
        inline void operator=(const complex& c);

        //- Assign zero
        inline void operator=(const Foam::zero);

        inline void operator+=(const complex& c);
        inline void operator-=(const complex& c);
        inline void operator*=(const complex& c);
        inline void operator/=(const complex& c);

        inline void operator=(const scalar s);
        inline void operator+=(const scalar s);
        inline void operator-=(const scalar s);
        inline void operator*=(const scalar s);
        inline void operator/=(const scalar s);

        inline bool operator==(const complex& c) const;
        inline bool operator!=(const complex& c) const;


    // Friend Functions

        friend scalar magSqr(const complex& c);
        friend scalar mag(const complex& c);
        friend complex sqr(const complex& c);
        friend const complex& min(const complex&, const complex&);
        friend const complex& max(const complex&, const complex&);

        friend complex limit(const complex&, const complex&);

        friend const complex& sum(const complex&);


    // Friend Operators

        friend complex operator+(const complex& c1, const complex& c2);
        friend complex operator-(const complex& c);
        friend complex operator-(const complex& c1, const complex& c2);
        friend complex operator*(const complex& c1, const complex& c2);
        friend complex operator/(const complex& c1, const complex& c2);

        friend complex operator*(const scalar s, const complex& c);
        friend complex operator*(const complex& c, const scalar s);
        friend complex operator/(const complex& c, const scalar s);
        friend complex operator/(const scalar s, const complex& c);


    // IOstream Operators

        friend Istream& operator>>(Istream& is, complex& c);
        friend Ostream& operator<<(Ostream& os, const complex& c);
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //


//- Return string representation of complex
word name(const complex& c);

//- Data associated with complex type are contiguous
template<>
inline bool contiguous<complex>() {return true;}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Complex conjugate
inline complex operator~(const complex& c);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
