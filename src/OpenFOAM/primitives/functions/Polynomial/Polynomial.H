/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::Polynomial

Description
    Polynomial templated on size (order):

        poly = sum(coeff_[i]*x^i)

    where 0 <= i <= n

    - integer powers, starting at zero
    - evaluate(x) to evaluate the poly for a given value
    - integrate(x1, x2) between two scalar values
    - integrate() to return a new, intergated coeff polynomial
      - increases the size (order)

SourceFiles
    Polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef Polynomial_H
#define Polynomial_H

#include "word.H"
#include "scalar.H"
#include "Ostream.H"
#include "VectorSpace.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<int PolySize>
class Polynomial;

// Forward declaration of friend functions
template<int PolySize>
Ostream& operator<<
(
    Ostream&,
    const Polynomial<PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<int PolySize>
class Polynomial
:
    public VectorSpace<Polynomial<PolySize>, scalar, PolySize>
{
private:

    // Private data

        //- Polynomial name
        word name_;


public:

    typedef VectorSpace<Polynomial<PolySize>, scalar, PolySize> polyType;

    typedef Polynomial<PolySize+1> intPolyType;

    //- Run-time type information
    TypeName("Polynomial")


    // Constructors

        //- Construct null
        Polynomial();

        //- Construct from name and Istream
        Polynomial(const word& name, Istream& is);

        //- Copy constructor
        Polynomial(const Polynomial& poly);

        //- Copy constructor with name
        Polynomial(const word& name, const Polynomial& poly);


    //- Destructor
    ~Polynomial();


    // Member Functions

        // Access

            //- Return const access to the polynomial name
            const word& name() const;


        // Evaluation

            //- Return polynomial value
            scalar evaluate(const scalar x) const;

            //- Return integrated polynomial coefficients
            //  argument becomes zeroth element (constant of integration)
            intPolyType integrate(const scalar intConstant = 0.0);

            //- Integrate between two values
            scalar integrateLimits(const scalar x1, const scalar x2) const;


    // Member operators

        void operator=(const Polynomial& poly);


    //- Ostream Operator
    friend Ostream& operator<< <PolySize>
    (
        Ostream&,
        const Polynomial&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePolynomial(PolySize)                                             \
                                                                             \
defineTemplateTypeNameAndDebugWithName                                       \
(                                                                            \
    Polynomial<PolySize>,                                                    \
    "Polynomial<"#PolySize">",                                               \
    0                                                                        \
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Polynomial.C"
#   include "PolynomialIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
