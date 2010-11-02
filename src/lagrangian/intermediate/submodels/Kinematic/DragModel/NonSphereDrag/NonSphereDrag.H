/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

Class
    Foam::NonSphereDrag

Description
    Drag model for non-spherical particles.

    Takes the form of

        24.0/Re*(1.0 + a_*pow(Re, b_)) + Re*c_/(Re + d_);

    Where a(phi), b(phi), c(phi) and d(phi) are model coefficients, with phi
    defined as:

              area of sphere with same volume as particle
        phi = -------------------------------------------
                       actual particle area

    Equation used is Eqn (11) of reference below - good to within 2 to 4 % of
    RMS values from experiment.

    H and L also give a simplified model with greater error compared to
    results from experiment - Eqn 12 - but since phi is presumed
    constant, it offers little benefit.

    Reference:
    @verbatim
        "Drag coefficient and terminal velocity of spherical and nonspherical
        particles"
        A. Haider and O. Levenspiel,
        Powder Technology
        Volume 58, Issue 1, May 1989, Pages 63-70
    @endverbatim


\*---------------------------------------------------------------------------*/

#ifndef NonSphereDrag_H
#define NonSphereDrag_H

#include "DragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class NonSphereDrag Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NonSphereDrag
:
    public DragModel<CloudType>
{
protected:

    // Protected Data

        //- Flag to indicate `Simple model'
        bool simpleModel_;

        //- Ratio of surface of sphere having same volume as particle to
        //  actual surface area of particle (0 < phi <= 1)
        scalar phi_;


        // Model coefficients

            scalar a_;

            scalar b_;

            scalar c_;

            scalar d_;


public:

    //- Runtime type information
    TypeName("NonSphereDrag");


    // Constructors

        //- Construct from dictionary
        NonSphereDrag(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NonSphereDrag(const NonSphereDrag<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DragModel<CloudType> > clone() const
        {
            return autoPtr<DragModel<CloudType> >
            (
                new NonSphereDrag<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NonSphereDrag();


    // Member Functions

        //- Return drag coefficient
        scalar Cd(const scalar Re) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NonSphereDrag.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
