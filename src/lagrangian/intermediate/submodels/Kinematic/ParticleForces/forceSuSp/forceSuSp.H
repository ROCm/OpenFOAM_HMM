/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::forceSuSp

Description
    Helper container for force Su and Sp terms

SourceFiles
    forceSuSpI.H

\*---------------------------------------------------------------------------*/

#ifndef forceSuSp_H
#define forceSuSp_H

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions
class forceSuSp;
inline forceSuSp operator+(const forceSuSp& susp1, const forceSuSp& susp2);
inline forceSuSp operator*(const forceSuSp& susp, const scalar s);

/*---------------------------------------------------------------------------*\
                          Class forceSuSp Declaration
\*---------------------------------------------------------------------------*/

class forceSuSp
:
    public Tuple2<vector, scalar>
{

public:

    //- Component labeling enumeration
    enum components
    {
        SU = 0,
        SP = 1
    };


    // Constructors

        //- Construct null
        inline forceSuSp();

        //- Construct given VectorSpace
        inline forceSuSp(const Tuple2<vector, scalar>& susp);

        //- Construct given two components
        inline forceSuSp(const vector& Su, const scalar Sp);

        //- Construct from Istream
        inline forceSuSp(Istream&);


    // Member Functions

        // Access

            //- Return const access to the explicit contribution
            inline const vector& Su() const;

            //- Return const access to the implicit contribution
            inline scalar Sp() const;


        // Edit

            //- Return reference to the explicit contribution
            inline vector& Su();

            //- Return reference to the explicit contribution
            inline scalar& Sp();


        // Operators

            //- Assignment
            inline void operator=(const forceSuSp& susp);

            //- Addition
            inline void operator+=(const forceSuSp& susp);

            //- Subtraction
            inline void operator-=(const forceSuSp& susp);


        // Friend operators

            //- Addition
            friend inline forceSuSp operator*
            (
                const forceSuSp& susp1,
                const forceSuSp& susp2
            );

            //- Multiplication
            friend inline forceSuSp operator*
            (
                const scalar s,
                const forceSuSp& susp
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "forceSuSpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
