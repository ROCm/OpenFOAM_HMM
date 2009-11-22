/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::Distribution

Description
    Accumulating histogram of component values.
    Specified bin resolution automatic generation of bins.

SourceFiles
    DistributionI.H
    Distribution.C
    DistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef Distribution_H
#define Distribution_H

#include "Map.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Distribution;

template<class Type>
Ostream& operator<<(Ostream&, const Distribution<Type>&);

/*---------------------------------------------------------------------------*\
                         Class Distribution Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Distribution
:
    public List< Map<label> >
{
    // Private data

        Type binWidth_;


    // Private Member Functions

        // //- Disallow default bitwise copy construct
        // Distribution(const Distribution<Type>&);

        // //- Disallow default bitwise assignment
        // void operator=(const Distribution<Type>&);


public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;


    // Static data members

        // static void write
        // (
        //     const fileName& file,
        //     const List<Pair<scalar> >& pairs
        // );


    // Constructors

        //- Construct null
        Distribution();

        //- Construct from separate binWidth for each component
        Distribution(const Type& binWidth);

        // //- Construct from single binWidth for each component
        // Distribution(const cmptType& binWidth);

        //- Construct as copy
        Distribution(const Distribution& d);


    //- Destructor
        ~Distribution();

    // Member Functions

        // Type totalEntries() const;

        // Type mean() const;

        // Type median();

        void add(const Type& valueToAdd);

        void insertMissingKeys();

        // List<Pair<scalar> > normalised();

        // List<Pair<scalar> > normalisedMinusMean();

        // List<Pair<scalar> > normalisedShifted(scalar shift Value);

        List< List < Pair<scalar> > > raw();

    // Access

        inline const Type& binWidth() const;

    // Member Operators

        void operator=(const Distribution<Type>&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const Distribution<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Distribution.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
