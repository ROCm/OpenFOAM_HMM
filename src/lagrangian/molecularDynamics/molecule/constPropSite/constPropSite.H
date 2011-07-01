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


Description

SourceFiles
    constPropSiteI.H
    constPropSite.C

\*---------------------------------------------------------------------------*/

#ifndef constPropSite_H
#define constPropSite_H

#include "vector.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class constPropSite;
Istream& operator>>(Istream&, constPropSite&);
Ostream& operator<<(Ostream&, const constPropSite&);


/*---------------------------------------------------------------------------*\
                         Class constPropSite Declaration
\*---------------------------------------------------------------------------*/

class constPropSite
{
    // Private data

        //-
        vector siteReferencePosition_;

        //-
        scalar siteMass_;

        //-
        scalar siteCharge_;

        //-
        label siteId_;

        //-
        word name_;

        //-
        bool pairPotentialSite_;

        //-
        bool electrostaticSite_;


public:

    // Constructors

        //- Construct null
        constPropSite();

        //- Construct from components
        constPropSite
        (
            const vector& siteReferencePosition,
            const scalar& siteMass,
            const scalar& siteCharge,
            const label& siteId,
            const word& name,
            const bool& pairPotentialSite,
            const bool& electrostaticSite
        );


    // Selectors

    // Destructor

        ~constPropSite();

    // Member Functions

        // Access

            //-
            inline const vector& siteReferencePosition() const;

            //-
            inline vector& siteReferencePosition();

            //-
            inline const scalar& siteMass() const;

            //-
            inline scalar& siteMass();

            //-
            inline const scalar& siteCharge() const;

            //-
            inline scalar& siteCharge();

            //-
            inline const label& siteId() const;

            //-
            inline label& siteId();

            //-
            inline const word& name() const;

            //-
            inline word& name();

            //-
            inline const bool& pairPotentialSite() const;

            //-
            inline bool& pairPotentialSite();

            //-
            inline const bool& electrostaticSite() const;

            //-
            inline bool& electrostaticSite();


    // Member Operators

        inline bool operator==(const constPropSite&) const;
        inline bool operator!=(const constPropSite&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, constPropSite&);
        friend Ostream& operator<<(Ostream&, const constPropSite&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constPropSiteI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
