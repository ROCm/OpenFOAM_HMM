/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::phasePairKey

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef phasePairKey_H
#define phasePairKey_H

#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class phasePairKey;

bool operator==(const phasePairKey& a, const phasePairKey& b);
bool operator!=(const phasePairKey& a, const phasePairKey& b);

Istream& operator>>(Istream& is, phasePairKey& key);
Ostream& operator<<(Ostream& os, const phasePairKey& key);

/*---------------------------------------------------------------------------*\
                         Class phasePairKey Declaration
\*---------------------------------------------------------------------------*/

class phasePairKey
:
    public Pair<word>
{
public:

        class hash
        :
            public Hash<phasePairKey>
        {
        public:

            // Constructors

                // Construct null
                hash();


            // Member operators

                // Generate a hash from a phase pair key
                label operator()(const phasePairKey& key) const;
        };


private:

    // Private data

        //- Flag to indicate whether ordering is important
        bool ordered_;

public:


    // Constructors

        //- Construct null
        phasePairKey();

        //- Construct from names and the ordering flag
        phasePairKey
        (
            const word& name1,
            const word& name2,
            const bool ordered
        );


    //- Destructor
    virtual ~phasePairKey();


    // Access

        //- Return the ordered flag
        bool ordered() const;


    // Friend Operators

        //- Test for equality
        friend bool operator==(const phasePairKey& a, const phasePairKey& b);

        //- Test for inequality
        friend bool operator!=(const phasePairKey& a, const phasePairKey& b);

        //- Read from Istream
        friend Istream& operator>>(Istream& is, phasePairKey& key);

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const phasePairKey& key);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
