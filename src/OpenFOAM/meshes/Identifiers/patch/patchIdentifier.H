/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::patchIdentifier

Description
    Identifies a patch by name, patch index and physical type

SourceFiles
    patchIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef patchIdentifier_H
#define patchIdentifier_H

#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                       Class patchIdentifier Declaration
\*---------------------------------------------------------------------------*/

class patchIdentifier
{
    // Private Data

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label index_;

        //- Optional physical type
        mutable word physicalType_;

        //- Optional groups to which the patch belongs
        wordList inGroups_;

public:

    // Static Member Functions

        //- Default patch name: "patch" or "patchN"
        static word defaultName(const label n = -1)
        {
            return
            (
                n < 0
              ? word("patch", false)
              : word("patch" + std::to_string(n), false)
            );
        }


    // Generated Methods

        //- Copy construct
        patchIdentifier(const patchIdentifier&) = default;

        //- Copy assignment
        patchIdentifier& operator=(const patchIdentifier&) = default;

        //- Destructor
        virtual ~patchIdentifier() = default;


    // Constructors

        //- Default construct, with index zero
        patchIdentifier();

        //- Construct from mandatory components
        patchIdentifier(const word& name, const label index);

        //- Construct from components
        patchIdentifier
        (
            const word& name,
            const label index,
            const word& physicalType,
            const wordList& inGroups = wordList()
        );

        //- Construct from dictionary
        patchIdentifier
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Copy construct, resetting the index
        patchIdentifier
        (
            const patchIdentifier& p,
            const label index
        );


    // Member Functions

        //- The patch name
        const word& name() const
        {
            return name_;
        }

        //- Modifiable patch name
        word& name()
        {
            return name_;
        }

        //- The (optional) physical type of the patch
        const word& physicalType() const
        {
            return physicalType_;
        }

        //- Modifiable (optional) physical type of the patch
        word& physicalType()
        {
            return physicalType_;
        }

        //- The index of this patch in the boundaryMesh
        label index() const
        {
            return index_;
        }

        //- Modifiable index of this patch in the boundaryMesh
        label& index()
        {
            return index_;
        }

        //- The (optional) groups that the patch belongs to
        const wordList& inGroups() const
        {
            return inGroups_;
        }

        //- Modifiable (optional) groups that the patch belongs to
        wordList& inGroups()
        {
            return inGroups_;
        }

        //- True if the patch is in named group
        bool inGroup(const word& name) const
        {
            return inGroups_.found(name);
        }

        //- Write (physicalType, inGroups) dictionary entries
        //- (without surrounding braces)
        void write(Ostream& os) const;
};


// Global Operators

//- Write (physicalType, inGroups) dictionary entries
//- (without surrounding braces)
Ostream& operator<<(Ostream& os, const patchIdentifier& p);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
