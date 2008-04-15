/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::PatchProperties

Description

SourceFiles
    PatchProperties.C

\*---------------------------------------------------------------------------*/

#ifndef PatchProperties_H
#define PatchProperties_H

// Foam header files.
#include "dictionary.H"
#include "Dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class  PatchProperties Declaration
\*---------------------------------------------------------------------------*/

class PatchProperties
:
    public Foam::Dictionary<PatchProperties>::link
{
    // Private data

        Foam::word patchName_;
        Foam::word patchType_;
        Foam::word physicalType_;
        Foam::label startFace_;
        Foam::label nFaces_;

        bool modified_;

public:

    // Constructors

        PatchProperties(const Foam::word& patchName);


    // Destructor

        ~PatchProperties();


    // Member Functions

        const Foam::word& patchName() const;
        void patchName(const Foam::word&);

        //- Keyword used by Dictionary
        const word& keyword() const
        {
            return patchName_;
        }

        const Foam::word& patchType() const;
        void patchType(const word&);

        const Foam::word& physicalType() const;
        void physicalType(const Foam::word&);

        Foam::label startFace() const;
        void startFace(Foam::label);

        Foam::label nFaces() const;
        void nFaces(Foam::label);

        bool modified() const;
        void modified(bool);

        // Persistence

            void load(const Foam::dictionary& patchDict);
            void save(DictionaryWriter& dictWriter);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
