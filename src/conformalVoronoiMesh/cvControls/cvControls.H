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
    Foam::cvControls

Description
    Controls for the conformalVoronoiMesh mesh generator

SourceFiles
    cvControls.C

\*---------------------------------------------------------------------------*/

#ifndef cvControls_H
#define cvControls_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;


/*---------------------------------------------------------------------------*\
                          Class cvControls Declaration
\*---------------------------------------------------------------------------*/

class cvControls
:
{
    // Private data

        // Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cvControls(const cvControls&);

        //- Disallow default bitwise assignment
        void operator=(const cvControls&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        cvControls();

        //- Construct from components
        cvControls(const dataType& data);

        //- Construct from Istream
        cvControls(Istream&);

        //- Construct as copy
        cvControls(const cvControls&);


    // Selectors

        //- Select null constructed
        static autoPtr<cvControls> New();


    //- Destructor
    ~cvControls();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const cvControls&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, cvControls&);
        friend Ostream& operator<<(Ostream&, const cvControls&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cvControlsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
