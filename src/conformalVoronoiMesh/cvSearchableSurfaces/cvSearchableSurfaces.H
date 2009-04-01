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
    Foam::cvSearchableSurfaces

Description
    The surface geometry to be meshed with a conformalVoronoiMesh, adding
    inside/outside queries and surfaceFeatures

SourceFiles
    cvSearchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef cvSearchableSurfaces_H
#define cvSearchableSurfaces_H

#include "searchableSurfaces.H"
#include "surfaceFeatures.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cvSearchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class cvSearchableSurfaces
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cvSearchableSurfaces(const cvSearchableSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const cvSearchableSurfaces&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        cvSearchableSurfaces();

        //- Construct from components
        cvSearchableSurfaces(const dataType& data);

        //- Construct from Istream
        cvSearchableSurfaces(Istream&);

        //- Construct as copy
        cvSearchableSurfaces(const cvSearchableSurfaces&);


    // Selectors

        //- Select null constructed
        static autoPtr<cvSearchableSurfaces> New();


    //- Destructor
    ~cvSearchableSurfaces();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const cvSearchableSurfaces&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, cvSearchableSurfaces&);
        friend Ostream& operator<<(Ostream&, const cvSearchableSurfaces&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cvSearchableSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
