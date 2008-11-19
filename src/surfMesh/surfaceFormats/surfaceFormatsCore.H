/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::fileFormats::surfaceFormatsCore

Description
    A collection of helper functions for reading/writing surface formats.

SourceFiles
    surfaceFormatsCore.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFormatsCore_H
#define surfaceFormatsCore_H

#include "surfPatchIdentifierList.H"
#include "surfGroupList.H"
#include "labelList.H"
#include "Map.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class IFstream;
class Time;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class surfaceFormatsCore Declaration
\*---------------------------------------------------------------------------*/

class surfaceFormatsCore
{
public:

    // Static Data

    //- Return the mesh sub-directory name (usually "meshedSurface")
    static word meshSubDir;

    //- The file extension corresponding to 'native' surface format
    static word nativeExt;

    // Static Member Functions

        //- Check if file extension corresponds to 'native' surface format
        static bool isNative(const word&);

        //- Read non-comment line
        static string getLineNoComment(IFstream&);

        //- Name of UnsortedMeshedSurface directory to use.
        static fileName findMeshInstance(const Time&, const word& subdirName);

        //- Name of UnsortedMeshedSurface directory to use.
        static fileName findMeshName(const Time&, const word& subdirName);

        //- Name of UnsortedMeshedSurface directory to use.
        static fileName findMeshInstance(const Time&);

        //- Name of UnsortedMeshedSurface directory to use.
        static fileName findMeshName(const Time&);

        //- Determine the sort order from the region list.
        //  Returns patch list and sets faceMap to indices within faceLst
        static surfGroupList sortedPatchRegions
        (
            const UList<label>& regionLst,
            const Map<word>& patchNames,
            labelList& faceMap
        );

        static bool checkSupport
        (
            const wordHashSet& available,
            const word& ext,
            const bool verbose,
            const word& functionName
        );

    // Constructors

        //- Construct null
        surfaceFormatsCore();

    // Destructor

        virtual ~surfaceFormatsCore();

    // Member Functions

    // Access

    // Edit

    // Read

    // Write

    // Member operators

    // Ostream Operator
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
