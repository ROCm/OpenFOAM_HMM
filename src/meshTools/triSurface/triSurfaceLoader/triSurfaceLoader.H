/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::triSurfaceLoader

Description
    Convenience class for loading single or multiple surface files
    from the constant/triSurface (or other) directory.

    Surfaces selection based on word, wordRe, wordReList.
    If multiple surfaces are selected, they are concatenated into a
    single surface with offset faces,points,regions.

SourceFiles
    triSurfaceLoader.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceLoader_H
#define triSurfaceLoader_H

#include "triSurface.H"
#include "wordReList.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurfaceLoader;
class Time;

/*---------------------------------------------------------------------------*\
                      Class triSurfaceLoader Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceLoader
{
    // Private data

        //- The directory to load from (eg, case/constant/triSurface)
        fileName directory_;

        //- All available files
        hashedWordList available_;

        //- Selected files
        hashedWordList selected_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        triSurfaceLoader(const triSurfaceLoader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceLoader&) = delete;

public:

    // Constructors

        //- Construct with directory name
        triSurfaceLoader(const fileName& directory);

        //- Construct with time. Selects "constant/triSurface" directory.
        triSurfaceLoader(const Time& runTime);


    //- Destructor
    ~triSurfaceLoader();


    // Member Functions

      // Access

        //- The directory being used
        inline const fileName& directory() const
        {
            return directory_;
        }

        //- The list of available files
        inline const hashedWordList& available() const
        {
            return available_;
        }

        //- The list of selected files
        inline const hashedWordList& selected() const
        {
            return selected_;
        }


      // Edit

        //- Read directory and populate the 'available' files.
        //  Automatically called on construction.
        label readDir();

        //- Populates 'selected' with all available files.
        label selectAll();

        //- Populates 'selected' with a subset of the available files.
        label select(const word& name);

        //- Populates 'selected' with a subset of the available files.
        label select(const wordRe& mat);

        //- Populates 'selected' with a subset of the available files.
        label select(const wordReList& matcher);

        //- Load a single file, or load and combine multiple selected files
        autoPtr<triSurface> load() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
