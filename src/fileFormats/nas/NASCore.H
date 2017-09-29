/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fileFormats::NASCore

Description
    Core routines used when reading/writing NASTRAN files.

SourceFiles
    NASCore.C

\*---------------------------------------------------------------------------*/

#ifndef NASCore_H
#define NASCore_H

#include "scalar.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::NASCore Declaration
\*---------------------------------------------------------------------------*/

class NASCore
{
public:

    // Public Member Functions

        //- Extract numbers from things like "-2.358-8" (same as "-2.358e-8")
        static scalar readNasScalar(const string& str);

        //- Extract numbers from things like "-2.358-8" (same as "-2.358e-8")
        // \deprecated use readNasScalar instead (deprecated Sep 2017)
        inline static scalar parseNASCoord(const string& str)
        {
            return readNasScalar(str);
        }


    // Constructors

        //- Construct null
        NASCore();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
