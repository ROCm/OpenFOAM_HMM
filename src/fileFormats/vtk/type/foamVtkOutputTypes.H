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

InNamespace
    Foam::foamVtkOutput

Description
    Enumerations and typedefs for VTK output.

SourceFiles
    foamVtkOutputTypes.H

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutputTypes_H
#define foamVtkOutputTypes_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace foamVtkOutput
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Enumerations

    //- The output format type for file contents.
    //  Upper bits for output type, lower bits for the format itself.
    enum class formatType
    {
        /** XML inline ASCII, using the asciiFormatter */
        INLINE_ASCII  = 0,
        /** XML inline base64, using the base64Formatter */
        INLINE_BASE64 = 0x01,
        /** XML append base64, using the appendBase64Formatter */
        APPEND_BASE64 = 0x11,
        /** XML append raw binary, using the appendRawFormatter */
        APPEND_BINARY = 0x12,
        /** Legacy ASCII, using the legacyAsciiFormatter */
        LEGACY_ASCII  = 0x20,
        /** Legacy raw binary, using the legacyRawFormatter */
        LEGACY_BINARY = 0x22,
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace foamVtkOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
