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
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Ostream;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                    Class fileFormats::NASCore Declaration
\*---------------------------------------------------------------------------*/

class NASCore
{
public:

        //- File field formats
        enum fieldFormat
        {
            SHORT,      //<! Short format (field width = 8)
            LONG,       //<! Long format (field width = 16)
            FREE        //<! Free format (comma-separated fields)
        };

        //- Selection names for the NASTRAN file field formats
        static const Enum<fieldFormat> fieldFormatNames;


    // Constructors

        //- Construct null
        NASCore() = default;


    // Public Static Member Functions

        //- Extract numbers from things like "-2.358-8" (same as "-2.358e-8")
        static scalar readNasScalar(const std::string& str);

        //- Extract numbers from things like "-2.358-8" (same as "-2.358e-8")
        // \deprecated use readNasScalar instead (deprecated Sep 2017)
        inline static scalar parseNASCoord(const std::string& str)
        {
            return readNasScalar(str);
        }

        //- A string::substr() to handle fixed-format and free-format NASTRAN.
        //  Returns the substr to the next comma (if found) or the given length
        //
        //  \param str The string to extract from
        //  \param pos On input, the position of the first character of the
        //      substring. On output, advances to the next position to use.
        //  \param len The fixed-format length to use if a comma is not found.
        static std::string nextNasField
        (
            const std::string& str,
            std::string::size_type& pos,
            std::string::size_type len
        );


        //- Set output stream precision and format flags
        static void setPrecision(Ostream& os, const fieldFormat format);


        //- Write initial keyword (eg, 'GRID' or 'GRID*') followed by the
        //- requisite number of spaces for the field-width
        static Ostream& writeKeyword
        (
            Ostream& os,
            const word& keyword,
            const fieldFormat format
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
