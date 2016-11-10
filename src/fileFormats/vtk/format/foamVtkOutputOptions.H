/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    foamVtkOutputOptions

Description
    Encapsulate combinations of output format options.

SourceFiles
    foamVtkOutputOptions.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutputOptions_H
#define foamVtkOutputOptions_H

#include "autoPtr.H"
#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class Ostream;

/*---------------------------------------------------------------------------*\
                    Class foamVtkOutputOptions Declaration
\*---------------------------------------------------------------------------*/

class foamVtkOutputOptions
{
    // Private data

        //- The supported output/format types
        enum foamVtkOptionTypes
        {
            ASCII  = 0x0000,    //!< ASCII formatting for data
            BINARY = 0x0001,    //!< Raw binary formatting for data
            BASE64 = 0x0002,    //!< Base64 encoding for data
            LEGACY = 0x0100,    //!< Legacy vtk file format
            APPEND = 0x0200     //!< XML append format
        };

        //- The output style tuning
        enum foamVtkStyleOptions
        {
            NONE   = 0x0000,    //!< Normal
            HEADER = 0x0001     //!< Emit xml header
        };


        //- The output format type
        unsigned short type_;

        //- The output style tuning
        unsigned short style_;


        //- ASCII write precision
        mutable unsigned precision_;


public:

    // Constructors

        //- Construct null - XML insitu ASCII format with default precision
        foamVtkOutputOptions();


    // Selectors

        //- Return new data formatter based on the writer options
        autoPtr<foamVtkFormatter> newFormatter(std::ostream&) const;


    // Member Functions

    // Access

        //- True if writer uses legacy file format
        inline bool legacy() const;

        //- True if writer uses XML file format (non-legacy)
        inline bool xml() const;

        //- True if output format uses an append mode
        inline bool append() const;

        //- True if output format does not use an append mode
        inline bool insitu() const;

        //- True if output format is ASCII
        inline bool ascii() const;


    // Edit

        //- Toggle ASCII mode on/off.
        //  In append mode, this switches between base64 and raw binary.
        //  In XML mode, this switches between ASCII and base64.
        //  In legacy mode, this switches between ASCII and binary.
        void ascii(bool);

        //- Toggle append mode on/off.
        void append(bool);

        //- Toggle legacy mode on/off.
        void legacy(bool);

        //- Set the write precision to be used for new ASCII formatters
        void precision(unsigned val) const;


    // Other

        //- Report information about the options
        Ostream& info(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkOutputOptionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
