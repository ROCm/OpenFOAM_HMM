/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::vtk::asciiFormatter

Description
    Inline ASCII output.
    Adds spaces between entries and a newline every 6 items
    (for consistency with what VTK itself outputs).

SourceFiles
    foamVtkAsciiFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkAsciiFormatter_H
#define foamVtkAsciiFormatter_H

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class asciiFormatter Declaration
\*---------------------------------------------------------------------------*/

class asciiFormatter
:
    public formatter
{
    // Private Data Members

        static const char* name_;
        static const outputOptions opts_;

        //- Track the current output position
        unsigned short pos_;


    // Private Member Functions

        //- Advance to next position, adding space or newline as needed
        inline void next();

        //- Finish an output line, adding newline as needed
        inline void done();


        //- Disallow default bitwise copy construct
        asciiFormatter(const asciiFormatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const asciiFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream, use default precision
        asciiFormatter(std::ostream& os);

        //- Construct and attach to an output stream, use specified precision
        asciiFormatter(std::ostream& os, unsigned precision);


    //- Destructor
    virtual ~asciiFormatter();


    // Member Functions

        //- The output is INLINE_ASCII.
        virtual const outputOptions& opts() const;

        //- Name for the XML output type ("ascii")
        virtual const char* name() const;

        //- Name for the XML append encoding - unused.
        //  Currently identical to name(), but do not rely on this.
        virtual const char* encoding() const;


        //- Write leading size - this is a no-op for ascii output
        virtual void writeSize(const uint64_t ignored);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- Write a newline if needed to finish a line of output.
        virtual void flush();

        //- The encoded length for ascii output is not applicable.
        virtual std::size_t encodedLength(std::size_t ignored) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
