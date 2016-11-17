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
    foamVtkAsciiFormatter

Description
    Inline ASCII binary output.
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

/*---------------------------------------------------------------------------*\
                    Class foamVtkAsciiFormatter Declaration
\*---------------------------------------------------------------------------*/

class foamVtkAsciiFormatter
:
    public foamVtkFormatter
{
    // Private Data Members

        static const char* name_;

        //- Track the current output position
        unsigned short pos_;


    // Private Member Functions

        //- Advance to next position, adding space or newline as required
        inline void next();


        //- Disallow default bitwise copy construct
        foamVtkAsciiFormatter(const foamVtkAsciiFormatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const foamVtkAsciiFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream, use default precision
        foamVtkAsciiFormatter(std::ostream&);

        //- Construct and attach to an output stream, use specified precision
        foamVtkAsciiFormatter(std::ostream&, unsigned precision);


    //- Destructor
    virtual ~foamVtkAsciiFormatter();


    // Member Functions

        //- Name for the XML output type ("ascii")
        //  The legacy output type is an uppercase version of this.
        virtual const char* name() const;

        //- Name for the XML append encoding - unused.
        //  Currently simply "ASCII", but this should not be relied upon.
        virtual const char* encoding() const;


        //- Write leading size - this is a no-op for ascii output
        virtual void writeSize(const uint64_t);

        virtual void write(const uint8_t);
        virtual void write(const label);
        virtual void write(const float);
        virtual void write(const double);
        virtual void flush();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
