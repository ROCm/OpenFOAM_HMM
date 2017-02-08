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
    foamVtkAppendRawFormatter

Description
    Appended raw binary output.

SourceFiles
    foamVtkAppendRawFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkAppendRawFormatter_H
#define foamVtkAppendRawFormatter_H

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class foamVtkAppendRawFormatter Declaration
\*---------------------------------------------------------------------------*/

class foamVtkAppendRawFormatter
:
    public foamVtkFormatter
{
    // Private Data Members

        static const char* name_;
        static const char* encoding_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamVtkAppendRawFormatter(const foamVtkAppendRawFormatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const foamVtkAppendRawFormatter&) = delete;


protected:

    // Protected Member Functions

        //- Write
        void write(const char* s, std::streamsize n);


public:

    // Constructors

        //- Construct and attach to an output stream
        foamVtkAppendRawFormatter(std::ostream& os);


    //- Destructor
    virtual ~foamVtkAppendRawFormatter();


    // Member Functions

        //- Output name for XML type ("append")
        virtual const char* name() const;

        //- Output name for append encoding type ("raw")
        virtual const char* encoding() const;


        //- Write leading size for binary output
        virtual void writeSize(const uint64_t nBytes);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- A no-op for this format
        virtual void flush();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
