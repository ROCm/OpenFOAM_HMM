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
    foamVtkBase64Formatter

Description
    Inline base-64 encoded binary output.
    Uses an output filter layer to write base-64 encoded content.

\*---------------------------------------------------------------------------*/

#ifndef foamVtkBase64Formatter_H
#define foamVtkBase64Formatter_H

#include "foamVtkBase64Layer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class foamVtkBase64Formatter Declaration
\*---------------------------------------------------------------------------*/

class foamVtkBase64Formatter
:
    public foamVtkBase64Layer
{
    // Private Data Members

        static const char* name_;
        static const char* encoding_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamVtkBase64Formatter(const foamVtkBase64Formatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const foamVtkBase64Formatter&) = delete;

public:

    // Constructors

        //- Construct and attach to an output stream
        foamVtkBase64Formatter(std::ostream& os);


    //- Destructor
    virtual ~foamVtkBase64Formatter();


    // Member Functions

        //- Name for the XML output type ("binary")
        //  The lowercase version of the Legacy output type.
        virtual const char* name() const;


        //- End the encoding sequence (padding the final characters with '=')
        //  and write a newline to the output if any encoding was done.
        virtual void flush();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
