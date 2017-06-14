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
    Foam::vtk::appendBase64Formatter

Description
    Appended base-64 encoded binary output.
    Uses an output filter layer to write base-64 encoded content.

SourceFiles
    foamVtkAppendBase64Formatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkAppendBase64Formatter_H
#define foamVtkAppendBase64Formatter_H

#include "foamVtkBase64Layer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class appendBase64Formatter Declaration
\*---------------------------------------------------------------------------*/

class appendBase64Formatter
:
    public foamVtkBase64Layer
{
    // Private Data Members

        static const char* name_;
        static const outputOptions opts_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        appendBase64Formatter(const appendBase64Formatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const appendBase64Formatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream
        appendBase64Formatter(std::ostream& os);


    //- Destructor
    virtual ~appendBase64Formatter();


    // Member Functions

        //- The output is APPEND_BASE64.
        virtual const outputOptions& opts() const;

        //- Output name for XML type ("append")
        virtual const char* name() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
