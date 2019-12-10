/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider <bgschaid@hfd-research.com>
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::expressions::exprDriverWriter

Description
    Registered input/output for an expressions::fvExprDriver

SourceFiles
    exprDriverWriter.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_exprDriverWriter_H
#define expressions_exprDriverWriter_H

#include "fvExprDriver.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprDriverWriter Declaration
\*---------------------------------------------------------------------------*/

class exprDriverWriter
:
    public regIOobject
{
    // Private Data

        //- The driver to read/write
        fvExprDriver& driver_;


    // Private Member Functions

        //- No null constructor
        exprDriverWriter() = delete;

        //- No copy construct
        exprDriverWriter(const exprDriverWriter&) = delete;

        //- No copy assignment
        void operator=(const exprDriverWriter&) = delete;


public:

    //- Runtime type information
    TypeNameNoDebug("exprDriverWriter");


    // Constructors

        //- Construct for named driver
        exprDriverWriter(const word& name, fvExprDriver& driver);


    //- Destructor
    virtual ~exprDriverWriter() = default;


    // Member Functions

        virtual bool readData(Istream& is);
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
