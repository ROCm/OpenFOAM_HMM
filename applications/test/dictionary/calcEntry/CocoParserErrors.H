/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::CocoParserErrors

Description
    Templated class to shadow the error handling for Coco/R parsers

\*---------------------------------------------------------------------------*/

#ifndef CocoParserErrors_H
#define CocoParserErrors_H

#include "error.H"
#include "wchar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CocoParserErrors Declaration
\*---------------------------------------------------------------------------*/

template<class BaseClass>
class CocoParserErrors
:
    public BaseClass
{
    // Private data

        //- The name issued in warnings and errors
        word   name_;

public:

    // Constructors

        //- Construct with given name
        CocoParserErrors(const word& name)
        :
            BaseClass(),
            name_(name)
        {}


    //- Destructor
    virtual ~CocoParserErrors()
    {}


    // Member functions

        //- Return the name issued for warnings
        virtual const word& name() const
        {
            return name_;
        }

        //- Return the name issued for warnings
        virtual word& name()
        {
            return name_;
        }


    // Error Handling

        //- Handle a general warning 'msg'
        virtual void Warning(const wchar_t* msg)
        {
            WarningIn(name_)
                << msg << endl;
        }

        //- Handle a general warning 'msg'
        virtual void Warning(int line, int col, const wchar_t* msg)
        {
            WarningIn(name_)
                <<"line " << line << " col " << col << ": "
                << msg << endl;
        }

        //- Handle general error 'msg' (eg, a semantic error)
        virtual void Error(int line, int col, const wchar_t* msg)
        {
            FatalErrorIn(name_)
                << "line " << line << " col " << col <<": " << msg << endl
                << exit(FatalError);
        }

        //- Handle general error 'msg' (eg, a semantic error)
        virtual void Error(const wchar_t* msg)
        {
            FatalErrorIn(name_)
                << msg << endl
                << exit(FatalError);
        }

        //- Handle a general exception 'msg'
        virtual void Exception(const wchar_t* msg)
        {
            this->Error(msg);
        }

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
