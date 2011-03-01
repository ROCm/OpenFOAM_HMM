/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::codeProperties

Description
    IOdictionary with an internal flag to explicitly track when a file changed.

SourceFiles
    codeProperties.C

\*---------------------------------------------------------------------------*/

#ifndef codeProperties_H
#define codeProperties_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class codeProperties Declaration
\*---------------------------------------------------------------------------*/

class codeProperties
:
    public IOdictionary
{
    // Private data

        //- Monitor file change
        mutable bool modified_;

public:

    // Declare name of the class and its debug switch
    ClassName("codeDict");

    // Constructors

        //- Construct from IOobject
        codeProperties(const IOobject&);


    // Member Functions

        //- A sticky version of regIOobject::modified()
        //  Must explicitly clear with setUnmodified()
        virtual bool modified() const;

        //- Make modification unsticky
        void setUnmodified() const;

        //- Read the dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
