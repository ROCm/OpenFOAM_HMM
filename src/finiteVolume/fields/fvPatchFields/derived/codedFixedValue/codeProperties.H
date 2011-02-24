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
    IOdictionary + flag whether file has changed.

SourceFiles
    codeProperties.C

\*---------------------------------------------------------------------------*/

#ifndef codeProperties_H
#define codeProperties_H

#include "MeshObject.H"
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

        //- File change
        mutable bool modified_;

public:

    // Declare name of the class and its debug switch
    ClassName("codeDict");

    // Constructors

        //- Construct from IOobject
        codeProperties(const IOobject&);


    // Member Functions

        bool modified() const
        {
            return modified_;
        }

        void setUnmodified() const
        {
            modified_ = false;
        }

        //- Read the dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
