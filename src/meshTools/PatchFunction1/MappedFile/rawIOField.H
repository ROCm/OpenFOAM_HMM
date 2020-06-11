/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::rawIOField

Description
    Like IOField but falls back to raw IFstream if no header found.
    Optionally reads average value. For use in MappedFile container.

SourceFiles
    rawIOField.C

\*---------------------------------------------------------------------------*/

#ifndef rawIOField_H
#define rawIOField_H

#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rawIOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rawIOField
:
    public regIOobject,
    public Field<Type>
{
    // Private Data

        //- The average of the field (Zero if not used)
        Type average_;


public:

    TypeName("rawField");


    // Constructors

        //- Default copy construct
        rawIOField(const rawIOField&) = default;

        //- Construct from IOobject
        explicit rawIOField(const IOobject& io, const bool readAverage);


    //- Destructor
    virtual ~rawIOField() = default;


    // Member Functions

        const Type& average() const
        {
            return average_;
        }

        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy or move assignment of entries
        using Field<Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rawIOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
