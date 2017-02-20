/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::FlatOutput

Description
    Simple output adapter for list output on a single line.
    The backend type must support a two-argument \c writeList() method.

\*---------------------------------------------------------------------------*/

#ifndef FlatOutput_H
#define FlatOutput_H

#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Container> class FlatOutput;
template<class Container>
Ostream& operator<<(Ostream& os, const FlatOutput<Container>& obj);


/*---------------------------------------------------------------------------*\
                         Class FlatOutput Declaration
\*---------------------------------------------------------------------------*/

template<class Container>
class FlatOutput
{
    const Container& ref_;
    const label len_;

public:

    //- Construct from components
    inline FlatOutput(const Container& obj, label len)
    :
        ref_(obj),
        len_(len)
    {}

    //- Ostream operator
    inline friend Ostream& operator<<
    (
        Ostream& os,
        const FlatOutput<Container>& wrapped
    )
    {
        return wrapped.ref_.writeList(os, wrapped.len_);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global flatOutput function
template<class Container>
FlatOutput<Container> flatOutput(const Container& obj, label len=0)
{
    return FlatOutput<Container>(obj, len);
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
