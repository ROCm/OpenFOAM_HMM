/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Application
    Test-parallel-broadcast

Description
    Test for various broadcast routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "vector.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

using namespace Foam;


// This is what our new scatter will look like inside
template<class Type>
void testBroadcastStream
(
    Type& value,
    const label comm = UPstream::worldComm
)
{
    Info<< nl << "is_contiguous:" << is_contiguous<Type>::value << endl;
    Pout<< "pre-broadcast: " << value << endl;

    if (is_contiguous<Type>::value)
    {
         UPstream::broadcast
         (
             reinterpret_cast<char*>(&value),
             sizeof(Type),
             comm,
             UPstream::masterNo()
         );
    }
    else
    {
        if (UPstream::master())
        {
            OPBstream toAll(UPstream::masterNo(), comm);
            toAll << value;
        }
        else
        {
            IPBstream fromMaster(UPstream::masterNo(), comm);
            fromMaster >> value;
        }
    }

    Pout<< "post-broadcast: " << value << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    #include "setRootCase.H"
    #include "createTime.H"

    {
        label value = -1;
        if (Pstream::master())
        {
            value = UPstream::nProcs();
        }
        testBroadcastStream(value);
    }

    {
        labelList values;
        if (Pstream::master())
        {
            values = identity(UPstream::nProcs());
        }
        testBroadcastStream(values);
    }

    {
        wordList values;
        if (Pstream::master())
        {
            values.resize(UPstream::nProcs());
            forAll(values, i)
            {
                values[i] = "value_" + Foam::name(i);
            }
        }
        testBroadcastStream(values);
    }

    {
        vector values(vector::uniform(-1));
        if (Pstream::master())
        {
            values = vector(1,2,3);
        }
        testBroadcastStream(values);
    }

    {
        FixedList<vector, 3> values(vector::uniform(-1));
        if (Pstream::master())
        {
            values = vector(1,2,3);
        }
        testBroadcastStream(values);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
