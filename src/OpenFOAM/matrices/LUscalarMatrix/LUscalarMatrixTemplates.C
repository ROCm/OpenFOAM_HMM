/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "LUscalarMatrix.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::LUscalarMatrix::solve
(
    List<Type>& x,
    const UList<Type>& source
) const
{
    // If x and source are different initialize x = source
    if (&x != &source)
    {
        x = source;
    }

    if (Pstream::parRun())
    {
        List<Type> X; // scratch space (on master)

        if (Pstream::master(comm_))
        {
            X.resize(m());

            SubList<Type>(X, x.size()) = x;

            for (const int slave : Pstream::subProcs(comm_))
            {
                IPstream::read
                (
                    Pstream::commsTypes::scheduled,
                    slave,
                    reinterpret_cast<char*>
                    (
                        &(X[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave+1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            OPstream::write
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                x.cdata_bytes(),
                x.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }

        if (Pstream::master(comm_))
        {
            LUBacksubstitute(*this, pivotIndices_, X);

            x = SubList<Type>(X, x.size());

            for (const int slave : Pstream::subProcs(comm_))
            {
                OPstream::write
                (
                    Pstream::commsTypes::scheduled,
                    slave,
                    reinterpret_cast<const char*>
                    (
                        &(X[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave+1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            IPstream::read
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                x.data_bytes(),
                x.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }
    }
    else
    {
        LUBacksubstitute(*this, pivotIndices_, x);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LUscalarMatrix::solve
(
    const UList<Type>& source
) const
{
    auto tx(tmp<Field<Type>>::New(m()));

    solve(tx.ref(), source);

    return tx;
}


// ************************************************************************* //
