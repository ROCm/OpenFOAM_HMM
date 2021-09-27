/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "globalIndex.H"
#include "PstreamBuffers.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void Foam::vtk::write
(
    vtk::formatter& fmt,
    const Type& val,
    const label n
)
{
    const direction nCmpt = pTraits<Type>::nComponents;

    for (label i=0; i < n; ++i)
    {
        for (direction cmpt=0; cmpt < nCmpt; ++cmpt)
        {
            fmt.write(component(val, cmpt));
        }
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values
)
{
    for (const Type& val : values)
    {
        vtk::write(fmt, val);
    }
}


template<class Type, unsigned N>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const FixedList<Type, N>& values
)
{
    for (const Type& val : values)
    {
        vtk::write(fmt, val);
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const labelUList& addressing
)
{
    for (const label idx : addressing)
    {
        vtk::write(fmt, values[idx]);
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const bitSet& selected
)
{
    for (const label idx : selected)
    {
        vtk::write(fmt, values[idx]);
    }
}


template<class Type>
void Foam::vtk::writeLists
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const UList<Type>& indirect,
    const labelUList& addressing
)
{
    vtk::writeList(fmt, values);
    vtk::writeList(fmt, indirect, addressing);
}


template<class Type>
void Foam::vtk::writeValueParallel
(
    vtk::formatter& fmt,
    const Type& val,
    const label count
)
{
    if (Pstream::master())
    {
        vtk::write(fmt, val, count);

        label subCount;
        Type subValue;

        // Receive each [size, value] tuple
        for (const int proci : Pstream::subProcs())
        {
            IPstream is(Pstream::commsTypes::blocking, proci);
            is >> subCount >> subValue;

            vtk::write(fmt, subValue, subCount);
        }
    }
    else
    {
        OPstream os
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        // Send [size, value] tuple
        os << count << val;
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values
)
{
    // List sizes
    const globalIndex sizes(values.size());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send to master
    if (!Pstream::master())
    {
        UOPstream os(Pstream::masterNo(), pBufs);
        if (is_contiguous<Type>::value)
        {
            os.write
            (
                reinterpret_cast<const char*>(values.cdata()),
                values.size_bytes()
            );
        }
        else
        {
            os << values;
        }
    }

    pBufs.finishedSends();

    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values);

        // Receive and write
        for (const int proci : Pstream::subProcs())
        {
            UIPstream is(proci, pBufs);

            {
                List<Type> recv(sizes.localSize(proci));

                if (is_contiguous<Type>::value)
                {
                    is.read
                    (
                        reinterpret_cast<char*>(recv.data()),
                        recv.size_bytes()
                    );
                }
                else
                {
                    is >> recv;
                }
                vtk::writeList(fmt, recv);
            }
        }
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const labelUList& addressing
)
{
    UIndirectList<Type> send(values, addressing);

    // List sizes
    const globalIndex sizes(send.size());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send to master
    if (!Pstream::master())
    {
        UOPstream os(Pstream::masterNo(), pBufs);
        os << send;
    }

    pBufs.finishedSends();

    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values, addressing);

        // Receive and write
        for (const int proci : Pstream::subProcs())
        {
            UIPstream is(proci, pBufs);

            {
                List<Type> recv;
                is >> recv;
                vtk::writeList(fmt, recv);
            }
        }
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values,
    const bitSet& selected
)
{
    List<Type> send;
    if (!Pstream::master())
    {
        send = subset(selected, values);
    }

    // List sizes.
    // NOTE okay to skip proc0 since we only need sizes (not offsets)
    const globalIndex sizes(send.size());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send to master
    if (!Pstream::master())
    {
        UOPstream os(Pstream::masterNo(), pBufs);
        if (is_contiguous<Type>::value)
        {
            os.write
            (
                reinterpret_cast<const char*>(send.cdata()),
                send.size_bytes()
            );
        }
        else
        {
            os << send;
        }
    }

    pBufs.finishedSends();

    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values, selected);

        // Receive and write
        for (const int proci : Pstream::subProcs())
        {
            UIPstream is(proci, pBufs);

            {
                List<Type> recv(sizes.localSize(proci));

                if (is_contiguous<Type>::value)
                {
                    is.read
                    (
                        reinterpret_cast<char*>(recv.data()),
                        recv.size_bytes()
                    );
                }
                else
                {
                    is >> recv;
                }
                vtk::writeList(fmt, recv);
            }
        }
    }
}


template<class Type>
void Foam::vtk::writeListsParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values1,
    const UList<Type>& values2
)
{
    // List sizes
    const globalIndex sizes1(values1.size());
    const globalIndex sizes2(values2.size());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send to master
    if (!Pstream::master())
    {
        UOPstream os(Pstream::masterNo(), pBufs);
        if (is_contiguous<Type>::value)
        {
            os.write
            (
                reinterpret_cast<const char*>(values1.cdata()),
                values1.size_bytes()
            );
            os.write
            (
                reinterpret_cast<const char*>(values2.cdata()),
                values2.size_bytes()
            );
        }
        else
        {
            os << values1 << values2;
        }
    }

    pBufs.finishedSends();

    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2);

        // Reserve max receive size
        DynamicList<Type> recv
        (
            max(sizes1.maxNonLocalSize(), sizes2.maxNonLocalSize())
        );

        // Receive and write
        for (const int proci : Pstream::subProcs())
        {
            UIPstream is(proci, pBufs);

            // values1
            {
                List<Type> recv(sizes1.localSize(proci));
                if (is_contiguous<Type>::value)
                {
                    is.read
                    (
                        reinterpret_cast<char*>(recv.data()),
                        recv.size_bytes()
                    );
                }
                else
                {
                    is >> recv;
                }
                vtk::writeList(fmt, recv);
            }

            // values2
            {
                List<Type> recv(sizes2.localSize(proci));
                if (is_contiguous<Type>::value)
                {
                    is.read
                    (
                        reinterpret_cast<char*>(recv.data()),
                        recv.size_bytes()
                    );
                }
                else
                {
                    is >> recv;
                }
                vtk::writeList(fmt, recv);
            }
        }
    }
}


template<class Type>
void Foam::vtk::writeListsParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values1,
    const UList<Type>& values2,
    const labelUList& addressing
)
{
    UIndirectList<Type> send2(values2, addressing);

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send to master
    if (!Pstream::master())
    {
        UOPstream os(Pstream::masterNo(), pBufs);
        os << values1 << send2;
    }

    pBufs.finishedSends();

    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2, addressing);

        // Receive and write
        for (const int proci : Pstream::subProcs())
        {
            UIPstream is(proci, pBufs);

            // values1
            {
                List<Type> recv;
                is >> recv;
                vtk::writeList(fmt, recv);
            }

            // values2 (send2)
            {
                List<Type> recv;
                is >> recv;
                vtk::writeList(fmt, recv);
            }
        }
    }
}


// ************************************************************************* //
