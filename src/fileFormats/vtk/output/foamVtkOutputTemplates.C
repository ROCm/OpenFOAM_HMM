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
#include "Pstream.H"
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
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }

    // Gather [count, value] tuples, including from master
    const List<std::pair<label, Type>> countValues
    (
        UPstream::listGatherValues<std::pair<label, Type>>
        (
            std::pair<label, Type>(count, val)
        )
    );

    if (Pstream::master())
    {
        for (const auto& countVal : countValues)
        {
            // Write [count, value] tuple
            vtk::write(fmt, countVal.first, countVal.second);
        }
    }
}


template<class Type>
void Foam::vtk::writeListParallel
(
    vtk::formatter& fmt,
    const UList<Type>& values
)
{
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }


    // Gather sizes - master information, offsets are irrelevant
    const globalIndex procAddr
    (
        UPstream::listGatherValues<label>(values.size()),
        globalIndex::SIZES
    );


    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values);

        // Receive and write
        DynamicList<Type> recvData(procAddr.maxNonLocalSize());

        for (const label proci : procAddr.subProcs())
        {
            recvData.resize_nocopy(procAddr.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);
        }
    }
    else
    {
        // Send
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            values.cdata_bytes(),
            values.size_bytes()
        );
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
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }


    List<Type> sendData;
    if (!Pstream::master())
    {
        sendData = UIndirectList<Type>(values, addressing);
    }

    // Gather sizes - master information, offsets are irrelevant
    const globalIndex procAddr
    (
        UPstream::listGatherValues<label>(sendData.size()),
        globalIndex::SIZES
    );


    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values, addressing);

        // Receive and write
        DynamicList<Type> recvData(procAddr.maxNonLocalSize());

        for (const label proci : procAddr.subProcs())
        {
            recvData.resize_nocopy(procAddr.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);
        }
    }
    else
    {
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            sendData.cdata_bytes(),
            sendData.size_bytes()
        );
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
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }


    List<Type> sendData;
    if (!Pstream::master())
    {
        sendData = subset(selected, values);
    }

    // Gather sizes - master information, offsets are irrelevant
    const globalIndex procAddr
    (
        UPstream::listGatherValues<label>(sendData.size()),
        globalIndex::SIZES
    );


    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values, selected);

        // Receive and write
        DynamicList<Type> recvData(procAddr.maxNonLocalSize());

        for (const label proci : procAddr.subProcs())
        {
            recvData.resize_nocopy(procAddr.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);
        }
    }
    else
    {
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            sendData.cdata_bytes(),
            sendData.size_bytes()
        );
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
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }


    // Gather sizes - master information and offsets are irrelevant
    const globalIndex procAddr1
    (
        UPstream::listGatherValues<label>(values1.size()),
        globalIndex::SIZES
    );
    const globalIndex procAddr2
    (
        UPstream::listGatherValues<label>(values2.size()),
        globalIndex::SIZES
    );


    if (Pstream::master())
    {
        // Write master data
        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2);

        // Receive and write
        DynamicList<Type> recvData
        (
            max(procAddr1.maxNonLocalSize(), procAddr2.maxNonLocalSize())
        );

        for (const label proci : procAddr1.subProcs())
        {
            // values1
            recvData.resize_nocopy(procAddr1.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);

            // values2
            recvData.resize_nocopy(procAddr2.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);
        }
    }
    else
    {
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            values1.cdata_bytes(),
            values1.size_bytes()
        );

        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            values2.cdata_bytes(),
            values2.size_bytes()
        );
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
    if (!is_contiguous<Type>::value)
    {
        // Non-contiguous data does not make sense
        FatalErrorInFunction
            << "Contiguous data only" << endl
            << Foam::exit(FatalError);
    }


    List<Type> sendData2;
    if (!Pstream::master())
    {
        sendData2 = UIndirectList<Type>(values2, addressing);
    }


    // Gather sizes - master information, offsets are irrelevant
    const globalIndex procAddr1
    (
        UPstream::listGatherValues<label>(values1.size()),
        globalIndex::SIZES
    );
    const globalIndex procAddr2
    (
        UPstream::listGatherValues<label>(sendData2.size()),
        globalIndex::SIZES
    );


    if (Pstream::master())
    {
        // Write master data

        vtk::writeList(fmt, values1);
        vtk::writeList(fmt, values2, addressing);

        // Receive and write
        DynamicList<Type> recvData
        (
            max(procAddr1.maxNonLocalSize(), procAddr2.maxNonLocalSize())
        );

        for (const label proci : procAddr1.subProcs())
        {
            // values1
            recvData.resize_nocopy(procAddr1.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);

            // values2
            recvData.resize_nocopy(procAddr2.localSize(proci));
            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                proci,
                recvData.data_bytes(),
                recvData.size_bytes()
            );
            vtk::writeList(fmt, recvData);
        }
    }
    else
    {
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            values1.cdata_bytes(),
            values1.size_bytes()
        );

        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            sendData2.cdata_bytes(),
            sendData2.size_bytes()
        );
    }
}


// ************************************************************************* //
