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

#include "ensightOutputCloud.H"
#include "ensightPTraits.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::ensightOutput::Detail::writeCloudFieldContent
(
    ensightFile& os,
    const UList<Type>& field,
    label count
)
{
    // Write master data
    for (Type val : field)          // <-- working on a copy!
    {
        if (mag(val) < 1e-90)       // approximately root(ROOTVSMALL)
        {
            val = Zero;
        }

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            os.write(component(val, cmpt));

            if (++count % 6 == 0)
            {
                os.newline();
            }
        }
    }

    return count;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os
)
{
    if (returnReduce(field.empty(), andOp<bool>()))
    {
        return false;
    }

    // Size information (offsets are irrelevant)
    globalIndex procAddr;
    if (Pstream::parRun())
    {
        procAddr.reset(UPstream::listGatherValues<label>(field.size()));
    }
    else
    {
        procAddr.reset(labelList(Foam::one{}, field.size()));
    }


    if (Pstream::master())
    {
        // 6 values per line
        label count = 0;

        // Write master data
        count = ensightOutput::Detail::writeCloudFieldContent
        (
            os,
            field,
            count
        );

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

            count = ensightOutput::Detail::writeCloudFieldContent
            (
                os,
                recvData,
                count
            );
        }

        // Add final newline if required
        if (count % 6)
        {
            os.newline();
        }
    }
    else
    {
        // Send
        UOPstream::write
        (
            UPstream::commsTypes::scheduled,
            Pstream::masterNo(),
            field.cdata_bytes(),
            field.size_bytes()
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writeCloudField
(
    const IOobject& io,
    const bool exists,
    autoPtr<ensightFile>& output
)
{
    if (exists)
    {
        // When exists == true, it exists globally,
        // but can still be missing on the local processor.
        // Handle this by READ_IF_PRESENT instead.

        IOobject fieldObj(io);
        fieldObj.readOpt(IOobject::READ_IF_PRESENT);

        IOField<Type> field(fieldObj);

        writeCloudField(field, output.ref());
    }

    return true;
}


// ************************************************************************* //
