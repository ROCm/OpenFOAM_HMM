/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "ensightOutput.H"
#include "ensightPTraits.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::Detail::copyComponent
(
    List<scalar>& cmptBuffer,
    const FieldContainer<Type>& input,
    const direction cmpt
)
{
    if (cmptBuffer.size() < input.size())
    {
        FatalErrorInFunction
            << "Component buffer too small: "
            << cmptBuffer.size() << " < " << input.size() << nl
            << exit(FatalError);
    }

    auto iter = cmptBuffer.begin();

    for (const Type& val : input)
    {
        *iter = component(val, cmpt);
        ++iter;
    }
}


template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::Detail::writeFieldContent
(
    ensightFile& os,
    const FieldContainer<Type>& fld,
    bool parallel
)
{
    // already checked prior to calling, but extra safety
    parallel = parallel && Pstream::parRun();

    // Size information (offsets are irrelevant)
    globalIndex procAddr;
    if (parallel)
    {
        procAddr.reset(UPstream::listGatherValues<label>(fld.size()));
    }
    else
    {
        // Master size
        procAddr.reset(labelList(Foam::one{}, fld.size()));
    }


    if (Pstream::master())
    {
        DynamicList<scalar> cmptBuffer(procAddr.maxSize());

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            // Write master data
            cmptBuffer.resize_nocopy(procAddr.localSize(0));
            copyComponent(cmptBuffer, fld, cmpt);
            os.writeList(cmptBuffer);

            // Receive and write
            for (const label proci : procAddr.subProcs())
            {
                cmptBuffer.resize_nocopy(procAddr.localSize(proci));
                UIPstream::read
                (
                    UPstream::commsTypes::scheduled,
                    proci,
                    cmptBuffer.data_bytes(),
                    cmptBuffer.size_bytes()
                );
                os.writeList(cmptBuffer);
            }
        }
    }
    else if (parallel)
    {
        // Send
        List<scalar> cmptBuffer(fld.size());

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            copyComponent(cmptBuffer, fld, cmpt);
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                Pstream::masterNo(),
                cmptBuffer.cdata_bytes(),
                cmptBuffer.size_bytes()
            );
        }
    }
}


template<template<typename> class FieldContainer>
bool Foam::ensightOutput::Detail::writeCoordinates
(
    ensightGeoFile& os,
    const label partId,
    const word& partName,
    const label nPoints,
    const FieldContainer<Foam::point>& fld,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    if (Pstream::master())
    {
        os.beginPart(partId, partName);
        os.beginCoordinates(nPoints);
    }

    ensightOutput::Detail::writeFieldContent(os, fld, parallel);

    return true;
}


template<template<typename> class FieldContainer, class Type>
bool Foam::ensightOutput::Detail::writeFieldComponents
(
    ensightFile& os,
    const char* key,
    const FieldContainer<Type>& fld,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks
    {
        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.writeKeyword(key);
    }

    ensightOutput::Detail::writeFieldContent(os, fld, parallel);

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceSubField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFaces::key(etype),
            SubField<Type>(fld, part.range(etype)),
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceLocalField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }

    bool validAddressing = (part.size() == part.faceOrder().size());

    if (parallel)
    {
        reduce(validAddressing, orOp<bool>());
    }

    if (!validAddressing)
    {
        FatalErrorInFunction
            << "Called without faceOrder having been set" << nl
            << exit(FatalError);
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFaces::key(etype),
            UIndirectList<Type>(fld, part.faceOrder(etype)),
            parallel
        );
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightCells& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightCells::key(etype),
            UIndirectList<Type>(fld, part.cellIds(etype)),
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFaces::key(etype),
            UIndirectList<Type>(fld, part.faceIds(etype)),
            parallel
        );
    }

    return true;
}


// ************************************************************************* //
