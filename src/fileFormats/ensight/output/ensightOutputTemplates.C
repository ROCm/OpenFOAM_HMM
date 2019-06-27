/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ensightOutput.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<typename> class FieldContainer, class Type>
bool Foam::ensightOutput::Detail::writeFieldComponents
(
    const char* key,
    const FieldContainer<Type>& fld,
    ensightFile& os,
    bool parallel
)
{
    // Preliminary checks
    // ~~~~~~~~~~~~~~~~~~

    parallel = parallel && Pstream::parRun();

    bool hasField = !fld.empty();

    if (parallel)
    {
        reduce(hasField, orOp<bool>());
    }

    // Nothing to write
    if (!hasField) return false;

    // End preliminary checks
    // ~~~~~~~~~~~~~~~~~~~~~~


    if (Pstream::master())
    {
        os.writeKeyword(key);

        if (!parallel)
        {
            // Serial output
            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                const label cmpt = ensightPTraits<Type>::componentOrder[d];

                os.writeList(fld.component(cmpt));
            }
        }
        else
        {
            // Parallel (master)

            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                const label cmpt = ensightPTraits<Type>::componentOrder[d];

                os.writeList(fld.component(cmpt));

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    scalarField received(fromSlave);
                    os.writeList(received);
                }
            }
        }
    }
    else if (parallel)
    {
        // Parallel (slaves)

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const label cmpt = ensightPTraits<Type>::componentOrder[d];

            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << fld.component(cmpt);
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os,
    bool parallel
)
{
    // Preliminary checks
    // ~~~~~~~~~~~~~~~~~~

    parallel = parallel && Pstream::parRun();

    bool hasGeom = (parallel ? part.total() : part.size());
    bool hasField = !fld.empty();

    if (parallel)
    {
        // Used 'pre-reduced' information for hasGeom
        reduce(hasField, orOp<bool>());
    }

    // Nothing to write
    if (!hasGeom || !hasField) return false;

    // End preliminary checks
    // ~~~~~~~~~~~~~~~~~~~~~~


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const ensightFaces::elemType what = ensightFaces::elemType(typei);

        writeFieldComponents
        (
            ensightFaces::key(what),
            Field<Type>(fld, part.faceIds(what)),
            os,
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceSubField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os,
    bool parallel
)
{
    // Preliminary checks
    // ~~~~~~~~~~~~~~~~~~

    parallel = parallel && Pstream::parRun();

    bool hasGeom = (parallel ? part.total() : part.size());
    bool hasField = !fld.empty();

    if (parallel)
    {
        // Used 'pre-reduced' information for hasGeom
        reduce(hasField, orOp<bool>());
    }

    // Nothing to write
    if (!hasGeom || !hasField) return false;

    // End preliminary checks
    // ~~~~~~~~~~~~~~~~~~~~~~


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }


    label start = 0; // The start of the sub-list
    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const ensightFaces::elemType what = ensightFaces::elemType(typei);

        const label size = part.faceIds(what).size();

        writeFieldComponents
        (
            ensightFaces::key(what),
            SubField<Type>(fld, size, start),
            os,
            parallel
        );

        start += size;  // Advance the start for next sub-list
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Serial::writePointField
(
    const Field<Type>& fld,
    const ensightFaces& part,
    ensightFile& os
)
{
    // Preliminary checks
    // ~~~~~~~~~~~~~~~~~~

    bool parallel = false && Pstream::parRun();

    bool hasGeom = (parallel ? part.total() : part.size());
    bool hasField = !fld.empty();

    if (parallel)
    {
        // Used 'pre-reduced' information for hasGeom
        reduce(hasField, orOp<bool>());
    }

    // Nothing to write
    if (!hasGeom || !hasField) return false;

    // End preliminary checks
    // ~~~~~~~~~~~~~~~~~~~~~~


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    ensightOutput::Detail::writeFieldComponents
    (
        "coordinates",
        fld,
        os,
        parallel
    );

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeCellField
(
    const Field<Type>& fld,
    const ensightCells& part,
    ensightFile& os,
    bool parallel
)
{
    // Preliminary checks
    // ~~~~~~~~~~~~~~~~~~

    parallel = parallel && Pstream::parRun();

    bool hasGeom = (parallel ? part.total() : part.size());
    bool hasField = !fld.empty();

    if (parallel)
    {
        // Used 'pre-reduced' information for hasGeom
        reduce(hasField, orOp<bool>());
    }

    // Nothing to write
    if (!hasGeom || !hasField) return false;

    // End preliminary checks
    // ~~~~~~~~~~~~~~~~~~~~~~


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const ensightCells::elemType what = ensightCells::elemType(typei);

        Detail::writeFieldComponents
        (
            ensightCells::key(what),
            Field<Type>(fld, part.cellIds(what)),
            os,
            parallel
        );
    }

    return true;
}


// ************************************************************************* //
