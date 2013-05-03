/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "globalIndex.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalIndex::gather
(
    const label comm,
    const labelList& procIDs,
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag
) const
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allFld.setSize(size());

        // Assign my local data
        SubList<Type>(allFld, fld.size(), 0).assign(fld);

        for (label i = 1; i < procIDs.size(); i++)
        {
            SubList<Type> procSlot
            (
                allFld,
                offsets_[i+1]-offsets_[i],
                offsets_[i]
            );

            if (contiguous<Type>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    procIDs[i],
                    reinterpret_cast<char*>(procSlot.begin()),
                    procSlot.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromSlave
                (
                    Pstream::scheduled,
                    procIDs[i],
                    0,
                    tag,
                    comm
                );
                fromSlave >> procSlot;
            }
        }
    }
    else
    {
        if (contiguous<Type>())
        {
            OPstream::write
            (
                Pstream::scheduled,
                procIDs[0],
                reinterpret_cast<const char*>(fld.begin()),
                fld.byteSize(),
                tag,
                comm
            );
        }
        else
        {
            OPstream toMaster
            (
                Pstream::scheduled,
                procIDs[0],
                0,
                tag,
                comm
            );
            toMaster << fld;
        }
    }
}


template<class Type>
void Foam::globalIndex::gather
(
    const label comm,
    const labelList& procIDs,
    List<Type>& fld,
    const int tag
) const
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        List<Type> allFld(size());

        // Assign my local data
        SubList<Type>(allFld, fld.size(), 0).assign(fld);

        for (label i = 1; i < procIDs.size(); i++)
        {
            SubList<Type> procSlot
            (
                allFld,
                offsets_[i+1]-offsets_[i],
                offsets_[i]
            );

            if (contiguous<Type>())
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    procIDs[i],
                    reinterpret_cast<char*>(procSlot.begin()),
                    procSlot.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromSlave
                (
                    Pstream::scheduled,
                    procIDs[i],
                    0,
                    tag,
                    comm
                );
                fromSlave >> procSlot;
            }
        }

        fld.transfer(allFld);
    }
    else
    {
        if (contiguous<Type>())
        {
            OPstream::write
            (
                Pstream::scheduled,
                procIDs[0],
                reinterpret_cast<const char*>(fld.begin()),
                fld.byteSize(),
                tag,
                comm
            );
        }
        else
        {
            OPstream toMaster
            (
                Pstream::scheduled,
                procIDs[0],
                0,
                tag,
                comm
            );
            toMaster << fld;
        }
    }
}


template<class Type>
void Foam::globalIndex::scatter
(
    const label comm,
    const labelList& procIDs,
    const UList<Type>& allFld,
    UList<Type>& fld,
    const int tag
) const
{
    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        fld.assign
        (
            SubList<Type>
            (
                allFld,
                offsets_[1]-offsets_[0]
            )
        );

        for (label i = 1; i < procIDs.size(); i++)
        {
            const SubList<Type> procSlot
            (
                allFld,
                offsets_[i+1]-offsets_[i],
                offsets_[i]
            );

            if (contiguous<Type>())
            {
                OPstream::write
                (
                    Pstream::scheduled,
                    procIDs[i],
                    reinterpret_cast<const char*>(procSlot.begin()),
                    procSlot.byteSize(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toSlave
                (
                    Pstream::scheduled,
                    procIDs[i],
                    0,
                    tag,
                    comm
                );
                toSlave << procSlot;
            }
        }
    }
    else
    {
        if (contiguous<Type>())
        {
            IPstream::read
            (
                Pstream::scheduled,
                procIDs[0],
                reinterpret_cast<char*>(fld.begin()),
                fld.byteSize(),
                tag,
                comm
            );
        }
        else
        {
            IPstream fromMaster
            (
                Pstream::scheduled,
                procIDs[0],
                0,
                tag,
                comm
            );
            fromMaster >> fld;
        }
    }
}


// ************************************************************************* //
