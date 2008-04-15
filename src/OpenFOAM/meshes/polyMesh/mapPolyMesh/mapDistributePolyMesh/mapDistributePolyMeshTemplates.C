/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Distribute list.
template<class T>
void Foam::mapDistributePolyMesh::distribute
(
    const label newSize,
    const labelListList& subMap,        // from subset to original
    const labelListList& constructMap,  // from subset to new
    List<T>& field
)
{
    // Send sub field to neighbour
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        if (domain != Pstream::myProcNo())
        {
            OPstream toNbr(Pstream::blocking, domain);
            toNbr << IndirectList<T>(field, subMap[domain])();
        }
    }

    // Subset myself
    List<T> subField(IndirectList<T>(field, subMap[Pstream::myProcNo()]));

    // Receive sub field from myself (subField)
    field.setSize(newSize);

    const labelList& map = constructMap[Pstream::myProcNo()];

    forAll(map, i)
    {
        field[map[i]] = subField[i];
    }


    // Receive sub field from neighbour
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        if (domain != Pstream::myProcNo())
        {
            IPstream fromNbr(Pstream::blocking, domain);
            List<T> subField(fromNbr);

            const labelList& map = constructMap[domain];

            forAll(map, i)
            {
                field[map[i]] = subField[i];
            }
        }
    }
}


// ************************************************************************* //
