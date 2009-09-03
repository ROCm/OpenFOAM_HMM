/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "DirectInteractionList.H"
#include "InteractionLists.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::DirectInteractionList<ParticleType>::buildDirectInteractionList
(
    bool pointPointListBuild
)
{
    Info<< "    Building list of direct interaction neighbours" << endl;

    const polyMesh& mesh(il_.mesh());

    List<DynamicList<label> > DirectInteractionList(mesh.nCells());

    if (pointPointListBuild)
    {
        Info<< "        Point-Point direct interaction list build." << endl;

        label pointJIndex;

        forAll (mesh.points(), pointIIndex)
        {
            for
            (
                pointJIndex = pointIIndex;
                pointJIndex != mesh.points().size();
                ++pointJIndex
            )
            {
                if (il_.testPointPointDistance(pointIIndex, pointJIndex))
                {
                    const labelList& ptICells
                    (
                        mesh.pointCells()[pointIIndex]
                    );

                    const labelList& ptJCells
                    (
                        mesh.pointCells()[pointJIndex]
                    );

                    forAll(ptICells, pIC)
                    {
                        const label cellI(ptICells[pIC]);

                        forAll(ptJCells, pJC)
                        {
                            const label cellJ(ptJCells[pJC]);

                            if (cellJ > cellI)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellI],
                                        cellJ
                                    )
                                 == -1
                                )
                                {
                                    DirectInteractionList[cellI].append(cellJ);
                                }
                            }

                            if (cellI > cellJ)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellJ],
                                        cellI
                                    )
                                 ==
                                    -1
                                )
                                {
                                    DirectInteractionList[cellJ].append(cellI);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        Info<< "        Point-Face, Edge-Edge direct interaction list build."
            << endl;

        forAll(mesh.points(), p)
        {
            forAll(mesh.faces(), f)
            {
                if (il_.testPointFaceDistance(p, f))
                {
                    const labelList& pCells(mesh.pointCells()[p]);

                    const label cellO(mesh.faceOwner()[f]);

                    forAll(pCells, pC)
                    {
                        const label cellI(pCells[pC]);

                        // cells are not added to their own DIL

                        if (cellO > cellI)
                        {
                            if
                            (
                                findIndex
                                (
                                    DirectInteractionList[cellI],
                                    cellO
                                )
                             ==
                                -1
                            )
                            {
                                DirectInteractionList[cellI].append(cellO);
                            }
                        }

                        if (cellI > cellO)
                        {
                            if
                            (
                                findIndex
                                (
                                    DirectInteractionList[cellO],
                                    cellI
                                )
                             ==
                                -1
                            )
                            {
                                DirectInteractionList[cellO].append(cellI);
                            }
                        }

                        if (mesh.isInternalFace(f))
                        {
                            // boundary faces will not have neighbour
                            // information

                            const label cellN(mesh.faceNeighbour()[f]);

                            if (cellN > cellI)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellI],
                                        cellN
                                    )
                                 ==
                                    -1
                                )
                                {
                                    DirectInteractionList[cellI].append(cellN);
                                }
                            }

                            if (cellI > cellN)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellN],
                                        cellI
                                    )
                                 ==
                                    -1
                                )
                                {
                                    DirectInteractionList[cellN].append(cellI);
                                }
                            }
                        }
                    }
                }
            }
        }

        label edgeJIndex;

        forAll(mesh.edges(), edgeIIndex)
        {
            const edge& eI(mesh.edges()[edgeIIndex]);

            for
            (
                edgeJIndex = edgeIIndex + 1;
                edgeJIndex != mesh.edges().size();
                ++edgeJIndex
            )
            {
                const edge& eJ(mesh.edges()[edgeJIndex]);

                if (il_.testEdgeEdgeDistance(eI, eJ))
                {
                    const labelList& eICells(mesh.edgeCells()[edgeIIndex]);

                    const labelList& eJCells(mesh.edgeCells()[edgeJIndex]);

                    forAll(eICells, eIC)
                    {
                        const label cellI(eICells[eIC]);

                        forAll(eJCells, eJC)
                        {
                            const label cellJ(eJCells[eJC]);

                            if (cellJ > cellI)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellI],
                                        cellJ
                                    )
                                 ==
                                    -1
                                )
                                {
                                    DirectInteractionList[cellI].append(cellJ);
                                }
                            }

                            if (cellI > cellJ)
                            {
                                if
                                (
                                    findIndex
                                    (
                                        DirectInteractionList[cellJ],
                                        cellI
                                    )
                                 ==
                                    -1
                                )
                                {
                                    DirectInteractionList[cellJ].append(cellI);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    forAll(DirectInteractionList, transDIL)
    {
        (*this)[transDIL].transfer
        (
            DirectInteractionList[transDIL].shrink()
        );
    }

    // sorting DILs

    forAll((*this), dIL)
    {
        sort((*this)[dIL]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::DirectInteractionList<ParticleType>::DirectInteractionList
(
    const InteractionLists<ParticleType>& il,
    bool pointPointListBuild
)
:
    labelListList(il.mesh().nCells()),
    il_(il)
{
    if ((*this).size() > 1)
    {
        buildDirectInteractionList(pointPointListBuild);
    }
    else if ((*this).size() == 1)
    {
        Info<< "    Single cell mesh, no direct interaction lists required."
            << endl;

        (*this)[0].setSize(0);
    }
}


template<class ParticleType>
Foam::DirectInteractionList<ParticleType>::DirectInteractionList
(
    const InteractionLists<ParticleType>& il
)
:
    labelListList(il.mesh().nCells()),
    il_(il)
{
    Info<< "    Read DirectInteractionList from disk not implemented" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::DirectInteractionList<ParticleType>::~DirectInteractionList()
{}


// ************************************************************************* //
