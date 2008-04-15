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

Class
    Foam::dataSchedule

Description
    Gets for every data item the source processor and the send processor.

    Determines:
    - which processor do I need to send to / receive from
    - what are the indices into the data I need to send / I am receiving.

    Example:
    - we have global sampling points
    - for every sampling point we have the cell/processor containing the point
      (sampleCells, sampleCellProcs)
    - for every sampling point we also have the destination cell/processor
      (receiveCells, receiveCellProcs)
    - sampleCellProcs, receiveCellProcs have to be the same on all processors!

    then we construct dataSchedule:
    @code
        dataSchedule mySchedule(sampleCellProcs, receiveCellProcs);

        const labelListList& sendOrder = mySchedule.sendOrder();
        const labelListList& receiveOrder = mySchedule.receiveOrder();

        // Do all remote data
        // ~~~~~~~~~~~~~~~~~~

        forAll(mySchedule, commI)
        {
            const labelPair& twoProcs = mySchedule[i];
            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am sender. Send to recvProc.

                // Labels of my cells to sample
                labelList cellLabels
                (
                    IndirectList<label>(sampleCells, sendOrder[recvProc])
                );

                OPstream toProc(recvProc);
                toProc<< IndirectList<vector>(fld, cellLabels)();
            }
            else
            {
                // I am receiver. Receive from sendProc.
                IPstream fromProc(sendProc);

                vectorField fromFld(fromProc);

                // Destination cells
                labelList cellLabels
                (
                    IndirectList<label>(receiveCells, receiveOrder[sendProc])
                );

                forAll(fromFld, i)
                {
                    fld[cellLabels[i]] = fromFld[i];
                }
            }
        }

        // Do local data
        // ~~~~~~~~~~~~~

        labelList sendCells
        (
            IndirectList<label>(sampleCells, sendOrder[Pstream::myProcNo()])
        );

        vectorField toMyself(IndirectList<vector>(fld, sendCells)());

        // Destination cells
        labelList recvCells
        (
            IndirectList<label>
            (
                patchFaces,
                receiveOrder[Pstream::myProcNo()]
            )
        );

        forAll(toMyself, i)
        {
            recvCells[i] = toMySelf[i];
        }
    @endcode


SourceFiles
    dataSchedule.C

\*---------------------------------------------------------------------------*/

#ifndef dataSchedule_H
#define dataSchedule_H

#include "commSchedule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dataSchedule Declaration
\*---------------------------------------------------------------------------*/

class dataSchedule
:
    public List<labelPair>
{
    // Private data

        //- Per processor to send to the order in which the data needs to
        //  be packed.
        labelListList sendOrder_;

        //- Per processor to receive from the order how to unpack the
        //  received data.
        labelListList receiveOrder_;


public:

    ClassName("dataSchedule");

    // Constructors

        //- Construct from components
        dataSchedule(const labelList& sendProcs, const labelList& receiveProcs);


    // Member Functions

        //- Send order per processor
        const labelListList& sendOrder() const
        {
            return sendOrder_;
        }

        //- Receive order per processor
        const labelListList& receiveOrder() const
        {
            return receiveOrder_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
