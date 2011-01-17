/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Class
    Foam::mapDistribute

Description
    Class containing processor-to-processor mapping information.

    We store mapping from the bits-to-send to the complete starting list
    (subXXXMap) and from the received bits to their location in the new
    list (constructXXXMap).

Note:
    Schedule is a list of processor pairs (one send, one receive. One of
    them will be myself) which forms a scheduled (i.e. non-buffered) exchange.
    See distribute on how to use it.
    Note2: number of items send on one processor have to equal the number
    of items received on the other processor.

    Constructors using compact numbering: all my own elements first
    (whether used or not) followed by used-only remote elements.
    So e.g 4 procs and on proc 1 the compact
    table will first have all globalIndex.localSize() elements from proc1
    followed by used-only elements of proc0, proc2, proc3.
    The constructed mapDistribute sends the local elements from and
    receives the remote elements into their compact position.
    compactMap[procI] is the position of elements from procI in the compact
    map. compactMap[myProcNo()] is empty since trivial addressing. The indices
    into compactMap[procI] are local, not global, indices.

SourceFiles
    mapDistribute.C
    mapDistributeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistribute_H
#define mapDistribute_H

#include "labelList.H"
#include "labelPair.H"
#include "Pstream.H"
#include "boolList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;
class globalIndex;
class PstreamBuffers;
class globalIndexAndTransform;

/*---------------------------------------------------------------------------*\
                           Class mapDistribute Declaration
\*---------------------------------------------------------------------------*/

class mapDistribute
{
    // Private data

        //- Size of reconstructed data
        label constructSize_;

        //- Maps from subsetted data back to original data
        labelListList subMap_;

        //- Maps from subsetted data to new reconstructed data
        labelListList constructMap_;

        // Optional transformation

            //- For every globalIndexAndTransform::transformPermutations
            //  gives the elements that need to be transformed
            labelListList transformElements_;

            //- Destination in constructMap for transformed elements
            labelList transformStart_;


        //- Schedule
        mutable autoPtr<List<labelPair> > schedulePtr_;


   // Private Member Functions

        static void checkReceivedSize
        (
            const label procI,
            const label expectedSize,
            const label receivedSize
        );

        void printLayout
        (
            const label localSize,
            List<Map<label> >& compactMap,
            Ostream& os
        ) const;

        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelList& elements,
            List<Map<label> >& compactMap
        ) const;
        void calcCompactAddressing
        (
            const globalIndex& globalNumbering,
            const labelListList& elements,
            List<Map<label> >& compactMap
        ) const;

        void exchangeAddressing
        (
            const globalIndex& globalNumbering,
            labelList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );
        void exchangeAddressing
        (
            const globalIndex& globalNumbering,
            labelListList& elements,
            List<Map<label> >& compactMap,
            labelList& compactStart
        );


        //- Helper function: copy transformElements without transformation
        template<class T>
        void copyElements(List<T>& field) const;

        template<class T, class CombineOp>
        void applyTransforms
        (
            const globalIndexAndTransform& globalTransforms,
            List<T>& field,
            const bool isPosition,
            const CombineOp& cop
        ) const;
        template<class T, class CombineOp>
        void applyInverseTransforms
        (
            const globalIndexAndTransform& globalTransforms,
            List<T>& field,
            const bool isPosition,
            const CombineOp& cop
        ) const;


public:

    // Declare name of the class and its debug switch
    ClassName("mapDistribute");


    // Constructors

        //- Construct from components
        mapDistribute
        (
            const label constructSize,
            const Xfer<labelListList>& subMap,
            const Xfer<labelListList>& constructMap
        );

        //- Construct from components
        mapDistribute
        (
            const label constructSize,
            const Xfer<labelListList>& subMap,
            const Xfer<labelListList>& constructMap,
            const Xfer<labelListList>& transformElements,
            const Xfer<labelList>& transformStart
        );

        //- Construct from reverse addressing: per data item the send
        //  processor and the receive processor. All processors get same data.
        mapDistribute
        (
            const labelList& sendProcs,
            const labelList& recvProcs
        );

        //- Construct from list of (possibly) remote elements in globalIndex
        //  numbering (or -1). Determines compact numbering (see above) and
        //  distribute map to get data into this ordering and renumbers the
        //  elements to be in compact numbering.
        mapDistribute
        (
            const globalIndex&,
            labelList& elements,
            List<Map<label> >& compactMap
        );

        //- Construct from list of (possibly) remote elements in globalIndex
        //  numbering (or -1). Determines compact numbering (see above) and
        //  distribute map to get data into this ordering and renumbers the
        //  elements to be in compact numbering.
        mapDistribute
        (
            const globalIndex&,
            labelList& untransformedElements,
            const globalIndexAndTransform&,
            const labelPairList& transformedElements,
            labelList& transformedIndices,
            List<Map<label> >& compactMap
        );

        //- Special variant that works with the info sorted into bins
        //  according to local indices. E.g. think cellCells where
        //  cellCells[localCellI] is a list of global cells
        mapDistribute
        (
            const globalIndex&,
            labelListList& cellCells,
            List<Map<label> >& compactMap
        );

        //- As above but with transformed elements as well
        mapDistribute
        (
            const globalIndex&,
            labelListList& cellCells,
            const globalIndexAndTransform& ,
            const List<labelPairList>& transformedElements,
            labelListList& transformedIndices,
            List<Map<label> >& compactMap
        );

        //- Construct copy
        mapDistribute(const mapDistribute&);


    // Member Functions

        // Access

            //- Constructed data size
            label constructSize() const
            {
                return constructSize_;
            }

            //- Constructed data size
            label& constructSize()
            {
                return constructSize_;
            }

            //- From subsetted data back to original data
            const labelListList& subMap() const
            {
                return subMap_;
            }

            //- From subsetted data back to original data
            labelListList& subMap()
            {
                return subMap_;
            }

            //- From subsetted data to new reconstructed data
            const labelListList& constructMap() const
            {
                return constructMap_;
            }

            //- From subsetted data to new reconstructed data
            labelListList& constructMap()
            {
                return constructMap_;
            }

            //- Calculate a schedule. See above.
            static List<labelPair> schedule
            (
                const labelListList& subMap,
                const labelListList& constructMap
            );

            //- Return a schedule. Demand driven. See above.
            const List<labelPair>& schedule() const;


        // Other

            //- Helper for construct from globalIndex. Renumbers element
            //  (in globalIndex numbering) into compact indices.
            static label renumber
            (
                const globalIndex&,
                const List<Map<label> >& compactMap,
                const label globalElement
            );

            //- Compact maps. Gets per field a bool whether it is used (locally)
            //  and works out itself what this side and sender side can remove
            //  from maps.
            void compact(const boolList& elemIsUsed);


            //- Distribute data. Note:schedule only used for Pstream::scheduled
            //  for now, all others just use send-to-all, receive-from-all.
            template<class T>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&
            );

            //- Distribute data. If multiple processors writing to same
            //  position adds contributions using cop.
            template<class T, class CombineOp>
            static void distribute
            (
                const Pstream::commsTypes commsType,
                const List<labelPair>& schedule,
                const label constructSize,
                const labelListList& subMap,
                const labelListList& constructMap,
                List<T>&,
                const CombineOp& cop,
                const T& nullValue
            );

            //- Distribute data using default commsType.
            template<class T>
            void distribute(List<T>& fld) const;

            //- Same but with transforms
            template<class T>
            void distribute
            (
                const globalIndexAndTransform&,
                List<T>& fld,
                const bool isPosition
            ) const;

            //- Reverse distribute data using default commsType.
            template<class T>
            void reverseDistribute(const label constructSize, List<T>&) const;

            //- Same but with transforms
            template<class T>
            void reverseDistribute
            (
                const globalIndexAndTransform&,
                const label constructSize,
                List<T>& fld,
                const bool isPosition
            ) const;

            //- Reverse distribute data using default commsType.
            //  Since constructSize might be larger than supplied size supply
            //  a nullValue
            template<class T>
            void reverseDistribute
            (
                const label constructSize,
                const T& nullValue,
                List<T>& fld
            ) const;

            //- Same but with transforms
            template<class T>
            void reverseDistribute
            (
                const globalIndexAndTransform&,
                const label constructSize,
                const T& nullValue,
                List<T>& fld,
                const bool isPosition
            ) const;

            //- Do all sends using PstreamBuffers
            template<class T>
            void send(PstreamBuffers&, const List<T>&) const;
            //- Do all receives using PstreamBuffers
            template<class T>
            void receive(PstreamBuffers&, List<T>&) const;

            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistribute::updateMesh(const mapPolyMesh&)"
                );
            }

    // Member Operators

        void operator=(const mapDistribute&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mapDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
