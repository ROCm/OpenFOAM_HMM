/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::Pstream

Description
    Inter-processor communications stream.

SourceFiles
    Pstream.C
    PstreamBroadcast.C
    PstreamGather.C
    PstreamCombineGather.C
    PstreamGatherList.C
    PstreamExchange.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Pstream_H
#define Foam_Pstream_H

#include "UPstream.H"
#include "DynamicList.H"

// Legacy
// #define Foam_Pstream_scatter_nobroadcast

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Pstream Declaration
\*---------------------------------------------------------------------------*/

class Pstream
:
    public UPstream
{
   // Private Static Functions

        //- Exchange contiguous data. Sends sendBufs, receives into recvBufs.
        //  Data provided and received as container.
        template<class Container, class T>
        static void exchangeContainer
        (
            const UList<Container>& sendBufs,
            const labelUList& recvSizes,  //!< Num of recv elements (not bytes)
            List<Container>& recvBufs,
            const int tag,
            const label comm,
            const bool wait               //!< Wait for requests to complete
        );

        //- Exchange contiguous data. Sends sendBufs, receives into recvBufs.
        //  Data provided and received as pointers.
        template<class T>
        static void exchangeBuf
        (
            const labelUList& sendSizes,  //!< Num of send elements (not bytes)
            const UList<const char*>& sendBufs,
            const labelUList& recvSizes,  //!< Num of recv elements (not bytes)
            List<char*>& recvBufs,
            const int tag,
            const label comm,
            const bool wait               //!< Wait for requests to complete
        );


protected:

    // Protected Data

        //- Allocated transfer buffer (can be used for send or receive)
        DynamicList<char> transferBuf_;


public:

    // Declare name of the class and its debug switch
    ClassName("Pstream");


    // Constructors

        //- Construct for given commsTypes, with optional buffer size
        explicit Pstream
        (
            const UPstream::commsTypes commsType,
            const label bufSize = 0
        )
        :
            UPstream(commsType)
        {
            if (bufSize > 0)
            {
                transferBuf_.setCapacity(bufSize + 2*sizeof(scalar) + 1);
            }
        }


    // Static Functions

    // Broadcast

        //- Broadcast buffer content to all processes in communicator.
        using UPstream::broadcast;

        //- Broadcast content (contiguous or non-contiguous)
        //- to all processes in communicator.
        template<class Type>
        static void broadcast
        (
            Type& value,
            const label comm = UPstream::worldComm
        );

        //- Broadcast multiple items to all processes in communicator.
        template<class Type, class... Args>
        static void broadcasts(const label comm, Type& arg1, Args&&... args);


    // Gather

        //- Gather (reduce) data, appyling \c bop to combine \c value
        //- from different processors. The basis for Foam::reduce().
        //  Uses the specified communication schedule.
        template<class T, class BinaryOp>
        static void gather
        (
            const List<commsStruct>& comms,
            T& value,
            const BinaryOp& bop,
            const int tag,
            const label comm
        );

        //- Gather (reduce) data, applying \c bop to combine \c value
        //- from different processors. The basis for Foam::reduce().
        //  Uses linear/tree communication.
        template<class T, class BinaryOp>
        static void gather
        (
            T& value,
            const BinaryOp& bop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Gather/combine data
    // Inplace combine values from processors.
    // (Uses construct from Istream instead of <<)

        //- Gather data, applying \c cop to inplace combine \c value
        //- from different processors.
        //  Uses the specified communication schedule.
        template<class T, class CombineOp>
        static void combineGather
        (
            const List<commsStruct>& comms,
            T& value,
            const CombineOp& cop,
            const int tag,
            const label comm
        );

        //- Gather data, applying \c cop to inplace combine \c value
        //- from different processors.
        //  Uses linear/tree communication.
        template<class T, class CombineOp>
        static void combineGather
        (
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Gather data, applying \c cop to inplace combine \c value
        //- from different processors.
        //- After completion all processors have the same data.
        //  Uses the specified communication schedule.
        //  Wraps combineGather/broadcast (may change in the future).
        template<class T, class CombineOp>
        static void combineAllGather
        (
            const List<commsStruct>& comms,
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Gather data, applying \c cop to inplace combine \c value
        //- from different processors.
        //- After completion all processors have the same data.
        //  Uses linear/tree communication.
        //  Wraps combineGather/broadcast (may change in the future).
        template<class T, class CombineOp>
        static void combineAllGather
        (
            T& value,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Combine variants working on whole List at a time.

        template<class T, class CombineOp>
        static void listCombineGather
        (
            const List<commsStruct>& comms,
            List<T>& values,
            const CombineOp& cop,
            const int tag,
            const label comm
        );

        //- Like above but switches between linear/tree communication
        template<class T, class CombineOp>
        static void listCombineGather
        (
            List<T>& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- After completion all processors have the same data.
        template<class T, class CombineOp>
        static void listCombineAllGather
        (
            const List<commsStruct>& comms,
            List<T>& values,
            const CombineOp& cop,
            const int tag,
            const label comm
        );

        //- After completion all processors have the same data.
        template<class T, class CombineOp>
        static void listCombineAllGather
        (
            List<T>& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Combine variants working on whole map at a time.
    // Container needs iterators, find() and insert methods defined.

        template<class Container, class CombineOp>
        static void mapCombineGather
        (
            const List<commsStruct>& comms,
            Container& values,
            const CombineOp& cop,
            const int tag,
            const label comm
        );

        //- Like above but switches between linear/tree communication
        template<class Container, class CombineOp>
        static void mapCombineGather
        (
            Container& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


        //- After completion all processors have the same data.
        template<class Container, class CombineOp>
        static void mapCombineAllGather
        (
            const List<commsStruct>& comms,
            Container& values,
            const CombineOp& cop,
            const int tag,
            const label comm
        );

        //- After completion all processors have the same data.
        template<class Container, class CombineOp>
        static void mapCombineAllGather
        (
            Container& values,
            const CombineOp& cop,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


    // Gather/scatter keeping the individual processor data separate.
    // The values is a List of size UPstream::nProcs() where
    // values[UPstream::myProcNo()] is the data for the current processor.

        //- Gather data, but keep individual values separate.
        //- Uses the specified communication schedule.
        template<class T>
        static void gatherList
        (
            const List<commsStruct>& comms,
            List<T>& values,
            const int tag,
            const label comm
        );

        //- Gather data, but keep individual values separate.
        //- Uses linear/tree communication.
        template<class T>
        static void gatherList
        (
            List<T>& values,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );

        //- Gather data, but keep individual values separate.
        //- Uses the specified communication schedule.
        //  After completion all processors have the same data.
        //  Wraps gatherList/scatterList (may change in the future).
        template<class T>
        static void allGatherList
        (
            const List<commsStruct>& comms,
            List<T>& values,
            const int tag,
            const label comm
        );

        //- Gather data, but keep individual values separate.
        //- Uses linear/tree communication.
        //  After completion all processors have the same data.
        //  Wraps gatherList/scatterList (may change in the future).
        template<class T>
        static void allGatherList
        (
            List<T>& values,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm
        );


        // Scatter

            //- Broadcast data: Distribute without modification.
            //  \note comms and tag parameters only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void scatter
            (
                const List<commsStruct>& comms,
                T& value,
                const int tag,
                const label comm
            );

            //- Broadcast data: Distribute without modification.
            //  \note tag parameter only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void scatter
            (
                T& value,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Broadcast data: Distribute without modification.
            //  \note tag parameter only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void combineScatter
            (
                const List<commsStruct>& comms,
                T& value,
                const int tag,
                const label comm
            );

            //- Broadcast data: Distribute without modification.
            //  \note tag parameter only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void combineScatter
            (
                T& value,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Broadcast data: Distribute without modification.
            //  \note comms and tag parameters only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void listCombineScatter
            (
                const List<commsStruct>& comms,
                List<T>& value,
                const int tag,
                const label comm
            );

            //- Broadcast data: Distribute without modification.
            //  \note comms and tag parameters only used when
            //  Foam_Pstream_scatter_nobroadcast is defined
            template<class T>
            static void listCombineScatter
            (
                List<T>& value,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Broadcast data: Distribute without modification.
            template<class Container>
            static void mapCombineScatter
            (
                const List<commsStruct>& comms,
                Container& values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class Container>
            static void mapCombineScatter
            (
                Container& values,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );


            //- Scatter data. Reverse of gatherList
            template<class T>
            static void scatterList
            (
                const List<commsStruct>& comms,
                List<T>& values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void scatterList
            (
                List<T>& values,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );


        // Exchange

            //- Helper: exchange sizes of sendData for specified
            //- set of send/receive processes.
            template<class Container>
            static void exchangeSizes
            (
                const labelUList& sendProcs,
                const labelUList& recvProcs,
                const Container& sendData,
                labelList& sizes,
                const label tag = UPstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Helper: exchange sizes of sendData.
            //- The sendData is the data per processor (in the communicator).
            //  Returns sizes of sendData on the sending processor.
            template<class Container>
            static void exchangeSizes
            (
                const Container& sendData,
                labelList& sizes,
                const label comm = UPstream::worldComm
            );


            //- Helper: exchange contiguous data.
            //- Sends sendData, receives into recvData.
            //  If wait=true will wait for all transfers to finish.
            template<class Container, class T>
            static void exchange
            (
                const UList<Container>& sendData,
                const labelUList& recvSizes,
                List<Container>& recvData,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm,
                const bool wait = true  //!< Wait for requests to complete
            );

            //- Exchange contiguous data.
            //- Sends sendData, receives into recvData.
            //- Determines sizes to receive.
            //  If wait=true will wait for all transfers to finish.
            template<class Container, class T>
            static void exchange
            (
                const UList<Container>& sendData,
                List<Container>& recvData,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm,
                const bool wait = true  //!< Wait for requests to complete
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PstreamBroadcast.C"
    #include "PstreamGather.C"
    #include "PstreamCombineGather.C"
    #include "PstreamGatherList.C"
    #include "PstreamExchange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
