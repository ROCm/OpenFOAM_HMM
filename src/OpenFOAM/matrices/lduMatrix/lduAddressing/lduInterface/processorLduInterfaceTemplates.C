/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "processorLduInterface.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions * * *  * * * * * * * * * * //

// template<class Type>
// Foam::List<Type>& Foam::processorLduInterface::setSendBuf(const label nElems)
// const
// {
//     if (!contiguous<Type>())
//     {
//         FatalErrorIn("processorLduInterface::setSendBuf(const label) const")
//             << "Cannot return the binary size of a list of "
//                "non-primitive elements"
//             << abort(FatalError);
//     }
// 
//     label nBytes = nElems*sizeof(Type);
//     sendBuf_.setSize(nBytes);
// 
//     return reinterpret_cast<List<Type>&>(sendBuf_);
// }
// 
// 
// template<class Type>
// Foam::List<Type>& Foam::processorLduInterface::setReceiveBuf
// (
//     const label nElems
// ) const
// {
//     if (!contiguous<Type>())
//     {
//         FatalErrorIn("processorLduInterface::setReceiveBuf(const label) const")
//             << "Cannot return the binary size of a list of "
//                "non-primitive elements"
//             << abort(FatalError);
//     }
// 
//     label nBytes = nElems*sizeof(Type);
// 
//     //receiveBuf_.setSize(nBytes, '\0');    // necessary because of expanding
//                                             // compression?
//     receiveBuf_.setSize(nBytes);
// 
//     return reinterpret_cast<List<Type>&>(receiveBuf_);
// }


template<class Type>
void Foam::processorLduInterface::send
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        OPstream::write
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<const char*>(f.begin()),
            f.byteSize()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        //setReceiveBuf<Type>(f.size());
        resizeBuf(receiveBuf_, f.size()*sizeof(Type));

        IPstream::read
        (
            commsType,
            neighbProcNo(),
            receiveBuf_.begin(),
            receiveBuf_.size()
        );

        //setSendBuf<Type>(f.size());
        resizeBuf(sendBuf_, f.byteSize());
        memcpy(sendBuf_.begin(), f.begin(), f.byteSize());

        OPstream::write
        (
            commsType,
            neighbProcNo(),
            sendBuf_.begin(),
            f.byteSize()
        );
    }
    else
    {
        FatalErrorIn("processorLduInterface::send")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}


template<class Type>
void Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
    {
        IPstream::read
        (
            commsType,
            neighbProcNo(),
            reinterpret_cast<char*>(f.begin()),
            f.byteSize()
        );
    }
    else if (commsType == Pstream::nonBlocking)
    {
        memcpy(f.begin(), receiveBuf_.begin(), f.byteSize());
    }
    else
    {
        FatalErrorIn("processorLduInterface::receive")
            << "Unsupported communications type " << commsType
            << exit(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::processorLduInterface::receive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<Field<Type> > tf(new Field<Type>(size));
    receive(commsType, tf());
    return tf;
}


template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const Pstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];
        //setSendBuf<float>(nFloats);
        resizeBuf(sendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(sendBuf_.begin());

        for (register label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        reinterpret_cast<Type&>(fArray[nm1]) = f[f.size() - 1];

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes
            );
        }
        else if (commsType == Pstream::nonBlocking)
        {
            //setReceiveBuf<float>(nFloats);
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                receiveBuf_.size()
            );

            OPstream::write
            (
                commsType,
                neighbProcNo(),
                sendBuf_.begin(),
                nBytes
            );
        }
        else
        {
            FatalErrorIn("processorLduInterface::compressedSend")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }
    }
    else
    {
        this->send(commsType, f);
    }
}

template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if (sizeof(scalar) != sizeof(float) && Pstream::floatTransfer && f.size())
    {
        static const label nCmpts = sizeof(Type)/sizeof(scalar);
        label nm1 = (f.size() - 1)*nCmpts;
        label nlast = sizeof(Type)/sizeof(float);
        label nFloats = nm1 + nlast;
        label nBytes = nFloats*sizeof(float);

        if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
        {
            //setReceiveBuf<float>(nFloats);
            resizeBuf(receiveBuf_, nBytes);

            IPstream::read
            (
                commsType,
                neighbProcNo(),
                receiveBuf_.begin(),
                nBytes
            );
        }
        else if (commsType != Pstream::nonBlocking)
        {
            FatalErrorIn("processorLduInterface::compressedReceive")
                << "Unsupported communications type " << commsType
                << exit(FatalError);
        }

        const float *fArray =
            reinterpret_cast<const float*>(receiveBuf_.begin());
        f[f.size() - 1] = reinterpret_cast<const Type&>(fArray[nm1]);
        scalar *sArray = reinterpret_cast<scalar*>(f.begin());
        const scalar *slast = &sArray[nm1];

        for (register label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }
    }
    else
    {
        this->receive<Type>(commsType, f);
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::processorLduInterface::compressedReceive
(
    const Pstream::commsTypes commsType,
    const label size
) const
{
    tmp<Field<Type> > tf(new Field<Type>(size));
    compressedReceive(commsType, tf());
    return tf;
}


// template<class Type>
// void Foam::processorLduInterface::compressedBufferSend
// (
//     const Pstream::commsTypes commsType
// ) const
// {
//     // Optionally inline compress sendBuf
//     if
//     (
//         sizeof(scalar) > sizeof(float)
//      && sendBuf_.size()
//      && Pstream::floatTransfer
//     )
//     {
//         const List<Type>& f = reinterpret_cast<const List<Type>&>(sendBuf_);
//         label fSize = f.size()/sizeof(Type);
// 
//         // Inplace compress
//         static const label nCmpts = sizeof(Type)/sizeof(scalar);
//         label nm1 = (fSize - 1)*nCmpts;
//         label nlast = sizeof(Type)/sizeof(float);
//         label nFloats = nm1 + nlast;
// 
//         const scalar *sArray = reinterpret_cast<const scalar*>(f.begin());
//         const scalar *slast = &sArray[nm1];
//         float *fArray = reinterpret_cast<float*>(sendBuf_.begin());
// 
//         for (register label i=0; i<nm1; i++)
//         {
//             fArray[i] = sArray[i] - slast[i%nCmpts];
//         }
// 
//         reinterpret_cast<Type&>(fArray[nm1]) = f[fSize - 1];
// 
//         // Trim
//         setSendBuf<float>(nFloats);
//     }
// 
//     // Send sendBuf
//     if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
//     {
//         OPstream::write
//         (
//             commsType,
//             neighbProcNo(),
//             sendBuf_.begin(),
//             sendBuf_.size()
//         );
//     }
//     else if (commsType == Pstream::nonBlocking)
//     {
//         setReceiveBuf<char>(sendBuf_.size());
// 
//         IPstream::read
//         (
//             commsType,
//             neighbProcNo(),
//             receiveBuf_.begin(),
//             receiveBuf_.size()
//         );
// 
//         OPstream::write
//         (
//             commsType,
//             neighbProcNo(),
//             sendBuf_.begin(),
//             sendBuf_.size()
//         );
//     }
//     else
//     {
//         FatalErrorIn("processorLduInterface::compressedBufferSend")
//             << "Unsupported communications type " << commsType
//             << exit(FatalError);
//     }
// }
// 
// 
// template<class Type>
// const Foam::List<Type>& Foam::processorLduInterface::compressedBufferReceive
// (
//     const Pstream::commsTypes commsType,
//     const label size
// ) const
// {
//     if (sizeof(scalar) > sizeof(float) && size && Pstream::floatTransfer)
//     {
//         static const label nCmpts = sizeof(Type)/sizeof(scalar);
//         label nm1 = (size - 1)*nCmpts;
//         label nlast = sizeof(Type)/sizeof(float);
//         label nFloats = nm1 + nlast;
// 
//         if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
//         {
//             setReceiveBuf<float>(nFloats);
// 
//             IPstream::read
//             (
//                 commsType,
//                 neighbProcNo(),
//                 receiveBuf_.begin(),
//                 receiveBuf_.size()
//             );
//         }
//         else if (commsType != Pstream::nonBlocking)
//         {
//             FatalErrorIn("processorLduInterface::compressedBufferReceive")
//                 << "Unsupported communications type " << commsType
//                 << exit(FatalError);
//         }
// 
//         // Inline expand
//         List<Type>& f = setReceiveBuf<Type>(size);
//         label fSize = f.size()/sizeof(Type);
// 
//         const float *fArray =
//             reinterpret_cast<const float*>(receiveBuf_.begin());
//         f[fSize - 1] = reinterpret_cast<const Type&>(fArray[nm1]);
//         scalar *sArray = reinterpret_cast<scalar*>(f.begin());
//         const scalar *slast = &sArray[nm1];
// 
//         for (register label i=0; i<nm1; i++)
//         {
//             sArray[i] = fArray[i] + slast[i%nCmpts];
//         }
//     }
//     else
//     {
//         if (commsType == Pstream::blocking || commsType == Pstream::scheduled)
//         {
//             setReceiveBuf<Type>(size);
// 
//             IPstream::read
//             (
//                 commsType,
//                 neighbProcNo(),
//                 receiveBuf_.begin(),
//                 receiveBuf_.size()
//             );
//         }
//         else if (commsType != Pstream::nonBlocking)
//         {
//             FatalErrorIn("processorLduInterface::compressedBufferReceive")
//                 << "Unsupported communications type " << commsType
//                 << exit(FatalError);
//         }
//     }
//     return reinterpret_cast<List<Type>&>(receiveBuf_);
// }


// ************************************************************************* //
