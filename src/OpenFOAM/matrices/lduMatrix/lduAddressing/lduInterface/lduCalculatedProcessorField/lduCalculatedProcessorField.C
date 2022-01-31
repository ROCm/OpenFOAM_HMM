/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "lduCalculatedProcessorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduInterface& interface,
    const Field<Type>& iF
)
:
    LduInterfaceField<Type>(interface),
    procInterface_(refCast<const lduPrimitiveProcessorInterface>(interface)),
    field_(iF),
    sendBuf_(procInterface_.faceCells().size()),
    receiveBuf_(procInterface_.faceCells().size()),
    scalarSendBuf_(procInterface_.faceCells().size()),
    scalarReceiveBuf_(procInterface_.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduCalculatedProcessorField<Type>& ptf
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf)),
    procInterface_(ptf.procInterface_),
    field_(ptf.field_),
    sendBuf_(procInterface_.faceCells().size()),
    receiveBuf_(procInterface_.faceCells().size()),
    scalarSendBuf_(procInterface_.faceCells().size()),
    scalarReceiveBuf_(procInterface_.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::lduCalculatedProcessorField<Type>::ready() const
{
    if
    (
        this->outstandingSendRequest_ >= 0
     && this->outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        if (!UPstream::finishedRequest(this->outstandingSendRequest_))
        {
            return false;
        }
    }
    this->outstandingSendRequest_ = -1;

    if
    (
        this->outstandingRecvRequest_ >= 0
     && this->outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        if (!UPstream::finishedRequest(this->outstandingRecvRequest_))
        {
            return false;
        }
    }
    this->outstandingRecvRequest_ = -1;

    return true;
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Bypass patchInternalField since uses fvPatch addressing
    const labelList& fc = lduAddr.patchAddr(patchId);

    scalarSendBuf_.setSize(fc.size());
    forAll(fc, i)
    {
        scalarSendBuf_[i] = psiInternal[fc[i]];
    }

    if (!this->ready())
    {
        FatalErrorInFunction
            << "On patch "
            << " outstanding request."
            << abort(FatalError);
    }


    scalarReceiveBuf_.setSize(scalarSendBuf_.size());
    outstandingRecvRequest_ = UPstream::nRequests();

    UIPstream::read
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarReceiveBuf_.data_bytes(),
        scalarReceiveBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    outstandingSendRequest_ = UPstream::nRequests();

    UOPstream::write
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarSendBuf_.cdata_bytes(),
        scalarSendBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = false;
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::addToInternalField
(
    solveScalarField& result,
    const bool add,
    const scalarField& coeffs,
    const solveScalarField& vals
) const
{
    const labelUList& faceCells = this->procInterface_.faceCells();

    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*vals[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*vals[elemI];
        }
    }
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        UPstream::waitRequest(outstandingRecvRequest_);
    }
    // Recv finished so assume sending finished as well.
    outstandingSendRequest_ = -1;
    outstandingRecvRequest_ = -1;

    // Consume straight from scalarReceiveBuf_. Note use of our own
    // helper to avoid using fvPatch addressing
    addToInternalField(result, !add, coeffs, scalarReceiveBuf_);

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = true;
}


// ************************************************************************* //
