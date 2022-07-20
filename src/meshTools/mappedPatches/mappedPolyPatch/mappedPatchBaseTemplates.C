/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

template<class Type>
void Foam::mappedPatchBase::distribute(List<Type>& lst) const
{
    const label myComm = getCommunicator();  // Get or create
    const label oldWarnComm(Pstream::warnComm);

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm(Pstream::worldComm);
            const auto& interp = AMI();

            Pstream::warnComm = myComm;
            Pstream::worldComm = myComm;

            if (sameWorld())
            {
                // lst is the other side's values
                lst = interp.interpolateToSource(Field<Type>(std::move(lst)));
            }
            else
            {
                // lst is my local data. Now the mapping in the AMI is
                // from my side to other side. Each processor contains either
                // faces from one side or from the other side.

                if (masterWorld())
                {
                    // I have lst.size() faces on my side, zero of the other
                    // side

                    tmp<Field<Type>> tmasterFld
                    (
                        interp.interpolateToSource(Field<Type>(0))
                    );
                    (void)interp.interpolateToTarget
                    (
                        Field<Type>(std::move(lst))
                    );

                    // We've received in our interpolateToSource the
                    // contribution from the other side
                    lst = tmasterFld;
                }
                else
                {
                    (void)interp.interpolateToSource
                    (
                        Field<Type>(std::move(lst))
                    );
                    tmp<Field<Type>> tmasterFld
                    (
                        interp.interpolateToTarget(Field<Type>(0))
                    );

                    // We've received in our interpolateToTarget the
                    // contribution from the other side
                    lst = tmasterFld;
                }
            }
            Pstream::worldComm = oldWorldComm;
            Pstream::warnComm = oldWarnComm;
            break;
        }
        default:
        {
            const auto& m = map();

            Pstream::warnComm = m.comm();
            m.distribute(lst);
            Pstream::warnComm = oldWarnComm;
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label myComm = getCommunicator();  // Get or create
    const label oldWarnComm(Pstream::warnComm);

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm(Pstream::worldComm);
            const auto& interp = AMI();
            Pstream::warnComm = myComm;
            Pstream::worldComm = myComm;
            lst = interp.interpolateToSource(Field<Type>(std::move(lst)), cop);
            Pstream::worldComm = oldWorldComm;
            Pstream::warnComm = oldWarnComm;
            break;
        }
        default:
        {
            // Force early construction of parallel data
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();

            Pstream::warnComm = myComm;
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                m.schedule(),
                m.constructSize(),
                m.subMap(),
                false,
                m.constructMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                myComm
            );
            Pstream::warnComm = oldWarnComm;
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute(List<Type>& lst) const
{
    const label myComm = getCommunicator();  // Get or create
    const label oldWarnComm(Pstream::warnComm);

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm(Pstream::worldComm);
            const auto& interp = AMI();
            Pstream::warnComm = myComm;
            Pstream::worldComm = myComm;
            lst = interp.interpolateToTarget(Field<Type>(std::move(lst)));
            Pstream::worldComm = oldWorldComm;
            Pstream::warnComm = oldWarnComm;
            break;
        }
        default:
        {
            // Force early construction of parallel data
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();

            Pstream::warnComm = m.comm();
            m.reverseDistribute(sampleSize(), lst);
            Pstream::warnComm = oldWarnComm;
            break;
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label myComm = getCommunicator();  // Get or create
    const label oldWarnComm(Pstream::warnComm);

    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm(Pstream::worldComm);
            const auto& interp = AMI();
            Pstream::warnComm = myComm;
            Pstream::worldComm = myComm;
            lst = interp.interpolateToTarget(Field<Type>(std::move(lst)), cop);
            Pstream::worldComm = oldWorldComm;
            Pstream::warnComm = oldWarnComm;
            break;
        }
        default:
        {
            (void)patch_.boundaryMesh().mesh().tetBasePtIs();
            const auto& m = map();
            const label cSize = sampleSize();
            Pstream::warnComm = myComm;
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                m.schedule(),
                cSize,
                m.constructMap(),
                false,
                m.subMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                myComm
            );
            Pstream::warnComm = oldWarnComm;
            break;
        }
    }
}


template<class Type>
bool Foam::mappedPatchBase::writeIOField
(
    const regIOobject& obj,
    dictionary& dict
)
{
    const auto* fldPtr = isA<IOField<Type>>(obj);
    if (fldPtr)
    {
        const auto& fld = *fldPtr;

        token tok;
        tok = new token::Compound<List<Type>>(fld);

        primitiveEntry* pePtr = new primitiveEntry
        (
            fld.name(),
            tokenList
            (
                one(),
                std::move(tok)
            )
        );

        dict.set(pePtr);
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool Foam::mappedPatchBase::constructIOField
(
    const word& name,
    token& tok,
    Istream& is,
    objectRegistry& obr
)
{
    const word tag = "List<" + word(pTraits<Type>::typeName) + '>';

    if (tok.isCompound() && tok.compoundToken().type() == tag)
    {
        IOField<Type>* fldPtr = obr.findObject<IOField<Type>>(name);
        if (fldPtr)
        {
            fldPtr->transfer
            (
                dynamicCast<token::Compound<List<Type>>>
                (
                    tok.transferCompoundToken(is)
                )
            );
        }
        else
        {
            IOField<Type>* fldPtr = new IOField<Type>
            (
                IOobject
                (
                    name,
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                label(0)
            );
            fldPtr->transfer
            (
                dynamicCast<token::Compound<List<Type>>>
                (
                    tok.transferCompoundToken(is)
                )
            );
            objectRegistry::store(fldPtr);
        }
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void Foam::mappedPatchBase::storeField
(
    objectRegistry& obr,
    const word& fieldName,
    const Field<Type>& values
)
{
    IOField<Type>* fldPtr = obr.findObject<IOField<Type>>(fieldName);
    if (fldPtr)
    {
        *fldPtr = values;
    }
    else
    {
        fldPtr = new IOField<Type>
        (
            IOobject
            (
                fieldName,
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            values
        );
        objectRegistry::store(fldPtr);
    }
}


// ************************************************************************* //
