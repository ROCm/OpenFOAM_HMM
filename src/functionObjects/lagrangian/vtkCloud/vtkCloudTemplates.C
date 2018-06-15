/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "IOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::functionObjects::vtkCloud::writeFields
(
    autoPtr<vtk::formatter>& format,
    const objectRegistry& obrTmp,
    const label nTotParcels
) const
{
    const int nCmpt(pTraits<Type>::nComponents);

    const bool useIntField =
        std::is_integral<typename pTraits<Type>::cmptType>();

    // Fields are not always on all processors (eg, multi-component parcels).
    // Thus need to resolve names between all processors.

    wordList fieldNames(obrTmp.names<IOField<Type>>());
    Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);

    // Sort to get identical order of fields on all processors
    Foam::sort(fieldNames);

    for (const word& fieldName : fieldNames)
    {
        const auto* fldPtr = obrTmp.lookupObjectPtr<IOField<Type>>(fieldName);

        if (Pstream::master())
        {
            if (useIntField)
            {
                const uint64_t payLoad(nTotParcels * nCmpt * sizeof(label));

                format().openDataArray<label, nCmpt>(fieldName)
                    .closeTag();

                format().writeSize(payLoad);

                if (fldPtr)
                {
                    // Data on master
                    const auto& fld = *fldPtr;

                    // Ensure consistent output width
                    for (const Type& val : fld)
                    {
                        for (int cmpt=0; cmpt < nCmpt; ++cmpt)
                        {
                            format().write(label(component(val, cmpt)));
                        }
                    }
                }

                // Slaves - recv
                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    Field<Type> recv(fromSlave);

                    for (const Type& val : recv)
                    {
                        for (int cmpt=0; cmpt < nCmpt; ++cmpt)
                        {
                            format().write(label(component(val, cmpt)));
                        }
                    }
                }
            }
            else
            {
                const uint64_t payLoad(nTotParcels * nCmpt * sizeof(float));

                format().openDataArray<float, nCmpt>(fieldName)
                    .closeTag();

                format().writeSize(payLoad);

                if (fldPtr)
                {
                    // Data on master
                    vtk::writeList(format(), *fldPtr);
                }

                // Slaves - recv
                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    Field<Type> recv(fromSlave);

                    vtk::writeList(format(), recv);
                }
            }

            format().flush();

            format()
                .endDataArray();
        }
        else
        {
            // Slaves - send

            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            if (fldPtr)
            {
                toMaster
                    << *fldPtr;
            }
            else
            {
                toMaster
                    << Field<Type>();
            }
        }
    }

    return fieldNames;
}


// ************************************************************************* //
