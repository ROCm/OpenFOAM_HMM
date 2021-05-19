/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    const direction nCmpt(pTraits<Type>::nComponents);

    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );

    // Other integral types (eg, bool etc) would need cast/convert to label.
    // Similarly for labelVector etc.


    // Fields are not always on all processors (eg, multi-component parcels).
    // Thus need to resolve names between all processors.

    wordList fieldNames(obrTmp.names<IOField<Type>>());
    Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);

    // Consistent order on all processors
    Foam::sort(fieldNames);

    for (const word& fieldName : fieldNames)
    {
        const List<Type>* fldPtr = obrTmp.findObject<IOField<Type>>(fieldName);
        const List<Type>& values = (fldPtr ? *fldPtr : List<Type>());

        if (Pstream::master())
        {
            if (std::is_same<label, typename pTraits<Type>::cmptType>::value)
            {
                const uint64_t payLoad =
                    vtk::sizeofData<label, nCmpt>(nTotParcels);

                format().beginDataArray<label, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
            else
            {
                const uint64_t payLoad =
                    vtk::sizeofData<float, nCmpt>(nTotParcels);

                format().beginDataArray<float, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }

        if (applyFilter_)
        {
            vtk::writeListParallel(format.ref(), values, parcelAddr_);
        }
        else
        {
            vtk::writeListParallel(format.ref(), values);
        }

        if (Pstream::master())
        {
            // Non-legacy
            format().flush();
            format().endDataArray();
        }
    }

    return fieldNames;
}


// ************************************************************************* //
