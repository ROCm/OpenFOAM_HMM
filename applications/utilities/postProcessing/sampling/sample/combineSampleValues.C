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


#include "combineSampleValues.H"
#include "Pstream.H"
#include "ListListOps.H"
#include "IndirectList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Combine values from all processors. Valid result only on master processor.
template<class T>
void combineSampleValues
(
    const PtrList<volFieldSampler<T> >& sampledFields,
    const labelListList& indexSets,
    PtrList<volFieldSampler<T> >& masterFields
)
{
    forAll(sampledFields, fieldI)
    {
        List<Field<T> > masterValues(indexSets.size());

        forAll(indexSets, setI)
        {
            // Collect data from all processors
            List<Field<T> > gatheredData(Pstream::nProcs());
            gatheredData[Pstream::myProcNo()] = sampledFields[fieldI][setI];
            Pstream::gatherList(gatheredData);

            if (Pstream::master())
            {
                Field<T> allData
                (
                    ListListOps::combine<Field<T> >
                    (
                        gatheredData,
                        Foam::accessOp<Field<T> >()
                    )
                );

                masterValues[setI] =
                    IndirectList<T>(allData, indexSets[setI])();
            }
        }
        masterFields.set
        (
            fieldI,
            new volFieldSampler<T>
            (
                masterValues,
                sampledFields[fieldI].name()
            )
        );
    }
}


// ************************************************************************* //

