/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    Test-gatherValues1

Description
    Test list gather functionality

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    const labelList localValues
    (
        identity(2 *(Pstream::myProcNo()+1), -5*Pstream::myProcNo())
    );

    // Test resize
    {
        globalIndex globIdx(localValues.size());

        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;

        globIdx.setLocalSize(4, 0);
        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;
        globIdx.setLocalSize(3, 0);
        Info<< "globIdx = " << flatOutput(globIdx.offsets()) << nl;
    }


    // Gather all values
    {
        const auto& sendData = localValues;

        // One-sided sizing!  master only
        const globalIndex allProcAddr
        (
            UPstream::listGatherValues<label>(sendData.size()),
            globalIndex::SIZES
        );

        Pout<< "listGather sizes: " << flatOutput(allProcAddr.sizes()) << nl;

        // Collect all values
        labelList allValues
        (
            allProcAddr.mpiGather(sendData)
        );

        Pout<< "all-data: " << allValues << endl;
    }

    {
        const labelList::subList& sendData =
        (
            Pstream::master()
          ? SubList<label>(localValues, 0)  // exclude
          : SubList<label>(localValues)
        );

        const labelList sendSizes
        (
            UPstream::listGatherValues<label>(sendData.size())
        );

        const label sendSize
        (
            UPstream::listScatterValues<label>(sendSizes)
        );

        const globalIndex subProcAddr(sendSizes, globalIndex::SIZES);

        Pout<< "listGather "
            << localValues.size() << " = " << flatOutput(sendSizes)
            << " offsets " << flatOutput(subProcAddr.offsets())
            << nl;

        label newLocalValue = 5 + UPstream::listScatterValues(sendSizes);

        Pout<< "listScattered: " << newLocalValue << nl;

        // Can also scatter a longer list
        Pout<< "listScatter off "
            << UPstream::listScatterValues(subProcAddr.offsets()) << nl;


        Pout<< endl << "local list [" << Pstream::myProcNo() << "] "
            << flatOutput(localValues) << nl;


        Pout<< endl << "local send [" << Pstream::myProcNo() << "] "
            << sendSize << nl;


        // Collect all off-processor values
        labelList allValues
        (
            subProcAddr.mpiGather(sendData)
        );

        Pout<< "off-proc: " << allValues << endl;

        if (Pstream::master())
        {
            Info<< "master: " << flatOutput(localValues) << nl;

            label proci = 0;
            for (const labelRange& range : subProcAddr)
            {
                Info<< proci << ": " << flatOutput(allValues.slice(range)) << nl;
                ++proci;
            }

            Info<< nl << "verify ranges" << nl;

            {
                globalIndex glob;
                Info<< "empty:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
            {
                globalIndex glob(labelList(Foam::one{}, 0), globalIndex::OFFSETS);
                Info<< "degenerate:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
            {
                globalIndex glob(labelList(Foam::one{}, 0), globalIndex::SIZES);
                Info<< "single:" << nl;
                for (const labelRange& range : glob)
                {
                    Info<< "    range: " << range << endl;
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
