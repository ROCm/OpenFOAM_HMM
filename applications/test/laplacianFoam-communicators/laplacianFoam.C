/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Determine parent rank
label processorToRank
(
    const label comm,
    const label procI
)
{
    const List<int>& parentRanks = UPstream::procID(comm);
    label parentComm = UPstream::parent(comm);

    if (parentComm == -1)
    {
        return findIndex(parentRanks, procI);
    }
    else
    {
        label parentRank = processorToRank(parentComm, procI);
        return findIndex(parentRanks, parentRank);
    }
}


// Determine rank in new communicator
label rank
(
    const label currentComm,
    const label currentRank,
    const label newComm
)
{
    // Determine actual processor (walk up the tree to the top)
    label procID = currentRank;
    label comm = currentComm;

    while (UPstream::parent(comm) != -1)
    {
        const List<int>& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = UPstream::parent(comm);
    }

    //Pout<< "Walked to top : comm:" << comm << " procID:" << procID
    //    << endl;

    // Find out what procID in the newComm is.
    return processorToRank(newComm, procID);
}


void setCommunicator(fvMesh& mesh, const label newComm)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // The current state is consistent with the mesh so check where the new
    // communicator is and adjust accordingly.

    forAll(pbm, patchI)
    {
        if (isA<processorPolyPatch>(pbm[patchI]))
        {
            processorPolyPatch& ppp = const_cast<processorPolyPatch&>
            (
                refCast
                <
                    const processorPolyPatch
                >(pbm[patchI])
            );

            label thisRank = rank(ppp.comm(), ppp.myProcNo(), newComm);
            label nbrRank = rank(ppp.comm(), ppp.neighbProcNo(), newComm);

            ppp.comm() = newComm;
            ppp.myProcNo() = thisRank;
            ppp.neighbProcNo() = nbrRank;
        }
    }
    mesh.polyMesh::comm() = newComm;
}


void printIt(Ostream& os)
{
    os<< "**BLA**" << endl;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    simpleControl simple(mesh);

    //const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;


    // Get a subset of processors
    labelList subProcs(2);
    subProcs[0] = 1;
    subProcs[1] = 2;


    const UPstream::communicator newComm
    (
        UPstream::worldComm,
        subProcs,
        true
    );


    Pout<< "procIDs world  :" << UPstream::procID(UPstream::worldComm) << endl;
    Pout<< "procIDs newComm:" << UPstream::procID(newComm) << endl;


    {
        label oldComm = mesh.comm();
        setCommunicator(mesh, newComm);
        Pout<< "** oldcomm:" << oldComm
            << "  newComm:" << mesh.comm() << endl;

        printIt(Info(mesh.comm()));

        setCommunicator(mesh, oldComm);
        Pout<< "** reset mesh to:" << mesh.comm() << endl;
    }


    //{
    //    Pout<< "World:" << UPstream::worldComm
    //        << " procID:" << 2
    //        << " subComm:" << newComm
    //        << " rank:" << rank(UPstream::worldComm, 2, newComm)
    //        << endl;
    //}
    //
    //
    //{
    //    setCommunicator(mesh, newComm);
    //    Pout<< "Mesh comm:" << mesh.comm() << endl;
    //    forAll(pbm, patchI)
    //    {
    //        if (isA<processorPolyPatch>(pbm[patchI]))
    //        {
    //            const processorPolyPatch& ppp =
    //            refCast
    //            <
    //                const processorPolyPatch
    //            >(pbm[patchI]);
    //
    //
    //            Pout<< "    " << ppp.name()
    //                << " myRank:" << ppp.myProcNo()
    //                << " nbrRank:" << ppp.neighbProcNo()
    //                << endl;
    //        }
    //    }
    //    setCommunicator(mesh, UPstream::worldComm);
    //}


    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix Teqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
            );


            {
                label oldWarn = UPstream::warnComm;
                UPstream::warnComm = newComm;

                label oldComm = mesh.comm();
                setCommunicator(mesh, newComm);
                Pout<< "** oldcomm:" << oldComm
                    << "  newComm:" << mesh.comm() << endl;

                if (Pstream::myProcNo(mesh.comm()) != -1)
                {
                    solve(Teqn);
                }

                setCommunicator(mesh, oldComm);
                Pout<< "** reset mesh to:" << mesh.comm() << endl;

                UPstream::warnComm = oldWarn;
            }
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
