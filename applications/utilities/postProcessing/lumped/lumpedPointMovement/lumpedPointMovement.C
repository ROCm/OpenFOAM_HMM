/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    lumpedPointMovement

Description
    Thus utility can be used to produce VTK files to visualize the response
    points/rotations and the corresponding movement of the building surfaces.
    Uses the tabulated responses from the specified file.
    Optionally, it can also be used to a dummy responder for the
    externalFileCoupler logic, which makes it useful as a debugging facility
    as well demonstrating how an external application could communicate
    with the lumpedPointMovement point-patch boundary condition.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"

#include "OFstream.H"

#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Visualize lumpedPoint movements or provide a slave responder "
        "for diagnostic purposes."
    );

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects
    argList::addOption
    (
        "max",
        "N",
        "maximum number of outputs"
    );
    argList::addOption
    (
        "span",
        "N",
        "increment each input by factor N (default: 1)"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "relaxation/scaling factor for movement (default: 1)"
    );
    argList::addBoolOption
    (
        "removeLock",
        "remove lock-file on termination of slave"
    );
    argList::addBoolOption
    (
        "slave",
        "invoke as a slave responder for testing"
    );
    argList::addArgument("responseFile");

    #include "setRootCase.H"

    const label maxOut =
        Foam::max(0, args.optionLookupOrDefault<label>("max", 0));

    const label span =
        Foam::max(1, args.optionLookupOrDefault<label>("span", 1));

    const scalar relax = args.optionLookupOrDefault<scalar>("scale", 1);

    const bool slave = args.optionFound("slave");
    const bool removeLock = args.optionFound("removeLock");

    #include "createTime.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New
    (
        runTime
    );

    if (!movement.valid())
    {
        Info<< "no valid movement given" << endl;
        return 1;
    }

    List<lumpedPointStateTuple> responseTable =
        lumpedPointTools::lumpedPointStates(args[1]);

    Info<< "Using response table with " << responseTable.size()
        << " entries" << endl;

    Info << "Increment input by " << span << nl;

    if (maxOut)
    {
        Info<< "Stopping after " << maxOut << " outputs" << endl;
    }

    if (slave)
    {
        Info<< "Running as slave responder" << endl;

        externalFileCoupler& coupler = movement().coupler();

        label count = 0;
        for (label index = 0; index < responseTable.size(); index += span)
        {
            Info<< args.executable() << ": waiting for master" << endl;

            // Wait for master, but stop if status=done was seen
            if (!coupler.waitForMaster())
            {
                Info<< args.executable()
                    << ": stopping status=done was detected" << endl;
                break;
            }

            lumpedPointState state = responseTable[index].second();
            state.relax(relax, movement().state0());

            // Generate input for OpenFOAM
            OFstream os(coupler.resolveFile(movement().inputName()));
            if
            (
                movement().inputFormat()
             == lumpedPointState::inputFormatType::PLAIN
            )
            {
                state.writePlain(os);
            }
            else
            {
                os.writeEntry("time", responseTable[index].first());
                state.writeDict(os);
            }

            Info<< args.executable()
                << ": updated to state " << index
                << " - switch to master"
                << endl;

            // Let OpenFOAM know that it can continue
            coupler.useMaster();

            if (maxOut && ++count >= maxOut)
            {
                Info<< args.executable()
                    <<": stopping after " << maxOut << " outputs" << endl;
                break;
            }
        }

        if (removeLock)
        {
            Info<< args.executable() <<": removing lock file" << endl;
            coupler.useSlave();  // This removes the lock-file
        }
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);

        #include "createNamedPolyMesh.H"

        const labelList patchLst = lumpedPointTools::lumpedPointPatchList(mesh);
        if (patchLst.empty())
        {
            Info<< "no patch list found" << endl;
            return 2;
        }

        pointIOField points0 = lumpedPointTools::points0Field(mesh);
        movement().setBoundBox(mesh, patchLst, points0);

        label index = 0;

        // Initial geometry
        movement().writeVTP("geom_init.vtp", mesh, patchLst, points0);

        forAll(responseTable, i)
        {
            const bool output = ((i % span) == 0);
            lumpedPointState state = responseTable[i].second();
            state.relax(relax, movement().state0());

            if (output)
            {
                Info<<"output [" << i << "/"
                    << responseTable.size() << "]" << endl;
            }
            else
            {
                continue;
            }

            // State/response = what comes back from FEM
            {
                const word outputName = Foam::name("state_%06d.vtp", index);
                Info<<"    " << outputName << endl;

                state.writeVTP(outputName, movement().axis());
            }

            {
                const word outputName = Foam::name("geom_%06d.vtp", index);
                Info<<"    " << outputName << endl;

                movement().writeVTP(outputName, state, mesh, patchLst, points0);
            }

            {
                ++index;

                bool canOutput = !maxOut || (index <= maxOut);
                if (!canOutput)
                {
                    Info<<"stopping output after "
                        << maxOut << " outputs" << endl;
                    break;
                }
            }
        }

    }

    Info<< args.executable() << ": End\n" << endl;

    return 0;
}


// ************************************************************************* //
