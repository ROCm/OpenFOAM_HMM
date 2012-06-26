/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    patchSummary

Description
    Writes fields and boundary condition info for each patch at each requested
    time instance.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "patchSummaryTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"
    argList::addBoolOption
    (
        "collate",
        "Combine similar patches"
    );
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool collate = args.optionFound("collate");


#   include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const IOobjectList fieldObjs(mesh, runTime.timeName());
        const wordList objNames = fieldObjs.names();

        PtrList<volScalarField> vsf(objNames.size());
        PtrList<volVectorField> vvf(objNames.size());
        PtrList<volSphericalTensorField> vsptf(objNames.size());
        PtrList<volSymmTensorField> vsytf(objNames.size());
        PtrList<volTensorField> vtf(objNames.size());

        Info<< "Valid fields:" << endl;

        forAll(objNames, objI)
        {
            IOobject obj
            (
                objNames[objI],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (obj.headerOk())
            {
                addToFieldList<scalar>(vsf, obj, objI, mesh);
                addToFieldList<vector>(vvf, obj, objI, mesh);
                addToFieldList<sphericalTensor>(vsptf, obj, objI, mesh);
                addToFieldList<symmTensor>(vsytf, obj, objI, mesh);
                addToFieldList<tensor>(vtf, obj, objI, mesh);
            }
        }

        Info<< endl;

        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        DynamicList<HashTable<word> > fieldToTypes(bm.size());
        DynamicList<DynamicList<label> > groupToPatches(bm.size());
        forAll(bm, patchI)
        {
            if (collate)
            {
                HashTable<word> fieldToType;
                collectFieldList<scalar>(vsf, patchI, fieldToType);
                collectFieldList<vector>(vvf, patchI, fieldToType);
                collectFieldList<sphericalTensor>(vsptf, patchI, fieldToType);
                collectFieldList<symmTensor>(vsytf, patchI, fieldToType);
                collectFieldList<tensor>(vtf, patchI, fieldToType);

                label groupI = findIndex(fieldToTypes, fieldToType);
                if (groupI == -1)
                {
                    DynamicList<label> group(1);
                    group.append(patchI);
                    groupToPatches.append(group);
                    fieldToTypes.append(fieldToType);
                }
                else
                {
                    groupToPatches[groupI].append(patchI);
                }
            }
            else
            {
                Info<< bm[patchI].type() << ": " << bm[patchI].name() << nl;
                outputFieldList<scalar>(vsf, patchI);
                outputFieldList<vector>(vvf, patchI);
                outputFieldList<sphericalTensor>(vsptf, patchI);
                outputFieldList<symmTensor>(vsytf, patchI);
                outputFieldList<tensor>(vtf, patchI);
                Info<< endl;
            }
        }


        bool hasGroups = false;
        forAll(groupToPatches, groupI)
        {
            const DynamicList<label>& patchIDs = groupToPatches[groupI];
            if (patchIDs.size() > 1)
            {
                hasGroups = true;
                break;
            }
        }


        if (hasGroups)
        {
            Info<< "Collated:" << endl;
            forAll(groupToPatches, groupI)
            {
                const DynamicList<label>& patchIDs = groupToPatches[groupI];

                if (patchIDs.size() > 1)
                {
                    Info<< '(' << bm[patchIDs[0]].name();
                    for (label i = 1; i < patchIDs.size(); i++)
                    {
                        Info<< ' ' << bm[patchIDs[i]].name();
                    }
                    Info<< ')' << endl;
                    outputFieldList<scalar>(vsf, patchIDs[0]);
                    outputFieldList<vector>(vvf, patchIDs[0]);
                    outputFieldList<sphericalTensor>(vsptf, patchIDs[0]);
                    outputFieldList<symmTensor>(vsytf, patchIDs[0]);
                    outputFieldList<tensor>(vtf, patchIDs[0]);
                    Info<< endl;
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
