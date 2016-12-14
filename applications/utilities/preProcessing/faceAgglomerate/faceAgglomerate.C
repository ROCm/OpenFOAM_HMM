/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    faceAgglomerate

Group
    grpPreProcessingUtilities

Description
    Agglomerate boundary faces using the pairPatchAgglomeration algorithm.

    It writes a map from the fine to coarse grid.

SeeAlso
    pairPatchAgglomeration.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "CompactListList.H"
#include "unitConversion.H"
#include "pairPatchAgglomeration.H"
#include "labelListIOList.H"
#include "syncTools.H"
#include "globalIndex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("viewFactorsDict");

    #include "setConstantMeshDictionaryIO.H"

    // Read control dictionary
    const IOdictionary agglomDict(dictIO);

    bool writeAgglom = readBool(agglomDict.lookup("writeFacesAgglomeration"));

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        boundary.size()
    );

    label nCoarseFaces = 0;

    forAllConstIter(dictionary, agglomDict, iter)
    {
        labelList patchids = boundary.findIndices(iter().keyword());
        forAll(patchids, i)
        {
            label patchi =  patchids[i];
            const polyPatch& pp = boundary[patchi];

            if (!pp.coupled())
            {
                Info << "\nAgglomerating patch : " << pp.name() << endl;

                pairPatchAgglomeration agglomObject
                (
                    pp.localFaces(),
                    pp.localPoints(),
                    agglomDict.subDict(pp.name())
                );

                agglomObject.agglomerate();

                finalAgglom[patchi] =
                    agglomObject.restrictTopBottomAddressing();

                if (finalAgglom[patchi].size())
                {
                    nCoarseFaces += max(finalAgglom[patchi] + 1);
                }
            }
        }
    }


    // All patches which are not agglomerated are identity for finalAgglom
    forAll(boundary, patchi)
    {
        if (finalAgglom[patchi].size() == 0)
        {
            finalAgglom[patchi] = identity(boundary[patchi].size());
        }
    }

    // Sync agglomeration across coupled patches
    labelList nbrAgglom(mesh.nFaces() - mesh.nInternalFaces(), -1);

    forAll(boundary, patchi)
    {
        const polyPatch& pp = boundary[patchi];
        if (pp.coupled())
        {
            finalAgglom[patchi] = identity(pp.size());
            forAll(pp, i)
            {
                const label agglomi = pp.start() - mesh.nInternalFaces() + i;
                nbrAgglom[agglomi] = finalAgglom[patchi][i];
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, nbrAgglom);
    forAll(boundary, patchi)
    {
        const polyPatch& pp = boundary[patchi];
        if (pp.coupled() && !refCast<const coupledPolyPatch>(pp).owner())
        {
            forAll(pp, i)
            {
                const label agglomi = pp.start() - mesh.nInternalFaces() + i;
                finalAgglom[patchi][i] = nbrAgglom[agglomi];
            }
        }
    }

    finalAgglom.write();

    if (writeAgglom)
    {
        globalIndex index(nCoarseFaces);
        volScalarField facesAgglomeration
        (
            IOobject
            (
                "facesAgglomeration",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("facesAgglomeration", dimless, 0)
        );

        volScalarField::Boundary& facesAgglomerationBf =
            facesAgglomeration.boundaryFieldRef();

        label coarsePatchIndex = 0;
        forAll(boundary, patchi)
        {
            const polyPatch& pp = boundary[patchi];
            if (pp.size() > 0)
            {
                fvPatchScalarField& bFacesAgglomeration =
                    facesAgglomerationBf[patchi];

                forAll(bFacesAgglomeration, j)
                {
                    bFacesAgglomeration[j] =
                        index.toGlobal
                        (
                            Pstream::myProcNo(),
                            finalAgglom[patchi][j] + coarsePatchIndex
                        );
                }

                coarsePatchIndex += max(finalAgglom[patchi]) + 1;
            }
        }

        Info<< "\nWriting facesAgglomeration" << endl;
        facesAgglomeration.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
