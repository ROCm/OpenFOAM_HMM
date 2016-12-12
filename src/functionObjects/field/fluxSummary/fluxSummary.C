/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "fluxSummary.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "Time.H"
#include "syncTools.H"
#include "meshTools.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"
#include "globalIndex.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluxSummary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fluxSummary,
        dictionary
    );
}
template<>
const char* NamedEnum
<
    functionObjects::fluxSummary::modeType,
    3
>::names[] =
{
    "faceZone",
    "faceZoneAndDirection",
    "cellZoneAndDirection"
};
}


const Foam::NamedEnum<Foam::functionObjects::fluxSummary::modeType, 3>
Foam::functionObjects::fluxSummary::modeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fluxSummary::initialiseFaceZone
(
    const word& faceZoneName,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<List<scalar>>& faceSign
) const
{
    label zonei = mesh_.faceZones().findZoneID(faceZoneName);

    if (zonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    faceZoneNames.append(faceZoneName);

    const faceZone& fZone = mesh_.faceZones()[zonei];

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<scalar> faceSigns(fZone.size());

    forAll(fZone, i)
    {
        label facei = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceID = facei;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(facei);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = facei - pp.start();
            }
            else
            {
                faceID = -1;
                facePatchID = -1;
            }
        }

        if (faceID >= 0)
        {
            // Orientation set by faceZone flip map
            if (fZone.flipMap()[facei])
            {
                faceSigns.append(-1);
            }
            else
            {
                faceSigns.append(1);
            }

            faceIDs.append(faceID);
            facePatchIDs.append(facePatchID);
        }
    }

    faceID.append(faceIDs);
    facePatchID.append(facePatchIDs);
    faceSign.append(faceSigns);
}


void Foam::functionObjects::fluxSummary::initialiseFaceZoneAndDirection
(
    const word& faceZoneName,
    const vector& dir,
    DynamicList<vector>& zoneRefDir,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<List<scalar>>& faceSign
) const
{
    vector refDir = dir/(mag(dir) + ROOTVSMALL);

    label zonei = mesh_.faceZones().findZoneID(faceZoneName);

    if (zonei == -1)
    {
         FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    faceZoneNames.append(faceZoneName);
    zoneRefDir.append(refDir);

    const faceZone& fZone = mesh_.faceZones()[zonei];

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<scalar> faceSigns(fZone.size());

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    vector n(Zero);

    forAll(fZone, i)
    {
        label facei = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceID = facei;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(facei);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = facei - pp.start();
            }
            else
            {
                faceID = -1;
                facePatchID = -1;
            }
        }

        if (faceID >= 0)
        {
            // orientation set by comparison with reference direction
            if (facePatchID != -1)
            {
                n = Sf.boundaryField()[facePatchID][faceID]
                   /(magSf.boundaryField()[facePatchID][faceID] + ROOTVSMALL);
            }
            else
            {
                n = Sf[faceID]/(magSf[faceID] + ROOTVSMALL);
            }

            if ((n & refDir) > tolerance_)
            {
                faceSigns.append(1);
            }
            else
            {
                faceSigns.append(-1);
            }

            faceIDs.append(faceID);
            facePatchIDs.append(facePatchID);
        }
    }

    faceID.append(faceIDs);
    facePatchID.append(facePatchIDs);
    faceSign.append(faceSigns);
}


void Foam::functionObjects::fluxSummary::initialiseCellZoneAndDirection
(
    const word& cellZoneName,
    const vector& dir,
    DynamicList<vector>& zoneRefDir,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<List<scalar>>& faceSign
) const
{
    vector refDir = dir/(mag(dir) + ROOTVSMALL);

    const label cellZonei = mesh_.cellZones().findZoneID(cellZoneName);

    if (cellZonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find cellZone " << cellZoneName
            << ". Valid zones are: " << mesh_.cellZones().names()
            << exit(FatalError);
    }

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    labelList cellAddr(mesh_.nCells(), -1);
    const labelList& cellIDs = mesh_.cellZones()[cellZonei];
    UIndirectList<label>(cellAddr, cellIDs) = identity(cellIDs.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct boundary values for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Collect faces
    DynamicList<label> faceIDs(floor(0.1*mesh_.nFaces()));
    DynamicList<label> facePatchIDs(faceIDs.size());
    DynamicList<label> faceLocalPatchIDs(faceIDs.size());
    DynamicList<scalar> faceSigns(faceIDs.size());

    // Internal faces
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if (((own != -1) && (nbr == -1)) || ((own == -1) && (nbr != -1)))
        {
            vector n = mesh_.faces()[facei].normal(mesh_.points());
            n /= mag(n) + ROOTVSMALL;

            if ((n & refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                faceSigns.append(1);
            }
            else if ((n & -refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                faceSigns.append(-1);
            }
        }
    }

    // Loop over boundary faces
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, localFacei)
        {
            const label facei = pp.start() + localFacei;
            const label own = cellAddr[mesh_.faceOwner()[facei]];
            const label nbr = nbrFaceCellAddr[facei - nInternalFaces];

            if ((own != -1) && (nbr == -1))
            {
                vector n = mesh_.faces()[facei].normal(mesh_.points());
                n /= mag(n) + ROOTVSMALL;

                if ((n & refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    faceSigns.append(1);
                }
                else if ((n & -refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    faceSigns.append(-1);
                }
            }
        }
    }

    // Convert into primitivePatch for convenience
    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), faceIDs),
        mesh_.points()
    );

    if (debug)
    {
        OBJstream os(mesh_.time().path()/"patch.obj");
        faceList faces(patch);
        os.write(faces, mesh_.points(), false);
    }


    // Data on all edges and faces
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());

    bool search = true;

    DebugInfo
        << "initialiseCellZoneAndDirection: "
        << "Starting walk to split patch into faceZones"
        << endl;

    globalIndex globalFaces(patch.size());

    label oldFaceID = 0;
    label regioni = 0;
    while (search)
    {
        DynamicList<label> changedEdges;
        DynamicList<patchEdgeFaceRegion> changedInfo;

        label seedFacei = labelMax;
        for (; oldFaceID < patch.size(); oldFaceID++)
        {
            if (allFaceInfo[oldFaceID].region() == -1)
            {
                seedFacei = globalFaces.toGlobal(oldFaceID);
                break;
            }
        }
        reduce(seedFacei, minOp<label>());

        if (seedFacei == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(seedFacei))
        {
            label localFacei = globalFaces.toLocal(seedFacei);
            const labelList& fEdges = patch.faceEdges()[localFacei];

            forAll(fEdges, i)
            {
                if (allEdgeInfo[fEdges[i]].region() != -1)
                {
                    WarningInFunction
                        << "Problem in edge face wave: attempted to assign a "
                        << "value to an edge that has already been visited. "
                        << "Edge info: " << allEdgeInfo[fEdges[i]]
                        << endl;
                }

                changedEdges.append(fEdges[i]);
                changedInfo.append(regioni);
            }
        }


        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchEdgeFaceRegion
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        if (debug)
        {
            label nCells = 0;
            forAll(allFaceInfo, facei)
            {
                if (allFaceInfo[facei].region() == regioni)
                {
                    nCells++;
                }
            }

            Info<< "*** region:" << regioni
                << "  found:" << returnReduce(nCells, sumOp<label>())
                << " faces" << endl;
        }

        regioni++;
    }

    // Collect the data per region
    label nRegion = regioni;

    List<DynamicList<label>> regionFaceIDs(nRegion);
    List<DynamicList<label>> regionFacePatchIDs(nRegion);
    List<DynamicList<scalar>> regionFaceSigns(nRegion);

    forAll(allFaceInfo, facei)
    {
        regioni = allFaceInfo[facei].region();

        regionFaceIDs[regioni].append(faceLocalPatchIDs[facei]);
        regionFacePatchIDs[regioni].append(facePatchIDs[facei]);
        regionFaceSigns[regioni].append(faceSigns[facei]);
    }

    // Transfer to persistent storage
    forAll(regionFaceIDs, regioni)
    {
        const word zoneName = cellZoneName + ":faceZone" + Foam::name(regioni);
        faceZoneNames.append(zoneName);
        zoneRefDir.append(refDir);
        faceID.append(regionFaceIDs[regioni]);
        facePatchID.append(regionFacePatchIDs[regioni]);
        faceSign.append(regionFaceSigns[regioni]);

        // Write OBj of faces to file
        if (debug)
        {
            OBJstream os(mesh_.time().path()/zoneName + ".obj");
            faceList faces(mesh_.faces(), regionFaceIDs[regioni]);
            os.write(faces, mesh_.points(), false);
        }
    }

    if (log)
    {
        Info<< type() << " " << name() << " output:" << nl
            << "    Created " << faceID.size()
            << " separate face zones from cell zone " << cellZoneName << nl;

        forAll(faceZoneNames, i)
        {
            label nFaces = returnReduce(faceID[i].size(), sumOp<label>());
            Info<< "    " << faceZoneNames[i] << ": "
                << nFaces << " faces" << nl;
        }

        Info<< endl;
    }
}


void Foam::functionObjects::fluxSummary::initialiseFaceArea()
{
    faceArea_.setSize(faceID_.size(), 0);

    const surfaceScalarField& magSf = mesh_.magSf();

    forAll(faceID_, zonei)
    {
        const labelList& faceIDs = faceID_[zonei];
        const labelList& facePatchIDs = facePatchID_[zonei];

        scalar sumMagSf = 0;

        forAll(faceIDs, i)
        {
            label facei = faceIDs[i];

            if (facePatchIDs[i] == -1)
            {
                sumMagSf += magSf[facei];
            }
            else
            {
                label patchi = facePatchIDs[i];
                sumMagSf += magSf.boundaryField()[patchi][facei];
            }
        }

        faceArea_[zonei] = returnReduce(sumMagSf, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fluxSummary::fluxSummary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    mode_(mdFaceZone),
    scaleFactor_(1),
    phiName_("phi"),
    faceZoneName_(),
    refDir_(),
    faceID_(),
    facePatchID_(),
    faceSign_(),
    faceArea_(),
    filePtrs_(),
    tolerance_(0.8)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fluxSummary::~fluxSummary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fluxSummary::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    mode_ = modeTypeNames_.read(dict.lookup("mode"));
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    scaleFactor_ = dict.lookupOrDefault<scalar>("scaleFactor", 1.0);
    tolerance_   = dict.lookupOrDefault<scalar>("tolerance", 0.8);

    // Initialise with capacity of 10 faceZones
    DynamicList<vector> refDir(10);
    DynamicList<word> faceZoneName(refDir.size());
    DynamicList<List<label>> faceID(refDir.size());
    DynamicList<List<label>> facePatchID(refDir.size());
    DynamicList<List<scalar>> faceSign(refDir.size());

    switch (mode_)
    {
        case mdFaceZone:
        {
            List<word> zones(dict.lookup("faceZones"));

            forAll(zones, i)
            {
                initialiseFaceZone
                (
                    zones[i],
                    faceZoneName,
                    faceID,
                    facePatchID,
                    faceSign
                );
            }
            break;
        }
        case mdFaceZoneAndDirection:
        {
            List<Tuple2<word, vector>>
                zoneAndDirection(dict.lookup("faceZoneAndDirection"));

            forAll(zoneAndDirection, i)
            {
                initialiseFaceZoneAndDirection
                (
                    zoneAndDirection[i].first(),
                    zoneAndDirection[i].second(),
                    refDir,
                    faceZoneName,
                    faceID,
                    facePatchID,
                    faceSign
                );
            }
            break;
        }
        case mdCellZoneAndDirection:
        {
            List<Tuple2<word, vector>>
                zoneAndDirection(dict.lookup("cellZoneAndDirection"));

            forAll(zoneAndDirection, i)
            {
                initialiseCellZoneAndDirection
                (
                    zoneAndDirection[i].first(),
                    zoneAndDirection[i].second(),
                    refDir,
                    faceZoneName,
                    faceID,
                    facePatchID,
                    faceSign
                );
            }
            break;
        }
        default:
        {
            FatalIOErrorInFunction(dict)
                << "unhandled enumeration " << modeTypeNames_[mode_]
                << abort(FatalIOError);
        }
    }

    faceZoneName_.transfer(faceZoneName);
    refDir_.transfer(refDir);
    faceID_.transfer(faceID);
    facePatchID_.transfer(facePatchID);
    faceSign_.transfer(faceSign);

    initialiseFaceArea();

    if (writeToFile())
    {
        filePtrs_.setSize(faceZoneName_.size());

        forAll(filePtrs_, filei)
        {
            const word& fzName = faceZoneName_[filei];
            filePtrs_.set(filei, createFile(fzName));
            writeFileHeader
            (
                fzName,
                faceArea_[filei],
                refDir_[filei],
                filePtrs_[filei]
            );
        }
    }

    Info<< type() << " " << name() << " output:" << nl;

    forAll(faceZoneName_, zonei)
    {
        const word& zoneName = faceZoneName_[zonei];
        scalar zoneArea = faceArea_[zonei];

        Info<< "    Zone: " << zoneName << ", area: " << zoneArea << nl;
    }

    Info<< endl;

    return true;
}


void Foam::functionObjects::fluxSummary::writeFileHeader
(
    const word& fzName,
    const scalar area,
    const vector& refDir,
    Ostream& os
) const
{
    writeHeader(os, "Flux summary");
    writeHeaderValue(os, "Face zone", fzName);
    writeHeaderValue(os, "Total area", area);

    switch (mode_)
    {
        case mdFaceZoneAndDirection:
        case mdCellZoneAndDirection:
        {
            writeHeaderValue(os, "Reference direction", refDir);
            break;
        }
        default:
        {}
    }

    writeHeaderValue(os, "Scale factor", scaleFactor_);

    writeCommented(os, "Time");
    os  << tab << "positive"
        << tab << "negative"
        << tab << "net"
        << tab << "absolute"
        << endl;
}


bool Foam::functionObjects::fluxSummary::execute()
{
    return true;
}


bool Foam::functionObjects::fluxSummary::write()
{
    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);

    word flowType;
    if (phi.dimensions() == dimVolume/dimTime)
    {
        flowType = "volumetric";
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        flowType = "mass";
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported flux field " << phi.name() << " with dimensions "
            << phi.dimensions() << ".  Expected either mass flow or volumetric "
            << "flow rate" << abort(FatalError);
    }

    Log << type() << " " << name() << ' ' << flowType << " write:" << nl;

    forAll(faceZoneName_, zonei)
    {
        const labelList& faceID = faceID_[zonei];
        const labelList& facePatchID = facePatchID_[zonei];
        const scalarList& faceSign = faceSign_[zonei];

        scalar phiPos = scalar(0);
        scalar phiNeg = scalar(0);
        scalar phif = scalar(0);

        forAll(faceID, i)
        {
            label facei = faceID[i];
            label patchi = facePatchID[i];

            if (patchi != -1)
            {
                phif = phi.boundaryField()[patchi][facei];
            }
            else
            {
                phif = phi[facei];
            }

            phif *= faceSign[i];

            if (phif > 0)
            {
                phiPos += phif;
            }
            else
            {
                phiNeg += phif;
            }
        }

        reduce(phiPos, sumOp<scalar>());
        reduce(phiNeg, sumOp<scalar>());

        phiPos *= scaleFactor_;
        phiNeg *= scaleFactor_;

        scalar netFlux = phiPos + phiNeg;
        scalar absoluteFlux = phiPos - phiNeg;

        Log << "    faceZone " << faceZoneName_[zonei] << ':' << nl
            << "        positive : " << phiPos << nl
            << "        negative : " << phiNeg << nl
            << "        net      : " << netFlux << nl
            << "        absolute : " << absoluteFlux
            << nl << endl;

        if (writeToFile())
        {
            filePtrs_[zonei]
                << time_.value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
