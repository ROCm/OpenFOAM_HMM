/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxSummary, 0);

    template<>
    const char* NamedEnum
    <
        fluxSummary::modeType,
        3
    >::names[] =
    {
        "faceZone",
        "faceZoneAndDirection",
        "cellZoneAndDirection"
    };
}


const Foam::NamedEnum<Foam::fluxSummary::modeType, 3>
Foam::fluxSummary::modeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fluxSummary::initialiseFaceZone
(
    const word& faceZoneName,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label> >& faceID,
    DynamicList<List<label> >& facePatchID,
    DynamicList<List<scalar> >& faceSign
) const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    label zoneI = mesh.faceZones().findZoneID(faceZoneName);

    if (zoneI == -1)
    {
        FatalErrorIn
        (
            "void  Foam::fluxSummary::initialiseFaceZone"
            "("
                "const word&, "
                "DynamicList<word>&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<scalar> >&"
            ") const"
        )
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh.faceZones().names()
            << exit(FatalError);
    }

    faceZoneNames.append(faceZoneName);

    const faceZone& fZone = mesh.faceZones()[zoneI];

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<scalar> faceSigns(fZone.size());

    forAll(fZone, i)
    {
        label faceI = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh.isInternalFace(faceI))
        {
            faceID = faceI;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(faceI);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = faceI - pp.start();
            }
            else
            {
                faceID = -1;
                facePatchID = -1;
            }
        }

        if (faceID >= 0)
        {
            // orientation set by faceZone flip map
            if (fZone.flipMap()[faceI])
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


void Foam::fluxSummary::initialiseFaceZoneAndDirection
(
    const word& faceZoneName,
    const vector& dir,
    DynamicList<vector>& zoneRefDir,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label> >& faceID,
    DynamicList<List<label> >& facePatchID,
    DynamicList<List<scalar> >& faceSign
) const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    vector refDir = dir/(mag(dir) + ROOTVSMALL);

    label zoneI = mesh.faceZones().findZoneID(faceZoneName);

    if (zoneI == -1)
    {
        FatalErrorIn
        (
            "void  Foam::fluxSummary::initialiseFaceZoneAndDirection"
            "("
                "const word&, "
                "const vector&, "
                "DynamicList<vector>&, "
                "DynamicList<word>&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<scalar> >&"
            ") const"
        )
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh.faceZones().names()
            << exit(FatalError);
    }

    faceZoneNames.append(faceZoneName);
    zoneRefDir.append(refDir);

    const faceZone& fZone = mesh.faceZones()[zoneI];

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<scalar> faceSigns(fZone.size());

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    vector n = vector::zero;

    forAll(fZone, i)
    {
        label faceI = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh.isInternalFace(faceI))
        {
            faceID = faceI;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(faceI);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = faceI - pp.start();
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


void Foam::fluxSummary::initialiseCellZoneAndDirection
(
    const word& cellZoneName,
    const vector& dir,
    DynamicList<vector>& zoneRefDir,
    DynamicList<word>& faceZoneNames,
    DynamicList<List<label> >& faceID,
    DynamicList<List<label> >& facePatchID,
    DynamicList<List<scalar> >& faceSign
) const
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    vector refDir = dir/(mag(dir) + ROOTVSMALL);

    const label cellZoneI = mesh.cellZones().findZoneID(cellZoneName);

    if (cellZoneI == -1)
    {
        FatalErrorIn
        (
            "void Foam::fluxSummary::initialiseCellZoneAndDirection"
            "("
                "const word&, "
                "const vector&, "
                "DynamicList<vector>&, "
                "DynamicList<word>&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<label> >&, "
                "DynamicList<List<scalar> >&"
            ") const"
        )
            << "Unable to find cellZone " << cellZoneName
            << ". Valid zones are: " << mesh.cellZones().names()
            << exit(FatalError);
    }

    const label nInternalFaces = mesh.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList cellAddr(mesh.nCells(), -1);
    const labelList& cellIDs = mesh.cellZones()[cellZoneI];
    UIndirectList<label>(cellAddr, cellIDs) = identity(cellIDs.size());
    labelList nbrFaceCellAddr(mesh.nFaces() - nInternalFaces, -1);

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label nbrFaceI = faceI - nInternalFaces;
                label own = mesh.faceOwner()[faceI];
                nbrFaceCellAddr[nbrFaceI] = cellAddr[own];
            }
        }
    }

    // correct boundary values for parallel running
    syncTools::swapBoundaryFaceList(mesh, nbrFaceCellAddr);

    // collect faces
    DynamicList<label> faceIDs(floor(0.1*mesh.nFaces()));
    DynamicList<label> facePatchIDs(faceIDs.size());
    DynamicList<label> faceLocalPatchIDs(faceIDs.size());
    DynamicList<scalar> faceSigns(faceIDs.size());

    // internal faces
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const label own = cellAddr[mesh.faceOwner()[faceI]];
        const label nbr = cellAddr[mesh.faceNeighbour()[faceI]];

        if (((own != -1) && (nbr == -1)) || ((own == -1) && (nbr != -1)))
        {
            vector n = mesh.faces()[faceI].normal(mesh.points());
            n /= mag(n) + ROOTVSMALL;

            if ((n & refDir) > tolerance_)
            {
                faceIDs.append(faceI);
                faceLocalPatchIDs.append(faceI);
                facePatchIDs.append(-1);
                faceSigns.append(1);
            }
            else if ((n & -refDir) > tolerance_)
            {
                faceIDs.append(faceI);
                faceLocalPatchIDs.append(faceI);
                facePatchIDs.append(-1);
                faceSigns.append(-1);
            }
        }
    }

    // loop of boundary faces
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        forAll(pp, localFaceI)
        {
            const label faceI = pp.start() + localFaceI;
            const label own = cellAddr[mesh.faceOwner()[faceI]];
            const label nbr = nbrFaceCellAddr[faceI - nInternalFaces];

            if ((own != -1) && (nbr == -1))
            {
                vector n = mesh.faces()[faceI].normal(mesh.points());
                n /= mag(n) + ROOTVSMALL;

                if ((n & refDir) > tolerance_)
                {
                    faceIDs.append(faceI);
                    faceLocalPatchIDs.append(localFaceI);
                    facePatchIDs.append(patchI);
                    faceSigns.append(1);
                }
                else if ((n & -refDir) > tolerance_)
                {
                    faceIDs.append(faceI);
                    faceLocalPatchIDs.append(localFaceI);
                    facePatchIDs.append(patchI);
                    faceSigns.append(-1);
                }
            }
        }
    }

    // convert into primitivePatch for convenience
    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh.faces(), faceIDs),
        mesh.points()
    );

    if (debug)
    {
        OBJstream os(mesh.time().path()/"patch.obj");
        faceList faces(patch);
        os.write(faces, mesh.points(), false);
    }


    // data on all edges and faces
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());

    bool search = true;

    if (debug)
    {
        Info<< "initialiseCellZoneAndDirection: "
            << "Starting walk to split patch into faceZones"
            << endl;
    }

    globalIndex globalFaces(patch.size());

    label oldFaceID = 0;
    label regionI = 0;
    while (search)
    {
        DynamicList<label> changedEdges;
        DynamicList<patchEdgeFaceRegion> changedInfo;

        label seedFaceI = labelMax;
        for (; oldFaceID < patch.size(); oldFaceID++)
        {
            if (allFaceInfo[oldFaceID].region() == -1)
            {
                seedFaceI = globalFaces.toGlobal(oldFaceID);
                break;
            }
        }
        reduce(seedFaceI, minOp<label>());

        if (seedFaceI == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(seedFaceI))
        {
            label localFaceI = globalFaces.toLocal(seedFaceI);
            const labelList& fEdges = patch.faceEdges()[localFaceI];

            forAll(fEdges, i)
            {
                if (allEdgeInfo[fEdges[i]].region() != -1)
                {
                    WarningIn
                    (
                        "void Foam::fluxSummary::initialiseCellZoneAndDirection"
                        "("
                            "const word&, "
                            "const vector&, "
                            "DynamicList<vector>&, "
                            "DynamicList<word>&, "
                            "DynamicList<List<label> >&, "
                            "DynamicList<List<label> >&, "
                            "DynamicList<List<scalar> >&"
                        ") const"
                    )   << "Problem in edge face wave: attempted to assign a "
                        << "value to an edge that has already been visited. "
                        << "Edge info: " << allEdgeInfo[fEdges[i]]
                        << endl;
                }

                changedEdges.append(fEdges[i]);
                changedInfo.append(regionI);
            }
        }


        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchEdgeFaceRegion
        > calc
        (
            mesh,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        label nCells = 0;
        forAll(allFaceInfo, faceI)
        {
            if (allFaceInfo[faceI].region() == regionI)
            {
                nCells++;
            }
        }

        if (debug)
        {
            Info<< "*** region:" << regionI
                << "  found:" << returnReduce(nCells, sumOp<label>())
                << " faces" << endl;
        }

        regionI++;
    }

    // collect the data per region
    label nRegion = regionI;

    List<DynamicList<label> > regionFaceIDs(nRegion);
    List<DynamicList<label> > regionFacePatchIDs(nRegion);
    List<DynamicList<scalar> > regionFaceSigns(nRegion);

    forAll(allFaceInfo, faceI)
    {
        regionI = allFaceInfo[faceI].region();

        regionFaceIDs[regionI].append(faceLocalPatchIDs[faceI]);
        regionFacePatchIDs[regionI].append(facePatchIDs[faceI]);
        regionFaceSigns[regionI].append(faceSigns[faceI]);
    }

    // transfer to persistent storage
    forAll(regionFaceIDs, regionI)
    {
        const word zoneName = cellZoneName + ":faceZone" + Foam::name(regionI);
        faceZoneNames.append(zoneName);
        zoneRefDir.append(refDir);
        faceID.append(regionFaceIDs[regionI]);
        facePatchID.append(regionFacePatchIDs[regionI]);
        faceSign.append(regionFaceSigns[regionI]);

        // write OBJ of faces to file
        if (debug)
        {
            OBJstream os(mesh.time().path()/zoneName + ".obj");
            faceList faces(mesh.faces(), regionFaceIDs[regionI]);
            os.write(faces, mesh.points(), false);
        }
    }

    if (log_)
    {
        Info<< type() << " " << name_ << " output:" << nl
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


void Foam::fluxSummary::initialiseFaceArea()
{
    faceArea_.setSize(faceID_.size(), 0);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(faceID_, zoneI)
    {
        const labelList& faceIDs = faceID_[zoneI];
        const labelList& facePatchIDs = facePatchID_[zoneI];

        scalar sumMagSf = 0;

        forAll(faceIDs, i)
        {
            label faceI = faceIDs[i];

            if (facePatchIDs[i] == -1)
            {
                sumMagSf += magSf[faceI];
            }
            else
            {
                label patchI = facePatchIDs[i];
                sumMagSf += magSf.boundaryField()[patchI][faceI];
            }
        }

        faceArea_[zoneI] = returnReduce(sumMagSf, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSummary::fluxSummary
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
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
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "fluxSummary::fluxSummary"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSummary::~fluxSummary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSummary::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_ = dict.lookupOrDefault<Switch>("log", true);

        mode_ = modeTypeNames_.read(dict.lookup("mode"));
        phiName_= dict.lookupOrDefault<word>("phiName", "phi");
        dict.readIfPresent("scaleFactor", scaleFactor_);
        dict.readIfPresent("tolerance", tolerance_);

        // initialise with capacity of 10 faceZones
        DynamicList<vector> refDir(10);
        DynamicList<word> faceZoneName(refDir.size());
        DynamicList<List<label> > faceID(refDir.size());
        DynamicList<List<label> > facePatchID(refDir.size());
        DynamicList<List<scalar> > faceSign(refDir.size());

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
                List<Tuple2<word, vector> >
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
                List<Tuple2<word, vector> >
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
                FatalIOErrorIn
                (
                    "void Foam::fluxSummary::read(const dictionary&)",
                    dict
                )
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

            forAll(filePtrs_, fileI)
            {
                const word& fzName = faceZoneName_[fileI];
                filePtrs_.set(fileI, createFile(fzName));
                writeFileHeader
                (
                    fzName,
                    faceArea_[fileI],
                    refDir_[fileI],
                    filePtrs_[fileI]
                );
            }
        }
    }
}


void Foam::fluxSummary::writeFileHeader
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


void Foam::fluxSummary::execute()
{
    // Do nothing - only valid on write
}


void Foam::fluxSummary::end()
{
    // Do nothing - only valid on write
}


void Foam::fluxSummary::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::fluxSummary::write()
{
    if (!active_)
    {
        return;
    }

    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    word flowType = "";
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
        FatalErrorIn("void Foam::fluxSummary::write()")
            << "Unsupported flux field " << phi.name() << " with dimensions "
            << phi.dimensions() << ".  Expected eithe mass flow or volumetric "
            << "flow rate" << abort(FatalError);
    }

    if (log_)
    {
        Info<< type() << " " << name_ << ' ' << flowType << " flux output:"
            << nl;
    }

    forAll(faceZoneName_, zoneI)
    {
        const labelList& faceID = faceID_[zoneI];
        const labelList& facePatchID = facePatchID_[zoneI];
        const scalarList& faceSign = faceSign_[zoneI];

        scalar phiPos = scalar(0);
        scalar phiNeg = scalar(0);
        scalar phif = scalar(0);

        forAll(faceID, i)
        {
            label faceI = faceID[i];
            label patchI = facePatchID[i];

            if (patchI != -1)
            {
                phif = phi.boundaryField()[patchI][faceI];
            }
            else
            {
                phif = phi[faceI];
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

        if (log_)
        {
            Info<< "    faceZone " << faceZoneName_[zoneI] << ':' << nl
                << "        positive : " << phiPos << nl
                << "        negative : " << phiNeg << nl
                << "        net      : " << netFlux << nl
                << "        absolute : " << absoluteFlux
                << nl << endl;
        }

        if (writeToFile())
        {
            filePtrs_[zoneI]
                << obr_.time().value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    if (log_) Info<< endl;
}



// ************************************************************************* //
