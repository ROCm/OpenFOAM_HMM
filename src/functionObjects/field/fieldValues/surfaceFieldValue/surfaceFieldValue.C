/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "surfaceFieldValue.H"
#include "fvMesh.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"
#include "mergePoints.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "meshedSurfRef.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(surfaceFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, surfaceFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, surfaceFieldValue, dictionary);
}
}
}


const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes
>
Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypeNames_
{
    { regionTypes::stFaceZone, "faceZone" },
    { regionTypes::stPatch, "patch" },
    { regionTypes::stSurface, "surface" },
    { regionTypes::stSampledSurface, "sampledSurface" },
};


const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType
>
Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_
{
    // Normal operations
    { operationType::opNone, "none" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opSum, "sum" },
    { operationType::opSumMag, "sumMag" },
    { operationType::opSumDirection, "sumDirection" },
    { operationType::opSumDirectionBalance, "sumDirectionBalance" },
    { operationType::opAverage, "average" },
    { operationType::opAreaAverage, "areaAverage" },
    { operationType::opAreaIntegrate, "areaIntegrate" },
    { operationType::opCoV, "CoV" },
    { operationType::opAreaNormalAverage, "areaNormalAverage" },
    { operationType::opAreaNormalIntegrate, "areaNormalIntegrate" },

    // Using weighting
    { operationType::opWeightedSum, "weightedSum" },
    { operationType::opWeightedAverage, "weightedAverage" },
    { operationType::opWeightedAreaAverage, "weightedAreaAverage" },
    { operationType::opWeightedAreaIntegrate, "weightedAreaIntegrate" },

    // Using absolute weighting
    { operationType::opAbsWeightedSum, "absWeightedSum" },
    { operationType::opAbsWeightedAverage, "absWeightedAverage" },
    { operationType::opAbsWeightedAreaAverage, "absWeightedAreaAverage" },
    { operationType::opAbsWeightedAreaIntegrate, "absWeightedAreaIntegrate" },
};

const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationType
>
Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationTypeNames_
{
    { postOperationType::postOpNone, "none" },
    { postOperationType::postOpSqrt, "sqrt" },
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::fieldValues::surfaceFieldValue::obr() const
{
    if (regionType_ == stSurface)
    {
        return mesh_.lookupObject<objectRegistry>(regionName_);
    }
    else
    {
        return mesh_;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces()
{
    const label zoneId = mesh_.faceZones().findZoneID(regionName_);

    if (zoneId < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown face zone name: " << regionName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneId];

    DynamicList<label> faceIds(fZone.size());
    DynamicList<label> facePatchIds(fZone.size());
    DynamicList<bool> faceFlip(fZone.size());

    forAll(fZone, i)
    {
        const label facei = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            faceIds.append(faceId);
            facePatchIds.append(facePatchId);
            faceFlip.append(fZone.flipMap()[i] ? true : false);
        }
    }

    faceId_.transfer(faceIds);
    facePatchId_.transfer(facePatchIds);
    faceFlip_.transfer(faceFlip);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << faceId_.size() << endl;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchFaces()
{
    const label patchid = mesh_.boundaryMesh().findPatchID(regionName_);

    if (patchid < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Unknown patch name: " << regionName_
            << ". Valid patch names are: "
            << mesh_.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh_.boundaryMesh()[patchid];

    label nFaces = pp.size();
    if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceFlip_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, facei)
    {
        faceId_[facei] = facei;
        facePatchId_[facei] = patchid;
        faceFlip_[facei] = false;
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    List<faceList> allFaces(Pstream::nProcs());
    List<pointField> allPoints(Pstream::nProcs());

    labelList globalFacesIs(faceId_);
    forAll(globalFacesIs, i)
    {
        if (facePatchId_[i] != -1)
        {
            const label patchi = facePatchId_[i];
            globalFacesIs[i] += mesh_.boundaryMesh()[patchi].start();
        }
    }

    // Add local faces and points to the all* lists
    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh_.faces(), globalFacesIs),
        mesh_.points()
    );
    allFaces[Pstream::myProcNo()] = pp.localFaces();
    allPoints[Pstream::myProcNo()] = pp.localPoints();

    Pstream::gatherList(allFaces);
    Pstream::gatherList(allPoints);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, proci)
    {
        nFaces += allFaces[proci].size();
        nPoints += allPoints[proci].size();
    }

    faces.setSize(nFaces);
    points.setSize(nPoints);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const faceList& fcs = allFaces[Pstream::myProcNo()];
        for (const face& f : fcs)
        {
            face& newF = faces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        for (const point& pt : pts)
        {
            points[nPoints++] = pt;
        }
    }

    // Other proc data follows
    forAll(allFaces, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            const faceList& fcs = allFaces[proci];
            for (const face& f : fcs)
            {
                face& newF = faces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[proci];
            for (const point& pt : pts)
            {
                points[nPoints++] = pt;
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,
        SMALL,
        false,
        oldToNew,
        newPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << points.size()
                << " down to " << newPoints.size() << " points" << endl;
        }

        points.transfer(newPoints);
        for (face& f : faces)
        {
            inplaceRenumber(oldToNew, f);
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::
combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        if (Pstream::parRun())
        {
            // Dimension as fraction of surface
            const scalar mergeDim = 1e-10*boundBox(s.points(), true).mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
    else if (surfacePtr_.valid())
    {
        const sampledSurface& s = surfacePtr_();

        if (Pstream::parRun())
        {
            // Dimension as fraction of mesh bounding box
            const scalar mergeDim = 1e-10*mesh_.bounds().mag();

            labelList pointsMap;

            PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch
                (
                    SubList<face>(s.faces(), s.faces().size()),
                    s.points()
                ),
                points,
                faces,
                pointsMap
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::totalArea() const
{
    scalar totalArea;

    if (regionType_ == stSurface)
    {
        const surfMesh& s = dynamicCast<const surfMesh>(obr());

        totalArea = gSum(s.magSf());
    }
    else if (surfacePtr_.valid())
    {
        totalArea = gSum(surfacePtr_().magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh_.magSf()));
    }

    return totalArea;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::usesSf() const
{
    // Only a few operations do not require the Sf field
    switch (operation_)
    {
        case opNone:
        case opMin:
        case opMax:
        case opSum:
        case opSumMag:
        case opAverage:
        case opSumDirection:
        case opSumDirectionBalance:
        {
            return false;
        }
        default:
        {
            return true;
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::initialise
(
    const dictionary& dict
)
{
    dict.lookup("name") >> regionName_;

    switch (regionType_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            surfacePtr_.clear();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            surfacePtr_.clear();
            break;
        }
        case stSurface:
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            nFaces_ = returnReduce(s.size(), sumOp<label>());

            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();
            surfacePtr_.clear();
            break;
        }
        case stSampledSurface:
        {
            faceId_.clear();
            facePatchId_.clear();
            faceFlip_.clear();

            surfacePtr_ = sampledSurface::New
            (
                name(),
                mesh_,
                dict.subDict("sampledSurfaceDict")
            );
            surfacePtr_().update();
            nFaces_ =
                returnReduce(surfacePtr_().faces().size(), sumOp<label>());

            break;
        }
        default:
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << int(regionType_) << "(" << regionName_ << "):"
                << nl << "    Unknown region type. Valid region types are:"
                << regionTypeNames_ << nl
                << exit(FatalError);
        }
    }

    if (nFaces_ == 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << regionTypeNames_[regionType_] << "(" << regionName_ << "):" << nl
            << "    Region has no faces" << exit(FatalError);
    }

    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    totalArea_ = totalArea();

    Info<< type() << " " << name() << ":" << nl
        << "    operation     = ";

    if (postOperation_ != postOpNone)
    {
        Info<< postOperationTypeNames_[postOperation_] << '('
            << operationTypeNames_[operation_] << ')'  << nl;
    }
    else
    {
        Info<< operationTypeNames_[operation_] << nl;
    }
    Info<< "    total faces   = " << nFaces_ << nl
        << "    total area    = " << totalArea_ << nl;


    weightFieldName_ = "none";
    if (usesWeight())
    {
        if (regionType_ == stSampledSurface)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot use weighted operation '"
                << operationTypeNames_[operation_]
                << "' for sampledSurface"
                << exit(FatalIOError);
        }

        if (dict.readIfPresent("weightField", weightFieldName_))
        {
            Info<< "    weight field  = " << weightFieldName_ << nl;
        }
        else
        {
            // Suggest possible alternative unweighted operation?
            FatalIOErrorInFunction(dict)
                << "The '" << operationTypeNames_[operation_]
                << "' operation is missing a weightField." << nl
                << "Either provide the weightField, "
                << "use weightField 'none' to suppress weighting," << nl
                << "or use a different operation."
                << exit(FatalIOError);
        }
    }

    // Backwards compatibility for v1612+ and older
    List<word> orientedFields;
    if (dict.readIfPresent("orientedFields", orientedFields))
    {
        WarningInFunction
            << "The 'orientedFields' option is deprecated.  These fields can "
            << "and have been added to the standard 'fields' list."
            << endl;

        fields_.append(orientedFields);
    }


    surfaceWriterPtr_.clear();
    if (writeFields_)
    {
        const word surfaceFormat(dict.lookup("surfaceFormat"));

        if (surfaceFormat != "none")
        {
            surfaceWriterPtr_.reset
            (
                surfaceWriter::New
                (
                    surfaceFormat,
                    dict.subOrEmptyDict("formatOptions").
                        subOrEmptyDict(surfaceFormat)
                ).ptr()
            );

            Info<< "    surfaceFormat = " << surfaceFormat << nl;
        }
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    if (operation_ != opNone)
    {
        writeCommented(os, "Region type : ");
        os << regionTypeNames_[regionType_] << " " << regionName_ << endl;

        writeHeaderValue(os, "Faces", nFaces_);
        writeHeaderValue(os, "Area", totalArea_);
        writeHeaderValue(os, "Scale factor", scaleFactor_);

        if (weightFieldName_ != "none")
        {
            writeHeaderValue(os, "Weight field", weightFieldName_);
        }

        writeCommented(os, "Time");
        if (writeArea_)
        {
            os  << tab << "Area";
        }

        for (const word& fieldName : fields_)
        {
            os  << tab << operationTypeNames_[operation_]
                << "(" << fieldName << ")";
        }

        os  << endl;
    }
}


template<>
Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<scalar>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            const vector n(dict_.lookup("direction"));
            return gSum(pos0(values*(Sf & n))*mag(values));
        }
        case opSumDirectionBalance:
        {
            const vector n(dict_.lookup("direction"));
            const scalarField nv(values*(Sf & n));

            return gSum(pos0(nv)*mag(values) - neg(nv)*mag(values));
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::vector
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opSumDirectionBalance:
        {
            vector n(dict_.lookup("direction"));
            n /= mag(n) + ROOTVSMALL;

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opAreaNormalAverage:
        {
            const scalar val = gSum(values & Sf)/gSum(mag(Sf));
            return vector(val, 0, 0);
        }
        case opAreaNormalIntegrate:
        {
            const scalar val = gSum(values & Sf);
            return vector(val, 0, 0);
        }
        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField
) const
{
    if (usesMag())
    {
        return mag(weightField);
    }
    else
    {
        // pass through
        return weightField;
    }
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField,
    const vectorField& Sf
) const
{
    // scalar * Area
    if (returnReduce(weightField.empty(), andOp<bool>()))
    {
        // No weight field - revert to unweighted form
        return mag(Sf);
    }
    else if (usesMag())
    {
        return mag(weightField * mag(Sf));
    }
    else
    {
        return (weightField * mag(Sf));
    }
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<vector>& weightField,
    const vectorField& Sf
) const
{
    // vector (dot) Area
    if (returnReduce(weightField.empty(), andOp<bool>()))
    {
        // No weight field - revert to unweighted form
        return mag(Sf);
    }
    else if (usesMag())
    {
        return mag(weightField & Sf);
    }
    else
    {
        return (weightField & Sf);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    regionType_(regionTypeNames_.lookup("regionType", dict)),
    operation_(operationTypeNames_.lookup("operation", dict)),
    postOperation_
    (
        postOperationTypeNames_.lookupOrDefault
        (
            "postOperation",
            dict,
            postOperationType::postOpNone
        )
    ),
    weightFieldName_("none"),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    regionType_(regionTypeNames_.lookup("regionType", dict)),
    operation_(operationTypeNames_.lookup("operation", dict)),
    postOperation_
    (
        postOperationTypeNames_.lookupOrDefault
        (
            "postOperation",
            dict,
            postOperationType::postOpNone
        )
    ),
    weightFieldName_("none"),
    writeArea_(dict.lookupOrDefault("writeArea", false)),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::~surfaceFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
    if (surfacePtr_.valid())
    {
        surfacePtr_().update();
    }

    if (operation_ != opNone)
    {
        fieldValue::write();

        if (Pstream::master())
        {
            writeTime(file());
        }
    }

    if (writeArea_)
    {
        totalArea_ = totalArea();
        Log << "    total area = " << totalArea_ << endl;

        if (operation_ != opNone && Pstream::master())
        {
            file() << tab << totalArea_;
        }
    }

    // Many operations use the Sf field
    vectorField Sf;
    if (usesSf())
    {
        if (regionType_ == stSurface)
        {
            const surfMesh& s = dynamicCast<const surfMesh>(obr());
            Sf = s.Sf();
        }
        else if (surfacePtr_.valid())
        {
            Sf = surfacePtr_().Sf();
        }
        else
        {
            Sf = filterField(mesh_.Sf());
        }
    }

    // Faces and points for surface format (if specified)
    faceList faces;
    pointField points;

    if (surfaceWriterPtr_.valid())
    {
        if (regionType_ == stSurface || surfacePtr_.valid())
        {
            combineSurfaceGeometry(faces, points);
        }
        else
        {
            combineMeshGeometry(faces, points);
        }
    }

    meshedSurfRef surfToWrite(points, faces);

    // Only a few weight types (scalar, vector)
    if (weightFieldName_ != "none")
    {
        if (validField<scalar>(weightFieldName_))
        {
            scalarField weightField
            (
                getFieldValues<scalar>(weightFieldName_, true)
            );

            // Process the fields
            writeAll(Sf, weightField, surfToWrite);
        }
        else if (validField<vector>(weightFieldName_))
        {
            vectorField weightField
            (
                getFieldValues<vector>(weightFieldName_, true)
            );

            // Process the fields
            writeAll(Sf, weightField, surfToWrite);
        }
        else
        {
            FatalErrorInFunction
                << "weightField " << weightFieldName_
                << " not found or an unsupported type"
                << abort(FatalError);
        }
    }
    else
    {
        // Default is a zero-size scalar weight field (ie, weight = 1)
        scalarField weightField;

        // Process the fields
        writeAll(Sf, weightField, surfToWrite);
    }

    if (operation_ != opNone && Pstream::master())
    {
        file() << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
