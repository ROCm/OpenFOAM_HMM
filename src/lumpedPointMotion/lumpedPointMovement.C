/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "lumpedPointMovement.H"
#include "lumpedPointIOMovement.H"
#include "demandDrivenData.H"
#include "linearInterpolationWeights.H"
#include "Fstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrMap.H"
#include "faceZoneMesh.H"
#include "PstreamReduceOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::lumpedPointMovement::outputFormatType
>
Foam::lumpedPointMovement::formatNames
({
    { outputFormatType::PLAIN, "plain" },
    { outputFormatType::DICTIONARY, "dictionary" },
});


const Foam::Enum
<
    Foam::lumpedPointMovement::scalingType
>
Foam::lumpedPointMovement::scalingNames
({
    { scalingType::LENGTH, "length" },
    { scalingType::FORCE, "force" },
    { scalingType::MOMENT, "moment" },
});


const Foam::word
Foam::lumpedPointMovement::canonicalName("lumpedPointMovement");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//! \cond fileScope
//- Space-separated vector value (ASCII)
static inline Ostream& putPlain(Ostream& os, const vector& val)
{
    os  << val.x() << ' ' << val.y() << ' ' << val.z();
    return os;
}


//! \cond fileScope
//- Space-separated vector value (ASCII)
static inline Ostream& putTime(Ostream& os, const Time& t)
{
    os  <<"Time index=" << t.timeIndex()
        << " value=" << t.timeOutputValue();

    return os;
}


//! \cond fileScope
//- Write list content with size, bracket, content, bracket one-per-line.
//  This makes for consistent for parsing, regardless of the list length.
template <class T>
static void writeList(Ostream& os, const string& header, const UList<T>& list)
{
    const label len = list.size();

    // Header string
    os  << header.c_str() << nl;

    // Write size and start delimiter
    os  << len << nl << token::BEGIN_LIST << nl;

    // Write contents
    for (label i=0; i < len; ++i)
    {
        os << list[i] << nl;
    }

    // Write end delimiter
    os << token::END_LIST << token::END_STATEMENT << nl << nl;
}
//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lumpedPointMovement::calcThresholds() const
{
    thresholdPtr_ = new scalarField(locations_);

    scalarField& thresh = *thresholdPtr_;

    for (label i=1; i < thresh.size(); ++i)
    {
        thresh[i-1] =
        (
            locations_[i-1]
          + (division_ * (locations_[i] - locations_[i-1]))
        );
    }
}


Foam::label Foam::lumpedPointMovement::threshold(scalar pos) const
{
    const scalarField& thresh = this->thresholds();
    const label n = thresh.size();

    // Clamp above and below
    //
    // ? may need special treatment to clamp values below the first threshold ?
    // ? special treatment when 'axis' is negative ?

    for (label i=0; i < n; ++i)
    {
        if (pos < thresh[i])
        {
            return i;
        }
    }

    // everything above goes into the last zone too
    return n-1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointMovement::lumpedPointMovement()
:
    axis_(0,0,1),
    locations_(),
    division_(0),
    relax_(1),
    interpolationScheme_(linearInterpolationWeights::typeName),
    ownerId_(-1),
    boundBox_(),
    centre_(),
    autoCentre_(true),
    forcesDict_(),
    coupler_(),
    inputName_("positions.in"),
    outputName_("forces.out"),
    logName_("movement.log"),
    inputFormat_(lumpedPointState::inputFormatType::DICTIONARY),
    outputFormat_(outputFormatType::DICTIONARY),
    scaleInput_(-1.0),
    scaleOutput_(-1.0),
    state0_(),
    state_(),
    thresholdPtr_(0),
    interpolatorPtr_()
{}


Foam::lumpedPointMovement::lumpedPointMovement
(
    const dictionary& dict,
    label ownerId
)
:
    axis_(0,0,1),
    locations_(),
    division_(0),
    relax_(1),
    interpolationScheme_
    (
        dict.lookupOrDefault<word>
        (
            "interpolationScheme",
            linearInterpolationWeights::typeName
        )
    ),
    ownerId_(ownerId),
    boundBox_(),
    centre_(),
    autoCentre_(true),
    forcesDict_(),
    coupler_(),
    inputName_("positions.in"),
    outputName_("forces.out"),
    logName_("movement.log"),
    scaleInput_(-1.0),
    scaleOutput_(-1.0),
    state0_(),
    state_(),
    thresholdPtr_(0),
    interpolatorPtr_()
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::lumpedPointMovement::~lumpedPointMovement()
{
    deleteDemandDrivenData(thresholdPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointMovement::readDict(const dictionary& dict)
{
    // assume the worst
    deleteDemandDrivenData(thresholdPtr_);

    dict.readEntry("axis", axis_);

    division_ = 0;
    if (dict.readIfPresent("division", division_))
    {
        if (division_ < 0)
        {
            division_ = 0;
        }
        else if (division_ > 1)
        {
            division_ = 1;
        }
    }

    dict.readIfPresent("relax", relax_);

    dict.readEntry("locations", locations_);

    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        interpolatorPtr_.clear();
    }

    autoCentre_ = !dict.readIfPresent("centre", centre_);

    // Initial state from initial locations, with zero rotation
    tmp<vectorField> pts = locations_ * axis_;
    state0_ = lumpedPointState(pts);
    state_  = state0_;

    forcesDict_.merge(dict.subOrEmptyDict("forces"));

    const dictionary& commDict = dict.subDict("communication");
    coupler_.readDict(commDict);

    // TODO: calcFrequency_  = dict.lookupOrDefault("calcFrequency", 1);

    commDict.readEntry("inputName", inputName_);
    commDict.readEntry("outputName", outputName_);
    commDict.readIfPresent("logName", logName_);

    inputFormat_ =
        lumpedPointState::formatNames.get("inputFormat", commDict);

    outputFormat_ =
        lumpedPointMovement::formatNames.get("outputFormat", commDict);

    scaleInput_  = -1;
    scaleOutput_ = -1;

    const dictionary* scaleDict = nullptr;

    if ((scaleDict = commDict.findDict("scaleInput")))
    {
        for (int i=0; i < scaleInput_.size(); ++i)
        {
            const word& key = scalingNames[scalingType(i)];

            if
            (
                scaleDict->readIfPresent(key, scaleInput_[i])
             && scaleInput_[i] > 0
            )
            {
                Info<<"Using input " << key << " multiplier: "
                    << scaleInput_[i] << nl;
            }
        }
    }

    if ((scaleDict = commDict.findDict("scaleOutput")))
    {
        for (int i=0; i < scaleOutput_.size(); ++i)
        {
            const word& key = scalingNames[scalingType(i)];

            if
            (
                scaleDict->readIfPresent(key, scaleOutput_[i])
             && scaleOutput_[i] > 0
            )
            {
                Info<<"Using output " << key << " multiplier: "
                    << scaleOutput_[i] << nl;
            }
        }
    }
}


void Foam::lumpedPointMovement::setBoundBox
(
    const polyMesh& mesh,
    const labelUList& patchIds,
    const pointField& points0
)
{
    boundBox_ = boundBox::invertedBox;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    for (const label patchId : patchIds)
    {
        boundBox_.add(points0, patches[patchId].meshPoints());
    }

    boundBox_.reduce();

    if (autoCentre_)
    {
        centre_ = boundBox_.centre();
        centre_ -= (centre_ & axis_) * axis_;
        if (lumpedPointIOMovement::debug)
        {
            Pout<<"autoCentre on " << centre_
                << " boundBox: " << boundBox_ << endl;
        }
    }
    else
    {
        // User-specified centre
        if (lumpedPointIOMovement::debug)
        {
            Pout<<"centre on " << centre_
                << " boundBox: " << boundBox_ << endl;
        }
    }
}


void Foam::lumpedPointMovement::setMapping
(
    const polyMesh& mesh,
    const labelUList& patchIds,
    const pointField& points0
)
{
    setBoundBox(mesh, patchIds, points0);

    faceZones_.clear();

    // Local storage location (boundary faces only)
    const label firstFace = mesh.nInternalFaces();
    labelList faceToZoneID(mesh.nFaces() - firstFace, -1);

    // Number of faces per zone
    labelList nFaces(thresholds().size(), Zero);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (const label patchId : patchIds)
    {
        const polyPatch& pp = patches[patchId];

        // Classify faces into respective zones according to their centres
        label meshFaceI = pp.start();

        forAll(pp, i)
        {
            const face& f = mesh.faces()[meshFaceI];
            label zoneId = this->threshold(f.centre(points0));

            // localize storage location
            faceToZoneID[meshFaceI - firstFace] = zoneId;
            ++nFaces[zoneId];

            ++meshFaceI;
        }
    }

    if (lumpedPointIOMovement::debug)
    {
        Pout<<"faces per zone:" << nFaces << endl;
    }

    faceZones_.setSize(nFaces.size());

    label nZonedFaces = 0;
    forAll(nFaces, zonei)
    {
        label n = nFaces[zonei];
        labelList& addr = faceZones_[zonei];
        addr.setSize(n);

        n = 0;
        forAll(faceToZoneID, faceI)
        {
            if (faceToZoneID[faceI] == zonei)
            {
                // from local storage to global mesh face
                label meshFaceI = faceI + firstFace;

                addr[n] = meshFaceI;
                ++n;
            }
        }

        addr.setSize(n);
        nZonedFaces += n;

        if (lumpedPointIOMovement::debug)
        {
            Pout<< "Created faceZone[" << zonei << "] "
                << returnReduce(n, sumOp<label>()) << " faces, "
                << thresholds()[zonei] << " threshold" << nl;
        }
    }
}


bool Foam::lumpedPointMovement::forcesAndMoments
(
    const polyMesh& pmesh,
    List<vector>& forces,
    List<vector>& moments
) const
{
    forces.setSize(faceZones_.size());
    moments.setSize(faceZones_.size());

    if (faceZones_.empty())
    {
        WarningInFunction
            << "Attempted to calculate forces without setMapping()"
            << nl << endl;

        forces.setSize(thresholds().size(), Zero);
        moments.setSize(thresholds().size(), Zero);
        return false;
    }

    // Initialize with zero
    forces = Zero;
    moments = Zero;

    // The mass centres
    const pointField& lumpedCentres = state().points();

    const polyBoundaryMesh& patches = pmesh.boundaryMesh();

    const word pName(forcesDict_.lookupOrDefault<word>("p", "p"));
    scalar pRef   = forcesDict_.lookupOrDefault<scalar>("pRef",   0.0);
    scalar rhoRef = forcesDict_.lookupOrDefault<scalar>("rhoRef", 1.0);


    // Calculated force per patch - cache
    PtrMap<vectorField> forceOnPatches;

    const volScalarField* pPtr = pmesh.findObject<volScalarField>(pName);

    // fvMesh and has pressure field
    if (isA<fvMesh>(pmesh) && pPtr)
    {
        const fvMesh& mesh = dynamicCast<const fvMesh>(pmesh);
        const volScalarField& p = *pPtr;

        // face centres (on patches)
        const surfaceVectorField::Boundary& patchCf =
            mesh.Cf().boundaryField();

        // face areas (on patches)
        const surfaceVectorField::Boundary& patchSf =
            mesh.Sf().boundaryField();

        // Pressure (on patches)
        const volScalarField::Boundary& patchPress =
            p.boundaryField();

        // rhoRef if the pressure field is dynamic, i.e. p/rho otherwise 1
        rhoRef = (p.dimensions() == dimPressure ? 1.0 : rhoRef);

        // Scale pRef by density for incompressible simulations
        pRef /= rhoRef;

        forAll(faceZones_, zonei)
        {
            const labelList& zn = faceZones_[zonei];

            forAll(zn, localFaceI)
            {
                const label meshFaceI = zn[localFaceI];

                // locate which patch this face is on,
                // and its offset within the patch
                const label patchI = patches.whichPatch(meshFaceI);
                if (patchI == -1)
                {
                    // Should not happen - could also warn
                    continue;
                }

                const label patchFaceI = meshFaceI - patches[patchI].start();

                if (!forceOnPatches.found(patchI))
                {
                    // Patch faces are +ve outwards,
                    // so the forces (exerted by fluid on solid)
                    // already have the correct sign
                    forceOnPatches.set
                    (
                        patchI,
                        (
                            rhoRef
                          * patchSf[patchI] * (patchPress[patchI] - pRef)
                        ).ptr()
                    );
                }

                const vectorField& forceOnPatch = *forceOnPatches[patchI];

                // Force from the patch-face into sum
                forces[zonei] += forceOnPatch[patchFaceI];

                // effective torque arm:
                // - translated into the lumped-points coordinate system
                //   prior to determining the distance
                const vector lever =
                (
                    patchCf[patchI][patchFaceI] - centre_
                  - lumpedCentres[zonei]
                );

                // Moment from the patch-face into sum
                moments[zonei] += lever ^ forceOnPatch[patchFaceI];
            }
        }
    }
    else
    {
        Info<<"no pressure" << endl;
    }

    Pstream::listCombineGather(forces, plusEqOp<vector>());
    Pstream::listCombineScatter(forces);

    Pstream::listCombineGather(moments, plusEqOp<vector>());
    Pstream::listCombineScatter(moments);

    return true;
}


// \cond fileScope
template<class T>
static inline T interpolatedValue
(
    const Foam::UList<T> input,
    const Foam::labelUList& indices,
    const Foam::scalarField& weights
)
{
    T result = weights[0] * input[indices[0]];
    for (Foam::label i=1; i < indices.size(); ++i)
    {
        result += weights[i] * input[indices[i]];
    }

    return result;
}
// \endcond


Foam::tmp<Foam::pointField>
Foam::lumpedPointMovement::displacePoints
(
    const pointField& points0,
    const labelList& pointLabels
) const
{
    return displacePoints(state(), points0, pointLabels);
}


Foam::tmp<Foam::pointField>
Foam::lumpedPointMovement::displacePoints
(
    const lumpedPointState& state,
    const pointField& points0,
    const labelList& pointLabels
) const
{
    // storage for the interpolation indices/weights
    labelList   indices;
    scalarField weights;
    const interpolationWeights& interp = interpolator();

    // The mass centres
    const pointField& lumpedCentres = state.points();

    // local-to-global transformation tensor
    const tensorField& localToGlobal = state.rotations();

    // could also verify the sizes (state vs original)

    tmp<pointField> tdisp(new pointField(pointLabels.size()));
    pointField& disp = tdisp.ref();

    forAll(pointLabels, ptI)
    {
        const point& p0 = points0[pointLabels[ptI]];

        // Interpolation factor based on the axis component
        // of the original point (location)
        scalar where = p0 & axis_;

        interp.valueWeights(where, indices, weights);

        vector origin    = interpolatedValue(lumpedCentres, indices, weights);
        tensor rotTensor = interpolatedValue(localToGlobal, indices, weights);

        if (indices.size() == 1)
        {
            // out-of-bounds, use corresponding end-point
            where = locations_[indices[0]];
        }

        // local point:
        // - translated into the lumped-points coordinate system
        // - relative to the interpolation (where) location
        const vector local = (p0 - (where * axis_)) - centre_;

        // local-to-global rotation and translation
        // - translate back out of the lumped-points coordinate system
        //
        // -> subtract p0 to return displacements
        disp[ptI] = (rotTensor & local) + origin + centre_ - p0;
    }

    return tdisp;
}


void Foam::lumpedPointMovement::writeDict(Ostream& os) const
{
    os.writeEntry("axis", axis_);
    os.writeEntry("locations", locations_);
    os.writeEntry("division", division_);
    // os.writeEntry("interpolationScheme", interpolationScheme_);
}


bool Foam::lumpedPointMovement::readState()
{
    lumpedPointState prev = state_;

    bool status = state_.readData
    (
        inputFormat_,
        coupler().resolveFile(inputName_)
    );

    state_.scalePoints(scaleInput_[scalingType::LENGTH]);

    state_.relax(relax_, prev);

    return status;
}


bool Foam::lumpedPointMovement::writeData
(
    Ostream& os,
    const UList<vector>& forces,
    const UList<vector>& moments,
    const outputFormatType fmt,
    const Time* timeinfo
) const
{
    const bool writeMoments = (moments.size() == forces.size());

    if (fmt == outputFormatType::PLAIN)
    {
        os  <<"########" << nl;
        if (timeinfo)
        {
            os <<"# ";
            putTime(os, *timeinfo) << nl;
        }
        os  <<"# size=" << this->size() << nl
            <<"# columns (points) (forces)";

        if (writeMoments)
        {
            os << " (moments)";
        }

        os << nl;

        bool report = false;
        scalar scaleLength = scaleOutput_[scalingType::LENGTH];
        scalar scaleForce  = scaleOutput_[scalingType::FORCE];
        scalar scaleMoment = scaleOutput_[scalingType::MOMENT];

        if (scaleLength > 0)
        {
            report = true;
        }
        else
        {
            scaleLength = 1.0;
        }

        if (scaleForce > 0)
        {
            report = true;
        }
        else
        {
            scaleForce = 1.0;
        }

        if (writeMoments)
        {
            if (scaleMoment > 0)
            {
                report = true;
            }
            else
            {
                scaleMoment = 1.0;
            }
        }

        if (report)
        {
            os  <<"# scaling points=" << scaleLength
                <<" forces=" << scaleForce;

            if (writeMoments)
            {
                os  <<" moments=" << scaleMoment;
            }

            os << nl;
        }

        os <<"########" << nl;

        forAll(locations_, i)
        {
            const vector pos = scaleLength * (locations_[i] * axis_);

            putPlain(os, pos) << ' ';

            if (i < forces.size())
            {
                const vector val(scaleForce * forces[i]);
                putPlain(os, val);
            }
            else
            {
                putPlain(os, vector::zero);
            }

            if (writeMoments)
            {
                os << ' ';
                if (i < moments.size())
                {
                    const vector val(scaleMoment * moments[i]);
                    putPlain(os, val);
                }
                else
                {
                    putPlain(os, vector::zero);
                }
            }

            os << nl;
        }
    }
    else
    {
        // Make it easier for external programs to parse
        // - exclude the usual OpenFOAM 'FoamFile' header
        // - ensure lists have consistent format

        os  <<"////////" << nl;
        if (timeinfo)
        {
            os <<"// ";
            putTime(os, *timeinfo) << nl;
        }
        os  << nl;

        writeList(os, "points", (locations_*axis_)());
        writeList(os, "forces", forces);

        if (writeMoments)
        {
            writeList(os, "moments", moments);
        }
    }

    return true;
}


bool Foam::lumpedPointMovement::writeData
(
    const UList<vector>& forces,
    const UList<vector>& moments,
    const Time* timeinfo
) const
{
    if (!Pstream::master())
    {
        return false;
    }

    // Regular output
    {
        const fileName output(coupler().resolveFile(outputName_));
        OFstream os(output, IOstream::ASCII);

        writeData(os, forces, moments, outputFormat_, timeinfo);
    }

    // Log output
    {
        const fileName output(coupler().resolveFile(logName_));

        OFstream os
        (
            output,
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true // append mode
        );

        writeData(os, forces, moments, outputFormatType::PLAIN, timeinfo);
    }

    return true;
}


const Foam::interpolationWeights&
Foam::lumpedPointMovement::interpolator() const
{
    if (!interpolatorPtr_.valid())
    {
        interpolatorPtr_ = interpolationWeights::New
        (
            interpolationScheme_,
            locations_
        );
    }

    return *interpolatorPtr_;
}


// ************************************************************************* //
