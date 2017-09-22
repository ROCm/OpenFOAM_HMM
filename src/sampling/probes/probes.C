/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "probes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        probes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh)
{
    DebugInfo<< "probes: resetting sample locations" << endl;

    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    processor_.setSize(size());
    processor_ = -1;

    forAll(*this, probei)
    {
        const vector& location = operator[](probei);

        const label celli = mesh.findCell(location);

        elementList_[probei] = celli;

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];
            const vector& cellCentre = mesh.cellCentres()[celli];
            scalar minDistance = GREAT;
            label minFaceID = -1;
            forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector dist = mesh.faceCentres()[facei] - cellCentre;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = facei;
                }
            }
            faceList_[probei] = minFaceID;
        }
        else
        {
            faceList_[probei] = -1;
        }

        if (debug && (elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "probes : found point " << location
                << " in cell " << elementList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(elementList_, probei)
    {
        const vector& location = operator[](probei);
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        processor_[probei] = (celli != -1 ? Pstream::myProcNo() : -1);

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());
        reduce(processor_[probei], maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[probei] != -1 && elementList_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::probes::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        DebugInfo
            << "Probing fields: " << currentFields << nl
            << "Probing locations: " << *this << nl
            << endl;


        fileName probeDir;
        fileName probeSubDir = name();

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/mesh_.name();
        }
        probeSubDir = "postProcessing"/probeSubDir/mesh_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = mesh_.time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = mesh_.time().path()/probeSubDir;
        }
        // Remove ".."
        probeDir.clean();

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                DebugInfo<< "close probe stream: " << iter()->name() << endl;

                delete probeFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(probeDir);

            OFstream* fPtr = new OFstream(probeDir/fieldName);

            OFstream& fout = *fPtr;

            DebugInfo<< "open probe stream: " << fout.name() << endl;

            probeFilePtrs_.insert(fieldName, fPtr);

            unsigned int w = IOstream::defaultPrecision() + 7;

            forAll(*this, probei)
            {
                fout<< "# Probe " << probei << ' ' << operator[](probei);

                if (processor_[probei] == -1)
                {
                    fout<< "  # Not Found";
                }
                fout<< endl;
            }

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Probe";

            forAll(*this, probei)
            {
                if (includeOutOfBounds_ || processor_[probei] != -1)
                {
                    fout<< ' ' << setw(w) << probei;
                }
            }
            fout<< endl;

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << endl;
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    stateFunctionObject(name, runTime),
    pointField(0),
    mesh_
    (
        refCast<const fvMesh>
        (
            runTime.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    loadFromFiles_(loadFromFiles),
    fieldSelection_(),
    fixedLocations_(true),
    interpolationScheme_("cell"),
    includeOutOfBounds_(true)
{
    if (readFields)
    {
        read(dict);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::probes::read(const dictionary& dict)
{
    dict.lookup("probeLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;

    dict.readIfPresent("fixedLocations", fixedLocations_);
    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored"
                << endl;
        }
    }
    dict.readIfPresent("includeOutOfBounds", includeOutOfBounds_);

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    prepare();

    return true;
}


bool Foam::probes::execute()
{
    return true;
}


bool Foam::probes::write()
{
    if (size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


void Foam::probes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "probes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        findElements(mesh_);
    }
    else
    {
        DebugInfo<< "probes: remapping sample locations" << endl;

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                if (celli != -1)
                {
                    label newCelli = reverseMap[celli];
                    if (newCelli == -1)
                    {
                        // cell removed
                    }
                    else if (newCelli < -1)
                    {
                        // cell merged
                        elems.append(-newCelli - 2);
                    }
                    else
                    {
                        // valid new cell
                        elems.append(newCelli);
                    }
                }
                else
                {
                    // Keep -1 elements so the size stays the same
                    elems.append(-1);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            forAll(faceList_, i)
            {
                label facei = faceList_[i];
                if (facei != -1)
                {
                    label newFacei = reverseMap[facei];
                    if (newFacei == -1)
                    {
                        // face removed
                    }
                    else if (newFacei < -1)
                    {
                        // face merged
                        elems.append(-newFacei - 2);
                    }
                    else
                    {
                        // valid new face
                        elems.append(newFacei);
                    }
                }
                else
                {
                    // Keep -1 elements
                    elems.append(-1);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::probes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        findElements(mesh_);
    }
}


// ************************************************************************* //
