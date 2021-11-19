/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "regionSizeDistribution.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(regionSizeDistribution, 0);
        addToRunTimeSelectionTable
        (
            functionObject,
            regionSizeDistribution,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::regionSizeDistribution::writeGraph
(
    const coordSet& coords,
    const word& valueName,
    const scalarField& values
) const
{
    const wordList valNames(1, valueName);

    fileName outputPath = baseTimeDir();
    mkDir(outputPath);

    OFstream str(outputPath/formatterPtr_().getFileName(coords, valNames));

    Log << "    Writing distribution of " << valueName << " to " << str.name()
        << endl;

    List<const scalarField*> valPtrs(1);
    valPtrs[0] = &values;
    formatterPtr_().write(coords, valNames, valPtrs, str);
}


void Foam::functionObjects::regionSizeDistribution::writeAlphaFields
(
    const regionSplit& regions,
    const Map<label>& patchRegions,
    const Map<scalar>& regionVolume,
    const volScalarField& alpha
) const
{
    const scalar maxDropletVol = 1.0/6.0*pow3(maxDiam_);

    // Split alpha field
    // ~~~~~~~~~~~~~~~~~
    // Split into
    //  - liquidCore            : region connected to inlet patches
    //  - per region a volume   : for all other regions
    //  - backgroundAlpha       : remaining alpha


    // Construct field
    volScalarField liquidCore
    (
        IOobject
        (
            scopedName(alphaName_ + "_liquidCore"),
            obr_.time().timeName(),
            obr_,
            IOobject::NO_READ
        ),
        alpha,
        fvPatchField<scalar>::calculatedType()
    );

    volScalarField backgroundAlpha
    (
        IOobject
        (
            scopedName(alphaName_ + "_background"),
            obr_.time().timeName(),
            obr_,
            IOobject::NO_READ
        ),
        alpha,
        fvPatchField<scalar>::calculatedType()
    );


    // Knock out any cell not in patchRegions
    forAll(liquidCore, celli)
    {
        const label regioni = regions[celli];
        if (patchRegions.found(regioni))
        {
            backgroundAlpha[celli] = 0;
        }
        else
        {
            liquidCore[celli] = 0;

            const scalar regionVol = regionVolume[regioni];
            if (regionVol < maxDropletVol)
            {
                backgroundAlpha[celli] = 0;
            }
        }
    }
    liquidCore.correctBoundaryConditions();
    backgroundAlpha.correctBoundaryConditions();

    if (log)
    {
        Info<< "    Volume of liquid-core = "
            << fvc::domainIntegrate(liquidCore).value()
            << nl
            << "    Volume of background  = "
            << fvc::domainIntegrate(backgroundAlpha).value()
            << endl;
    }

    Log << "    Writing liquid-core field to " << liquidCore.name() << endl;
    liquidCore.write();

    Log<< "    Writing background field to " << backgroundAlpha.name() << endl;
    backgroundAlpha.write();
}


Foam::Map<Foam::label>
Foam::functionObjects::regionSizeDistribution::findPatchRegions
(
    const regionSplit& regions
) const
{
    // Mark all regions starting at patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Count number of patch faces (just for initial sizing)
    const labelHashSet patchIDs(mesh_.boundaryMesh().patchSet(patchNames_));

    label nPatchFaces = 0;
    for (const label patchi : patchIDs)
    {
        nPatchFaces += mesh_.boundaryMesh()[patchi].size();
    }


    Map<label> patchRegions(nPatchFaces);
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        // Collect all regions on the patch
        const labelList& faceCells = pp.faceCells();

        for (const label celli : faceCells)
        {
            patchRegions.insert
            (
                regions[celli],
                Pstream::myProcNo()     // dummy value
            );
        }
    }


    // Make sure all the processors have the same set of regions
    Pstream::mapCombineGather(patchRegions, minEqOp<label>());
    Pstream::mapCombineScatter(patchRegions);

    return patchRegions;
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::regionSizeDistribution::divide
(
    const scalarField& num,
    const scalarField& denom
)
{
    auto tresult = tmp<scalarField>::New(num.size());
    auto& result = tresult.ref();

    forAll(denom, i)
    {
        if (denom[i] != 0)
        {
            result[i] = num[i]/denom[i];
        }
        else
        {
            result[i] = 0.0;
        }
    }
    return tresult;
}


void Foam::functionObjects::regionSizeDistribution::writeGraphs
(
    const word& fieldName,              // name of field
    const labelList& indices,           // index of bin for each region
    const scalarField& sortedField,     // per region field data
    const scalarField& binCount,        // per bin number of regions
    const coordSet& coords              // graph data for bins
) const
{
    if (Pstream::master())
    {
        // Calculate per-bin average
        scalarField binSum(nBins_, Zero);
        forAll(sortedField, i)
        {
            binSum[indices[i]] += sortedField[i];
        }

        scalarField binAvg(divide(binSum, binCount));

        // Per bin deviation
        scalarField binSqrSum(nBins_, Zero);
        forAll(sortedField, i)
        {
            binSqrSum[indices[i]] += Foam::sqr(sortedField[i]);
        }
        scalarField binDev
        (
            sqrt(divide(binSqrSum, binCount) - Foam::sqr(binAvg))
        );

        // Write average
        writeGraph(coords, fieldName + "_sum", binSum);
        // Write average
        writeGraph(coords, fieldName + "_avg", binAvg);
        // Write deviation
        writeGraph(coords, fieldName + "_dev", binDev);
    }
}


void Foam::functionObjects::regionSizeDistribution::writeGraphs
(
    const word& fieldName,              // name of field
    const scalarField& cellField,       // per cell field data
    const regionSplit& regions,         // per cell the region(=droplet)
    const labelList& sortedRegions,     // valid regions in sorted order
    const scalarField& sortedNormalisation,

    const labelList& indices,           // per region index of bin
    const scalarField& binCount,        // per bin number of regions
    const coordSet& coords              // graph data for bins
) const
{
    // Sum on a per-region basis. Parallel reduced.
    Map<scalar> regionField(regionSum(regions, cellField));

    // Extract in region order
    scalarField sortedField
    (
        sortedNormalisation
      * extractData
        (
            sortedRegions,
            regionField
        )
    );

    writeGraphs
    (
        fieldName,      // name of field
        indices,        // index of bin for each region
        sortedField,    // per region field data
        binCount,       // per bin number of regions
        coords          // graph data for bins
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionSizeDistribution::regionSizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    alphaName_(dict.get<word>("field")),
    patchNames_(dict.get<wordRes>("patches")),
    isoPlanes_(dict.getOrDefault("isoPlanes", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionSizeDistribution::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readEntry("nBins", nBins_);
    dict.readEntry("field", alphaName_);
    dict.readEntry("threshold", threshold_);
    dict.readEntry("maxDiameter", maxDiam_);
    minDiam_ = 0.0;
    dict.readIfPresent("minDiameter", minDiam_);
    dict.readEntry("patches", patchNames_);
    dict.readEntry("fields", fields_);

    const word format(dict.get<word>("setFormat"));
    formatterPtr_ = writer<scalar>::New(format);

    if (dict.found(coordinateSystem::typeName_()))
    {
        csysPtr_.reset
        (
            coordinateSystem::New(obr_, dict, coordinateSystem::typeName_())
        );

        Info<< "Transforming all vectorFields with coordinate system "
            << csysPtr_->name() << endl;
    }
    else
    {
        csysPtr_.clear();
    }

    if (isoPlanes_)
    {
         dict.readEntry("origin", origin_);
         dict.readEntry("direction", direction_);
         dict.readEntry("maxD", maxDiameter_);
         dict.readEntry("nDownstreamBins", nDownstreamBins_);
         dict.readEntry("maxDownstream", maxDownstream_);
         direction_.normalise();
    }

    return true;
}


bool Foam::functionObjects::regionSizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::regionSizeDistribution::write()
{
    Log << type() << " " << name() << " write:" << nl;

    autoPtr<volScalarField> alphaPtr;
    if (obr_.foundObject<volScalarField>(alphaName_))
    {
        Log << "    Looking up field " << alphaName_ << endl;
    }
    else
    {
        Info<< "    Reading field " << alphaName_ << endl;
        alphaPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    alphaName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );
    }


    const volScalarField& alpha =
    (
         alphaPtr
       ? *alphaPtr
       : obr_.lookupObject<volScalarField>(alphaName_)
    );

    Log << "    Volume of alpha          = "
        << fvc::domainIntegrate(alpha).value()
        << endl;

    const scalar meshVol = gSum(mesh_.V());
    const scalar maxDropletVol = 1.0/6.0*pow3(maxDiam_);
    const scalar delta = (maxDiam_-minDiam_)/nBins_;

    Log << "    Mesh volume              = " << meshVol << nl
        << "    Maximum droplet diameter = " << maxDiam_ << nl
        << "    Maximum droplet volume   = " << maxDropletVol
        << endl;


    // Determine blocked faces
    boolList blockedFace(mesh_.nFaces(), false);
    label nBlocked = 0;

    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            scalar ownVal = alpha[mesh_.faceOwner()[facei]];
            scalar neiVal = alpha[mesh_.faceNeighbour()[facei]];

            if
            (
                (ownVal < threshold_ && neiVal > threshold_)
             || (ownVal > threshold_ && neiVal < threshold_)
            )
            {
                blockedFace[facei] = true;
                nBlocked++;
            }
        }

        // Block coupled faces
        forAll(alpha.boundaryField(), patchi)
        {
            const fvPatchScalarField& fvp = alpha.boundaryField()[patchi];
            if (fvp.coupled())
            {
                tmp<scalarField> townFld(fvp.patchInternalField());
                const scalarField& ownFld = townFld();

                tmp<scalarField> tnbrFld(fvp.patchNeighbourField());
                const scalarField& nbrFld = tnbrFld();

                label start = fvp.patch().patch().start();

                forAll(ownFld, i)
                {
                    scalar ownVal = ownFld[i];
                    scalar neiVal = nbrFld[i];

                    if
                    (
                        (ownVal < threshold_ && neiVal > threshold_)
                     || (ownVal > threshold_ && neiVal < threshold_)
                    )
                    {
                        blockedFace[start+i] = true;
                        nBlocked++;
                    }
                }
            }
        }
    }


    regionSplit regions(mesh_, blockedFace);

    Log << "    Determined " << regions.nRegions()
        << " disconnected regions" << endl;


    if (debug)
    {
        volScalarField region
        (
            IOobject
            (
                "region",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        Info<< "    Dumping region as volScalarField to " << region.name()
            << endl;

        forAll(regions, celli)
        {
            region[celli] = regions[celli];
        }
        region.correctBoundaryConditions();
        region.write();
    }


    // Determine regions connected to supplied patches
    Map<label> patchRegions(findPatchRegions(regions));


    // Sum all regions
    const scalarField alphaVol(alpha.primitiveField()*mesh_.V());
    Map<scalar> allRegionVolume(regionSum(regions, mesh_.V()));
    Map<scalar> allRegionAlphaVolume(regionSum(regions, alphaVol));
    Map<label> allRegionNumCells
    (
        regionSum
        (
            regions,
            labelField(mesh_.nCells(), 1.0)
        )
    );

    if (debug)
    {
        Info<< "    " << token::TAB << "Region"
            << token::TAB << "Volume(mesh)"
            << token::TAB << "Volume(" << alpha.name() << "):"
            << token::TAB << "nCells"
            << nl;
        scalar meshSumVol = 0.0;
        scalar alphaSumVol = 0.0;
        label nCells = 0;

        auto vIter = allRegionVolume.cbegin();
        auto aIter = allRegionAlphaVolume.cbegin();
        auto numIter = allRegionNumCells.cbegin();
        for
        (
            ;
            vIter.good() && aIter.good();
            ++vIter, ++aIter, ++numIter
        )
        {
            Info<< "    " << token::TAB << vIter.key()
                << token::TAB << vIter()
                << token::TAB << aIter()
                << token::TAB << numIter()
                << nl;

            meshSumVol += vIter();
            alphaSumVol += aIter();
            nCells += numIter();
        }
        Info<< "    " << token::TAB << "Total:"
            << token::TAB << meshSumVol
            << token::TAB << alphaSumVol
            << token::TAB << nCells
            << endl;
    }


    if (log)
    {
        Info<< "    Patch connected regions (liquid core):" << nl;
        Info<< token::TAB << "    Region"
            << token::TAB << "Volume(mesh)"
            << token::TAB << "Volume(" << alpha.name() << "):"
            << nl;

        forAllConstIters(patchRegions, iter)
        {
            const label regioni = iter.key();
            Info<< "    " << token::TAB << regioni
                << token::TAB << allRegionVolume[regioni]
                << token::TAB << allRegionAlphaVolume[regioni] << nl;

        }
        Info<< endl;
    }

    if (log)
    {
        Info<< "    Background regions:" << nl;
        Info<< "    " << token::TAB << "Region"
            << token::TAB << "Volume(mesh)"
            << token::TAB << "Volume(" << alpha.name() << "):"
            << nl;

        auto vIter = allRegionVolume.cbegin();
        auto aIter = allRegionAlphaVolume.cbegin();

        for
        (
            ;
            vIter.good() && aIter.good();
            ++vIter, ++aIter
        )
        {
            if
            (
               !patchRegions.found(vIter.key())
             && vIter() >= maxDropletVol
            )
            {
                Info<< "    " << token::TAB << vIter.key()
                    << token::TAB << vIter()
                    << token::TAB << aIter() << nl;
            }
        }
        Info<< endl;
    }



    // Split alpha field
    // ~~~~~~~~~~~~~~~~~
    // Split into
    //  - liquidCore            : region connected to inlet patches
    //  - per region a volume   : for all other regions
    //  - backgroundAlpha       : remaining alpha
    writeAlphaFields(regions, patchRegions, allRegionVolume, alpha);


    // Extract droplet-only allRegionVolume, i.e. delete liquid core
    // (patchRegions) and background regions from maps.
    // Note that we have to use mesh volume (allRegionVolume) and not
    // allRegionAlphaVolume since background might not have alpha in it.
    // Deleting regions where the volume-alpha-weighted is lower than
    // threshold
    forAllIters(allRegionVolume, vIter)
    {
        const label regioni = vIter.key();
        if
        (
            patchRegions.found(regioni)
         || vIter() >= maxDropletVol
         || (allRegionAlphaVolume[regioni]/vIter() < threshold_)
        )
        {
            allRegionVolume.erase(vIter);
            allRegionAlphaVolume.erase(regioni);
            allRegionNumCells.erase(regioni);
        }
    }

    if (allRegionVolume.size())
    {
        // Construct mids of bins for plotting
        pointField xBin(nBins_);

        scalar x = 0.5*delta;
        forAll(xBin, i)
        {
            xBin[i] = point(x, 0, 0);
            x += delta;
        }

        const coordSet coords("diameter", "x", xBin, mag(xBin));


        // Get in region order the alpha*volume and diameter
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        const labelList sortedRegions = allRegionAlphaVolume.sortedToc();

        scalarField sortedVols
        (
            extractData
            (
                sortedRegions,
                allRegionAlphaVolume
            )
        );

        vectorField centroids(sortedVols.size(), Zero);

        // Check if downstream bins are calculated
        if (isoPlanes_)
        {
            vectorField alphaDistance
            (
                (alpha.primitiveField()*mesh_.V())
               *(mesh_.C().primitiveField() - origin_)()
            );

            Map<vector> allRegionAlphaDistance
            (
                regionSum
                (
                    regions,
                    alphaDistance
                )
            );

            // 2. centroid
            vectorField sortedMoment
            (
                extractData
                (
                    sortedRegions,
                    allRegionAlphaDistance
                )
            );

            centroids = sortedMoment/sortedVols + origin_;

            // Bin according to centroid
            scalarField distToPlane((centroids - origin_) & direction_);

            vectorField radialDistToOrigin
            (
                (centroids - origin_) - (distToPlane*direction_)
            );

            const scalar deltaX = maxDownstream_/nDownstreamBins_;
            labelList downstreamIndices(distToPlane.size(), -1);
            forAll(distToPlane, i)
            {
                if
                (
                    (mag(radialDistToOrigin[i]) < maxDiameter_)
                 && (distToPlane[i] < maxDownstream_)
                )
                {
                    downstreamIndices[i] = distToPlane[i]/deltaX;
                }
            }

            scalarField binDownCount(nDownstreamBins_, Zero);
            forAll(distToPlane, i)
            {
                if (downstreamIndices[i] != -1)
                {
                    binDownCount[downstreamIndices[i]] += 1.0;
                }
            }

            // Write
            if (Pstream::master())
            {
                // Construct mids of bins for plotting
                pointField xBin(nDownstreamBins_);

                scalar x = 0.5*deltaX;
                forAll(xBin, i)
                {
                    xBin[i] = point(x, 0, 0);
                    x += deltaX;
                }

                const coordSet coords("distance", "x", xBin, mag(xBin));
                writeGraph(coords, "isoPlanes", binDownCount);
            }

            // Write to log
            if (log)
            {
                Info<< "    Iso-planes Bins:" << nl
                    << "    " << token::TAB << "Bin"
                    << token::TAB << "Min distance"
                    << token::TAB << "Count:"
                    << nl;

                scalar delta = 0.0;
                forAll(binDownCount, bini)
                {
                    Info<< "    " << token::TAB << bini
                        << token::TAB << delta
                        << token::TAB << binDownCount[bini] << nl;
                    delta += deltaX;
                }
                Info<< endl;

            }
        }

        // Calculate the diameters
        scalarField sortedDiameters(sortedVols.size());
        forAll(sortedDiameters, i)
        {
            sortedDiameters[i] = Foam::cbrt
            (
                sortedVols[i]
               *6/constant::mathematical::pi
            );
        }

        // Determine the bin index for all the diameters
        labelList indices(sortedDiameters.size());
        forAll(sortedDiameters, i)
        {
            indices[i] = (sortedDiameters[i]-minDiam_)/delta;
        }

        // Calculate the counts per diameter bin
        scalarField binCount(nBins_, Zero);
        forAll(sortedDiameters, i)
        {
            binCount[indices[i]] += 1.0;
        }

        // Write counts
        if (Pstream::master())
        {
            writeGraph(coords, "count", binCount);
        }

        // Write to log
        if (log)
        {
            Info<< "    Bins:" << nl
                << "    " << token::TAB << "Bin"
                << token::TAB << "Min diameter"
                << token::TAB << "Count:"
                << nl;

            scalar diam = 0.0;
            forAll(binCount, bini)
            {
                Info<< "    " << token::TAB << bini
                    << token::TAB << diam
                    << token::TAB << binCount[bini] << nl;

                diam += delta;
            }

            Info<< endl;
        }


        // Write average and deviation of droplet volume.
        writeGraphs
        (
            "volume",           // name of field
            indices,            // per region the bin index
            sortedVols,         // per region field data
            binCount,           // per bin number of regions
            coords              // graph data for bins
        );

        // Collect some more field
        {
            wordList scalarNames(obr_.names(volScalarField::typeName));

            const labelList selected(fields_.matching(scalarNames));

            for (const label fieldi : selected)
            {
                const word& fldName = scalarNames[fieldi];

                Log << "    Scalar field " << fldName << endl;

                const scalarField& fld = obr_.lookupObject
                <
                    volScalarField
                >(fldName).primitiveField();

                writeGraphs
                (
                    fldName,            // name of field
                    alphaVol*fld,       // per cell field data

                    regions,            // per cell the region(=droplet)
                    sortedRegions,      // valid regions in sorted order
                    1.0/sortedVols,     // per region normalisation

                    indices,            // index of bin for each region
                    binCount,           // per bin number of regions
                    coords              // graph data for bins
                );
            }
        }
        {
            wordList vectorNames(obr_.names(volVectorField::typeName));

            const labelList selected(fields_.matching(vectorNames));

            for (const label fieldi : selected)
            {
                const word& fldName = vectorNames[fieldi];

                Log << "    Vector field " << fldName << endl;

                vectorField fld = obr_.lookupObject
                <
                    volVectorField
                >(fldName).primitiveField();

                if (csysPtr_)
                {
                    Log << "Transforming vector field " << fldName
                        << " with coordinate system "
                        << csysPtr_->name()
                        << endl;

                    fld = csysPtr_->localVector(fld);
                }


                // Components

                for (direction cmp = 0; cmp < vector::nComponents; cmp++)
                {
                    writeGraphs
                    (
                        fldName + vector::componentNames[cmp],
                        alphaVol*fld.component(cmp),// per cell field data

                        regions,        // per cell the region(=droplet)
                        sortedRegions,  // valid regions in sorted order
                        1.0/sortedVols, // per region normalisation

                        indices,        // index of bin for each region
                        binCount,       // per bin number of regions
                        coords          // graph data for bins
                    );
                }

                // Magnitude
                writeGraphs
                (
                    fldName + "mag",    // name of field
                    alphaVol*mag(fld),  // per cell field data

                    regions,            // per cell the region(=droplet)
                    sortedRegions,      // valid regions in sorted order
                    1.0/sortedVols,     // per region normalisation

                    indices,            // index of bin for each region
                    binCount,           // per bin number of regions
                    coords              // graph data for bins
                );
            }
        }
    }

    return true;
}


// ************************************************************************* //
