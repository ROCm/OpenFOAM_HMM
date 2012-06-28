/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "regionSizeDistribution.H"
#include "volFields.H"
#include "regionSplit.H"
#include "fvcVolumeIntegrate.H"
#include "Histogram.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSizeDistribution, 0);

    //- plus op for FixedList<scalar>
    template<class T, unsigned Size>
    class ListPlusEqOp
    {
        public:
        void operator()
        (
            FixedList<T, Size>& x,
            const FixedList<T, Size>& y
        ) const
        {
            forAll(x, i)
            {
                x[i] += y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSizeDistribution::regionSizeDistribution
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    alphaName_(dict.lookup("field")),
    patchNames_(dict.lookup("patches"))
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "regionSizeDistribution::regionSizeDistribution"
            "(const objectRegistry&, const dictionary&)"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSizeDistribution::~regionSizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSizeDistribution::read(const dictionary& dict)
{
    if (active_)
    {
        dict.lookup("field") >> alphaName_;
        dict.lookup("patches") >> patchNames_;
        dict.lookup("threshold") >> threshold_;
        dict.lookup("volFraction") >> volFraction_;
        dict.lookup("nBins") >> nBins_;

        word format(dict.lookup("setFormat"));
        formatterPtr_ = writer<scalar>::New(format);
    }
}


void Foam::regionSizeDistribution::execute()
{
    // Do nothing - only valid on write
}


void Foam::regionSizeDistribution::end()
{
    // Do nothing - only valid on write
}


void Foam::regionSizeDistribution::write()
{
    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        autoPtr<volScalarField> alphaPtr;
        if (obr_.foundObject<volScalarField>(alphaName_))
        {
            Info<< "Looking up field " << alphaName_ << endl;
        }
        else
        {
            Info<< "Reading field " << alphaName_ << endl;
            alphaPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        alphaName_,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }


        const volScalarField& alpha =
        (
             alphaPtr.valid()
           ? alphaPtr()
           : obr_.lookupObject<volScalarField>(alphaName_)
        );

        Info<< "Volume of alpha = "
            << fvc::domainIntegrate(alpha).value()
            << endl;

        const scalar meshVol = gSum(mesh.V());
        Info<< "Mesh volume = " << meshVol << endl;
        Info<< "Background region volume limit = " << volFraction_*meshVol
            << endl;


        // Determine blocked faces
        boolList blockedFace(mesh.nFaces(), false);
        label nBlocked = 0;

        {
            for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                scalar ownVal = alpha[mesh.faceOwner()[faceI]];
                scalar neiVal = alpha[mesh.faceNeighbour()[faceI]];

                if
                (
                    (ownVal < threshold_ && neiVal > threshold_)
                 || (ownVal > threshold_ && neiVal < threshold_)
                )
                {
                    blockedFace[faceI] = true;
                    nBlocked++;
                }
            }

            // Block coupled faces
            forAll(alpha.boundaryField(), patchI)
            {
                const fvPatchScalarField& fvp = alpha.boundaryField()[patchI];
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


        regionSplit regions(mesh, blockedFace);

        Info<< "Determined " << regions.nRegions() << " disconnected regions"
            << endl;


        if (debug)
        {
            volScalarField region
            (
                IOobject
                (
                    "region",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0)
            );
            Info<< "Dumping region as volScalarField to " << region.name()
                << endl;

            forAll(regions, cellI)
            {
                region[cellI] = regions[cellI];
            }
            region.correctBoundaryConditions();
            region.write();
        }


        // Sum all regions
        Map<Pair<scalar> > regionVolume(regions.nRegions()/Pstream::nProcs());
        forAll(alpha, cellI)
        {
            scalar cellVol = mesh.V()[cellI];
            scalar alphaVol = alpha[cellI]*cellVol;

            label regionI = regions[cellI];

            Map<Pair<scalar> >::iterator fnd = regionVolume.find(regionI);
            if (fnd == regionVolume.end())
            {
                regionVolume.insert
                (
                    regionI,
                    Pair<scalar>(cellVol, alphaVol)
                );
            }
            else
            {
                fnd().first() += cellVol;
                fnd().second() += alphaVol;
            }
        }
        Pstream::mapCombineGather(regionVolume, ListPlusEqOp<scalar, 2>());
        Pstream::mapCombineScatter(regionVolume);


        if (debug)
        {
            Info<< token::TAB << "Region"
                << token::TAB << "Volume(mesh)"
                << token::TAB << "Volume(" << alpha.name() << "):"
                << endl;
            scalar meshSumVol = 0.0;
            scalar alphaSumVol = 0.0;

            forAllConstIter(Map<Pair<scalar> >, regionVolume, iter)
            {
                Info<< token::TAB << iter.key()
                    << token::TAB << iter().first()
                    << token::TAB << iter().second() << endl;

                meshSumVol += iter().first();
                alphaSumVol += iter().second();
            }
            Info<< token::TAB << "Total:"
                << token::TAB << meshSumVol
                << token::TAB << alphaSumVol << endl;
            Info<< endl;
        }



        // Mark all regions starting at patches
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Count number of patch faces (just for initial sizing)
        label nPatchFaces = 0;
        forAll(patchNames_, i)
        {
            const word& pName = patchNames_[i];
            label patchI = mesh.boundaryMesh().findPatchID(pName);
            if (patchI == -1)
            {
                WarningIn("regionSizeDistribution::write()")
                    << "Cannot find patch " << pName << ". Valid patches are "
                    << mesh.boundaryMesh().names()
                    << endl;
            }
            else
            {
                nPatchFaces += mesh.boundaryMesh()[patchI].size();
            }
        }

        Map<label> keepRegions(nPatchFaces);
        forAll(patchNames_, i)
        {
            const word& pName = patchNames_[i];

            label patchI = mesh.boundaryMesh().findPatchID(pName);
            if (patchI != -1)
            {
                const polyPatch& pp = mesh.boundaryMesh()[patchI];

                // Collect all regions on the patch
                const labelList& faceCells = pp.faceCells();

                forAll(faceCells, i)
                {
                    keepRegions.insert
                    (
                        regions[faceCells[i]],
                        Pstream::myProcNo()
                    );
                }
            }
        }


        // Make sure all the processors have the same set of regions
        Pstream::mapCombineGather(keepRegions, minEqOp<label>());
        Pstream::mapCombineScatter(keepRegions);

        Info<< "Patch connected regions (liquid core):" << endl;
        forAllConstIter(Map<label>, keepRegions, iter)
        {
            label regionI = iter.key();
            Pair<scalar>& vols = regionVolume[regionI];
            Info<< token::TAB << iter.key()
                << token::TAB << vols.first()
                << token::TAB << vols.second() << endl;

        }
        Info<< endl;

        Info<< "Background regions:" << endl;
        forAllConstIter(Map<Pair<scalar> >, regionVolume, iter)
        {
            if
            (
               !keepRegions.found(iter.key())
             && iter().first() >= volFraction_*meshVol
            )
            {
                Info<< token::TAB << iter.key()
                    << token::TAB << iter().first()
                    << token::TAB << iter().second() << endl;
            }
        }
        Info<< endl;


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
                alphaName_ + "_liquidCore",
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
                alphaName_ + "_background",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            alpha,
            fvPatchField<scalar>::calculatedType()
        );


        // Knock out any cell not in keepRegions
        forAll(liquidCore, cellI)
        {
            label regionI = regions[cellI];
            if (keepRegions.found(regionI))
            {
                backgroundAlpha[cellI] = 0;
            }
            else
            {
                liquidCore[cellI] = 0;

                scalar regionVol = regionVolume[regionI].first();
                if (regionVol < volFraction_*meshVol)
                {
                    backgroundAlpha[cellI] = 0;
                }
            }
        }
        liquidCore.correctBoundaryConditions();
        backgroundAlpha.correctBoundaryConditions();

        Info<< "Volume of liquid-core = "
            << fvc::domainIntegrate(liquidCore).value()
            << endl;

        Info<< "Writing liquid-core field to " << liquidCore.name() << endl;
        liquidCore.write();

        Info<< "Volume of background = "
            << fvc::domainIntegrate(backgroundAlpha).value()
            << endl;

        Info<< "Writing background field to " << backgroundAlpha.name() << endl;
        backgroundAlpha.write();



        // Collect histogram
        if (Pstream::master())
        {
            DynamicList<scalar> diameters(regionVolume.size());
            forAllConstIter(Map<Pair<scalar> >, regionVolume, iter)
            {
                if (!keepRegions.found(iter.key()))
                {
                    if (iter().first() < volFraction_*meshVol)
                    {
                        scalar v = iter().second();
                      //scalar diam = Foam::cbrt(v*6/mathematicalConstant::pi);
                        scalar diam =
                            Foam::cbrt(v*6/constant::mathematical::pi);
                        diameters.append(diam);
                    }
                }
            }

            if (diameters.size())
            {
                scalar maxDiam = max(diameters);
                scalar minDiam = 0.0;

                Info<< "Maximum diameter:" << maxDiam << endl;

                Histogram<List<scalar> > bins
                (
                    minDiam,
                    maxDiam,
                    nBins_,
                    diameters
                );

                /* 1.7.x
                scalarField xBin(nBins_);

                scalar dx = (maxDiam-minDiam)/nBins_;
                scalar x = 0.5*dx;
                forAll(bins.counts(), i)
                {
                    xBin[i] = x;
                    x += dx;
                }

                scalarField normalisedCount(bins.counts().size());
                forAll(bins.counts(), i)
                {
                    normalisedCount[i] = 1.0*bins.counts()[i];
                }

                const coordSet coords
                (
                    "diameter",
                    "x",
                    xBin
                );
                */

                pointField xBin(nBins_);
                scalar dx = (maxDiam - minDiam)/nBins_;
                scalar x = 0.5*dx;
                forAll(bins.counts(), i)
                {
                    xBin[i] = point(x, 0, 0);
                    x += dx;
                }

                scalarField normalisedCount(bins.counts().size());
                forAll(bins.counts(), i)
                {
                    normalisedCount[i] = 1.0*bins.counts()[i];
                }

                const coordSet coords
                (
                    "diameter",
                    "x",
                    xBin,
                    mag(xBin)
                );
                const wordList valNames(1, "count");


                fileName outputPath;
                if (Pstream::parRun())
                {
                    outputPath = mesh.time().path()/".."/name_;
                }
                else
                {
                    outputPath = mesh.time().path()/name_;
                }

                if (mesh.name() != fvMesh::defaultRegion)
                {
                    outputPath = outputPath/mesh.name();
                }

                mkDir(outputPath/mesh.time().timeName());
                OFstream str
                (
                    outputPath
                  / mesh.time().timeName()
                  / formatterPtr_().getFileName(coords, valNames)
                );
                Info<< "Writing distribution to " << str.name() << endl;

                List<const scalarField*> valPtrs(1);
                valPtrs[0] = &normalisedCount;
                formatterPtr_().write(coords, valNames, valPtrs, str);
            }
        }
    }
}


// ************************************************************************* //
