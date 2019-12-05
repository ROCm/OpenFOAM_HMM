/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "unitConversion.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "medialAxisMeshMover.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const label n,
    const scalar totalOverFirst
) const
{
    if (n <= 1)
    {
        return 1.0;
    }

    const label maxIters = 20;
    const scalar tol = 1e-8;

    if (mag(n-totalOverFirst) < tol)
    {
        return 1.0;
    }

    // Calculate the bounds of the solution
    scalar minR;
    scalar maxR;

    if (totalOverFirst < n)
    {
        minR = 0.0;
        maxR = pow(totalOverFirst/n, 1/(n-1));
    }
    else
    {
        minR = pow(totalOverFirst/n, 1/(n-1));
        maxR = totalOverFirst/(n - 1);
    }

    // Starting guess
    scalar r = 0.5*(minR + maxR);

    for (label i = 0; i < maxIters; ++i)
    {
        const scalar prevr = r;

        const scalar fx = pow(r, n) - totalOverFirst*r - (1 - totalOverFirst);
        const scalar dfx = n*pow(r, n - 1) - totalOverFirst;
        r -= fx/dfx;

        if (mag(r - prevr) < tol)
        {
            break;
        }
    }
    return r;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh,
    const bool dryRun
)
:
    dict_(dict),
    numLayers_(boundaryMesh.size(), -1),
    relativeSizes_(meshRefinement::get<bool>(dict, "relativeSizes", dryRun)),
    layerSpec_(ILLEGAL),
    firstLayerThickness_(boundaryMesh.size(), -123),
    finalLayerThickness_(boundaryMesh.size(), -123),
    thickness_(boundaryMesh.size(), -123),
    expansionRatio_(boundaryMesh.size(), -123),
    minThickness_
    (
        boundaryMesh.size(),
        meshRefinement::get<scalar>(dict, "minThickness", dryRun)
    ),
    featureAngle_(meshRefinement::get<scalar>(dict, "featureAngle", dryRun)),
    mergePatchFacesAngle_
    (
        dict.lookupOrDefault<scalar>
        (
            "mergePatchFacesAngle",
            featureAngle_
        )
    ),
    concaveAngle_
    (
        dict.lookupOrDefault("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(meshRefinement::get<label>(dict, "nGrow", dryRun)),
    maxFaceThicknessRatio_
    (
        meshRefinement::get<scalar>(dict, "maxFaceThicknessRatio", dryRun)
    ),
    nBufferCellsNoExtrude_
    (
        meshRefinement::get<label>(dict, "nBufferCellsNoExtrude", dryRun)
    ),
    nLayerIter_(meshRefinement::get<label>(dict, "nLayerIter", dryRun)),
    nRelaxedIter_(labelMax),
    additionalReporting_(dict.lookupOrDefault("additionalReporting", false)),
    meshShrinker_
    (
        dict.lookupOrDefault
        (
            "meshShrinker",
            medialAxisMeshMover::typeName
        )
    ),
    dryRun_(dryRun)
{
    // Detect layer specification mode

    label nSpec = 0;

    bool haveFirst = dict.found("firstLayerThickness");
    if (haveFirst)
    {
        firstLayerThickness_ = scalarField
        (
            boundaryMesh.size(),
            dict.get<scalar>("firstLayerThickness")
        );
        nSpec++;
    }
    bool haveFinal = dict.found("finalLayerThickness");
    if (haveFinal)
    {
        finalLayerThickness_ = scalarField
        (
            boundaryMesh.size(),
            dict.get<scalar>("finalLayerThickness")
        );
        nSpec++;
    }
    bool haveTotal = dict.found("thickness");
    if (haveTotal)
    {
        thickness_ = scalarField
        (
            boundaryMesh.size(),
            dict.get<scalar>("thickness")
        );
        nSpec++;
    }
    bool haveExp = dict.found("expansionRatio");
    if (haveExp)
    {
        expansionRatio_ = scalarField
        (
            boundaryMesh.size(),
            dict.get<scalar>("expansionRatio")
        );
        nSpec++;
    }


    if (haveFirst && haveTotal)
    {
        layerSpec_ = FIRST_AND_TOTAL;
        Info<< "Layer thickness specified as first layer and overall thickness."
            << endl;
    }
    else if (haveFirst && haveExp)
    {
        layerSpec_ = FIRST_AND_EXPANSION;
        Info<< "Layer thickness specified as first layer and expansion ratio."
            << endl;
    }
    else if (haveFinal && haveTotal)
    {
        layerSpec_ = FINAL_AND_TOTAL;
        Info<< "Layer thickness specified as final layer and overall thickness."
            << endl;
    }
    else if (haveFinal && haveExp)
    {
        layerSpec_ = FINAL_AND_EXPANSION;
        Info<< "Layer thickness specified as final layer and expansion ratio."
            << endl;
    }
    else if (haveTotal && haveExp)
    {
        layerSpec_ = TOTAL_AND_EXPANSION;
        Info<< "Layer thickness specified as overall thickness"
            << " and expansion ratio." << endl;
    }


    if (layerSpec_ == ILLEGAL || nSpec != 2)
    {
        FatalIOErrorInFunction(dict)
            << "Over- or underspecified layer thickness."
            << " Please specify" << nl
            << "    first layer thickness ('firstLayerThickness')"
            << " and overall thickness ('thickness') or" << nl
            << "    first layer thickness ('firstLayerThickness')"
            << " and expansion ratio ('expansionRatio') or" << nl
            << "    final layer thickness ('finalLayerThickness')"
            << " and expansion ratio ('expansionRatio') or" << nl
            << "    final layer thickness ('finalLayerThickness')"
            << " and overall thickness ('thickness') or" << nl
            << "    overall thickness ('thickness')"
            << " and expansion ratio ('expansionRatio'"
            << exit(FatalIOError);
    }


    dict.readIfPresent("nRelaxedIter", nRelaxedIter_);

    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Layer iterations should be >= 0." << nl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalIOError);
    }


    const dictionary& layersDict = meshRefinement::subDict
    (
        dict,
        "layers",
        dryRun
    );

    for (const entry& dEntry : layersDict)
    {
        if (dEntry.isDict())
        {
            const keyType& key = dEntry.keyword();
            const dictionary& layerDict = dEntry.dict();

            const labelHashSet patchIDs
            (
                boundaryMesh.patchSet(List<wordRe>(1, wordRe(key)))
            );

            if (patchIDs.size() == 0)
            {
                IOWarningInFunction(layersDict)
                    << "Layer specification for " << key
                    << " does not match any patch." << endl
                    << "Valid patches are " << boundaryMesh.names() << endl;
            }
            else
            {
                for (const label patchi : patchIDs)
                {
                    numLayers_[patchi] =
                        layerDict.get<label>("nSurfaceLayers");

                    switch (layerSpec_)
                    {
                        case FIRST_AND_TOTAL:
                            layerDict.readIfPresent
                            (
                                "firstLayerThickness",
                                firstLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                        break;

                        case FIRST_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "firstLayerThickness",
                                firstLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        case FINAL_AND_TOTAL:
                            layerDict.readIfPresent
                            (
                                "finalLayerThickness",
                                finalLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                        break;

                        case FINAL_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "finalLayerThickness",
                                finalLayerThickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        case TOTAL_AND_EXPANSION:
                            layerDict.readIfPresent
                            (
                                "thickness",
                                thickness_[patchi]
                            );
                            layerDict.readIfPresent
                            (
                                "expansionRatio",
                                expansionRatio_[patchi]
                            );
                        break;

                        default:
                            FatalIOErrorInFunction(dict)
                                << "problem." << exit(FatalIOError);
                        break;
                    }

                    layerDict.readIfPresent
                    (
                        "minThickness",
                        minThickness_[patchi]
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerThickness
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_TOTAL:
        case FINAL_AND_TOTAL:
        case TOTAL_AND_EXPANSION:
        {
            return totalThickness;
        }
        break;

        case FIRST_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < SMALL)
            {
                return firstLayerThickess * nLayers;
            }
            else
            {
                return firstLayerThickess
                   *(1.0 - pow(expansionRatio, nLayers))
                   /(1.0 - expansionRatio);
            }
        }
        break;

        case FINAL_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < SMALL)
            {
                return finalLayerThickess * nLayers;
            }
            else
            {
                scalar invExpansion = 1.0 / expansionRatio;
                return finalLayerThickess
                   *(1.0 - pow(invExpansion, nLayers))
                   /(1.0 - invExpansion);
            }
        }
        break;

        default:
        {
            FatalErrorInFunction
                << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_EXPANSION:
        case FINAL_AND_EXPANSION:
        case TOTAL_AND_EXPANSION:
        {
            return expansionRatio;
        }
        break;

        case FIRST_AND_TOTAL:
        {
            return layerExpansionRatio
            (
                nLayers,
                totalThickness/firstLayerThickess
            );
        }
        break;

        case FINAL_AND_TOTAL:
        {
            return
                1.0
               /layerExpansionRatio
                (
                    nLayers,
                    totalThickness/finalLayerThickess
                );
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::firstLayerThickness
(
    const label nLayers,
    const scalar firstLayerThickess,
    const scalar finalLayerThickess,
    const scalar totalThickness,
    const scalar expansionRatio
) const
{
    switch (layerSpec_)
    {
        case FIRST_AND_EXPANSION:
        case FIRST_AND_TOTAL:
        {
            return firstLayerThickess;
        }

        case FINAL_AND_EXPANSION:
        {
            return finalLayerThickess*pow(1.0/expansionRatio, nLayers-1);
        }
        break;

        case FINAL_AND_TOTAL:
        {
            scalar r = layerExpansionRatio
            (
                nLayers,
                firstLayerThickess,
                finalLayerThickess,
                totalThickness,
                expansionRatio
            );
            return finalLayerThickess/pow(r, nLayers-1);
        }
        break;

        case TOTAL_AND_EXPANSION:
        {
            scalar r = finalLayerThicknessRatio
            (
                nLayers,
                expansionRatio
            );
            scalar finalThickness = r*totalThickness;
            return finalThickness/pow(expansionRatio, nLayers-1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::finalLayerThicknessRatio
(
    const label nLayers,
    const scalar expansionRatio
) const
{
    if (nLayers > 0)
    {
        if (mag(expansionRatio-1) < SMALL)
        {
            return 1.0/nLayers;
        }
        else
        {
            return
                pow(expansionRatio, nLayers - 1)
               *(1.0 - expansionRatio)
               /(1.0 - pow(expansionRatio, nLayers));
        }
    }
    else
    {
        return 0.0;
    }
}


// ************************************************************************* //
