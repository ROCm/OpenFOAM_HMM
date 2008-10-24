/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "CentredFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedCentredStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::CentredFitData<Polynomial>::CentredFitData
(
    const fvMesh& mesh,
    const extendedCentredStencil& stencil,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    MeshObject<fvMesh, CentredFitData<Polynomial> >(mesh),
    stencil_(stencil),
    linearLimitFactor_(linearLimitFactor),
    centralWeight_(centralWeight),
#   ifdef SPHERICAL_GEOMETRY
    dim_(2),
#   else
    dim_(mesh.nGeometricD()),
#   endif
    minSize_(Polynomial::nTerms(dim_)),
    coeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing CentredFitData<Polynomial>" << endl;
    }

    // Check input
    if (linearLimitFactor > 1)
    {
        FatalErrorIn("CentredFitData<Polynomial>::CentredFitData")
            << "linearLimitFactor requested = " << linearLimitFactor
            << " should not be less than one"
            << exit(FatalError);
    }

    calcFit();

    if (debug)
    {
        Info<< "CentredFitData<Polynomial>::CentredFitData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::CentredFitData<Polynomial>::findFaceDirs
(
    vector& idir,        // value changed in return
    vector& jdir,        // value changed in return
    vector& kdir,        // value changed in return
    const fvMesh& mesh,
    const label facei
)
{
    idir = mesh.faceAreas()[facei];
    idir /= mag(idir);

#   ifndef SPHERICAL_GEOMETRY
    if (mesh.nGeometricD() <= 2) // find the normal direction
    {
        if (mesh.directions()[0] == -1)
        {
            kdir = vector(1, 0, 0);
        }
        else if (mesh.directions()[1] == -1)
        {
            kdir = vector(0, 1, 0);
        }
        else
        {
            kdir = vector(0, 0, 1);
        }
    }
    else // 3D so find a direction in the plane of the face
    {
        const face& f = mesh.faces()[facei];
        kdir = mesh.points()[f[0]] - mesh.faceCentres()[facei];
    }
#   else
    // Spherical geometry so kdir is the radial direction
    kdir = mesh.faceCentres()[facei];
#   endif

    if (mesh.nGeometricD() == 3)
    {
        // Remove the idir component from kdir and normalise
        kdir -= (idir & kdir)*idir;

        scalar magk = mag(kdir);

        if (magk < SMALL)
        {
            FatalErrorIn("findFaceDirs") << " calculated kdir = zero"
                << exit(FatalError);
        }
        else
        {
            kdir /= magk;
        }
    }

    jdir = kdir ^ idir;
}


template<class Polynomial>
void Foam::CentredFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the cell/face centres in stencil order.
    // Centred face stencils no good for triangles of tets.
    // Need bigger stencils
    List<List<point> > stencilPoints(mesh.nFaces());
    stencil_.collectData(mesh.C(), stencilPoints);

    // find the fit coefficients for every face in the mesh

    const surfaceScalarField& w = this->mesh().surfaceInterpolation::weights();

    for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        calcFit(stencilPoints[facei], w[facei], facei);
    }

    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().patch().start();

            forAll(pw, i)
            {
                calcFit(stencilPoints[facei], pw[i], facei);
                facei++;
            }
        }
    }
}


template<class Polynomial>
Foam::label Foam::CentredFitData<Polynomial>::calcFit
(
    const List<point>& C,
    const scalar wLin,
    const label facei
)
{
    vector idir(1,0,0);
    vector jdir(0,1,0);
    vector kdir(0,0,1);
    findFaceDirs(idir, jdir, kdir, this->mesh(), facei);

    // Setup the point weights
    scalarList wts(C.size(), scalar(1));
    wts[0] = centralWeight_;
    wts[1] = centralWeight_;

    // Reference point
    point p0 = this->mesh().faceCentres()[facei];

    // p0 -> p vector in the face-local coordinate system
    vector d;

    // Local coordinate scaling
    scalar scale = 1;

    // Matrix of the polynomial components
    scalarRectangularMatrix B(C.size(), minSize_, scalar(0));

    for(label ip = 0; ip < C.size(); ip++)
    {
        const point& p = C[ip];

        d.x() = (p - p0)&idir;
        d.y() = (p - p0)&jdir;
#       ifndef SPHERICAL_GEOMETRY
        d.z() = (p - p0)&kdir;
#       else
        d.z() = mag(p) - mag(p0);
#       endif

        if (ip == 0)
        {
            scale = cmptMax(cmptMag((d)));
        }

        // Scale the radius vector
        d /= scale;

        Polynomial::addCoeffs
        (
            B[ip],
            d,
            wts[ip],
            dim_
        );
    }

    // Set the fit
    label stencilSize = C.size();
    coeffs_[facei].setSize(stencilSize);
    scalarList singVals(minSize_);
    label nSVDzeros = 0;

    bool goodFit = false;
    for(int iIt = 0; iIt < 10 && !goodFit; iIt++)
    {
        SVD svd(B, SMALL);

        scalar fit0 = wts[0]*svd.VSinvUt()[0][0];
        scalar fit1 = wts[1]*svd.VSinvUt()[0][1];

        goodFit =
            (mag(fit0 - wLin) < linearLimitFactor_*wLin)
         && (mag(fit1 - (1 - wLin)) < linearLimitFactor_*(1 - wLin));

        //scalar w0Err = fit0/wLin;
        //scalar w1Err = fit1/(1 - wLin);

        //goodFit =
        //    (w0Err > linearLimitFactor_ && w0Err < (1 + linearLimitFactor_))
        // && (w1Err > linearLimitFactor_ && w1Err < (1 + linearLimitFactor_));

        if (goodFit)
        {
            coeffs_[facei][0] = fit0;
            coeffs_[facei][1] = fit1;

            for(label i=2; i<stencilSize; i++)
            {
                coeffs_[facei][i] = wts[i]*svd.VSinvUt()[0][i];
            }

            singVals = svd.S();
            nSVDzeros = svd.nZeros();
        }
        else // (not good fit so increase weight in the centre and for linear)
        {
            wts[0] *= 10;
            wts[1] *= 10;

            for(label j = 0; j < B.m(); j++)
            {
                B[0][j] *= 10;
                B[1][j] *= 10;
            }
        }
    }

    // static const scalar L = 0.1;
    // static const scalar R = 0.2;

    // static const scalar beta = 1.0/(R - L);
    // static const scalar alpha = R*beta;

    if (goodFit)
    {
        // scalar limiter =
        // max
        // (
        //     min
        //     (
        //         min(alpha - beta*mag(coeffs_[facei][0] - wLin)/wLin, 1),
        //         min(alpha - beta*mag(coeffs_[facei][1] - (1 - wLin))
        //       /(1 - wLin), 1)
        //     ), 0
        // );

        //Info<< wLin << " " << coeffs_[facei][0]
        //    << " " << (1 - wLin) << " " << coeffs_[facei][1] << endl;

        // Remove the uncorrected linear ocefficients
        coeffs_[facei][0] -= wLin;
        coeffs_[facei][1] -= 1 - wLin;

        // if (limiter < 0.99)
        // {
        //     for(label i = 0; i < stencilSize; i++)
        //     {
        //         coeffs_[facei][i] *= limiter;
        //     }
        // }
    }
    else
    {
        if (debug)
        {
            WarningIn
            (
                "CentredFitData<Polynomial>::calcFit"
                "(const List<point>& C, const label facei"
            )   << "Could not fit face " << facei
                << ", reverting to linear." << nl
                << "    Weights "
                << coeffs_[facei][0] << " " << wLin << nl
                << "    Linear weights "
                << coeffs_[facei][1] << " " << 1 - wLin << endl;
        }

        coeffs_[facei] = 0;
    }

    return minSize_ - nSVDzeros;
}


template<class Polynomial>
bool Foam::CentredFitData<Polynomial>::movePoints()
{
    calcFit();
    return true;
}


// ************************************************************************* //
