/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "interpolationPoint.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolationPoint<Type>::interpolationPoint
(
    const GeometricField<Type, fvPatchField, volMesh>& psi
)
:
    interpolation<Type>(psi),
    psip_(volPointInterpolation::New(psi.mesh()).interpolate(psi))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//void Foam::interpolationPoint<Type>::calcWeights
//(
//    const vector& position,
//    const label cellI,
//    const label faceI,
//    scalarField& weights
//) const
//{
//    const polyMesh& mesh = this->pMesh_;
//    const pointField& points = mesh.points();
//
//
//    const scalar eps = 0.00000001;
//
//
//    // Addressing - face vertices to local points
//    const labelList& toGlobal = mesh.cellPoints()[cellI];
//    Map<label> toLocal(2*toGlobal.size());
//    forAll(toGlobal, i)
//    {
//        toLocal.insert(toGlobal[i], i);
//    }
//
//    // Initialise weights
//    weights.setSize(toGlobal.size());
//    weights = 0.0;
//
//    // Point-to-vertex vectors and distances
//    scalarField dist(toGlobal.size());
//    vectorField uVec(toGlobal.size());
//    forAll(toGlobal, pid)
//    {
//        const point& pt = points[toGlobal[pid]];//-cc;
//        uVec[pid] = pt-position;
//        dist[pid] = mag(uVec[pid]);
//
//        // Special case: point is close to vertex
//        if (dist[pid] < eps)
//        {
//            weights[pid] = 1.0;
//            return;
//        }
//    }
//
//    // Project onto unit sphere
//    uVec /= dist;
//
//
//    // Loop over all triangles of all polygons of cell to compute weights
//    DynamicList<scalar> alpha(100);
//    DynamicList<scalar> theta(100);
//
//    const cell& cFaces = mesh.cells()[cellI];
//
//    forAll(cFaces, iter)
//    {
//        label faceI = cFaces[iter];
//        const face& f = mesh.faces()[faceI];
//
//        Pout<< "face:" << faceI << " at:"
//            << pointField(mesh.points(), f)
//            << endl;
//
//        vector v(point::zero);
//        forAll(f, j)
//        {
//            label jPlus1 = f.fcIndex(j);
//            const point& uj = points[f[j]];//-cc;
//            const point& ujPlus1 = points[f[jPlus1]];//-cc;
//            Pout<< "    uj:" << uj << " ujPlus1:" << ujPlus1 << endl;
//
//            vector temp = uj ^ ujPlus1;
//            temp /= mag(temp);
//
//            scalar l = mag(uj-ujPlus1);
//            scalar angle = 2.0*Foam::asin(l/2.0);
//
//            v += 0.5*angle*temp;
//        }
//
//        scalar vNorm = mag(v);
//        v /= vNorm;
//
//        // Make sure v points towards the polygon
//        if ((v&points[f[0]]) < 0)
//        {
//            v = -v;
//        }
//
//        Pout<< "    v:" << v << endl;
//
//        // angles between edges
//        forAll(f, j)
//        {
//            label jPlus1 = f.fcIndex(j);
//            const point& uj = points[f[j]];//-cc;
//            const point& ujPlus1 = points[f[jPlus1]];//-cc;
//            Pout<< "    uj:" << uj << " ujPlus1:" << ujPlus1 << endl;
//
//            vector n0 = uj ^ v;
//            n0 /= mag(n0);
//            vector n1 = ujPlus1 ^ v;
//            n1 /= mag(n1);
//
//            scalar l = mag(n0-n1);
//            Pout<< "    l:" << l << endl;
//            alpha(j) = 2.0*Foam::asin(l/2.0);
//
//            vector temp = n0 ^ n1;
//            if ((temp&v) < 0.0)
//            {
//                alpha(j) = -alpha(j);
//            }
//
//            l = mag(uj-v);
//            Pout<< "    l:" << l << endl;
//            theta(j) = 2.0*Foam::asin(l/2.0);
//        }
//
//
//        bool outlierFlag = false;
//        forAll(f, j)
//        {
//            if (mag(theta(j)) < eps)
//            {
//                outlierFlag = true;
//
//                label pid = toLocal[f[j]];
//                weights[pid] += vNorm / dist[pid];
//                break;
//            }
//        }
//
//        if (outlierFlag)
//        {
//            continue;
//        }
//
//        scalar sum = 0.0;
//        forAll(f, j)
//        {
//            label jMin1 = f.rcIndex(j);
//            sum +=
//                1.0
//              / Foam::tan(theta(j))
//              * (Foam::tan(alpha(j)/2.0)+Foam::tan(alpha(jMin1)/2.0));
//        }
//
//        // The special case when x lies on the polygon, handle it using 2D mvc.
//        // In the 2D case, alpha = theta
//        if (mag(sum) < eps)
//        {
//            weights = 0.0;
//
//            // recompute theta, the theta computed previously are not robust
//            forAll(f, j)
//            {
//                label jPlus1 = f.fcIndex(j);
//                const point& uj = points[f[j]];//-cc;
//                const point& ujPlus1 = points[f[jPlus1]];//-cc;
//                scalar l = mag(uj-ujPlus1);
//                theta(j) = 2.0*Foam::asin(l/2.0);
//            }
//
//            scalar sumWeight = 0;
//            forAll(f, j)
//            {
//                label pid = toLocal[f[j]];
//                label jMin1 = f.rcIndex(j);
//                weights[pid] =
//                    1.0
//                  / dist[pid]
//                  * (Foam::tan(theta(jMin1)/2.0)+Foam::tan(theta(j)/2.0));
//                sumWeight += weights[pid];
//            }
//
//            if (sumWeight < eps)
//            {
//                return;
//            }
//            weights /= sumWeight;
//            return;
//        }
//
//
//        // Normal 3D case
//        forAll(f, j)
//        {
//            label pid = toLocal[f[j]];
//            label jMin1 = f.rcIndex(j);
//            weights[pid] +=
//                vNorm
//              / sum
//              / dist[pid]
//              / Foam::sin(theta(j))
//              * (Foam::tan(alpha(j)/2.0)+Foam::tan(alpha(jMin1)/2.0));
//        }
//    }
//
//    // normalise weights
//    scalar sumWeight = sum(weights);
//
//    if (mag(sumWeight) < eps)
//    {
//        return;
//    }
//    weights /= sumWeight;
//}


// ************************************************************************* //
