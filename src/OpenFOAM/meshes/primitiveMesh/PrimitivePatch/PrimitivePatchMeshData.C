/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "PrimitivePatch.H"
#include "Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcMeshData() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMeshData() : "
               "calculating mesh data in PrimitivePatch"
            << endl;
    }

    if (meshPointsPtr_ || localFacesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "meshPointsPtr_ or localFacesPtr_ already allocated"
            << abort(FatalError);
    }

    // Create a map for marking points.  Estimated size is 4 times the
    // number of faces in the patch
    Map<label> markedPoints(4*this->size());


    // Important:
    // ~~~~~~~~~~
    // In <= 1.5 the meshPoints would be in increasing order but this gives
    // problems in processor point synchronisation where we have to find out
    // how the opposite side would have allocated points.

    ////- 1.5 code:
    //// if the point is used, set the mark to 1
    //forAll(*this, facei)
    //{
    //    const face_type& curPoints = this->operator[](facei);
    //
    //    forAll(curPoints, pointi)
    //    {
    //        markedPoints.insert(curPoints[pointi], -1);
    //    }
    //}
    //
    //// Create the storage and store the meshPoints.  Mesh points are
    //// the ones marked by the usage loop above
    //meshPointsPtr_.reset(new labelList(markedPoints.toc()));
    //auto& pointPatch = *meshPointsPtr_;
    //
    //// Sort the list to preserve compatibility with the old ordering
    //sort(pointPatch);
    //
    //// For every point in map give it its label in mesh points
    //forAll(pointPatch, pointi)
    //{
    //    markedPoints.find(pointPatch[pointi])() = pointi;
    //}

    //- Unsorted version:
    DynamicList<label> meshPoints(2*this->size());
    for (const face_type& f : *this)
    {
        for (const label pointi : f)
        {
            if (markedPoints.insert(pointi, meshPoints.size()))
            {
                meshPoints.append(pointi);
            }
        }
    }
    // Transfer to straight list (reuses storage)
    meshPointsPtr_.reset(new labelList(meshPoints, true));

    // Create local faces. Deep-copy original faces to retain additional
    // data (e.g. region number of labelledTri)
    // The vertices will be overwritten later
    localFacesPtr_.reset(new List<face_type>(*this));
    auto& locFaces = *localFacesPtr_;

    for (face_type& f : locFaces)
    {
        for (label& pointi : f)
        {
            pointi = *(markedPoints.cfind(pointi));
        }
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMeshData() : "
               "finished calculating mesh data in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcMeshPointMap() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMeshPointMap() : "
               "calculating mesh point map in PrimitivePatch"
            << endl;
    }

    if (meshPointMapPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "meshPointMapPtr_ already allocated"
            << abort(FatalError);
    }

    const labelList& mp = meshPoints();

    meshPointMapPtr_.reset(new Map<label>(2*mp.size()));
    auto& mpMap = *meshPointMapPtr_;

    forAll(mp, i)
    {
        mpMap.insert(mp[i], i);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMeshPointMap() : "
               "finished calculating mesh point map in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcLocalPoints() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcLocalPoints() : "
               "calculating localPoints in PrimitivePatch"
            << endl;
    }

    if (localPointsPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "localPointsPtr_ already allocated"
            << abort(FatalError);
    }

    const labelList& meshPts = meshPoints();

    localPointsPtr_.reset(new Field<point_type>(meshPts.size()));
    auto& locPts = *localPointsPtr_;

    forAll(meshPts, pointi)
    {
        locPts[pointi] = points_[meshPts[pointi]];
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
            << "calcLocalPoints() : "
            << "finished calculating localPoints in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcPointNormals() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcPointNormals() : "
               "calculating pointNormals in PrimitivePatch"
            << endl;
    }

    if (pointNormalsPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "pointNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    const auto& faceUnitNormals = faceNormals();

    const labelListList& pf = pointFaces();

    pointNormalsPtr_.reset(new Field<point_type>(meshPoints().size(), Zero));
    auto& n = *pointNormalsPtr_;

    forAll(pf, pointi)
    {
        point_type& curNormal = n[pointi];

        const labelList& curFaces = pf[pointi];

        for (const label facei : curFaces)
        {
            curNormal += faceUnitNormals[facei];
        }

        curNormal.normalise();
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcPointNormals() : "
               "finished calculating pointNormals in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcFaceCentres() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceCentres() : "
               "calculating faceCentres in PrimitivePatch"
            << endl;
    }

    if (faceCentresPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "faceCentresPtr_ already allocated"
            << abort(FatalError);
    }

    faceCentresPtr_.reset(new Field<point_type>(this->size()));
    auto& c = *faceCentresPtr_;

    forAll(c, facei)
    {
        c[facei] = this->operator[](facei).centre(points_);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceCentres() : "
               "finished calculating faceCentres in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcMagFaceAreas() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMagFaceAreas() : "
               "calculating magFaceAreas in PrimitivePatch"
            << endl;
    }

    if (magFaceAreasPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "magFaceAreasPtr_ already allocated"
            << abort(FatalError);
    }

    magFaceAreasPtr_.reset(new Field<scalar>(this->size()));
    auto& a = *magFaceAreasPtr_;

    forAll(a, facei)
    {
        a[facei] = this->operator[](facei).mag(points_);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcMagFaceAreas() : "
               "finished calculating magFaceAreas in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcFaceAreas() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceAreas() : "
               "calculating faceAreas in PrimitivePatch"
            << endl;
    }

    if (faceAreasPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "faceAreasPtr_ already allocated"
            << abort(FatalError);
    }

    faceAreasPtr_.reset(new Field<point_type>(this->size()));
    auto& n = *faceAreasPtr_;

    forAll(n, facei)
    {
        n[facei] = this->operator[](facei).areaNormal(points_);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceAreas() : "
               "finished calculating faceAreas in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcFaceNormals() const
{
    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceNormals() : "
               "calculating faceNormals in PrimitivePatch"
            << endl;
    }

    if (faceNormalsPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "faceNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    faceNormalsPtr_.reset(new Field<point_type>(this->size()));
    auto& n = *faceNormalsPtr_;

    forAll(n, facei)
    {
        n[facei] = this->operator[](facei).unitNormal(points_);
    }

    if (debug)
    {
        Pout<< "PrimitivePatch<FaceList, PointField>::"
               "calcFaceNormals() : "
               "finished calculating faceNormals in PrimitivePatch"
            << endl;
    }
}


// ************************************************************************* //
