/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007-2009 OpenCFD Ltd.
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

#include "CV2D.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV2D::writePoints(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing points to " << fName << nl << endl;
    OFstream str(fName);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || vit->internalOrBoundaryPoint())
        {
            meshTools::writeOBJ(str, toPoint3D(vit->point()));
        }
    }
}


void Foam::CV2D::writeTriangles(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing triangles to " << fName << nl << endl;
    OFstream str(fName);

    labelList vertexMap(number_of_vertices());
    label verti = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || !vit->farPoint())
        {
            vertexMap[vit->index()] = verti++;
            meshTools::writeOBJ(str, toPoint3D(vit->point()));
        }
    }

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            !internalOnly
         || (
                fit->vertex(0)->internalOrBoundaryPoint()
             || fit->vertex(1)->internalOrBoundaryPoint()
             || fit->vertex(2)->internalOrBoundaryPoint()
            )
        )
        {
            str << "f " << vertexMap[fit->vertex(0)->index()] + 1
                << ' ' << vertexMap[fit->vertex(1)->index()] + 1
                << ' ' << vertexMap[fit->vertex(2)->index()] + 1
                << nl;
        }
    }
}


void Foam::CV2D::writeFaces(const fileName& fName, bool internalOnly) const
{
    Info<< "Writing dual faces to " << fName << nl << endl;
    OFstream str(fName);

    label dualVerti = 0;

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            !internalOnly
         || (
                fit->vertex(0)->internalOrBoundaryPoint()
             || fit->vertex(1)->internalOrBoundaryPoint()
             || fit->vertex(2)->internalOrBoundaryPoint()
            )
        )
        {
            fit->faceIndex() = dualVerti++;
            meshTools::writeOBJ(str, toPoint3D(circumcenter(fit)));
        }
        else
        {
            fit->faceIndex() = -1;
        }
    }

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || vit->internalOrBoundaryPoint())
        {
            Face_circulator fcStart = incident_faces(vit);
            Face_circulator fc = fcStart;

            str<< 'f';

            do
            {
                if (!is_infinite(fc))
                {
                    if (fc->faceIndex() < 0)
                    {
                        FatalErrorIn
                        (
                            "Foam::CV2D::writeFaces"
                            "(const fileName& fName, bool internalOnly)"
                        )<< "Dual face uses vertex defined by a triangle"
                            " defined by an external point"
                            << exit(FatalError);
                    }

                    str<< ' ' << fc->faceIndex() + 1;
                }
            } while (++fc != fcStart);

            str<< nl;
        }
    }
}


void Foam::CV2D::calcDual(point2DField& dualPoints, faceList& dualFaces) const
{
    // Dual points stored in triangle order.
    dualPoints.setSize(number_of_faces());
    label dualVerti = 0;

    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        if
        (
            fit->vertex(0)->internalOrBoundaryPoint()
         || fit->vertex(1)->internalOrBoundaryPoint()
         || fit->vertex(2)->internalOrBoundaryPoint()
        )
        {
            fit->faceIndex() = dualVerti;
            dualPoints[dualVerti++] = toPoint2D(circumcenter(fit));
        }
        else
        {
            fit->faceIndex() = -1;
        }
    }

    dualPoints.setSize(dualVerti);


    // Create dual faces
    // ~~~~~~~~~~~~~~~~~

    dualFaces.setSize(number_of_vertices());
    label dualFacei = 0;
    labelList faceVerts(maxNvert);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            Face_circulator fcStart = incident_faces(vit);
            Face_circulator fc = fcStart;
            label verti = 0;

            do
            {
                if (!is_infinite(fc))
                {
                    if (fc->faceIndex() < 0)
                    {
                        FatalErrorIn
                        (
                            "Foam::CV2D::calcDual"
                            "(point2DField& dualPoints, faceList& dualFaces)"
                        )<< "Dual face uses vertex defined by a triangle"
                            " defined by an external point"
                            << exit(FatalError);
                    }

                    // Look up the index of the triangle
                    faceVerts[verti++] = fc->faceIndex();
                }
            } while (++fc != fcStart);

            if (faceVerts.size() > 2)
            {
                dualFaces[dualFacei++] =
                    face(labelList::subList(faceVerts, verti));
            }
            else
            {
                Info<< "From triangle point:" << vit->index()
                    << " coord:" << toPoint2D(vit->point())
                    << " generated illegal dualFace:" << faceVerts
                    << endl;
            }
        }
    }

    dualFaces.setSize(dualFacei);
}


void Foam::CV2D::writePatch(const fileName& fName) const
{
    point2DField dual2DPoints;
    faceList dualFaces;

    calcDual(dual2DPoints, dualFaces);
    pointField dualPoints(dual2DPoints.size());
    forAll(dualPoints, ip)
    {
        dualPoints[ip] = toPoint3D(dual2DPoints[ip]);
    }

    // Dump as primitive patch to be read by extrudeMesh.
    OFstream str(fName);

    Info<< "Writing patch to be used with extrudeMesh to file " << fName
        << endl;

    str << dualPoints << nl << dualFaces << nl;
}


// ************************************************************************* //
