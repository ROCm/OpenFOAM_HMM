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

#include "CV3D.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV3D::writePoints(const fileName& fName, bool internalOnly) const
{
    Info << "Writing points to " << fName << nl << endl;
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
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


void Foam::CV3D::writeDual(const fileName& fName) const
{
    Info << "Writing dual points and faces to " << fName << nl << endl;
    OFstream str(fName);

    label dualVerti = 0;

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVerti;
            meshTools::writeOBJ(str, topoint(dual(cit)));
            dualVerti++;
        }
        else
        {
            cit->cellIndex() = -1;
        }
    }

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        if
        (
            eit->first->vertex(eit->second)->internalOrBoundaryPoint()
         || eit->first->vertex(eit->third)->internalOrBoundaryPoint()
        )
        {
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc = ccStart;

            str<< 'f';

            do
            {
                if (!is_infinite(cc))
                {
                    if (cc->cellIndex() < 0)
                    {
                        FatalErrorIn
                        (
                            "Foam::CV3D::writeDual(const fileName& fName)"
                        )<< "Dual face uses circumcenter defined by a Delaunay"
                            " tetrahedron with no internal or boundary points."
                            << exit(FatalError);
                    }
        
                    str<< ' ' << cc->cellIndex() + 1;
                }
            } while (++cc != ccStart);

            str<< nl;
        }
    }
}


void Foam::CV3D::writeTriangles(const fileName& fName, bool internalOnly) const
{
    Info << "Writing triangles to " << fName << nl << endl;
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
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle& c(fit->first);
        
        const int& oppositeVertex(fit->second);

        List<label> facetIndices(3,-1);

        bool writeFacet = true;

        for(label i = 0, k = 0;i < 4; i++)
        {            
            if(i != oppositeVertex)
            {
                if(!internalOnly || !c->vertex(i)->farPoint())
                {
                    facetIndices[k] = i;
                    k++;
                }
                else
                {
                    writeFacet = false;
                }
            }
        }
    
        if(writeFacet)
        {
                str << "f "
                       << vertexMap[c->vertex(facetIndices[0])->index()] + 1
                << ' ' << vertexMap[c->vertex(facetIndices[1])->index()] + 1
                << ' ' << vertexMap[c->vertex(facetIndices[2])->index()] + 1
                << nl;
        }
    }
}

void Foam::CV3D::writeMesh(const Time& runTime)
{
    Info<< nl << "Writing polyMesh." << endl;

    Info << nl << "Temporary hack to produce boundary" << endl;

    // for
    // (
    //     Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    //     vit != finite_vertices_end();
    //     ++vit
    // )
    // {
    //     if
    //     (
    //         mag(topoint(vit->point()).x()) > 1.5
    //      || mag(topoint(vit->point()).y()) > 1.5
    //      || mag(topoint(vit->point()).z()) > 1.5   
    //     )
    //     {
    // 	    vit->type() = Vb::FAR_POINT;
    //     }
    // }

    markNearBoundaryPoints();

    pointField points(0);
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
    wordList patchNames(0);
    labelList patchSizes(0);
    labelList patchStarts(0);

    calcDualMesh
    (
        points,
        faces,
        owner,
        neighbour,
        patchNames,
        patchSizes,
        patchStarts
    );

    IOobject io
    (
        Foam::polyMesh::defaultRegion,
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    polyMesh pMesh
    (
        io,
        points,
        faces,
        owner,
        neighbour
    );

    List<polyPatch*> patches(patchStarts.size());

    forAll (patches, p)
    {
        patches[p] = new polyPatch
        (
            patchNames[p],
            patchSizes[p],
            patchStarts[p],
            p,
            pMesh.boundaryMesh()
        );
    }

    pMesh.addPatches(patches);

    if (!pMesh.write())
    {
        FatalErrorIn("CV3D::writeMesh(const Time& runTime)")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// ************************************************************************* //
