#include "checkGeometry.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

Foam::label Foam::checkGeometry
(
    const polyMesh& mesh,
    bool checkPointNearness,
    bool checkCellDeterminant
)
{
    label noFailedChecks = 0;

    Info<< "\nChecking geometry..." << endl;

    boundBox bb(mesh.points());

    Pout<< "    Domain bounding box: "
           << bb.min() << " " << bb.max() << endl;

    // Get a small relative length from the bounding box
    const boundBox& globalBb = mesh.globalData().bb();

    if (Pstream::parRun())
    {
        Info<< "    Overall domain bounding box: "
            << globalBb.min() << " " << globalBb.max() << endl;
    }


    // Min length
    scalar minDistSqr = magSqr(1e-6*(globalBb.max() - globalBb.min()));

    
    if (mesh.checkClosedBoundary(true)) noFailedChecks++;

    {
        cellSet cells(mesh, "nonClosedCells", mesh.nCells()/100+1);
        cellSet aspectCells(mesh, "highAspectRatioCells", mesh.nCells()/100+1);
        if (mesh.checkClosedCells(true, &cells, &aspectCells))
        {
            noFailedChecks++;

            if (cells.size() > 0)
            {
                Pout<< "  <<Writing " << cells.size()
                    << " non closed cells to set " << cells.name() << endl;
                cells.write();
            }
        }
        if (aspectCells.size() > 0)
        {
            Pout<< "  <<Writing " << aspectCells.size()
                << " cells with high aspect ratio to set "
                << aspectCells.name() << endl;
            aspectCells.write();
        }
    }

    {
        faceSet faces(mesh, "zeroAreaFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAreas(true, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Pout<< "  <<Writing " << faces.size()
                    << " zero area faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    {
        cellSet cells(mesh, "zeroVolumeCells", mesh.nCells()/100 + 1);
        if (mesh.checkCellVolumes(true, &cells))
        {
            noFailedChecks++;

            if (cells.size() > 0)
            {
                Pout<< "  <<Writing " << cells.size()
                    << " zero volume cells to set " << cells.name() << endl;
                cells.write();
            }
        }
    }

    {
        faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceOrthogonality(true, &faces))
        {
            noFailedChecks++;
        }

        if (faces.size() > 0)
        {
            Pout<< "  <<Writing " << faces.size()
                << " non-orthogonal faces to set " << faces.name() << endl;
            faces.write();
        }
    }


    {
        faceSet faces(mesh, "wrongOrientedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFacePyramids(true, -SMALL, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Pout<< "  <<Writing " << faces.size()
                    << " faces with incorrect orientation to set "
                    << faces.name() << endl;
                faces.write();
            }
        }
    }

    {
        faceSet faces(mesh, "skewFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceSkewness(true, &faces))
        {
            noFailedChecks++;

            if (faces.size() > 0)
            {
                Pout<< "  <<Writing " << faces.size()
                    << " skew faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    if (checkPointNearness)
    {
        // Note use of nPoints since don't want edge construction.
        pointSet points(mesh, "shortEdges", mesh.nPoints()/1000 + 1);
        if (mesh.checkEdgeLength(true, minDistSqr, &points))
        {
            //noFailedChecks++;

            if (points.size() > 0)
            {
                Pout<< "  <<Writing " << points.size()
                    << " points on short edges to set " << points.name()
                    << endl;
                points.write();
            }
        }

        label nEdgeClose = points.size();

        if (mesh.checkPointNearness(false, minDistSqr, &points))
        {
            //noFailedChecks++;

            if (points.size() > nEdgeClose)
            {
                pointSet nearPoints(mesh, "nearPoints", points);
                Pout<< "  <<Writing " << nearPoints.size()
                    << " near (closer than " << Foam::sqrt(minDistSqr)
                    << " apart) points to set " << nearPoints.name() << endl;
                nearPoints.write();
            }
        }
    }

    {
        faceSet faces(mesh, "concaveFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAngles(true, 10, &faces))
        {
            //noFailedChecks++;

            if (faces.size() > 0)
            {
                Pout<< "  <<Writing " << faces.size()
                    << " faces with concave angles to set " << faces.name()
                    << endl;
                faces.write();
            }
        }
    }

    {
        faceSet faces(mesh, "warpedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceFlatness(true, 0.8, &faces))
        {
            //noFailedChecks++;

            if (faces.size() > 0)
            {
                Pout<< "  <<Writing " << faces.size()
                    << " warped faces to set " << faces.name() << endl;
                faces.write();
            }
        }
    }

    if (checkCellDeterminant)
    {
        cellSet cells(mesh, "underdeterminedCells", mesh.nCells()/100);
        if (mesh.checkCellDeterminant(true, &cells))
        {
            noFailedChecks++;

            Pout<< "  <<Writing " << cells.size()
                << " under-determines cells to set " << cells.name() << endl;
            cells.write();
        }
    }
    

    return noFailedChecks;
}
