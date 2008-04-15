/*
 */

package FoamX.Post.Mesh;

import javax.media.j3d.*;
import javax.vecmath.Point3f;

public class PrimitivePatch
{

    public static final Point3f[] points = 
    {
        new Point3f(0.0f, 0.0f, 0.0f),
        new Point3f(1.0f, 0.0f, 0.0f),
        new Point3f(1.0f, 1.0f, 0.0f),
        new Point3f(0.5f, 2.0f, 0.0f),
        new Point3f(0.0f, 1.0f, 0.0f),

        new Point3f(0.0f, 0.0f, 1.0f),
        new Point3f(1.0f, 0.0f, 1.0f),
        new Point3f(1.0f, 1.0f, 1.0f),
        new Point3f(0.5f, 2.0f, 1.0f),
        new Point3f(0.0f, 1.0f, 1.0f),

        new Point3f(0.5f, 0.8f, 0.0f),   // 10
        new Point3f(0.5f, 0.8f, 1.0f)    // 11

    };

    // Array of localFaces(but in 'sufaceMesh' pointLabels),
    // one for each patch
    static final int[][][] localFacesList = 
    {

        // Patch 0
        {
            { 0, 1, 2,10},      // Face 0
            {10, 2, 3, 4},      // Face 1
            { 0,10, 4},         // Face 2

            { 5, 6, 7,11},      // Face 3
            {11, 7, 8, 9},      // Face 4
            { 5,11, 9}          // Face 5
        },
        // Patch 1
        {
            {0, 5, 9, 6},       // Face 0
            {4, 9, 8, 3},       // Face 1
            {3, 8, 7, 2},       // Face 2
            {2, 7, 6, 1},       // Face 3
            {1, 6, 5, 0}        // Face 4
        }
    };

    static final int[][] localTriFacesList = 
    {

        // Patch 0
        {
             0,10, 1,
             1,10, 2,
             2,10, 3,
             3,10, 4,
             4,10, 0,

             5, 6,11,
            11, 6, 7,
            11, 7, 8,
            11, 8, 9,
            11, 9, 5
        },
        // Patch 1
        {
            0, 5, 9,
            0, 9, 4,
            4, 9, 8,
            4, 8, 3,
            3, 8, 7,
            3, 7, 2,
            2, 7, 6,
            2, 6, 1,
            1, 6, 5,
            1, 5, 0
        }
    };

    // For each triangle give index in localFacesList
    static final int[][] triToPolyList = 
    {

        // Patch 0
        {
            0,
            0,
            1,
            1,
            2,

            3,
            3,
            4,
            4,
            5
        },
        // Patch 1
        {
            0,
            0,
            1,
            1,
            2,
            2,
            3,
            3,
            4,
            4,
            5,
            5
        }
    };

    static final int[][] featureLinesList =
    {

        // Patch 0
        {
            0, 1, 1, 2, 2, 3, 3, 4, 4, 0,
            5, 6, 6, 7, 7, 8, 8, 9, 9, 5,
            3, 10                               // Loose "crease"
        },
        // Patch 1
        {
            0, 1, 1, 2, 2, 3, 3, 4, 4, 0,
            5, 6, 6, 7, 7, 8, 8, 9, 9, 5,
            0, 5, 4, 9, 3, 8, 2, 7, 1, 6
        }
    };


    int patchNo_;

    public PrimitivePatch(int patchNo)
    {
        patchNo_ = patchNo;
    }


    public int getPatchNo()
    {
        return patchNo_;
    }

    //- Points
    public int getNumPoints()
    {
        return points.length;
    }

    //- Points
    public Point3f[] getPoints()
    {
        return points;
    }

    //- List of faces
    public int[][] getLocalFaces()
    {
        return localFacesList[patchNo_];
    }

    //- List of triangulated faces
    public int[] getLocalTriFaces()
    {
        return localTriFacesList[patchNo_];
    }

    //- List of triToPoly
    public int[] getTriToPoly()
    {
        return triToPolyList[patchNo_];
    }

    //- List of labels. Every two are edge. Unsorted/unconnected.
    public int[] getFeatureLines()
    {
        return featureLinesList[patchNo_];
    }
}
