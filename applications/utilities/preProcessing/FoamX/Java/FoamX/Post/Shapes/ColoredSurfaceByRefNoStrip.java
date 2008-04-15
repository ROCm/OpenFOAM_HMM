/*
 */
package FoamX.Post.Shapes;

import com.sun.j3d.utils.geometry.GeometryInfo;
import com.sun.j3d.utils.geometry.NormalGenerator;
import com.sun.j3d.utils.geometry.Stripifier;
import com.sun.j3d.utils.geometry.Triangulator;
import com.sun.j3d.utils.picking.PickIntersection;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickTool;

import javax.media.j3d.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import org.j3d.util.interpolator.ColorInterpolator;

import FoamX.Post.Util.OrbitPickBehavior;
import FoamX.Post.InteractiveNodes.ColoredSurfaceByRefNode;
import FoamX.Post.Util.BinaryFoamFile;

public class ColoredSurfaceByRef extends Shape3D implements PickableShape
{
    double creaseAngle_;

    public ColoredSurfaceByRef
    (
        ColoredSurfaceByRefNode node,
        double creaseAngle,
        String pointsFileName,
        String facesFileName,
        String valuesFileName,
        float min,
        float max
    )
    {
        // Store reference to input data
        creaseAngle_ = creaseAngle;

        System.out.println("Loading files ...");
        float[] coords = BinaryFoamFile.readFloatPointArray(pointsFileName);
        int[] triFaces = BinaryFoamFile.readIntArray(facesFileName);
        float[] values = BinaryFoamFile.readFloatArray(valuesFileName);


        // Auto colour
        min = values[0];
        max = values[0];
        for(int i = 0; i < values.length; i++)
        {
            min = Math.min(min, values[i]);
            max = Math.max(max, values[i]);
        }
        System.out.println("min=" + min + "  max=" + max);

//        for(int i = 0; i < coords.length; i++)
//        {
//            System.out.println("i=" + i + " Coord=" + coords[i]);
//        }
//        for(int i = 0; i < triFaces.length; i++)
//        {
//            System.out.println("i=" + i + " TriFaces=" + triFaces[i]);
//        }
//        for(int i = 0; i < values.length; i++)
//        {
//            System.out.println("i=" + i + " Values=" + values[i]);
//        }

        System.out.println("Calculating smooth normals");
        // Storage for vertex normals
        float[] normals = new float[coords.length];
        for(int i = 0; i < normals.length; i++)
        {
            normals[i] = 0;
        }

        System.out.println("Calculating face normals");
        int triIndex = 0;
        for(int faceI = 0; faceI < triFaces.length; faceI += 3)
        {
            int v0offset = 3*triFaces[faceI];
            int v1offset = 3*triFaces[faceI + 1];
            int v2offset = 3*triFaces[faceI + 2];

            float edge1x = coords[v1offset] - coords[v0offset];
            float edge1y = coords[v1offset + 1] - coords[v0offset + 1];
            float edge1z = coords[v1offset + 2] - coords[v0offset + 2];

            float edge2x = coords[v2offset] - coords[v0offset];
            float edge2y = coords[v2offset + 1] - coords[v0offset + 1];
            float edge2z = coords[v2offset + 2] - coords[v0offset + 2];

            Vector3f normal = new Vector3f();
            normal.cross
            (
                new Vector3f(edge1x, edge1y, edge1z),
                new Vector3f(edge2x, edge2y, edge2z)
            );

            normal.normalize();

            normals[v0offset]   += normal.x;
            normals[v0offset+1] += normal.y;
            normals[v0offset+2] += normal.z;

            normals[v1offset]   += normal.x;
            normals[v1offset+1] += normal.y;
            normals[v1offset+2] += normal.z;

            normals[v2offset]   += normal.x;
            normals[v2offset+1] += normal.y;
            normals[v2offset+2] += normal.z;
        }

        System.out.println("Normalizing normals");
        // Normalize normals
        for(int i = 0; i < normals.length; i += 3)
        {
            float x = normals[i];
            float y = normals[i+1];
            float z = normals[i+2];

            float mag = (float)Math.sqrt(x*x + y*y + z*z);
            normals[i] /= mag;
            normals[i+1] /= mag;
            normals[i+2] /= mag;
        }
            

        System.out.println("Converting values to colours ...");
        float[] colors = new float[3*values.length];
        valuesToVertexColors(values, colors, min, max);

        System.out.println("Unindexing geometry ...");
        float[] triCoords = new float[3*triFaces.length];
        int vertI = 0;
        for(int i = 0; i < triFaces.length; i++)
        {
            triCoords[vertI++] = coords[3*triFaces[i]];       // X
            triCoords[vertI++] = coords[3*triFaces[i] + 1];   // Y
            triCoords[vertI++] = coords[3*triFaces[i] + 2];   // Z
        }

        System.out.println("Unindexing normals ...");
        float[] triNormals = new float[3*triFaces.length];
        vertI = 0;
        for(int i = 0; i < triFaces.length; i++)
        {
            triNormals[vertI++] = normals[3*triFaces[i]];       // X
            triNormals[vertI++] = normals[3*triFaces[i] + 1];   // Y
            triNormals[vertI++] = normals[3*triFaces[i] + 2];   // Z
        }

        System.out.println("Unindexing colours ...");
        float[] triColors = new float[3*triFaces.length];
        vertI = 0;
        for(int i = 0; i < triFaces.length; i++)
        {
            triColors[vertI++] = colors[3*triFaces[i]];       // R
            triColors[vertI++] = colors[3*triFaces[i] + 1];   // G
            triColors[vertI++] = colors[3*triFaces[i] + 2];   // B
        }

        TriangleArray array = new TriangleArray
        (
            triFaces.length,
            TriangleArray.COORDINATES
          | TriangleArray.COLOR_3
          | TriangleArray.NORMALS
          | TriangleArray.BY_REFERENCE
        );
        array.setCoordRefFloat(triCoords);
        array.setColorRefFloat(triColors);
        array.setNormalRefFloat(triNormals);

//        TriangleArray array = new TriangleArray
//        (
//            triFaces.length,
//            TriangleArray.COORDINATES
//          | TriangleArray.COLOR_3
//        );
//        array.setCoordinates(0, triCoords);
//        array.setColors(0, triColors);

        this.setGeometry(array);
        //PickTool.setCapabilities(this, PickTool.INTERSECT_FULL);
        //this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        //this.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        //array.setCapability(IndexedGeometryArray.ALLOW_COLOR_INDEX_WRITE);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
        System.out.println("Done making shape...");
    }

    // Convert values into colors
    static void valuesToVertexColors
    (
        float[] values, float[] colors, float min, float max
    )
    {
        if (colors.length != 3*values.length)
        {
            System.out.println
            (
                "Lengths different:" + values.length + " " + colors.length
            );
        }
        if (values.length == 0)
        {
            return;
        }

        ColorInterpolator ci = new ColorInterpolator(3);
        ci.addRGBKeyFrame(min, 0.0f, 0.0f, 1.0f, 1.0f);
        float mid = (float)(0.5*min+0.5*max);
        ci.addRGBKeyFrame(mid, 0.0f, 1.0f, 0.0f, 1.0f);
        ci.addRGBKeyFrame(max, 1.0f, 0.0f, 0.0f, 1.0f);

        int colorComponent = 0;
        for(int i = 0; i < values.length; i++)
        {
            float[] rgbColor = ci.floatRGBValue(values[i]);
            // Copy RGB
            colors[colorComponent++] = rgbColor[0];
            colors[colorComponent++] = rgbColor[1];
            colors[colorComponent++] = rgbColor[2];
        }
    }


    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on ColoredSurfaceByRef ");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUnpicked on ColoredSurfaceByRef ");
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        System.out.println("remove of ColoredSurfaceByRef");

        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }

    //--------------------------------------------------------------------------
}
