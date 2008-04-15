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

import FoamX.Post.InteractiveNodes.ColoredSurfaceVecNode;
import FoamX.Post.Util.BinaryFoamFile;
import FoamX.Post.Util.Colorizer;

public class ColoredSurfaceVec extends Shape3D implements PickableShape
{
    // reads triangle list.
    public ColoredSurfaceVec
    (
        final ColoredSurfaceVecNode node,
        final double creaseAngle,
        final boolean stripify,
        final String pointsFileName,
        final String facesFileName,
        final String valuesFileName,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        this
        (
            node,
            false,
            false,
            stripify,
            creaseAngle,
            pointsFileName,
            null,
            facesFileName,
            valuesFileName,
            colorizer,
            min,
            max
        );
    }
            

    // reads triangle strip.
    public ColoredSurfaceVec
    (
        final ColoredSurfaceVecNode node,
        final boolean readTriStrip,
        final boolean readNormals,
        final boolean stripify,
        final double creaseAngle,
        final String pointsFileName,
        final String normalsFileName,
        final String trianglesFileName,
        final String valuesFileName,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        int[][] triStrips = null;
        int[] triFaces = null;
        float[] normals = null;

        System.out.println("Loading files ...");
        float[] coords = BinaryFoamFile.readFloatPointArray(pointsFileName);

        if (readNormals)
        {
            normals = BinaryFoamFile.readFloatPointArray(normalsFileName);
        }

        if (readTriStrip)
        {
            triStrips = BinaryFoamFile.readIntIntArray(trianglesFileName);
        }
        else
        {
            triFaces = BinaryFoamFile.readIntArray(trianglesFileName);
        }

        float[] values = BinaryFoamFile.readFloatArray(valuesFileName);

        Geometry array = createGeometry
            (
                node,
                readTriStrip,
                stripify,
                creaseAngle,
                coords,
                normals,
                triFaces,
                triStrips,
                values,
                colorizer,
                min,
                max
            );
        this.setGeometry(array);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }


    // Creates coloured geometry from values
    public ColoredSurfaceVec
    (
        final ColoredSurfaceVecNode node,
        final boolean stripInput,
        final boolean stripify,
        final double creaseAngle,
        final float[] points,
        final float[] normals,
        final int[] triFaces,
        final int[][] triStrips,
        final float[] values,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        Geometry array = createGeometry
            (
                node,
                stripInput,
                stripify,
                creaseAngle,
                points,
                normals,
                triFaces,
                triStrips,
                values,
                colorizer,
                min,
                max
            );
        this.setGeometry(array);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }


    // Creates coloured geometry from values
    private Geometry createGeometry
    (
        final ColoredSurfaceVecNode node,
        final boolean stripInput,
        final boolean stripify,
        final double creaseAngle,
        final float[] coords,
        final float[] normals,
        final int[] triFaces,
        final int[][] triStrips,
        final float[] values,
        final Colorizer colorizer,
        float min,
        float max
    )
    {

        if ((min == 0.0f) && (max == 0.0f))
        {
            // Auto colour
            if (values.length > 0)
            {
                min = values[0];
                max = values[0];
                for(int i = 0; i < values.length; i++)
                {
                    min = Math.min(min, values[i]);
                    max = Math.max(max, values[i]);
                }
            }
        }
        System.out.println("min=" + min + "  max=" + max);


        float[] localNormals = normals;
        if (normals == null)
        {
            System.out.println("Calculating vertex normals");
            localNormals = new float[coords.length];
            if (stripInput)
            {
                calcStripNormals(triStrips, coords, localNormals);
            }
            else
            {
                calcNormals(triFaces, coords, localNormals);
            }
        }


        System.out.println("Converting values to colours ...");
        float[] colors = new float[3*values.length];
        colorizer.color(min, max, values, colors);

        float[] triCoords;
        float[] triNormals;
        float[] triColors;

        GeometryArray array;

        if (stripInput)
        {
            int[] triLengths = new int[triStrips.length];
            int totalNvertices = stripLength(triStrips, triLengths);

            triCoords = new float[3*totalNvertices];
            triNormals = new float[3*totalNvertices];
            triColors = new float[3*totalNvertices];
    
            System.out.println("Unindexing ...");


            if (colors.length != coords.length)
            {
                System.out.println("***Cannot unindex face based data");
            }

            unindexStrip
            (
                triStrips, coords, localNormals, colors,
                triCoords, triNormals, triColors
            );

            System.out.println("Creating TriangleStripArray ...");

            array = new TriangleStripArray
            (
                totalNvertices,
                TriangleArray.COORDINATES
              | TriangleArray.COLOR_3
              | TriangleArray.NORMALS
              | TriangleArray.BY_REFERENCE,
                triLengths
            );
        }
        else
        {
            if (stripify)
            {
                // Use geometryInfo to triangulate

                System.out.println("Stripifying ...");

                GeometryInfo geometryInfo =
                    new GeometryInfo(GeometryInfo.TRIANGLE_ARRAY);
                geometryInfo.setCoordinateIndices(triFaces);
                geometryInfo.setCoordinates(coords);
                Stripifier st = new Stripifier();
                st.stripify(geometryInfo);

                int[] triLengths = geometryInfo.getStripCounts();
                int[] triIndices = geometryInfo.getCoordinateIndices();

                int totalNvertices = 0;
                for(int i = 0; i < triLengths.length; i++)
                {
                    totalNvertices += triLengths[i];
                }

                triCoords = new float[3*totalNvertices];
                triNormals = new float[3*totalNvertices];
                triColors = new float[3*totalNvertices];

                System.out.println("Unindexing ...");

                unindex
                (
                    triIndices, coords, localNormals, colors,
                    triCoords, triNormals, triColors
                );

                System.out.println("Creating TriangleStripArray ...");

                array = new TriangleStripArray
                (
                    totalNvertices,
                    TriangleArray.COORDINATES
                  | TriangleArray.COLOR_3
                  | TriangleArray.NORMALS
                  | TriangleArray.BY_REFERENCE,
                    triLengths
                );
            }
            else
            {
                triCoords = new float[3*triFaces.length];
                triNormals = new float[3*triFaces.length];
                triColors = new float[3*triFaces.length];

                System.out.println("Unindexing ...");

                unindex
                (
                    triFaces, coords, localNormals, colors,
                    triCoords, triNormals, triColors
                );

                array = new TriangleArray
                (
                    triFaces.length,
                    TriangleArray.COORDINATES
                  | TriangleArray.COLOR_3
                  | TriangleArray.NORMALS
                  | TriangleArray.BY_REFERENCE
                );
            }
        }
        array.setCoordRefFloat(triCoords);
        array.setColorRefFloat(triColors);
        array.setNormalRefFloat(triNormals);

        System.out.println("Done making shape...");
        return array;
    }

    //--------------------------------------------------------------------------

    // Calculate normals
    static void calcStripNormals
    (
        final int[][] triStrips,
        final float[] coords,
        float[] normals
    )
    {
        for(int i = 0; i < normals.length; i++)
        {
            normals[i] = 0;
        }

        System.out.println("Calculating face normals");
        Vector3f firstNormal = null;
        for(int stripI = 0; stripI < triStrips.length; stripI++)
        {
            int[] strip = triStrips[stripI];

            int v0offset = 3*strip[0];
            int v1offset = 3*strip[1];

            float sign = 1;


            for(int vertI = 2; vertI < strip.length; vertI++)
            {
                int v2offset = 3*strip[vertI];

                // Get edge from vertex1 to vertex0
                float edge1x = coords[v1offset] - coords[v0offset];
                float edge1y = coords[v1offset + 1] - coords[v0offset + 1];
                float edge1z = coords[v1offset + 2] - coords[v0offset + 2];

                // Get edge from vertex2 to vertex0
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

                if (vertI == 2)
                {
                    // First triangle of current strip
                    if (stripI == 0)
                    {
                        // If first strip store normal
                        firstNormal = normal;
                    }
                    else
                    {
                        // If not first strip compare and flip
                        // so orientation is equal
                        if (firstNormal.dot(normal) < 0)
                        {
                            sign = -1.0f;
                        }
                    }
                }

                normal.x *= sign;
                normal.y *= sign;
                normal.z *= sign;

                normals[v0offset]   += normal.x;
                normals[v0offset+1] += normal.y;
                normals[v0offset+2] += normal.z;

                normals[v1offset]   += normal.x;
                normals[v1offset+1] += normal.y;
                normals[v1offset+2] += normal.z;

                normals[v2offset]   += normal.x;
                normals[v2offset+1] += normal.y;
                normals[v2offset+2] += normal.z;

                // Go to next triangle
                v0offset = v1offset;
                v1offset = v2offset;
                sign = -sign;
            }
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
    }


    //--------------------------------------------------------------------------

    // Calculate normals
    static void calcNormals
    (
        final int[] triFaces,
        final float[] coords,
        float[] normals
    )
    {
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
    }

    //--------------------------------------------------------------------------

    // Determine number of vertices in strips
    static int stripLength(final int[][] triStrips, int[] triLengths)
    {
        int totalNvertices = 0;
        for(int stripI = 0; stripI < triStrips.length; stripI++)
        {
            int len = triStrips[stripI].length;
            triLengths[stripI] = len;
            totalNvertices += len;
        }
        return totalNvertices;
    }

    //--------------------------------------------------------------------------


    // Unindex
    static void unindexStrip
    (
        final int[][] triStrips,
        final float[] coords,
        final float[] normals,
        final float[] colors,

        float[] triCoords,
        float[] triNormals,
        float[] triColors
    )
    {
        int coordI = 0;
        for(int stripI = 0; stripI < triStrips.length; stripI++)
        {
            int[] strip = triStrips[stripI];
            int len = strip.length;
            for(int vertI = 0; vertI < len; vertI++)
            {
                triCoords[coordI] = coords[3*strip[vertI]];
                triNormals[coordI] = normals[3*strip[vertI]];
                triColors[coordI] = colors[3*strip[vertI]];
                coordI++;

                triCoords[coordI] = coords[3*strip[vertI]+1];
                triNormals[coordI] = normals[3*strip[vertI]+1];
                triColors[coordI] = colors[3*strip[vertI]+1];
                coordI++;

                triCoords[coordI] = coords[3*strip[vertI]+2];
                triNormals[coordI] = normals[3*strip[vertI]+2];
                triColors[coordI] = colors[3*strip[vertI]+2];
                coordI++;
            }
        }
    }

    //--------------------------------------------------------------------------

    // Unindex strip in Java3D format
    static void unindexJava3DStrip
    (
        final int[] triLengths,
        final int[] triIndices,
        final float[] coords,
        final float[] normals,
        final float[] colors,

        float[] triCoords,
        float[] triNormals,
        float[] triColors
    )
    {
        int coordI = 0;
        for(int vertI = 0; vertI < triIndices.length; vertI++)
        {
            int indexx = 3*triIndices[vertI];
            int indexy = indexx + 1;
            int indexz = indexx + 2;

            triCoords[coordI] = coords[indexx];
            triNormals[coordI] = normals[indexx];
            triColors[coordI] = colors[indexx];
            coordI++;

            triCoords[coordI] = coords[indexy];
            triNormals[coordI] = normals[indexy];
            triColors[coordI] = colors[indexy];
            coordI++;

            triCoords[coordI] = coords[indexz];
            triNormals[coordI] = normals[indexz];
            triColors[coordI] = colors[indexz];
            coordI++;
        }
    }

    //--------------------------------------------------------------------------

    // Unindex
    static void unindex
    (
        final int[] triFaces,
        final float[] coords,
        final float[] normals,
        final float[] colors,

        float[] triCoords,
        float[] triNormals,
        float[] triColors
    )
    {
        int vertI = 0;

System.out.println
(
    "unindex: colors:" + colors.length + "  triFaces:" + triFaces.length
);

        if (colors.length == triFaces.length)
        {
            // Element wise colors
            for(int i = 0; i < triFaces.length; i++)
            {
                int indexx = 3*triFaces[i];
                int indexy = indexx + 1;
                int indexz = indexx + 2;

                int faceI = i / 3;

                triCoords[vertI] = coords[indexx];      // X
                triNormals[vertI] = normals[indexx];    // X
                triColors[vertI] = colors[3*faceI];     // R
                vertI++;

                triCoords[vertI] = coords[indexy];      // Y
                triNormals[vertI] = normals[indexy];    // Y
                triColors[vertI] = colors[3*faceI+1];   // G
                vertI++;

                triCoords[vertI] = coords[indexz];      // Z
                triNormals[vertI] = normals[indexz];    // Z
                triColors[vertI] = colors[3*faceI+2];   // B
                vertI++;
            }
        }
        else
        {
            // Vertex wise colors
            for(int i = 0; i < triFaces.length; i++)
            {
                int indexx = 3*triFaces[i];
                int indexy = indexx + 1;
                int indexz = indexx + 2;

                triCoords[vertI] = coords[indexx];       // X
                triNormals[vertI] = normals[indexx];     // X
                triColors[vertI] = colors[indexx];       // R
                vertI++;

                triCoords[vertI] = coords[indexy];       // Y
                triNormals[vertI] = normals[indexy];     // Y
                triColors[vertI] = colors[indexy];       // G
                vertI++;

                triCoords[vertI] = coords[indexz];       // Z
                triNormals[vertI] = normals[indexz];     // Z
                triColors[vertI] = colors[indexz];       // B
                vertI++;
            }
        }
    }


    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on ColoredSurfaceVec ");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUnpicked on ColoredSurfaceVec ");
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        System.out.println("remove of ColoredSurfaceVec");

        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }

    //--------------------------------------------------------------------------
}
