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

import FoamX.Post.Mesh.PrimitivePatch;
import FoamX.Post.InteractiveNodes.PatchFacesNode;
import FoamX.Post.Util.OrbitPickBehavior;

public class PatchFaces extends Shape3D implements PickableShape
{
    static final int UNPICKEDCOLOR = 0;
    static final int PICKEDCOLOR = 1;


    PrimitivePatch patch_;
    Color3f color_;
    double creaseAngle_;

    public PatchFaces
    (
        PatchFacesNode PatchFacesNode,
        PrimitivePatch patch,
        Color3f color,
        double creaseAngle
    )
    {
        // Store reference to input data
        patch_ = patch;
        color_ = color;
        creaseAngle_ = creaseAngle;


        Point3f[] coords = patch_.getPoints();
        int[] triFaces = patch_.getLocalTriFaces();

        GeometryInfo geometryInfo =
            new GeometryInfo(GeometryInfo.TRIANGLE_ARRAY);
        geometryInfo.setCoordinateIndices(triFaces);
        geometryInfo.setCoordinates(coords);


        Color3f[] colors = new Color3f[2];
        colors[UNPICKEDCOLOR] = color_;
        colors[PICKEDCOLOR] = new Color3f(0.5f, 1.0f, 0.2f);
        geometryInfo.setColors(colors);

        //- Start off with all triangles unpickedcolor
        int[] colorIndices = new int[triFaces.length];
        for(int i = 0; i < triFaces.length; i++)
        {
            colorIndices[i] = UNPICKEDCOLOR;
        }
        geometryInfo.setColorIndices(colorIndices);


        //- Generate normals
        NormalGenerator ng = new NormalGenerator(creaseAngle_);
        ng.generateNormals(geometryInfo);

        //- Make shape and set various
        IndexedGeometryArray array = geometryInfo.getIndexedGeometryArray();

        this.setGeometry(array);
        PickTool.setCapabilities(this, PickTool.INTERSECT_FULL);
        //this.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
        this.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
        array.setCapability(IndexedGeometryArray.ALLOW_COLOR_INDEX_WRITE);
    }

    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on patch " + patch_.getPatchNo());

        IndexedGeometryArray array = (IndexedGeometryArray)this.getGeometry();

        int nearest = OrbitPickBehavior.getNearestIndex(pickResult);
        PickIntersection pi = pickResult.getIntersection(nearest);

        int[] indices = pi.getPrimitiveVertexIndices();

        int triIndex = indices[0] / 3;  // 3 labels per triangle
        int faceIndex = patch_.getTriToPoly()[triIndex];

        System.out.print
        (
            "-- triangle " + triIndex + " on face " + faceIndex
          + " consisting of local points "
        );

        for(int j = 0; j < indices.length; j++)
        {
            System.out.print(indices[j] + " ");
        }
        System.out.println();

        // Color all triangles of face
        //Hack: loop over all triangles instead just of triangles for face
        int[] triFaces = patch_.getLocalTriFaces();
        for(int i = 0; i < triFaces.length; i++)
        {
            triIndex = i / 3;
            if (patch_.getTriToPoly()[triIndex] == faceIndex)
            {
                array.setColorIndex(i, PICKEDCOLOR);
            }
        }
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUnpicked on patch " + patch_.getPatchNo());

        IndexedGeometryArray array = (IndexedGeometryArray)this.getGeometry();

        int nearest = OrbitPickBehavior.getNearestIndex(pickResult);
        PickIntersection pi = pickResult.getIntersection(nearest);

        int[] indices = pi.getPrimitiveVertexIndices();

        int triIndex = indices[0] / 3;  // 3 labels per triangle
        int faceIndex = patch_.getTriToPoly()[triIndex];

        System.out.print
        (
            "-- triangle " + triIndex + " on face " + faceIndex
          + " consisting of local points "
        );

        for(int j = 0; j < indices.length; j++)
        {
            System.out.print(indices[j] + " ");
        }
        System.out.println();

        // Color all triangles of face
        //Hack: loop over all triangles instead just of triangles for face
        int[] triFaces = patch_.getLocalTriFaces();
        for(int i = 0; i < triFaces.length; i++)
        {
            triIndex = i / 3;
            if (patch_.getTriToPoly()[triIndex] == faceIndex)
            {
                array.setColorIndex(i, UNPICKEDCOLOR);
            }
        }
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        //System.out.println("remove:" + patch_.getPatchNo());
        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }

    //--------------------------------------------------------------------------
}
