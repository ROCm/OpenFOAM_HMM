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
import FoamX.Post.InteractiveNodes.PatchPointsNode;
import FoamX.Post.Util.OrbitPickBehavior;

public class PatchPoints extends Shape3D implements PickableShape
{
    static final int UNPICKEDCOLOR = 0;
    static final int PICKEDCOLOR = 1;


    PrimitivePatch patch_;
    Color3f color_;

    public PatchPoints
    (
        PatchPointsNode PatchPointsNode,
        PrimitivePatch patch,
        Color3f color
    )
    {
        // Store reference to input data
        patch_ = patch;
        color_ = color;

        // Color table
        Color3f[] colors = new Color3f[2];
        colors[UNPICKEDCOLOR] = color_;
        colors[PICKEDCOLOR] = new Color3f(0.5f, 1.0f, 0.2f);


        Point3f[] coords = patch_.getPoints();
        int[] coordIndices = new int[coords.length];
        for(int i = 0; i < coordIndices.length; i++)
        {
            coordIndices[i] = i;
        }

        IndexedPointArray array = new IndexedPointArray
        (
            coords.length,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3,
            coords.length
        );

        array.setCoordinates(0, coords);
        array.setCoordinateIndices(0, coordIndices);

        int[] colorIndices = new int[coords.length];
        for(int i = 0; i < colorIndices.length; i++)
        {
            colorIndices[i] = UNPICKEDCOLOR;
        }
        array.setColors(0, colors);
        array.setColorIndices(0, colorIndices);


        this.setGeometry(array);
        PickTool.setCapabilities(this, PickTool.INTERSECT_FULL);

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

        IndexedPointArray array = (IndexedPointArray)this.getGeometry();

        int nearest = OrbitPickBehavior.getNearestIndex(pickResult);
        PickIntersection pi = pickResult.getIntersection(nearest);

        System.out.println("pi=" + pi);
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUnpicked on patch " + patch_.getPatchNo());

        IndexedPointArray array = (IndexedPointArray)this.getGeometry();

        int nearest = OrbitPickBehavior.getNearestIndex(pickResult);
        PickIntersection pi = pickResult.getIntersection(nearest);

        System.out.println("pi=" + pi);
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
