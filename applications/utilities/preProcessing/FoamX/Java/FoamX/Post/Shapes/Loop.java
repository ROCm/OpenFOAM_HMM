/*
 */

package FoamX.Post.Shapes;


import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.InteractiveNodes.LoopNode;
import FoamX.Post.Mesh.PrimitivePatch;

public class Loop extends Shape3D implements PickableShape
{
    static final int UNPICKEDCOLOR = 0;
    static final int PICKEDCOLOR = 1;

    static final Color3f[] colors =
    {
        new Color3f(0.0f, 1.0f, 0.0f),
        new Color3f(1.0f, 0.0f, 1.0f)
    };

    FeatureLines featureLines_; 
    PrimitivePatch patch_;

    //--------------------------------------------------------------------------
    
    public Loop
    (
        LoopNode loopNode,
        PrimitivePatch patch,
        FeatureLines featureLines
    )
    {
        // Store reference to original data
        featureLines_ = featureLines;
        patch_ = patch;

        int[] lines = featureLines_.getSelectedEdges();
        Point3f[] coords = patch_.getPoints();

        //
        // Make lineArray from selected edges
        //

for(int i =0; i < coords.length; i++)
{
    System.out.println("Loop:coord:" + coords[i]);
}
for(int i =0; i < lines.length; i++)
{
    System.out.println("Loop:lines:" + lines[i]);
}


        IndexedLineArray loopEdges = new IndexedLineArray
        (
            coords.length,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3,
            lines.length
        );

        loopEdges.setCoordinates(0, coords);
        loopEdges.setCoordinateIndices(0, lines);

        loopEdges.setColors(0, colors);
        for(int i = 0; i < lines.length; i++)
        {
            System.out.println("i=" + i);
            loopEdges.setColorIndex(i, UNPICKEDCOLOR);
        }

        //- Make it pickable
        loopEdges.setCapability(IndexedLineArray.ALLOW_COORDINATE_READ);
        loopEdges.setCapability(IndexedLineArray.ALLOW_FORMAT_READ);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COUNT_READ);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COORDINATE_INDEX_READ);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COLOR_WRITE);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COLOR_INDEX_WRITE);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COLOR_READ);
        loopEdges.setCapability(IndexedLineArray.ALLOW_COLOR_INDEX_READ);

        this.setGeometry(loopEdges);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }


    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on Loop " + patch_.getPatchNo());
        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

        for(int i = 0; i < pickResult.numIntersections(); i++)
        {
            PickIntersection pi = pickResult.getIntersection(i);
            int[] indices = pi.getPrimitiveVertexIndices();

            System.out.print("-- line consisting of local points ");
            for(int j = 0; j < indices.length; j++)
            {
                //int vertexLabel = vertexLabels_[indices[j]];
                //Point3f coord = new Point3f();
                //lines.getCoordinate(j, coord);
                System.out.print(indices[j] + " ");

                lines.setColorIndex(indices[j], PICKEDCOLOR);
            }
            System.out.println();
        }
        System.out.println();
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUNPicked on Loop " + patch_.getPatchNo());
        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

        for(int i = 0; i < pickResult.numIntersections(); i++)
        {
            PickIntersection pi = pickResult.getIntersection(i);
            int[] indices = pi.getPrimitiveVertexIndices();

            System.out.print("-- line consisting of local points ");
            for(int j = 0; j < indices.length; j++)
            {
                System.out.print(indices[j] + " ");

                lines.setColorIndex(indices[j], UNPICKEDCOLOR);
            }
            System.out.println();
        }
        System.out.println();
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

} // end of class Loop
