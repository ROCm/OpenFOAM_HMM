/*
 */

package FoamX.Post.Shapes;


import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.InteractiveNodes.MeshLinesNode;

public class MeshLines extends Shape3D implements PickableShape
{
    static final int UNPICKEDCOLOR = 0;
    static final int PICKEDCOLOR = 1;

    static final Color3f[] colors =
    {
        new Color3f(1.0f, 1.0f, 0.0f),
        new Color3f(1.0f, 0.0f, 1.0f)
    };

    private MeshLinesNode meshLinesNode_;

    //--------------------------------------------------------------------------
    
    public MeshLines
    (
        MeshLinesNode meshLinesNode, float[] coords, int[] edgeList
    )
    {
        // Store reference to original patch
        meshLinesNode_ = meshLinesNode;

        IndexedLineArray lines = new IndexedLineArray
        (
            coords.length / 3,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3,
            //GeometryArray.COORDINATES,
            edgeList.length
        );
            
        lines.setCoordinates(0, coords);
        lines.setCoordinateIndices(0, edgeList);

        lines.setColors(0, colors);
        for(int i = 0; i < edgeList.length; i++)
        {
            lines.setColorIndex(i, UNPICKEDCOLOR);
        }

        //- Make it pickable
        lines.setCapability(IndexedLineArray.ALLOW_COORDINATE_READ);
        lines.setCapability(IndexedLineArray.ALLOW_COORDINATE_WRITE);
        lines.setCapability(IndexedLineArray.ALLOW_FORMAT_READ);
        lines.setCapability(IndexedLineArray.ALLOW_COUNT_READ);
        lines.setCapability(IndexedLineArray.ALLOW_COORDINATE_INDEX_READ);
        lines.setCapability(IndexedLineArray.ALLOW_COLOR_WRITE);
        lines.setCapability(IndexedLineArray.ALLOW_COLOR_INDEX_WRITE);
        lines.setCapability(IndexedLineArray.ALLOW_COLOR_READ);
        lines.setCapability(IndexedLineArray.ALLOW_COLOR_INDEX_READ);

        this.setGeometry(lines);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }

    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
        
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\nPicked");

        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

//        System.out.println("-- consisting of vertices ");
//        for(int j = 0; j < lines.getVertexCount(); j++)
//        {
//            Point3f coord = new Point3f();
//            lines.getCoordinate(j, coord);
//            System.out.print(coord + " ");
//        }
//        System.out.println();


        for(int i = 0; i < pickResult.numIntersections(); i++)
        {
            PickIntersection pi = pickResult.getIntersection(i);
            int[] indices = pi.getPrimitiveVertexIndices();

            System.out.print("-- line consisting of local points ");
            for(int j = 0; j < indices.length; j++)
            {
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
        System.out.println("\nUnpicked");

        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

//        System.out.println("-- consisting of vertices ");
//        for(int j = 0; j < lines.getVertexCount(); j++)
//        {
//            Point3f coord = new Point3f();
//            lines.getCoordinate(j, coord);
//           System.out.print(coord + " ");
//        }
//        System.out.println();


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
        //System.out.println("remove");
        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }

} // end of class MeshLines
