/*
 */

package FoamX.Post.Shapes;


import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.InteractiveNodes.FeatureLinesNode;
import FoamX.Post.Mesh.PrimitivePatch;

public class FeatureLines extends Shape3D implements PickableShape
{
    static final int UNPICKEDCOLOR = 0;
    static final int PICKEDCOLOR = 1;

    static final Color3f[] colors =
    {
        new Color3f(1.0f, 1.0f, 0.0f),
        new Color3f(1.0f, 0.0f, 1.0f)
    };


    PrimitivePatch patch_;

    //--------------------------------------------------------------------------
    
    public FeatureLines(FeatureLinesNode featureLinesNode, PrimitivePatch patch)
    {
        // Store reference to original patch
        patch_ = patch;

        int[] edgeList = patch_.getFeatureLines();
        Point3f[] coords = patch_.getPoints();

        IndexedLineArray lines = new IndexedLineArray
        (
            coords.length,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3,
            edgeList.length
        );
            
        lines.setCoordinates(0, coords);
        lines.setCoordinateIndices(0, edgeList);

        lines.setColors(0, colors);
        for(int i = 0; i < coords.length; i++)
        {
            lines.setColorIndex(i, UNPICKEDCOLOR);
        }

        //- Make it pickable
        lines.setCapability(IndexedLineArray.ALLOW_COORDINATE_READ);
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
    
    public PrimitivePatch getPatch()
    {
        return patch_;
    }

    //--------------------------------------------------------------------------
    
    public int numSelectedEdges()
    {
        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

        int[] colorIndices = new int[lines.getIndexCount()];
        lines.getColorIndices(0, colorIndices);

        //
        // Collect edges with both vertices same color
        //

        int i = 0;
        int nEdges = 0;
        while (i < colorIndices.length)
        {
            if
            (
                (colorIndices[i] == PICKEDCOLOR)
             && (colorIndices[i+1] == PICKEDCOLOR)
            )
            {
                nEdges++;
            }

            i += 2;
        }

        return nEdges;
    }

        
    //--------------------------------------------------------------------------
    
    public int[] getSelectedEdges()
    {
        IndexedLineArray lines = (IndexedLineArray)this.getGeometry();

        int[] edgeList = patch_.getFeatureLines();

        int[] loopEdges = new int[lines.getIndexCount()];

        int[] colorIndices = new int[lines.getIndexCount()];
        lines.getColorIndices(0, colorIndices);

//        System.out.println("***Colors");
//        for (int j = 0; j < colorIndices.length; j++)
//        {
//            System.out.print(colorIndices[j] + " ");
//        }
//        System.out.println();


        //
        // Collect edges with both vertices same color
        //

        int i = 0;
        int nIndices = 0;
        while (i < edgeList.length)
        {
            int start = edgeList[i];
            int end =  edgeList[i+1];

            System.out.println("start:" + start + "  end:" + end);
            System.out.println("color1:" + lines.getColorIndex(i));
            System.out.println("color2:" + lines.getColorIndex(i+1));

            if
            (
                (colorIndices[i] == PICKEDCOLOR)
             && (colorIndices[i+1] == PICKEDCOLOR)
            )
            {
                loopEdges[nIndices++] = start;
                loopEdges[nIndices++] = end;
            }

            i += 2;
        }


        System.out.println("Found " + nIndices/2 + " selected edges.");


        // Compact

        int[] newEdges = new int[nIndices];
        for (int j = 0; j < nIndices; j++)
        {
            newEdges[j] = loopEdges[j];
        }

        return newEdges;
    }
        

    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on patch " + patch_.getPatchNo());

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
        System.out.println("\n\nUnpicked on patch " + patch_.getPatchNo());

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
                //int vertexLabel = vertexLabels_[indices[j]];
                //Point3f coord = new Point3f();
                //lines.getCoordinate(j, coord);
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

} // end of class FeatureLines
