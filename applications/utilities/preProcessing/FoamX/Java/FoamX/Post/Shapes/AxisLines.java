/*
 * AxisLines.java
 */


package FoamX.Post.Shapes;

import java.awt.*;
//import java.awt.Font;

import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.geometry.Text2D;
import com.sun.j3d.utils.picking.PickResult;

import javax.media.j3d.*;
import javax.vecmath.*;

import FoamX.Post.InteractiveNodes.AxisNode;

public class AxisLines extends Shape3D implements PickableShape
{
    public AxisLines(AxisNode axisNode, Color3f color)
    {
        this.setGeometry(createGeometry(color));

        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }

    //--------------------------------------------------------------------------


    private Geometry createGeometry(Color3f color)
    {
        // create line for X axis
        IndexedLineArray axisLines = new IndexedLineArray
        (
            18,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3,
            30
        );

        axisLines.setCoordinate( 0, new Point3f( 0.0f, 0.0f, 0.0f));
        axisLines.setCoordinate( 1, new Point3f( 1.0f, 0.0f, 0.0f));
        axisLines.setCoordinate( 2, new Point3f( 0.9f, 0.1f, 0.1f));
        axisLines.setCoordinate( 3, new Point3f( 0.9f,-0.1f, 0.1f));
        axisLines.setCoordinate( 4, new Point3f( 0.9f, 0.1f,-0.1f));
        axisLines.setCoordinate( 5, new Point3f( 0.9f,-0.1f,-0.1f));
        axisLines.setCoordinate( 6, new Point3f( 0.0f, 0.0f, 0.0f));
        axisLines.setCoordinate( 7, new Point3f( 0.0f, 1.0f, 0.0f));
        axisLines.setCoordinate( 8, new Point3f( 0.1f, 0.9f, 0.1f));
        axisLines.setCoordinate( 9, new Point3f(-0.1f, 0.9f, 0.1f));
        axisLines.setCoordinate(10, new Point3f( 0.1f, 0.9f,-0.1f));
        axisLines.setCoordinate(11, new Point3f(-0.1f, 0.9f,-0.1f));
        axisLines.setCoordinate(12, new Point3f( 0.0f, 0.0f, 0.0f));
        axisLines.setCoordinate(13, new Point3f( 0.0f, 0.0f, 1.0f));
        axisLines.setCoordinate(14, new Point3f( 0.1f, 0.1f, 0.9f));
        axisLines.setCoordinate(15, new Point3f(-0.1f, 0.1f, 0.9f));
        axisLines.setCoordinate(16, new Point3f( 0.1f,-0.1f, 0.9f));
        axisLines.setCoordinate(17, new Point3f(-0.1f,-0.1f, 0.9f));

        axisLines.setCoordinateIndex( 0, 0);
        axisLines.setCoordinateIndex( 1, 1);
        axisLines.setCoordinateIndex( 2, 2);
        axisLines.setCoordinateIndex( 3, 1);
        axisLines.setCoordinateIndex( 4, 3);
        axisLines.setCoordinateIndex( 5, 1);
        axisLines.setCoordinateIndex( 6, 4);
        axisLines.setCoordinateIndex( 7, 1);
        axisLines.setCoordinateIndex( 8, 5);
        axisLines.setCoordinateIndex( 9, 1);
        axisLines.setCoordinateIndex(10, 6);
        axisLines.setCoordinateIndex(11, 7);
        axisLines.setCoordinateIndex(12, 8);
        axisLines.setCoordinateIndex(13, 7);
        axisLines.setCoordinateIndex(14, 9);
        axisLines.setCoordinateIndex(15, 7);
        axisLines.setCoordinateIndex(16,10);
        axisLines.setCoordinateIndex(17, 7);
        axisLines.setCoordinateIndex(18,11);
        axisLines.setCoordinateIndex(19, 7);
        axisLines.setCoordinateIndex(20,12);
        axisLines.setCoordinateIndex(21,13);
        axisLines.setCoordinateIndex(22,14);
        axisLines.setCoordinateIndex(23,13);
        axisLines.setCoordinateIndex(24,15);
        axisLines.setCoordinateIndex(25,13);
        axisLines.setCoordinateIndex(26,16);
        axisLines.setCoordinateIndex(27,13);
        axisLines.setCoordinateIndex(28,17);
        axisLines.setCoordinateIndex(29,13);

        axisLines.setColor(0, color);
        for(int i = 0; i < 30; i++)
        {
            axisLines.setColorIndex(i, 0);
        }

        axisLines.setCapability(IndexedLineArray.ALLOW_COORDINATE_READ);
        axisLines.setCapability(IndexedLineArray.ALLOW_FORMAT_READ);
        axisLines.setCapability(IndexedLineArray.ALLOW_COUNT_READ);
        axisLines.setCapability(IndexedLineArray.ALLOW_COORDINATE_INDEX_READ);

        return axisLines;

    } // end of createGeometry()

    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("ERROR: AxisLines Picked:");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("ERROR: AxisLines unPicked:");
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        System.out.println("remove axis:");

        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }
    //--------------------------------------------------------------------------



} // end of class AxisLines
