/*
 * Arrows.java
 */


package FoamX.Post.Shapes;

import java.awt.*;

import com.sun.j3d.utils.picking.PickResult;

import javax.media.j3d.*;
import javax.vecmath.*;

import FoamX.Post.Util.Colorizer;

public class Arrows extends Shape3D implements PickableShape
{
    public Arrows
    (
        final float[] sampleCoords,
        final float[] vecValues,
        final float sizeScale,
        final float[] scalValues,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        this.setGeometry
        (
            createGeometry
            (
                sampleCoords,
                vecValues,
                sizeScale,
                scalValues,
                colorizer,
                min,
                max
            )
        );

        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }


    //--------------------------------------------------------------------------

    /** creates set of arrows, one per samplePoint. 
     * Arrow built of single line
     */
    static public Geometry createGeometry
    (
        final float[] sampleCoords,
        final float[] vecValues,
        final float sizeScale,
        final float[] scalValues,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        // create line for X axis
        LineArray arrowLines = new LineArray
        (
            2*sampleCoords.length,
            GeometryArray.COORDINATES | GeometryArray.COLOR_3
        );

        int vertI = 0;
        int floatI = 0;

        float[] start = new float[3];
        float[] end = new float[3];

        int nVecs = sampleCoords.length / 3;
        for(int i = 0; i < nVecs; i++)
        {
            // Get coords and (vector) value
            start[0] = sampleCoords[floatI];
            float u = vecValues[floatI++];

            start[1] = sampleCoords[floatI];
            float v = vecValues[floatI++];

            start[2] = sampleCoords[floatI];
            float w = vecValues[floatI++];

            // Get color from scalar value

            // Head of arrow
            end[0] = start[0] + sizeScale*u;
            end[1] = start[1] + sizeScale*v;
            end[2] = start[2] + sizeScale*w;

            //System.out.println
            //(
            //"i:" + i + " Start:" + start[0] + " " + start[1] + " " + start[2] + 
            //" End:" + end[0] + " " + end[1] + " " + end[2]
            //);

            arrowLines.setCoordinate(vertI++, start);
            arrowLines.setCoordinate(vertI++, end);
        }

        System.out.println("Converting values to colours ...");
        float[] colors = new float[3*scalValues.length];
        colorizer.color(min, max, scalValues, colors);

        float[] arrowColors = new float[2*colors.length];
        int compI = 0;
        int destI = 0;
        for(int i = 0; i < nVecs; i++)
        {
            float r = colors[compI++];
            float g = colors[compI++];
            float b = colors[compI++];

            // start-vertex color
            arrowColors[destI] = r;
            arrowColors[destI+1] = g;
            arrowColors[destI+2] = b;

            // end-vertex color
            arrowColors[destI+3] = r;
            arrowColors[destI+4] = g;
            arrowColors[destI+5] = b;

            destI += 6;
        }
        arrowLines.setColors(0, arrowColors);

        return arrowLines;

    } // end of createGeometry()

    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("ERROR: Arrows Picked:");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("ERROR: Arrows unPicked:");
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        //System.out.println("remove arrows:");
        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }
    //--------------------------------------------------------------------------



} // end of class Arrows
