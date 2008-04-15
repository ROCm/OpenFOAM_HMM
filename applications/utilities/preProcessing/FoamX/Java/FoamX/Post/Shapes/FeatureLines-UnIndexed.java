/*
 */

//import java.awt.*;
//import java.awt.Font;

//import com.sun.j3d.utils.applet.MainFrame;
//import com.sun.j3d.utils.geometry.Text2D;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

public class FeatureLines extends Shape3D implements InteractiveNode
{

    static final Point3d[] coords = 
    {
        new Point3d(0,   0, 0),
        new Point3d(1,   0, 0),
        new Point3d(1,   1, 0),
        new Point3d(0.5, 2, 0),
        new Point3d(0,   1, 0),

        new Point3d(0,   0, 1),
        new Point3d(1,   0, 1),
        new Point3d(1,   1, 1),
        new Point3d(0.5, 2, 1),
        new Point3d(0,   1, 1)
    };

    static final Color3f UNPICKEDCOLOR = new Color3f(1.0f, 1.0f, 0.0f);
    static final Color3f PICKEDCOLOR = new Color3f(1.0f, 0.0f, 1.0f);


    int patch_;

    int[] vertexLabels_;
    Color3f[] vertexColors_;


    public FeatureLines(int patch, float angle)
    {
        patch_ = patch;

        LineArray axisLines = null;

        // Just hardcoded for now.

        if (patch_ == 0)
        {
            vertexLabels_ = new int[10];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 4;
            vertexLabels_[8] = 4;
            vertexLabels_[9] = 0;

            vertexColors_ =
                fillColors
                (
                    10,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                10,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[0]);
            axisLines.setCoordinate(1, coords[1]);
            axisLines.setCoordinate(2, coords[1]);
            axisLines.setCoordinate(3, coords[2]);
            axisLines.setCoordinate(4, coords[2]);
            axisLines.setCoordinate(5, coords[3]);
            axisLines.setCoordinate(6, coords[3]);
            axisLines.setCoordinate(7, coords[4]);
            axisLines.setCoordinate(8, coords[4]);
            axisLines.setCoordinate(9, coords[0]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 1)
        {
            vertexLabels_ = new int[10];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 4;
            vertexLabels_[8] = 4;
            vertexLabels_[9] = 0;

            vertexColors_ =
                fillColors
                (
                    10,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                10,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[5]);
            axisLines.setCoordinate(1, coords[6]);
            axisLines.setCoordinate(2, coords[6]);
            axisLines.setCoordinate(3, coords[7]);
            axisLines.setCoordinate(4, coords[7]);
            axisLines.setCoordinate(5, coords[8]);
            axisLines.setCoordinate(6, coords[8]);
            axisLines.setCoordinate(7, coords[9]);
            axisLines.setCoordinate(8, coords[9]);
            axisLines.setCoordinate(9, coords[5]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 2)
        {
            vertexLabels_ = new int[8];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 0;

            vertexColors_ =
                fillColors
                (
                    8,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                8,
                GeometryArray.COLOR_3 | GeometryArray.COORDINATES
            );

            axisLines.setCoordinate(0, coords[0]);
            axisLines.setCoordinate(1, coords[5]);
            axisLines.setCoordinate(2, coords[5]);
            axisLines.setCoordinate(3, coords[9]);
            axisLines.setCoordinate(4, coords[9]);
            axisLines.setCoordinate(5, coords[4]);
            axisLines.setCoordinate(6, coords[4]);
            axisLines.setCoordinate(7, coords[0]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 3)
        {
            vertexLabels_ = new int[8];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 0;

            vertexColors_ =
                fillColors
                (
                    8,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                8,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[4]);
            axisLines.setCoordinate(1, coords[9]);
            axisLines.setCoordinate(2, coords[9]);
            axisLines.setCoordinate(3, coords[8]);
            axisLines.setCoordinate(4, coords[8]);
            axisLines.setCoordinate(5, coords[3]);
            axisLines.setCoordinate(6, coords[3]);
            axisLines.setCoordinate(7, coords[4]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 4)
        {
            vertexLabels_ = new int[8];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 0;

            vertexColors_ =
                fillColors
                (
                    8,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                8,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[3]);
            axisLines.setCoordinate(1, coords[8]);
            axisLines.setCoordinate(2, coords[8]);
            axisLines.setCoordinate(3, coords[7]);
            axisLines.setCoordinate(4, coords[7]);
            axisLines.setCoordinate(5, coords[2]);
            axisLines.setCoordinate(6, coords[2]);
            axisLines.setCoordinate(7, coords[3]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 5)
        {
            vertexLabels_ = new int[8];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 0;

            vertexColors_ =
                fillColors
                (
                    8,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                8,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[2]);
            axisLines.setCoordinate(1, coords[7]);
            axisLines.setCoordinate(2, coords[7]);
            axisLines.setCoordinate(3, coords[6]);
            axisLines.setCoordinate(4, coords[6]);
            axisLines.setCoordinate(5, coords[1]);
            axisLines.setCoordinate(6, coords[1]);
            axisLines.setCoordinate(7, coords[2]);

            axisLines.setColors(0, vertexColors_);
        }
        else if (patch_ == 6)
        {
            vertexLabels_ = new int[8];
            vertexLabels_[0] = 0;
            vertexLabels_[1] = 1;
            vertexLabels_[2] = 1;
            vertexLabels_[3] = 2;
            vertexLabels_[4] = 2;
            vertexLabels_[5] = 3;
            vertexLabels_[6] = 3;
            vertexLabels_[7] = 0;

            vertexColors_ =
                fillColors
                (
                    8,
                    UNPICKEDCOLOR
                );

            axisLines = new LineArray
            (
                8,
                GeometryArray.COORDINATES | GeometryArray.COLOR_3
            );

            axisLines.setCoordinate(0, coords[1]);
            axisLines.setCoordinate(1, coords[6]);
            axisLines.setCoordinate(2, coords[6]);
            axisLines.setCoordinate(3, coords[5]);
            axisLines.setCoordinate(4, coords[5]);
            axisLines.setCoordinate(5, coords[0]);
            axisLines.setCoordinate(6, coords[0]);
            axisLines.setCoordinate(7, coords[1]);

            axisLines.setColors(0, vertexColors_);
        }
        else
        {
            System.out.println("**\n** Unknown patch " + patch_ + "**\n**");
            vertexLabels_ = null;
            vertexColors_ = null;
        }
    
        axisLines.setCapability(LineArray.ALLOW_COORDINATE_READ);
        axisLines.setCapability(LineArray.ALLOW_FORMAT_READ);
        axisLines.setCapability(LineArray.ALLOW_COUNT_READ);
        axisLines.setCapability(LineArray.ALLOW_COLOR_READ);
        axisLines.setCapability(LineArray.ALLOW_COLOR_WRITE);

        this.setGeometry(axisLines);
        this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }


    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPatch " + patch_);

        LineArray lines = (LineArray)this.getGeometry();

        System.out.println("-- consisting of vertices ");

        for(int j = 0; j < lines.getVertexCount(); j++)
        {
            Point3f coord = new Point3f();
            lines.getCoordinate(j, coord);
           System.out.print(coord + " ");
        }
        System.out.println();


        for(int i = 0; i < pickResult.numIntersections(); i++)
        {
            PickIntersection pi = pickResult.getIntersection(i);
            int[] indices = pi.getPrimitiveVertexIndices();

            System.out.print("-- line consisting of local points ");
            for(int j = 0; j < indices.length; j++)
            {
                int vertexLabel = vertexLabels_[indices[j]];
                //Point3f coord = new Point3f();
                //lines.getCoordinate(j, coord);
                System.out.print(vertexLabel + "(" + indices[j] + ") ");

                lines.setColor(indices[j], PICKEDCOLOR);
            }
            System.out.println();
        }
        System.out.println();
    }


    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nPatch " + patch_);

        LineArray lines = (LineArray)this.getGeometry();

        System.out.println("-- consisting of vertices ");

        for(int j = 0; j < lines.getVertexCount(); j++)
        {
            Point3f coord = new Point3f();
            lines.getCoordinate(j, coord);
           System.out.print(coord + " ");
        }
        System.out.println();


        for(int i = 0; i < pickResult.numIntersections(); i++)
        {
            PickIntersection pi = pickResult.getIntersection(i);
            int[] indices = pi.getPrimitiveVertexIndices();

            System.out.print("-- line consisting of local points ");
            for(int j = 0; j < indices.length; j++)
            {
                int vertexLabel = vertexLabels_[indices[j]];
                //Point3f coord = new Point3f();
                //lines.getCoordinate(j, coord);
                System.out.print(vertexLabel + "(" + indices[j] + ") ");

                lines.setColor(indices[j], UNPICKEDCOLOR);
            }
            System.out.println();
        }
        System.out.println();
    }

    // Utility function to create color array
    static Color3f[] fillColors(int length, Color3f color)
    {
        Color3f[] colors = new Color3f[length];

        for(int i = 0; i < length; i++)
        {
            colors[i] = color;
        }
        return colors;
    }

} // end of class FeatureLines
