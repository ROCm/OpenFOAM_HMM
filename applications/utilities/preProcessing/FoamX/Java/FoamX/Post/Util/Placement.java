/*
 * Placement.java
 */


package FoamX.Post.Util;

import java.awt.*;
import javax.media.j3d.*;
import javax.vecmath.*;

/**
 * Abstracts positioning of rectangle inside 2D window. 
 * Specify corner and spacing from this corner and nearest corner on component.
 *
 * getTopLeft will return coordinates of top-left maintaining spacing to corner.
 */

public class Placement
{
    public static final int TOPLEFT = 1;
    public static final int TOPRIGHT = 2;
    public static final int BOTTOMLEFT = 3;
    public static final int BOTTOMRIGHT = 4;

    Dimension spacing_;
    int corner_;

    public Placement(Dimension spacing, int corner)
    {
        spacing_ = new Dimension(spacing);

        if ((corner < TOPLEFT) || (corner > BOTTOMRIGHT))
        {
            System.out.println("ERROR: Placement incorrect corner:" + corner);
            return;
        }
        corner_ = corner;
    }

    //--------------------------------------------------------------------------

    public Dimension getSpacing()
    {
        return spacing_;
    }

    //--------------------------------------------------------------------------

    public int getCorner()
    {
        return corner_;
    }

    //--------------------------------------------------------------------------

    // Calculates topleft of item inside component.
    public Dimension getTopLeft(Component window, Dimension item)
    {
        int windowWidth = window.getWidth();
        int windowHeight = window.getHeight();

        int itemWidth = item.width;
        int itemHeight = item.height;

        if (corner_ == TOPLEFT)
        {
            return spacing_;
        }
        else if (corner_ == TOPRIGHT)
        {
            return new Dimension
            (
                windowWidth - spacing_.width - itemWidth,
                spacing_.height
            );
        }
        else if (corner_ == BOTTOMLEFT)
        {
            return new Dimension
            (
                spacing_.width,
                windowHeight - spacing_.height - itemHeight
            );
        }
        else
        {
            // Bottomright
            return new Dimension
            (
                windowWidth - spacing_.width - itemWidth,
                windowHeight - spacing_.height - itemHeight
            );
        }
    }

    //--------------------------------------------------------------------------

} // end of class Placement
