package FoamX.Post.Shapes2D;

//import com.xith.java3d.overlay.*;
import org.j3d.geom.overlay.*;

import javax.media.j3d.*;
import javax.vecmath.*;
import java.awt.*;
import java.awt.image.*;
import java.util.*;

import FoamX.Post.Util.Colorizer;

public class ColorBar2D extends OverlayBase
{
    int width_;
    int height_;
    Colorizer colorizer_;
    float min_;
    float max_;

    public ColorBar2D
    (
        final Canvas3D c,
        final int x, final int y,
        final int width, final int height,
        final Colorizer colorizer,
        final float min, final float max
    )
    {
        super(c, new java.awt.Rectangle(x, y, width, height));
        //this.setBackgroundMode(BACKGROUND_COPY);    // Copy background
        width_ = width;
        height_ = height;
        colorizer_ = colorizer;
        min_ = min;
        max_ = max;
    }

   /**
    * This is where the actualy drawing of the window takes place.  Override
    * this to alter the contents of what is show in the window.
    */

    public void paint(Graphics2D g)
    {
        colorizer_.paint(g, width_, height_, min_, max_);
    }

    //--------------------------------------------------------------------------
}
