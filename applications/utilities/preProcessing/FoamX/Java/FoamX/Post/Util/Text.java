/*
 * Text.java
 */


package FoamX.Post.Util;

import java.awt.*;
import java.awt.Color.*;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.*;

import javax.media.j3d.*;
import javax.vecmath.*;


public class Text
{
    Font smallFont_ = new java.awt.Font("Courier", 0, 8);


    // creates a simple Raster text label (similar to Text2D)
    public Shape3D createLabel
    (
        String szText, float x, float y, float z, Color color
    )
    {
	BufferedImage bufferedImage =
            new BufferedImage(25, 14, BufferedImage.TYPE_INT_RGB);
	Graphics g = bufferedImage.getGraphics();
//        g.setComposite
//        (
//            AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f)
//        );
        g.setColor(color);
        g.setFont(smallFont_);
	g.drawString(szText, 2, 12);

	ImageComponent2D imageComponent2D =
            new ImageComponent2D(ImageComponent2D.FORMAT_RGB, bufferedImage);

	// create the Raster for the image
	javax.media.j3d.Raster renderRaster =
            new javax.media.j3d.Raster
            (
                new Point3f(x, y, z),
		javax.media.j3d.Raster.RASTER_COLOR,
		0, 0,
		bufferedImage.getWidth(),
		bufferedImage.getHeight(),
		imageComponent2D,
		null
            );

	return new Shape3D(renderRaster);
    }


} // end of class Text
