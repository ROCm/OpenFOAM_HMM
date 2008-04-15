/*
 */
package FoamX.Post.Util;

//import com.sun.j3d.utils.geometry.GeometryInfo;
//import com.sun.j3d.utils.geometry.NormalGenerator;
//import com.sun.j3d.utils.geometry.Stripifier;
//import com.sun.j3d.utils.geometry.Triangulator;
//import com.sun.j3d.utils.picking.PickIntersection;
//import com.sun.j3d.utils.picking.PickResult;
//import com.sun.j3d.utils.picking.PickTool;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.Graphics2D;

import org.j3d.util.interpolator.ColorInterpolator;

import javax.media.j3d.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import FoamX.Post.PostApp;

public class Colorizer
{
    PostApp pp_;
    float min_;
    float max_;
    Color3f minColor_;
    Color3f maxColor_;
    Color3f[] colorMap_;

    ColorInterpolator colorInterpolator_;
    float[] keys_;
    float[] colors_;

    Font smallFont_;
    

    // n points according to eye curve
    public Colorizer(PostApp pp, int nKeys)
    {
        pp_ = pp;
        keys_ = new float[nKeys];
        colors_ = new float[3*nKeys];
        smallFont_ = new Font("Courier", 0, 8);

        // Eye curve
        float curve = 1.4f;
        float bias = 1.0f;
        float rfact = 0.5f * bias;

        float scale = 1.0f / (nKeys-1);
        int colorI = 0;
        for(int keyI = 0; keyI < nKeys; keyI++)
        {
            // compute s in [0,1]
            float s = (float) keyI * scale;
            float t = curve * (s - rfact);   // t in [curve*-0.5,curve*0.5)

            keys_[keyI] = s;
            colors_[colorI++] = (float) (0.5 + 0.5 * Math.atan( 7.0*t ) / 1.57);
            colors_[colorI++] = (float) (0.5 + 0.5 * (2 * Math.exp(-7*t*t) - 1));
            colors_[colorI++] = (float) (0.5 + 0.5 * Math.atan( -7.0*t ) / 1.57);
        }

        setColorInterpolator();
    }

    //--------------------------------------------------------------------------

    // simplest one:
    public Colorizer(PostApp pp)
    {
        pp_ = pp;

        keys_ = new float[4];
        colors_ = new float[3*4];

        // min = blue (001)
        keys_[0] = 0.0f;
        colors_[0] = 0.0f;
        colors_[1] = 0.0f;
        colors_[2] = 1.0f;

        // green (010)
        keys_[1] = 0.5f;
        colors_[3] = 0.0f;
        colors_[4] = 1.0f;
        colors_[5] = 0.0f;

        // yellow (110)
        keys_[2] = 0.75f;
        colors_[6] = 1.0f;
        colors_[7] = 1.0f;
        colors_[8] = 0.0f;

        // max = red (100)
        keys_[3] = 1.0f;
        colors_[9] = 1.0f;
        colors_[10] = 0.0f;
        colors_[11] = 0.0f;

        setColorInterpolator();
    }

    //--------------------------------------------------------------------------


    // Reads in key/rbg pairs from file.
    // Note: key values should be [0..1]
    public Colorizer(PostApp pp, File file)
    {
        pp_ = pp;
        System.out.println("Colorizer read");
    }

    //--------------------------------------------------------------------------

    public void write(PostApp pp, File file)
    {
        pp_ = pp;
        System.out.println("Colorizer write");
    }

    //--------------------------------------------------------------------------

    private void setColorInterpolator()
    {
        //System.out.println("Constructing colorizer:");

        // Fill ColorInterPolator
        colorInterpolator_ =
            new ColorInterpolator
            (
                keys_.length,
                ColorInterpolator.RGB_SPACE
            );

        int colorI = 0;
        for(int keyI = 0; keyI < keys_.length; keyI++)
        {
            float R = colors_[colorI++];
            float G = colors_[colorI++];
            float B = colors_[colorI++];
            //System.out.println
            //(
            //    "    key: " + keys_[keyI] + " RGB: "
            //  + R + " " + G + " " + B
            //);
            colorInterpolator_.addRGBKeyFrame
            (
                keys_[keyI],
                R,
                G,
                B,
                1.0f
            );
        }
    }

    //--------------------------------------------------------------------------

    // Convert values into colors
    public void color
    (
        final float min,
        final float max,
        final float[] values,
        float[] colors
    )
    {
        if (3*values.length != colors.length)
        {
            System.out.println("ERROR: colors array not big enough");
            System.out.println("    values:" + values.length);
            System.out.println("    colors:" + colors.length);
            return;
        }

        int colorI = 0;

        float scale = 1.0f / (max - min);
        for(int i = 0; i < values.length; i++)
        {
            // Scale to [0,1]
            float scaledValue = (values[i] - min) * scale;
            float[] rgbColor = colorInterpolator_.floatRGBValue(scaledValue);

            // Copy RGB
            colors[colorI++] = rgbColor[0];
            colors[colorI++] = rgbColor[1];
            colors[colorI++] = rgbColor[2];
        }
    }

    //--------------------------------------------------------------------------

    // Make 2D colorbar of dimensions 1x1
    // No text for now.
    public TriangleArray colorBar(final float min, final float max)
    {
        if (keys_.length == 0)
        {
            return null;
        }
        
        int nTri = 2*(keys_.length - 1);
        int nVert = 3*nTri;

        Point3f[] coords = new Point3f[nVert];
        float[] vertColors = new float[3*nVert];

        // Lower coordinates
        Point3f lowerLeft = new Point3f(0.0f, 0.0f, 0.0f);
        Point3f lowerRight = new Point3f(1.0f, 0.0f, 0.0f);

        int vertI = 0;
        int colorI = 0;
        for(int keyI = 0; keyI < keys_.length - 1; keyI++)
        {
            // Upper coordinates
            Point3f upperLeft =
                new Point3f(lowerLeft.x, keys_[keyI+1], lowerLeft.z);
            Point3f upperRight =
                new Point3f(lowerRight.x, keys_[keyI+1], lowerRight.z);

            // Triangle 1

            // copy coordinates + RGB
            coords[vertI++] = lowerLeft;
            vertColors[colorI++] = colors_[3*keyI];
            vertColors[colorI++] = colors_[3*keyI+1];
            vertColors[colorI++] = colors_[3*keyI+2];
            
            coords[vertI++] = lowerRight;
            vertColors[colorI++] = colors_[3*keyI];
            vertColors[colorI++] = colors_[3*keyI+1];
            vertColors[colorI++] = colors_[3*keyI+2];
            
            coords[vertI++] = upperRight;
            vertColors[colorI++] = colors_[3*(keyI+1)];
            vertColors[colorI++] = colors_[3*(keyI+1)+1];
            vertColors[colorI++] = colors_[3*(keyI+1)+2];

            // Triangle 2
            coords[vertI++] = upperRight;
            vertColors[colorI++] = colors_[3*(keyI+1)];
            vertColors[colorI++] = colors_[3*(keyI+1)+1];
            vertColors[colorI++] = colors_[3*(keyI+1)+2];

            coords[vertI++] = upperLeft;
            vertColors[colorI++] = colors_[3*(keyI+1)];
            vertColors[colorI++] = colors_[3*(keyI+1)+1];
            vertColors[colorI++] = colors_[3*(keyI+1)+2];

            coords[vertI++] = lowerLeft;
            vertColors[colorI++] = colors_[3*keyI];
            vertColors[colorI++] = colors_[3*keyI+1];
            vertColors[colorI++] = colors_[3*keyI+2];

            // Move one up
            lowerLeft = upperLeft;
            lowerRight = upperRight;
        }

        TriangleArray array = new TriangleArray
        (
            3*nTri,
            TriangleArray.COORDINATES
          | TriangleArray.COLOR_3
        );
        array.setCoordinates(0, coords);
        array.setColors(0, vertColors);

        return array;
    }

    //--------------------------------------------------------------------------

    public void paint
    (
        Graphics2D g,
        final int width,
        final int height,
        final float min,
        final float max
    )
    {
        // Make table of values from min to max.
        float[] values = new float[height];
        float val = min;
        float valInc = (max - min) / height;
        for(int y = 0; y < height; y++)
        {
            values[y] = val;
            val += valInc;
        }

        // Convert into colors
        float[] colors = new float[3*height];
        color(min, max, values, colors);

        // Paint into Graphics2D. Note pixel coords start at 0,0 top left

        g.setRenderingHint
        (
            RenderingHints.KEY_TEXT_ANTIALIASING,
            RenderingHints.VALUE_TEXT_ANTIALIAS_OFF
        );
        g.setRenderingHint
        (
            RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_OFF
        );
        g.setRenderingHint
        (
            RenderingHints.KEY_FRACTIONALMETRICS,
            RenderingHints.VALUE_FRACTIONALMETRICS_OFF
        );
        g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1.0f));
        
        // Draw coloured lines starting at right half

        int barX = width / 2;
        int colorI = 0;
        for(int pixelY = height -1 ; pixelY >= 0; pixelY--)
        {
            Color drawColor =
                new Color
                (
                    colors[colorI++], colors[colorI++], colors[colorI++]
                );

            g.setPaint(drawColor);
            g.drawLine(barX, pixelY, width, pixelY);
        }

        // Draw ticks and put numbers
        int tickX = barX - 5;


        FontMetrics fm = g.getFontMetrics();
        int textHeight = fm.getHeight();

        //String text = "abc";
        //int textWidth = fm.stringWidth(text);
        //Rectangle2D dims =
        //    fm.getStringBounds(text.toCharArray(), 0, 2, g);

        if (g.getFont() != smallFont_)
        {
            g.setFont(smallFont_);
        }
        g.setPaint(pp_.getForegroundColor());

        // min
        g.drawLine(tickX, height - 1, barX, height - 1);
        Float floatY = new Float(min);
        int screenY = height - 1 - fm.getLeading();
        g.drawString(floatY.toString(), 0, screenY);

        // max
        g.drawLine(tickX, 0, barX, 0);
        floatY = new Float(max);
        screenY = textHeight;
        g.drawString(floatY.toString(), 1, screenY);

//        int nTicks = 5;
//        float range = nTicks*(max - min);
//        for(int i = 0; i < nTicks; )
//
//        g.setPaint(pp_.getForeground());
//        for(int keyI = 0; keyI < keys_.length; keyI++)
//        {
//            Float floatY =
//                new Float(height * (1.0f - (keys_[keyI] - min) / (max - min)));
//            int y = floatY.intValue();
//            g.drawLine(0, y, width, y);
//            g.drawString(floatY.toString(), 0, y);
//        }
    }

    //--------------------------------------------------------------------------
}
