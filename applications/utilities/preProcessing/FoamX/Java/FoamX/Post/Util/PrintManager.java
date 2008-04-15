package FoamX.Post.Util;

import java.io.*;
import com.sun.image.codec.jpeg.*;
import sun.awt.image.codec.*;
import com.sun.j3d.utils.universe.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.event.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import FoamX.Post.Util.OffScreenCanvas3D;
import FoamX.Post.PostApp;
import FoamX.App;

public class PrintManager
{
    PostApp postApp_;
    SimpleUniverse universe_;
    Switch printSwitch_;
    Dimension currentDim_;
    OffScreenCanvas3D offScreenCanvas_;
    BranchGroup switchGroup_;

    // Switch to control rendering
    public PrintManager(PostApp postApp, Switch printSwitch)
    {
        postApp_ = postApp;
        universe_ = postApp_.getUniverse();
        printSwitch_ = printSwitch;
        currentDim_ = null;
        offScreenCanvas_ = null;
        switchGroup_ = null;

        // Make sure we can modify switch and its children
        printSwitch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        printSwitch_.setCapability(Switch.ALLOW_CHILDREN_EXTEND);
        printSwitch_.setCapability(Switch.ALLOW_CHILDREN_WRITE);
        //printSwitch_.setWhichChild(Switch.CHILD_NONE);
        printSwitch_.setWhichChild(Switch.CHILD_ALL);
    }


    public void print(Dimension wantedDim, String fileName)
    {
        Dimension windowDim = postApp_.getSize();
        if (wantedDim.width == 0)
        {
            if (wantedDim.height == 0)
            {
                wantedDim = windowDim;
                System.out.println("Using screen size : " + wantedDim);
            }
            else
            {
                Float wantedWidth =
                    new Float
                    (
                        wantedDim.height
                      * (float)windowDim.width
                      / windowDim.height
                    );
                wantedDim =
                    new Dimension(wantedWidth.intValue(), wantedDim.height);
                System.out.println("Using specified height : " + wantedDim);
            }
        }
        else
        {
            if (wantedDim.height == 0)
            {
                Float wantedHeight =
                    new Float
                    (
                        wantedDim.width
                      * (float)windowDim.height
                      / windowDim.width
                    );
                wantedDim =
                    new Dimension(wantedDim.width, wantedHeight.intValue());
                System.out.println("Using specified width : " + wantedDim);
            }
            else
            {
                System.out.println
                (
                    "Using specified width and height : " + wantedDim
                );
            }
        }


        // Check if different resolution
        if
        (
            (currentDim_ != null)
         && !currentDim_.equals(wantedDim)
        )
        {
            // Disconnect current since wrong resolution
            System.out.println("Different resolution:" + currentDim_);
            printSwitch_.removeChild(switchGroup_);
            switchGroup_ = null;

            // signal below to create new dim
            currentDim_ = null;
        }

        // Create new print raster
        if (currentDim_ == null)
        {
            currentDim_ = wantedDim;
            
            System.out.println("Creating raster:" + currentDim_);

            BufferedImage bImage =
                new BufferedImage
                (
                    currentDim_.width,
                    currentDim_.height,
                    BufferedImage.TYPE_INT_ARGB
                );

            ImageComponent2D buffer =
                new ImageComponent2D
                (
                    ImageComponent.FORMAT_RGBA,
                    bImage
                );
            buffer.setCapability(ImageComponent2D.ALLOW_IMAGE_READ);

            Raster drawRaster =
                new Raster
                (
                    new Point3f(1.0f, 1.0f, 1.0f),
                    Raster.RASTER_COLOR,
                    0, 0, 
                    currentDim_.width, currentDim_.height,
                    buffer,
                    null
                );

            drawRaster.setCapability(Raster.ALLOW_IMAGE_WRITE);

            // Make shape out of it
            Shape3D shape = new Shape3D(drawRaster);

            // Hang under bg under switch
            switchGroup_ = new BranchGroup();
            switchGroup_.setCapability(BranchGroup.ALLOW_DETACH);
            switchGroup_.addChild(shape);

            printSwitch_.addChild(switchGroup_);

            GraphicsConfiguration config =
                SimpleUniverse.getPreferredConfiguration();
            offScreenCanvas_ = new OffScreenCanvas3D(config, currentDim_, true);

            // set the offscreen to match the onscreen
            Screen3D sOff = offScreenCanvas_.getScreen3D();

            sOff.setSize(currentDim_);

            Screen3D sOn = universe_.getCanvas().getScreen3D();
            double pixelWidth = sOn.getPhysicalScreenWidth();
            double pixelHeight = sOn.getPhysicalScreenHeight();

            pixelHeight = pixelWidth * currentDim_.height / currentDim_.width;
//            pixelHeight = pixelWidth;

            System.out.println
            (
                "Setting physscreen to " + pixelWidth + " " + pixelHeight
            );

            sOff.setPhysicalScreenWidth(pixelWidth);
            sOff.setPhysicalScreenHeight(pixelHeight);
//            if (sOff.getPhysicalScreenWidth() < 1E-8)
//            {
//                System.out.println
//                (
//                    "Warning: physicalScreenWidth not set; using 0.5"
//                );
//                sOff.setPhysicalScreenWidth(0.5);
//            }
//            if (sOff.getPhysicalScreenHeight() < 1E-8)
//            {
//                System.out.println
//                (
//                    "Warning: physicalScreenHeight not set; using 1.0"
//                );
//                sOff.setPhysicalScreenHeight(1.0);
//            }
//            sOff.setSize(currentDim_);
        }


        //
        // Now valid raster and hung under switch in scene. All set to start
        // printing.

        offScreenCanvas_.setFileName(fileName);
        // Set switch on
        printSwitch_.setWhichChild(Switch.CHILD_ALL);
        // Add to screen view so it gets rendered with correct view.
        universe_.getViewer().getView().addCanvas3D(offScreenCanvas_);
        // Print
        offScreenCanvas_.print();

        int sleep = App.getOptions().getProperty("FoamX.Sleep", "1000");

        try
        {
            while (offScreenCanvas_.isPrinting())
            {
                System.out.println("Waiting a bit more");
                Thread.currentThread().sleep(sleep);     
            }
            System.out.println("Finished printing");
        }
        catch (InterruptedException ie)
        {
            System.out.println("Interrupted exception");
        }

        // Disable switch
        printSwitch_.setWhichChild(Switch.CHILD_NONE);

        try
        {
            Thread.currentThread().sleep(sleep);     
        }
        catch (InterruptedException ie)
        {
            System.out.println("Interrupted exception");
        }

        // Remove canvas from view
        universe_.getViewer().getView().removeCanvas3D(offScreenCanvas_);

        try
        {
            Thread.currentThread().sleep(sleep);     
        }
        catch (InterruptedException ie)
        {
            System.out.println("Interrupted exception");
        }
    }
}
