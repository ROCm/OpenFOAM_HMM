/*
 *
 */

package FoamX.Post.InteractiveNodes;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import java.awt.*;
import java.awt.event.*;

import FoamX.Post.Util.Colorizer;
import FoamX.Post.Util.Placement;
import FoamX.Post.Shapes2D.ColorBar2D;

public class ColorBar2DNode
    implements InteractiveNode, ComponentListener
{
    String myName_;
    Canvas3D canvas_;
    Placement placement_;
    Dimension dimension_;
    Colorizer colorizer_;
    float min_, max_;

    Dimension topLeft_;
    Switch switch_;
    BranchGroup branchGroup_;
    BranchGroup switchGroup_;
    ColorBar2D colorBar2D_;
    boolean isShowing_;

    public ColorBar2DNode
    (
        final String myName,
        final Canvas3D canvas,
        final Placement placement,
        final Dimension dimension,
        final Colorizer ci,
        final float min,
        final float max
    )
    {
        myName_ = myName;
        canvas_ = canvas;
        placement_ = placement;
        dimension_ = dimension;
        colorizer_ = ci;
        min_ = min;
        max_ = max;

        topLeft_ =  placement_.getTopLeft(canvas_, dimension_);

//System.out.println("dimension_:" + dimension_);
//System.out.println("placement_:" + placement_.getSpacing() + " " + placement_.getCorner());

        isShowing_ = false;

        switch_ = new Switch();
        switch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        switch_.setCapability(Switch.ALLOW_CHILDREN_EXTEND);
        switch_.setCapability(Switch.ALLOW_CHILDREN_WRITE);
        switch_.setCapability(Node.ALLOW_BOUNDS_READ);

        switchGroup_ = new BranchGroup();
        switchGroup_.setCapability(BranchGroup.ALLOW_DETACH);
        switchGroup_.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        switch_.addChild(switchGroup_);

        branchGroup_ = new BranchGroup();
        branchGroup_.setCapability(Node.ALLOW_BOUNDS_READ);

        colorBar2D_ =
            new ColorBar2D
            (
                canvas_, topLeft_.width, topLeft_.height,
                dimension_.width, dimension_.height,
                colorizer_,
                min_, max_
            );
        colorBar2D_.repaint();
        colorBar2D_.setVisible(true);
        switchGroup_.addChild(colorBar2D_.getRoot());

        branchGroup_.addChild(switch_);

        canvas_.addComponentListener(this);
    }

    //--------------------------------------------------------------------------

    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    //--------------------------------------------------------------------------

    public ColorBar2D getColorBar2D()
    {
        return colorBar2D_;
    }

    //--------------------------------------------------------------------------

    public void move(Dimension newTopLeft)
    {
        if (newTopLeft.equals(topLeft_))
        {
            return;
        }
        topLeft_ = newTopLeft;


        boolean isLive = this.isShowing();
        if (isLive) 
        {
            this.hide();
        }

        //switch_.removeAllChildren();
        colorBar2D_ = null;
        switchGroup_.removeAllChildren();

        colorBar2D_ =
            new ColorBar2D
            (
                canvas_, topLeft_.width, topLeft_.height,
                dimension_.width, dimension_.height,
                colorizer_,
                min_, max_
            );

        colorBar2D_.repaint();
        colorBar2D_.setVisible(true);

        switchGroup_.addChild(colorBar2D_.getRoot());
        //switch_.addChild(switchGroup_);

        if (isLive)
        {
            this.show();
        }
    }


    //--------------------------------------------------------------------------
    //---- InteractiveNode Interface
    //--------------------------------------------------------------------------

    public String name()
    {
        return myName_;
    }

    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("Picked ColorBar2DNode");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("UnPicked ColorBar2DNode");
    }

    //--------------------------------------------------------------------------

    // Include in display
    public void show()
    {
        //System.out.println("show:" + this.name());
        switch_.setWhichChild(Switch.CHILD_ALL);
        isShowing_ = true;
    }

    //--------------------------------------------------------------------------

    // Exclude from display
    public void hide()
    {
        //System.out.println("hide:" + this.name());
        switch_.setWhichChild(Switch.CHILD_NONE);
        isShowing_ = false;
    }

    //--------------------------------------------------------------------------

    // Check status
    public boolean isShowing()
    {
        return isShowing_;
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        //System.out.println("remove:" + this.name());

        hide();
        colorBar2D_ = null;
        switch_ = null;
    }


    //--------------------------------------------------------------------------
    //---- ComponentListener Interface
    //--------------------------------------------------------------------------

    public void componentMoved(ComponentEvent e)
    {
        System.out.println("componentMoved:" + this.name());
    }

    public void componentResized(ComponentEvent e)
    {
        System.out.println("componentResized:" + canvas_.getWidth() + " " + canvas_.getHeight());

        if (placement_.getCorner() != Placement.TOPLEFT)
        {
            Dimension topLeft = placement_.getTopLeft(canvas_, dimension_);
            this.move(topLeft);
        }
    }

    public void componentHidden(ComponentEvent e)
    {
        System.out.println("componentHidden:" + this.name());
    }

    public void componentShown(ComponentEvent e)
    {
        System.out.println("componentShown:" + this.name());
    }

} // end of class ColorBar2DNode
