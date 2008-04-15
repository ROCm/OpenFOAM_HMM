/*
 *
 */

package FoamX.Post.InteractiveNodes;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.Util.Colorizer;
import FoamX.Post.Shapes.ColorBar;

public class ColorBarNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    ColorBar colorBar_;
    boolean isShowing_;

    public ColorBarNode
    (
        final String myName,
        final Colorizer ci,
        final float min,
        final float max
    )
    {
        myName_ = myName;

        isShowing_ = false;

        switch_ = new Switch();
        switch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        switch_.setCapability(Node.ALLOW_BOUNDS_READ);

        branchGroup_ = new BranchGroup();
        branchGroup_.setCapability(Node.ALLOW_BOUNDS_READ);

        colorBar_ = new ColorBar(this, ci, min, max);

        switch_.addChild(colorBar_);

        branchGroup_.addChild(switch_);
    }

    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    public ColorBar getColorBar()
    {
        return colorBar_;
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
        colorBar_.picked(pickResult);
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        colorBar_.unPicked(pickResult);
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
        colorBar_.remove();
        colorBar_ = null;
        switch_ = null;
    }

} // end of class ColorBarNode
