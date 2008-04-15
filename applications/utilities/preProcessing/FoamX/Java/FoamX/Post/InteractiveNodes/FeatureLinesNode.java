/*
 *
 */

package FoamX.Post.InteractiveNodes;

//import java.awt.*;
//import java.awt.Font;

//import com.sun.j3d.utils.applet.MainFrame;
//import com.sun.j3d.utils.geometry.Text2D;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.Shapes.FeatureLines;
import FoamX.Post.Mesh.PrimitivePatch;

public class FeatureLinesNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    FeatureLines featureLines_;
    boolean isShowing_;

    public FeatureLinesNode(String myName, PrimitivePatch patch)
    {
        myName_ = myName;

        isShowing_ = false;

        switch_ = new Switch();
        switch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        switch_.setCapability(Node.ALLOW_BOUNDS_READ);

        branchGroup_ = new BranchGroup();
        branchGroup_.setCapability(Node.ALLOW_BOUNDS_READ);

        featureLines_ = new FeatureLines(this, patch);

        switch_.addChild(featureLines_);

        branchGroup_.addChild(switch_);
    }

    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    public FeatureLines getFeatureLines()
    {
        return featureLines_;
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
        featureLines_.picked(pickResult);
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        featureLines_.unPicked(pickResult);
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
        featureLines_.remove();
        featureLines_ = null;
        switch_ = null;
    }

} // end of class FeatureLinesNode
