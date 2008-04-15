/*
 *      LightNode.java
 */
package FoamX.Post.InteractiveNodes;

import java.awt.*;
import java.awt.Font;

import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.geometry.Text2D;
import com.sun.j3d.utils.picking.PickResult;


import javax.media.j3d.*;
import javax.vecmath.*;

import FoamX.Post.Util.Text;

public class AmbientLightNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    AmbientLight light_;
    boolean isShowing_;

    public AmbientLightNode
    (
        final String myName,
        final Color3f color,
        final BoundingSphere bounds
    )
    {
        myName_ = myName;

        isShowing_ = false;

        switch_ = new Switch();
        switch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        switch_.setCapability(Node.ALLOW_BOUNDS_READ);
        
        branchGroup_ = new BranchGroup();
        branchGroup_.setCapability(Node.ALLOW_BOUNDS_READ);

        // Create light
        light_ = new AmbientLight(color);
        light_.setInfluencingBounds(bounds);

        switch_.addChild(light_);

        branchGroup_.addChild(switch_);
    }

    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    //--------------------------------------------------------------------------
    //---- InteractiveNode Interface
    //--------------------------------------------------------------------------


    public String name()
    {
        return myName_;
    }


    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("ERROR: LightNode Picked:");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("ERROR: LightNode unPicked:");
    }

    //--------------------------------------------------------------------------

    // Include in display
    public void show()
    {
        //System.out.println("LightNode show:" + this.name());
        switch_.setWhichChild(Switch.CHILD_ALL);
        isShowing_ = true;
    }

    //--------------------------------------------------------------------------

    // Exclude from display
    public void hide()
    {
        //System.out.println("LightNode hide:" + this.name());
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
        //System.out.println("LightNode remove:" + this.name());

        hide();
        light_ = null;
        switch_ = null;
    }
} // end of class LightNode
