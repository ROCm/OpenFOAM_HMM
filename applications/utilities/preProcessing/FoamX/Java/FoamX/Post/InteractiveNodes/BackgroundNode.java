/*
 *      BackgroundNode.java
 */
package FoamX.Post.InteractiveNodes;

//import java.awt.*;
//import java.awt.Font;

//import com.sun.j3d.utils.applet.MainFrame;
//import com.sun.j3d.utils.geometry.Text2D;
import com.sun.j3d.utils.picking.PickResult;


import javax.media.j3d.*;
import javax.vecmath.*;

public class BackgroundNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    Background bg_;
    boolean isShowing_;

    public BackgroundNode
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

        // Create background
        bg_ = new Background(color);
        bg_.setApplicationBounds(bounds);

        switch_.addChild(bg_);

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
        System.out.println("ERROR: BackgroundNode Picked:");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("ERROR: BackGroundNode unPicked:");
    }

    //--------------------------------------------------------------------------

    // Include in display
    public void show()
    {
        //System.out.println("BackgroundNode show:" + this.name());
        switch_.setWhichChild(Switch.CHILD_ALL);
        isShowing_ = true;
    }

    //--------------------------------------------------------------------------

    // Exclude from display
    public void hide()
    {
        //System.out.println("BackgroundNode hide:" + this.name());
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
        //System.out.println("BackgroundNode remove:" + this.name());
        hide();
        bg_ = null;
        switch_ = null;
    }
} // end of class BackgroundNode
