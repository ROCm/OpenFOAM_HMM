/*
 */

package FoamX.Post.InteractiveNodes;

import com.sun.j3d.utils.picking.PickResult;

import javax.media.j3d.*;
import javax.vecmath.Color3f;

import FoamX.Post.Shapes.PatchFaces;
import FoamX.Post.Mesh.PrimitivePatch;


public class PatchFacesNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    PatchFaces patchFaces_;
    boolean isShowing_;


    public PatchFacesNode
    (
        String myName,
        PrimitivePatch patch,
        Color3f color,
        double creaseAngle
    )
    {
        myName_ = myName;

        isShowing_ = false;

        switch_ = new Switch();
        switch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        switch_.setCapability(Node.ALLOW_BOUNDS_READ);

        branchGroup_ = new BranchGroup();
        branchGroup_.setCapability(Node.ALLOW_BOUNDS_READ);

        patchFaces_ = new PatchFaces
        (
            this,
            patch,
            color,
            creaseAngle
        );

        switch_.addChild(patchFaces_);

        branchGroup_.addChild(switch_);
    }


    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    public PatchFaces getPatchFaces()
    {
        return patchFaces_;
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
        patchFaces_.picked(pickResult);
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        patchFaces_.unPicked(pickResult);
    }

    //--------------------------------------------------------------------------

    // Include in display
    public void show()
    {
        //System.out.println("Node show:" + this.name());
        switch_.setWhichChild(Switch.CHILD_ALL);
        isShowing_ = true;
    }

    //--------------------------------------------------------------------------

    // Exclude from display
    public void hide()
    {
        //System.out.println("Node hide:" + this.name());
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
        //System.out.println("Node remove:" + this.name());
        hide();
        patchFaces_.remove();
        patchFaces_ = null;
        switch_ = null;
    }

    //--------------------------------------------------------------------------
}
