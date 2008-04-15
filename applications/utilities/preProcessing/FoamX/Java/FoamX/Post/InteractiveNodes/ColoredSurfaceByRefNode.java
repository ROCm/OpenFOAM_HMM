/*
 *
 */

package FoamX.Post.InteractiveNodes;


import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickIntersection;

import FoamX.Post.Shapes.ColoredSurfaceByRef;
//import FoamX.Post.Mesh.PrimitivePatch;
import FoamX.Post.Util.Colorizer;

public class ColoredSurfaceByRefNode implements InteractiveNode
{
    String myName_;
    Switch switch_;
    BranchGroup branchGroup_;
    ColoredSurfaceByRef coloredSurfaceByRef_;
    boolean isShowing_;

    //- From unstripped file
    public ColoredSurfaceByRefNode
    (
        final String myName,
        final double creaseAngle,
        final boolean stripify,
        final String pointsFileName,
        final String facesFileName,
        final String valuesFileName,
        final Colorizer colorizer,
        final float min,
        final float max
    )
    {
        this
        (
            myName, false, false, stripify, creaseAngle,
            pointsFileName, null, facesFileName, valuesFileName,
            colorizer, min, max
        );
    }

    //- From file stripped or not stripped
    public ColoredSurfaceByRefNode
    (
        final String myName,
        final boolean readTriStrip,
        final boolean readNormals,
        final boolean stripify,
        final double creaseAngle,
        final String pointsFileName,
        final String normalsFileName,
        final String facesFileName,
        final String valuesFileName,
        final Colorizer colorizer,
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

        coloredSurfaceByRef_ =
            new ColoredSurfaceByRef
            (
                this,
                readTriStrip,
                readNormals,
                stripify,
                creaseAngle,
                pointsFileName,
                normalsFileName,
                facesFileName,
                valuesFileName,
                colorizer,
                min,
                max
            );

        switch_.addChild(coloredSurfaceByRef_);

        branchGroup_.addChild(switch_);
    }


    // From values
    public ColoredSurfaceByRefNode
    (
        final String myName,
        final boolean readTriStrip,
        final boolean stripify,
        final double creaseAngle,
        final float[] points,
        final float[] normals,
        final int[] triFaces,
        final float[] values,
        final Colorizer colorizer,
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

        coloredSurfaceByRef_ =
            new ColoredSurfaceByRef
            (
                this,
                readTriStrip,
                stripify,
                creaseAngle,
                points,
                normals,
                triFaces,
                null,
                values,
                colorizer,
                min,
                max
            );

        switch_.addChild(coloredSurfaceByRef_);

        branchGroup_.addChild(switch_);
    }

    public BranchGroup getRoot()
    {
        return branchGroup_;
    }

    public ColoredSurfaceByRef getColoredSurfaceByRef()
    {
        return coloredSurfaceByRef_;
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
        coloredSurfaceByRef_.picked(pickResult);
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        coloredSurfaceByRef_.unPicked(pickResult);
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
        coloredSurfaceByRef_.remove();
        coloredSurfaceByRef_ = null;
        switch_ = null;
    }

} // end of class ColoredSurfaceByRefNode
