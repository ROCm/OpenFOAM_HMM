/*
 */
package FoamX.Post.Shapes;

//import com.sun.j3d.utils.geometry.GeometryInfo;
//import com.sun.j3d.utils.geometry.NormalGenerator;
//import com.sun.j3d.utils.geometry.Stripifier;
//import com.sun.j3d.utils.geometry.Triangulator;
//import com.sun.j3d.utils.picking.PickIntersection;
import com.sun.j3d.utils.picking.PickResult;
//import com.sun.j3d.utils.picking.PickTool;

import javax.media.j3d.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import FoamX.Post.InteractiveNodes.ColorBarNode;
import FoamX.Post.Util.Colorizer;

public class ColorBar extends Shape3D implements PickableShape
{
    // Create coloured triangles.
    public ColorBar
    (
        final ColorBarNode node,
        final Colorizer ci,
        final float min,
        final float max
    )
    {
        this.setGeometry(ci.colorBar(min, max));
        this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
    }

    //--------------------------------------------------------------------------
    //---- PickableShape Interface
    //--------------------------------------------------------------------------

    // Callback for if picked.
    public void picked(PickResult pickResult)
    {
        System.out.println("\n\nPicked on ColorBar ");
    }

    //--------------------------------------------------------------------------

    // Callback for if unpicked.
    public void unPicked(PickResult pickResult)
    {
        System.out.println("\n\nUnpicked on ColorBar ");
    }

    //--------------------------------------------------------------------------

    // Completely remove
    public void remove()
    {
        System.out.println("remove of ColorBar");

        if (this.numGeometries() != 0)
        {
            this.removeGeometry(0);
        }
    }

    //--------------------------------------------------------------------------
}
