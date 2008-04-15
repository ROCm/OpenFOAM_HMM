/*
 *	PostApp.java
 *
 */
package FoamX.Post;

import java.applet.Applet;
import java.awt.BorderLayout;
//import java.awt.Color.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.Graphics;
import java.awt.GraphicsConfiguration;
import javax.swing.*;
import java.util.*;

import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.behaviors.vp.*;
import com.sun.j3d.utils.behaviors.mouse.*;
import com.sun.j3d.utils.behaviors.keyboard.*;
import com.sun.j3d.utils.picking.PickTool;
import com.sun.j3d.utils.universe.*;
import com.sun.j3d.utils.geometry.*;
import com.sun.j3d.utils.behaviors.keyboard.KeyNavigatorBehavior;
import com.sun.j3d.utils.behaviors.mouse.*;

import javax.media.j3d.*;
import javax.vecmath.*;
import javax.vecmath.Color3f;

import FoamXServer.CasePostServer.ICasePostServer;

import FoamX.Post.InteractiveNodes.*;
import FoamX.Post.Util.Colorizer;
import FoamX.Post.Util.KeyBehavior;
import FoamX.Post.Util.OrbitPickBehavior;
import FoamX.Post.Util.myOrbitBehavior;
import FoamX.Post.Util.myMouseScale;
//import FoamX.Post.Util.myMouseZoom;
//import FoamX.Post.Util.myViewOrbitBehavior;
import FoamX.Post.Util.PickDragBehavior;
//import FoamX.Post.Util.PrintManager;
import FoamX.Post.Util.Placement;
import FoamX.Post.Table.NodeTable;
import FoamX.Post.Util.OffScreenCanvas3D;

public class PostApp extends Applet implements MouseBehaviorCallback
{
    private static final int OFF_SCREEN_SCALE_X = 10;
    private static final int OFF_SCREEN_SCALE_Y = 10;

    ICasePostServer casePostServer_;
    String caseRoot_;
    String caseName_;

    NodeTable nodeTable_;
    PostWindow postWindow_;

    OrbitPickBehavior pickBeh_;
    KeyBehavior key_;
    //myOrbitBehavior orbit_;
    OrbitBehavior orbit_;
    myMouseScale ms_;

    InteractiveNodeList shapes_;
    BranchGroup scene_;
    SimpleUniverse u_ = null;
    Colorizer colorizer_;
//    PrintManager printManager_;
    Switch printSwitch_;

    Color bgColor_;
    Color3f bgColor3f_;
    Color fgColor_;
    Color3f fgColor3f_;

    private OffScreenCanvas3D offScreenCanvas3D_;

    public BranchGroup createSceneGraph(SimpleUniverse u)
    {

        // Create the root of the branch graph
        BranchGroup objRoot = new BranchGroup();

        // Switch for printer canvas rendering
        printSwitch_ = new Switch();
        printSwitch_.setCapability(Switch.ALLOW_SWITCH_WRITE);
        printSwitch_.setCapability(Switch.ALLOW_CHILDREN_EXTEND);
        printSwitch_.setCapability(Switch.ALLOW_CHILDREN_WRITE);
        printSwitch_.setWhichChild(Switch.CHILD_NONE);
        //printSwitch_.setCapability(Shape3D.ALLOW_BOUNDS_READ);
        objRoot.addChild(printSwitch_);


//        // Create a Transformgroup to scale all objects so they
//        // appear in the scene.
//        TransformGroup objScale = new TransformGroup();
//        Transform3D t3d = new Transform3D();
//        t3d.setScale(0.3);
//        objScale.setTransform(t3d);
//        objRoot.addChild(objScale);


//	    // Create the TransformGroup node and initialize it to the
//	    // identity. Enable the TRANSFORM_WRITE capability so that
//	    // our behavior code can modify it at run time. Add it to
//	    // the root of the subgraph.
//        TransformGroup objTrans = new TransformGroup();
//        objTrans.setCapability(TransformGroup.ENABLE_PICK_REPORTING);
//        objRoot.addChild(objTrans);

	// Axis
        AxisNode axis = new AxisNode("axis", new Color3f(1.0f, 1.0f, 0.0f));
        objRoot.addChild(axis.getRoot());
        shapes_.addElement(axis);
        shapes_.show("axis");

//        // Add colorBar
//        ColorBar2DNode colorBar = new ColorBar2DNode
//        (
//            "ColorBar_0",
//            getUniverse().getCanvas(),
//            new Placement(new Dimension(10, 10), Placement.TOPLEFT),
//            new Dimension(50, 128),
//            getColorizer(),
//            -1e-5f,
//            2e-5f
//        );
//        colorBar.show();
//        TransformGroup viewTrans =
//            getUniverse().getViewingPlatform()
//                .getViewPlatformTransform();
//        viewTrans.addChild(colorBar.getRoot());

//        // Patch0 featurelines
//        PrimitivePatch p0 = new PrimitivePatch(0);
//        FeatureLinesNode fl0 = new FeatureLinesNode("featureLines_0", p0);
//
//        objRoot.addChild(fl0.getRoot());
//        shapes_.addElement(fl0);
//
//
//        // Patch1 featurelines
//        PrimitivePatch p1 = new PrimitivePatch(1);
//        FeatureLinesNode fl1 = new FeatureLinesNode("featureLines_1", p1);
//
//        objRoot.addChild(fl1.getRoot());
//        shapes_.addElement(fl1);
//
//
//        // Patch0 faces
//        PatchFacesNode pf0 = new PatchFacesNode
//        (
//            "patchFaces_0",
//            p0,
//            new Color3f(0.8f, 0.4f, 0.0f),
//            0       // Math.PI / 4.5
//        );
//        PatchFaces pFaces = pf0.getPatchFaces();
//        pFaces.setAppearance(createMaterialAppearance());
//
//        objRoot.addChild(pf0.getRoot());
//        shapes_.addElement(pf0);
//
//
//        // Patch1 faces
//        PatchFacesNode pf1 = new PatchFacesNode
//        (
//            "patchFaces_1",
//            p1,
//            new Color3f(0.5f, 0.3f, 0.2f),
//            0       //Math.PI / 4.5
//        );
//        pFaces = pf1.getPatchFaces();
//        pFaces.setAppearance(createMaterialAppearance());
//
//        objRoot.addChild(pf1.getRoot());
//        shapes_.addElement(pf1);
//
//
//        // Patch0 points
//        PatchPointsNode pp0 = new PatchPointsNode
//        (
//            "patchPoints_0",
//            p0,
//            new Color3f(0.8f, 0.4f, 0.0f)
//        );
//        objRoot.addChild(pp0.getRoot());
//        shapes_.addElement(pp0);



        BoundingSphere bounds =
            new BoundingSphere(new Point3d(0.0,0.0,0.0), 100.0);

        // Set up the background
        BackgroundNode whiteBGNode =
            new BackgroundNode
            (
                "whiteBG",
                new Color3f(1.0f, 1.0f, 1.0f),
                bounds
            );
        objRoot.addChild(whiteBGNode.getRoot());
        shapes_.addElement(whiteBGNode);
//        Background bg = new Background(new Color3f(1.0f, 1.0f, 1.0f));
//        bg.setApplicationBounds(bounds);
//        objRoot.addChild(bg);

        // Set up the ambient light
        AmbientLightNode ambientLightNode =
            new AmbientLightNode
            (
                "ambient_0",
                new Color3f(1.0f, 1.0f, 1.0f),
                bounds
            );
        objRoot.addChild(ambientLightNode.getRoot());
        ambientLightNode.show();
        shapes_.addElement(ambientLightNode);

        // Set up the ambient lights
        DirectionalLightNode light0 =
            new DirectionalLightNode
            (
                "dirLight_0",
                new Color3f(1.0f, 1.0f, 1.0f),
                new Vector3f(1.0f, 1.0f, 1.0f),
                bounds
            );
        objRoot.addChild(light0.getRoot());
        //light0.show();
        shapes_.addElement(light0);
        
        DirectionalLightNode light1 =
            new DirectionalLightNode
            (
                "dirLight_1",
                new Color3f(1.0f, 1.0f, 1.0f),
                new Vector3f(-1.0f, -1.0f, -1.0f),
                bounds
            );
        objRoot.addChild(light1.getRoot());
        //light1.show();
        shapes_.addElement(light1);


        // Have Java 3D perform optimizations on this scene graph.
//        objRoot.compile();

        objRoot.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        objRoot.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        objRoot.setCapability(BranchGroup.ALLOW_BOUNDS_READ);

        return objRoot;
    }


    //--------------------------------------------------------------------------

    public PostApp
    (
        ICasePostServer casePostServer,
        String caseRoot,
        String caseName
    )
    {
        casePostServer_ = casePostServer;
        caseRoot_ = caseRoot;
        caseName_ = caseName;

        postWindow_ = new PostWindow
        (
            null,
            casePostServer_,
            this,
            caseRoot_,
            caseName_
        );
        MainFrame mf = new MainFrame(this, 300, 600);
    }

    //--------------------------------------------------------------------------

    public void init() {
        shapes_ = new InteractiveNodeList();

        colorizer_ = new Colorizer(this);

        bgColor_ = Color.black;
        bgColor3f_ = new Color3f(0.0f, 0.0f, 0.0f);
        fgColor_ = Color.white;
        fgColor3f_ = new Color3f(1.0f, 1.0f, 1.0f);

        nodeTable_ = new NodeTable(shapes_);
        nodeTable_.pack();
        nodeTable_.setVisible(true);

        postWindow_.setVisible(true);

        setLayout(new BorderLayout());
        GraphicsConfiguration config =
            SimpleUniverse.getPreferredConfiguration();

        Canvas3D c = new Canvas3D(config);
        add("Center", c);

        u_ = new SimpleUniverse(c);

        // Create a simple scene and attach it to the virtual universe
        scene_ = createSceneGraph(u_);

//        // PrintManager singleton
//        printManager_  = new PrintManager(this, printSwitch_);


        ViewingPlatform viewingPlatform = u_.getViewingPlatform();

        // add mouse behaviors to the ViewingPlatform
        View view = u_.getViewer().getView();

        //view.setProjectionPolicy(View.PERSPECTIVE_PROJECTION);
        //view.setProjectionPolicy(View.PARALLEL_PROJECTION);


        // This will move the ViewPlatform back a bit so the
        // objects in the scene can be viewed.
        viewingPlatform.setNominalViewingTransform();

        BoundingSphere bounds = (BoundingSphere)scene_.getBounds();
        System.out.println("bounds=" + bounds);

//        // Now create the simple picking behavior
//        pickBeh_ = new OrbitPickBehavior(c, scene_, bounds);
//        pickBeh_.setEnable(true);

        //scene_.compile();

        // Attach orbiting behaviour to mouse
        orbit_ = new OrbitBehavior
        (
            c,
//            myOrbitBehavior.DISABLE_ZOOM
//            OrbitBehavior.REVERSE_TRANSLATE
//          | OrbitBehavior.REVERSE_ZOOM
//          | OrbitBehavior.PROPORTIONAL_ZOOM
            OrbitBehavior.PROPORTIONAL_ZOOM
        );
        orbit_.setSchedulingBounds(bounds);
        orbit_.setRotationCenter(new Point3d(0.0, 0.0, 0.0));
        viewingPlatform.setViewPlatformBehavior(orbit_);


        TransformGroup viewTrans =
        viewingPlatform.getViewPlatformTransform();

        // Add keyboard navigation
        key_ = new KeyBehavior(this, viewTrans);
        key_.setSchedulingBounds(bounds);
        key_.setEnable(true);
        scene_.addChild(key_);

//        // Mouse scaling
//        ms_ = new myMouseScale(this, viewTrans);
//        ms_.setSchedulingBounds(bounds);
//        ms_.setEnable(true);
//        scene_.addChild(ms_);

//        Point3d eye = new Point3d(10.0, 0.0, 0.0);
//        Point3d center = new Point3d(0.0, 0.0, 0.0);
//        Vector3d up = new Vector3d(0.0, 0.0, 1.0);
        Transform3D look = new Transform3D();

//        look.ortho
//        (
//           -0.0,    // left
//            1.0,    // right
//           -0.0,    // bottom
//            1.0,    // top
//            0.1,    // near
//            1.1     // far
//        );
//        look.lookAt(eye, center, up);
//        look.invert();

        look.set(20.0f, new Vector3f(0.0f, 0.0f, 10.0f));
//        look.set(10.0f);

        //view.setFieldOfView(0.1);
        //System.out.println("View field:" + view.getFieldOfView());

        viewTrans.setTransform(look);

        setBounds(bounds);

    	u_.addBranchGraph(scene_);


        // Printing, Offscreen rendering.

	// Create the off-screen Canvas3D object
	offScreenCanvas3D_ = new OffScreenCanvas3D(config, true);
	// Set the off-screen size based on a scale factor times the
	// on-screen size
	Screen3D sOn = c.getScreen3D();
	Screen3D sOff = offScreenCanvas3D_.getScreen3D();
	Dimension dim = sOn.getSize();
	dim.width *= OFF_SCREEN_SCALE_X;
	dim.height *= OFF_SCREEN_SCALE_Y;
	sOff.setSize(dim);
	sOff.setPhysicalScreenWidth(sOn.getPhysicalScreenWidth() *
				    OFF_SCREEN_SCALE_X);
	sOff.setPhysicalScreenHeight(sOn.getPhysicalScreenHeight() *
				     OFF_SCREEN_SCALE_Y);
	// attach the offscreen canvas to the view
	u_.getViewer().getView().addCanvas3D(offScreenCanvas3D_);
    }

    //--------------------------------------------------------------------------

    public void destroy() {
        postWindow_.closeCase(PostWindow.KILL_SERVER);
        postWindow_.setVisible(false);


    	u_.removeAllLocales();

        casePostServer_ = null;

        nodeTable_.setVisible(false);
        nodeTable_ = null;

        postWindow_.setVisible(false);
        postWindow_ = null;

        pickBeh_ = null;
        shapes_ = null;
        scene_ = null;
        u_ = null;
        colorizer_ = null;
    }

    //--------------------------------------------------------------------------

    public InteractiveNodeList getShapes()
    {
        return shapes_;
    }

    //--------------------------------------------------------------------------

    public BranchGroup getScene()
    {
        return scene_;
    }

    //--------------------------------------------------------------------------

    public SimpleUniverse getUniverse()
    {
        return u_;
    }

    //--------------------------------------------------------------------------

    public Colorizer getColorizer()
    {
        return colorizer_;
    }

    //--------------------------------------------------------------------------

//    public PrintManager getPrintManager()
//    {
//        return printManager_;
//    }

    //--------------------------------------------------------------------------

    public OffScreenCanvas3D getOffScreenCanvas3D()
    {
        return offScreenCanvas3D_;
    }

    //--------------------------------------------------------------------------

    public Color getBackgroundColor()
    {
        return bgColor_;
    }

//    //--------------------------------------------------------------------------
//
//    public void setBackgroundColor(col)
//    {
//        bgColor_ = col;
//    }

    //--------------------------------------------------------------------------

    public Color3f getBackgroundColor3f()
    {
        return bgColor3f_;
    }

    //--------------------------------------------------------------------------
//
//    public void setBackgroundColor3f(col3f)
//    {
//        bgColor3f_ = col3f;
//    }
    //--------------------------------------------------------------------------

    public Color getForegroundColor()
    {
        return fgColor_;
    }

//    //--------------------------------------------------------------------------
//
//    public void setForegroundColor(col)
//    {
//        fgColor_ = col;
//    }

    //--------------------------------------------------------------------------

    public Color3f getForegroundColor3f()
    {
        return fgColor3f_;
    }

    //--------------------------------------------------------------------------
//
//    public void setForegroundColor3f(col3f)
//    {
//        fgColor3f_ = col3f;
//    }

    //--------------------------------------------------------------------------

    public PostWindow getPostWindow()
    {
        return postWindow_;
    }

    //--------------------------------------------------------------------------

    public NodeTable getNodeTable()
    {
        return nodeTable_;
    }

    //--------------------------------------------------------------------------

    static public Appearance createMaterialAppearance(){

        Appearance materialAppear = new Appearance();
        materialAppear.setCapability
        (
            Appearance.ALLOW_RENDERING_ATTRIBUTES_READ
        );

        PolygonAttributes polyAttrib = new PolygonAttributes();
        polyAttrib.setCullFace(PolygonAttributes.CULL_NONE);
        polyAttrib.setBackFaceNormalFlip(true); // bidirectional lighting
        materialAppear.setPolygonAttributes(polyAttrib);

        Material material = new Material();
        material.setDiffuseColor(new Color3f(1.0f, 0.0f, 0.0f));
        materialAppear.setMaterial(material);

        //TransparencyAttributes transp = new TransparencyAttributes();
        //transp.setTransparencyMode(TransparencyAttributes.FASTEST);
        //transp.setTransparency(0.5f);
        //materialAppear.setTransparencyAttributes(transp);

        return materialAppear;
    }

    //--------------------------------------------------------------------------

    static public Appearance createWireFrameAppearance(){

        Appearance materialAppear = new Appearance();
        PolygonAttributes polyAttrib = new PolygonAttributes();
        polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_LINE);
        materialAppear.setPolygonAttributes(polyAttrib);
        ColoringAttributes redColoring = new ColoringAttributes();
        redColoring.setColor(1.0f, 0.0f, 0.0f);
        materialAppear.setColoringAttributes(redColoring);

        return materialAppear;
    }

    //--------------------------------------------------------------------------

    static public Appearance createDotAppearance(){

        Appearance materialAppear = new Appearance();
        PolygonAttributes polyAttrib = new PolygonAttributes();
        polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_POINT);
        materialAppear.setPolygonAttributes(polyAttrib);
        ColoringAttributes redColoring = new ColoringAttributes();
        redColoring.setColor(1.0f, 0.0f, 0.0f);
        materialAppear.setColoringAttributes(redColoring);

        return materialAppear;
    }

    //--------------------------------------------------------------------------

    static public Transform3D viewAll(View view, BoundingSphere bounds)
    {
        Transform3D viewTrans = new Transform3D();

        Point3d center = new Point3d();
        bounds.getCenter(center);
        double radius = bounds.getRadius();

        double eyeDist = 1.2 * radius / Math.tan(view.getFieldOfView() / 2.0);

        Vector3d up = new Vector3d(0, 1, 0);
        Point3d eyePos = new Point3d(center);
        eyePos.z += eyeDist;


        System.out.println("eyePos : " + eyePos);
        System.out.println("center : " + center);
        System.out.println("up     : " + up);

        viewTrans.lookAt(eyePos, center, up);
        viewTrans.invert();

        if ( view.getBackClipDistance() < eyeDist )
        {
            view.setBackClipDistance(eyeDist);
            view.setFrontClipDistance(eyeDist / 3000);
        }

        return viewTrans;
    }

    //--------------------------------------------------------------------------

    static int getJava3DVersion(View view)
    {
        try
        {
            view.setTransparencySortingPolicy
            (
                View.TRANSPARENCY_SORT_GEOMETRY //Not in java3d-1.2
            );
            //AppearanceHelper.setJava3DVersion(AppearanceHelper.JAVA3D_1_3);

            System.out.println("You are running Java3D version 1.3x");
            return 3;
        }
        catch(NoSuchMethodError e )
        {
            System.out.println("You are running Java3D version 1.2x");
        }
        return 2;
    }


    //--------------------------------------------------------------------------

    static void printConfig(Canvas3D canvas)
    {
        System.out.println("\n Info about Canvas3D / Hardware \n");
        Map props = canvas.queryProperties();
        Iterator it = props.entrySet().iterator();
        while(it.hasNext())
        {
            Map.Entry entry = (Map.Entry) it.next();
            System.out.println(entry.getKey() + " = " + entry.getValue());
        }
        System.out.println("The End Info\n");
    }

    //--------------------------------------------------------------------------

    // Updates bounds on all behaviors.
    public void setBounds(BoundingSphere bounds)
    {
        System.out.println("setBounds : bounds=" + bounds);

        //pickBeh_.setBounds(bounds);
        key_.setSchedulingBounds(bounds);
        orbit_.setSchedulingBounds(bounds);
//        ms_.setSchedulingBounds(bounds);

        // Get data from bounds
        Point3d center = new Point3d();
        bounds.getCenter(center);
        double radius = bounds.getRadius();

        double front = 0;
        double back = 0;

        View view = u_.getViewer().getView();

        if (view.getProjectionPolicy() == View.PARALLEL_PROJECTION)
        {
            front = center.z - 2*radius;
            back = center.z + 2*radius;

            System.out.println("Updating clipplanes : front:" + front + "  back:" + back);
            view.setFrontClipDistance(front);
            view.setBackClipDistance(back);
        }
        else
        {
            double eyeDist = 1.2 * radius / Math.tan(view.getFieldOfView() / 2.0);
            System.out.println("eyeDist:" + eyeDist);
            //??

//        // Get inverse of current view transform
//        Transform3D viewTransform = new Transform3D();
//        TransformGroup viewTransGroup =
//        u_.getViewingPlatform().getViewPlatformTransform();
//        viewTransGroup.getTransform(viewTransform);
//        
//        Transform3D inverseViewTransform = new Transform3D();
//        inverseViewTransform.invert(viewTransform);
//
//        // Invert center
//        Point3d result = new Point3d();
//        inverseViewTransform.transform(center, result);
//
//        // Distance from screen to center.
//        double z = -result.z;
//
//        front = z - radius;
//        back = z + radius;
//
//        System.out.println("z:" + z + "  front:" + front + "  back:" + back);

        }

        front = view.getFrontClipDistance();
        back = view.getBackClipDistance();
        System.out.println("Front clip:" + front + " Back clip:" + back);
    }

    //--------------------------------------------------------------------------
    //---- InteractiveNode Interface
    //--------------------------------------------------------------------------
    public void transformChanged(int type, Transform3D zoomTransform)
    {
//        Vector3d translation = new Vector3d();
//        zoomTransform.get(translation);
//
//        double zTrans = translation.z;
//        System.out.println("zTrans:" + zTrans);
//
//        double halfFov = Math.atan(zTrans);
//        System.out.println("Fov   :" + 2*halfFov);
//
//        View view = u_.getViewer().getView();
//        view.setFieldOfView(2 * halfFov);
//
//        setBounds((BoundingSphere)scene_.getBounds());
//
//        System.out.println();

        System.out.println("zoomTransform:" + zoomTransform);
    }

}



