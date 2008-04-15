/*
 *	@(#)myOrbitBehavior.java 1.9 01/02/28 16:26:56
 *
 * Copyright (c) 1996-2001 Sun Microsystems, Inc. All Rights Reserved.
 *
 * Sun grants you ("Licensee") a non-exclusive, royalty free, license to use,
 * modify and redistribute this software in source and binary code form,
 * provided that i) this copyright notice and license appear on all copies of
 * the software; and ii) Licensee does not utilize the software in a manner
 * which is disparaging to Sun.
 *
 * This software is provided "AS IS," without a warranty of any kind. ALL
 * EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY
 * IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR
 * NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL NOT BE
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
 * OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER
 * CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE USE OF
 * OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line control of
 * aircraft, air traffic, aircraft navigation or aircraft communications; or in
 * the design, construction, operation or maintenance of any nuclear
 * facility. Licensee represents and warrants that it will not use or
 * redistribute the Software for such purposes.
 */

package FoamX.Post.Util;

import com.sun.j3d.utils.behaviors.vp.*;

import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.KeyEvent;
import java.awt.AWTEvent;
import java.awt.Component;
import java.awt.Cursor;
import javax.swing.SwingUtilities;

import javax.media.j3d.WakeupOnAWTEvent;
import javax.media.j3d.WakeupOnElapsedFrames;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.Transform3D;
import javax.media.j3d.View;
import javax.media.j3d.Canvas3D;

import javax.vecmath.Vector3d;
import javax.vecmath.Point3d;
import javax.vecmath.Matrix3d;

import com.sun.j3d.utils.universe.ViewingPlatform;

import com.sun.j3d.internal.J3dUtilsI18N;

/**
 * Moves the View around a point of interest when the mouse is dragged with
 * a mouse button pressed.  Includes rotation, zoom, and translation
 * actions.
 * <p>
 * This behavior must be added to the ViewingPlatform
 * using the <code>ViewingPlatform.addViewPlatformBehavior()</code> method.
 * <p>
 * The rotate action rotates the ViewPlatform around the point of interest
 * when the mouse is moved with the main mouse button pressed.  The
 * rotation is in the direction of the mouse movement, with a default
 * rotation of 0.01 radians for each pixel of mouse movement.
 * <p>
 * The zoom action moves the ViewPlatform closer to or further from the
 * point of interest when the mouse is moved with the middle mouse button
 * pressed (or Alt-main mouse button on systems without a middle mouse button).
 * The default zoom action is to translate the ViewPlatform 0.01 units for each
 * pixel of mouse movement.  Moving the mouse up moves the ViewPlatform closer,
 * moving the mouse down moves the ViewPlatform further away.
 * <p>
 * By default, the zoom action allows the ViewPlatform to move through
 * the center of rotation to orbit at a negative radius.
 * The <code>STOP_ZOOM</code> constructor flag will stop the ViewPlatform at
 * a minimum radius from the center.  The default minimum radius is 0.0
 * and can be set using the <code>setMinRadius</code> method.
 * <p>
 * The <code>PROPORTIONAL_ZOOM</code> constructor flag changes the zoom action
 * to move the ViewPlatform proportional to its distance from the center
 * of rotation.  For this mode, the default action is to move the ViewPlatform
 * by 1% of its distance from the center of rotation for each pixel of
 * mouse movement.
 * <p>
 * The translate action translates the ViewPlatform when the mouse is moved
 * with the right mouse button pressed (Shift-main mouse button on systems
 * without a right mouse button).  The translation is in the direction of the
 * mouse movement, with a default translation of 0.01 units for each pixel
 * of mouse movement.
 * <p>
 * The sensitivity of the actions can be scaled using the
 * <code>set</code><i>Action</i><code>Factor()</code> methods which scale
 * the default movement by the factor. The rotate and translate actions
 * have separate factors for x and y.
 * <p>
 * The actions can be reversed using the <code>REVERSE_</code><i>ACTION</i>
 * constructor flags.  The default action moves the ViewPlatform around the
 * objects in the scene.  The <code>REVERSE_</code><i>ACTION</i> flags can
 * make the objects in the scene appear to be moving in the direction
 * of the mouse movement.
 * <p>
 * The actions can be disabled by either using the
 * <code>DISABLE_</code><i>ACTION</i> constructor flags or the
 * <code>set</code><i>Action</i><code>Enable</code> methods.
 * <p>
 * The default center of rotation is (0, 0, 0) and can be set using the
 * <code>setRotationCenter()</code> method.
 *
 * @since Java 3D 1.2.1
 */
public class myOrbitBehavior extends ViewPlatformAWTBehavior {

    private Transform3D velocityTransform;
    private Transform3D longditudeTransform;
    private Transform3D rollTransform;
    private Transform3D latitudeTransform;
    private Transform3D rotateTransform = new Transform3D();

    // needed for integrateTransforms but don't want to new every time
    private Transform3D temp1 = new Transform3D();
    private Transform3D temp2 = new Transform3D();
    private Transform3D translation = new Transform3D();
    private Vector3d transVector = new Vector3d();
    private Vector3d distanceVector = new Vector3d();
    private Vector3d centerVector = new Vector3d();
    private Vector3d invertCenterVector = new Vector3d();

    private double longditude = 0.0;
    private double latitude = 0.0;
    private double rollAngle = 0.0;
    private double distanceFromCenter = 20.0;

    private final double MAX_MOUSE_ANGLE = Math.toRadians( 3 );
    private final double ZOOM_FACTOR = 1.0;
    private Point3d rotationCenter = new Point3d();
    private Vector3d centerToView = new Vector3d();
    private Matrix3d rotMatrix = new Matrix3d();
    private Vector3d oldPos = new Vector3d();
    private Vector3d currentPos = new Vector3d();
    
    private int mouseX = 0;
    private int mouseY = 0;

    private Canvas3D canvas;

    private double rotXMul;
    private double rotYMul;
    private double transXMul;
    private double transYMul;
    private double zoomMul;

    private double rotXFactor = 1.0;
    private double rotYFactor = 1.0;
    private double transXFactor = 1.0;
    private double transYFactor = 1.0;
    private double zoomFactor = 1.0;

    private double xtrans = 0.0;
    private double ytrans = 0.0;

    private boolean zoomEnabled = true;
    private boolean rotateEnabled = true;
    private boolean translateEnabled = true;
    private boolean reverseRotate = false;
    private boolean reverseTrans = false;
    private boolean reverseZoom = false;
    private boolean stopZoom = false;
    private boolean proportionalZoom = false;
    private double minRadius = 0.0;
    private int leftButton = ROTATE;
    private int rightButton = TRANSLATE;
    private int middleButton = ZOOM;
    
    /**
     * Constructor flag to reverse the rotate behavior
     */
    public static final int REVERSE_ROTATE = 0x010;

    /**
     * Constructor flag to reverse the translate behavior
     */
    public static final int REVERSE_TRANSLATE = 0x020;

    /**
     * Constructor flag to reverse the zoom behavior
     */
    public static final int REVERSE_ZOOM = 0x040;

    /**
     * Constructor flag to reverse all the behaviors
     */
    public static final int REVERSE_ALL = (REVERSE_ROTATE | REVERSE_TRANSLATE |
					   REVERSE_ZOOM);

    /**
     * Constructor flag that indicates zoom should stop when it reaches
     * the minimum orbit radius set by setMinRadius().  The minimus
     * radius default is 0.0.
     */
    public static final int STOP_ZOOM = 0x100;
    
    /**
     * Constructor flag to disable rotate
     */
    public static final int DISABLE_ROTATE = 0x200;

    /**
     * Constructor flag to disable translate
     */
    public static final int DISABLE_TRANSLATE = 0x400;

    /**
     * Constructor flag to disable zoom
     */
    public static final int DISABLE_ZOOM = 0x800;

    /**
     * Constructor flag to use proportional zoom, which determines
     * how much you zoom based on view's distance from the center of
     * rotation.  The percentage of distance that the viewer zooms
     * is determined by the zoom factor.
     */
    public static final int PROPORTIONAL_ZOOM = 0x1000;
    
    /**
     * Used to set the fuction for a mouse button to Rotate 
     */
    private static final int ROTATE = 0;

    /**
     * Used to set the function for a mouse button to Translate
     */
    private static final int TRANSLATE = 1;

    /**
     * Used to set the function for a mouse button to Zoom
     */
    private static final int ZOOM = 2;

    private static final double NOMINAL_ZOOM_FACTOR = .01;
    private static final double NOMINAL_PZOOM_FACTOR = 1.0;
    private static final double NOMINAL_ROT_FACTOR = .01;
    private static final double NOMINAL_TRANS_FACTOR = .01;

    /**
     * Creates a new myOrbitBehavior
     *
     * @param c The Canvas3D to add the behavior to
     */
    public myOrbitBehavior(Canvas3D c) { 
	this(c, 0 );
    }

    /**
     * Creates a new myOrbitBehavior
     *
     * @param c The Canvas3D to add the behavior to
     * @param flags The option flags
     */
    public myOrbitBehavior(Canvas3D c, int flags) {
	super(c, MOUSE_LISTENER | MOUSE_MOTION_LISTENER | flags );
	canvas = c;
        
	if ((flags & DISABLE_ROTATE) != 0) rotateEnabled = false;
	if ((flags & DISABLE_ZOOM) != 0) zoomEnabled = false;
	if ((flags & DISABLE_TRANSLATE) != 0) translateEnabled = false;
        if ((flags & REVERSE_TRANSLATE) != 0) reverseTrans = true;	
        if ((flags & REVERSE_ROTATE) != 0)  reverseRotate = true;
        if ((flags & REVERSE_ZOOM) != 0) reverseZoom = true;
        if ((flags & STOP_ZOOM) != 0) stopZoom = true;
	if ((flags & PROPORTIONAL_ZOOM) !=0) {
	    proportionalZoom = true;
	}

	rotXMul = NOMINAL_ROT_FACTOR * rotXFactor;
	rotYMul = NOMINAL_ROT_FACTOR * rotYFactor;
	transXMul = NOMINAL_TRANS_FACTOR * transXFactor;
	transYMul = NOMINAL_TRANS_FACTOR * transYFactor;
	if (proportionalZoom) {
	    zoomMul = NOMINAL_PZOOM_FACTOR * zoomFactor;
	}
	else {
	    zoomMul = NOMINAL_ZOOM_FACTOR * zoomFactor;
	}
        
        rollTransform = new Transform3D();
        latitudeTransform = new Transform3D();
        longditudeTransform = new Transform3D();
        targetTransform = new Transform3D();
        velocityTransform = new Transform3D();
    }

    protected synchronized void processAWTEvents( final java.awt.AWTEvent[] events ) {
        motion = false;
        for(int i=0; i<events.length; i++)
            if (events[i] instanceof MouseEvent) 
                processMouseEvent( (MouseEvent)events[i] );
    }

    protected void processMouseEvent( final MouseEvent evt ) {
        
        if (evt.getID()==MouseEvent.MOUSE_PRESSED) {
            mouseX = evt.getX();
            mouseY = evt.getY();
            motion=true;
        } else if (evt.getID()==MouseEvent.MOUSE_DRAGGED) {	    
            int xchange = evt.getX() - mouseX;
            int ychange = evt.getY() - mouseY;
	    // rotate
	    if (rotate(evt)) {
		if (reverseRotate) {
		    longditude -= xchange * rotXMul;
		    latitude -= ychange * rotYMul;
		}
		else {
		    longditude += xchange * rotXMul;
		    latitude += ychange * rotYMul;
		}
	    }
	    // translate
	    else if (translate(evt)) {
		if (reverseTrans) {
		    xtrans -= xchange * transXMul;
		    ytrans += ychange * transYMul;
		}
		else {
		    xtrans += xchange * transXMul;
		    ytrans -= ychange * transYMul;
		}
            }
	    // zoom
	    else if (zoom(evt)) {
		if (proportionalZoom) {
		    if (reverseZoom) {
			if ((distanceFromCenter -
			     (zoomMul*ychange*distanceFromCenter/100.0)) >
			    minRadius) {
			    distanceFromCenter -= (zoomMul*ychange*
						   distanceFromCenter/100.0);
			}
			else {
			    distanceFromCenter = minRadius;
			}			    
		    }
		    else {
			if ((distanceFromCenter +
			     (zoomMul*ychange*distanceFromCenter/100.0))
			     > minRadius) {
			    distanceFromCenter += (zoomMul*ychange*
						   distanceFromCenter/100.0);
			}
			else {
			    distanceFromCenter = minRadius;
			}
		    }
		}
		else {
		    if (stopZoom) {
			if (reverseZoom) {
			    if ((distanceFromCenter - ychange*zoomMul) > minRadius) {
				distanceFromCenter -= ychange*zoomMul;
			    }
			    else {
				distanceFromCenter = minRadius;
			    }
			}
			else {
			    if ((distanceFromCenter + ychange*zoomMul) > minRadius) {
				distanceFromCenter += ychange * zoomMul;
			    }
			    else {
				distanceFromCenter = minRadius;
			    }
			}
		    }
		    else {
			if (reverseZoom) {
			    distanceFromCenter -= ychange*zoomMul;
			}
			else {
			    distanceFromCenter += ychange*zoomMul;
			}
		    }
		}
	    }
            mouseX = evt.getX();
            mouseY = evt.getY();
	    motion = true;
	} else if (evt.getID()==MouseEvent.MOUSE_RELEASED ) {
        }
        
    }
    
    /**
     * Sets the ViewingPlatform for this behavior.  This method is
     * called by the ViewingPlatform.
     * If a sub-calls overrides this method, it must call
     * super.setViewingPlatform(vp).
     * NOTE: Applications should <i>not</i> call this method.    
     */
    public void setViewingPlatform(ViewingPlatform vp) {
        super.setViewingPlatform( vp );
	if (vp!=null) {
            resetViewPosition();
	    integrateTransforms();
	}
    }
    
    /**
     * Reset the orientation and distance of this behavior to the current
     * values in the ViewPlatform Transform Group
     */
    private void resetViewPosition() {
        targetTG.getTransform( targetTransform );

        targetTransform.get( rotMatrix, transVector );
        centerToView.sub( transVector, rotationCenter );
        distanceFromCenter = centerToView.length();
        
        System.out.println
        (
            "resetViewPosition : distanceFromCenter:"
          + distanceFromCenter
        );

        targetTransform.get( rotMatrix );
        rotateTransform.set( rotMatrix );
    }
    
    protected synchronized void integrateTransforms() {
            targetTransform.get( oldPos );
            
            targetTG.getTransform( targetTransform );
            
            targetTransform.get( currentPos );
            
            // Check if the transform has been changed by another
            // behavior
            if ( !currentPos.equals(oldPos))
                resetViewPosition();
            
	    longditudeTransform.rotY( longditude );            
	    latitudeTransform.rotX( latitude );

	    rotateTransform.mul(rotateTransform, latitudeTransform);
	    rotateTransform.mul(rotateTransform, longditudeTransform);
//	    rotateTransform.mul(longditudeTransform, rotateTransform);

//	    distanceVector.z = distanceFromCenter;
//  	    temp1.set(distanceVector);

  	    temp1.set(distanceFromCenter);
 	    temp1.mul(rotateTransform, temp1);
	    
	    // want to look at rotationCenter
	    transVector.x = rotationCenter.x + xtrans;
	    transVector.y = rotationCenter.y + ytrans;
	    transVector.z = rotationCenter.z;
	    translation.set(transVector);
	    targetTransform.mul(temp1, translation);

	    // handle rotationCenter
	    temp1.set(centerVector);
	    temp1.mul(targetTransform);
           
            invertCenterVector.x = -centerVector.x;
            invertCenterVector.y = -centerVector.y;
            invertCenterVector.z = -centerVector.z;
            
	    temp2.set(invertCenterVector);
	    targetTransform.mul(temp1, temp2);

            // renormalize for security
            double scale = targetTransform.getScale();
            targetTransform.normalizeCP();
            targetTransform.set(scale);

	    targetTG.setTransform(targetTransform);

	    // reset yaw and pitch angles
	    longditude = 0.0;
	    latitude = 0.0;        
    }

    /**
     * Sets the center around which the View rotates.
     * The default is (0,0,0).
     * @param center The Point3d to set the center of rotation to
     */
    public synchronized void setRotationCenter(Point3d center) {
	rotationCenter.x = center.x;
	rotationCenter.y = center.y;
	rotationCenter.z = center.z;
	centerVector.set(rotationCenter);
    }

    /**
     * Places the value of the center around which the View rotates
     * into the Point3d.
     * @param center The Point3d
     */
    public void getRotationCenter(Point3d center) {
	center.x = rotationCenter.x;
	center.y = rotationCenter.y;
	center.z = rotationCenter.z;
    }

    // TODO
    // Need to add key factors for Rotate, Translate and Zoom
    // Method calls should just update MAX_KEY_ANGLE, KEY_TRANSLATE and
    // KEY_ZOOM
    //
    // Methods also need to correctly set sign of variables depending on
    // the Reverse settings.
    
    /**
     * Sets the rotation x and y factors.  The factors are used to determine
     * how many radians to rotate the view for each pixel of mouse movement.
     * The view is rotated factor * 0.01 radians for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param xfactor The x movement multiplier
     * @param yfactor The y movement multiplier
     **/   
    public synchronized void setRotFactors(double xfactor, double yfactor) {
	rotXFactor = xfactor;
	rotYFactor = yfactor;
	rotXMul = NOMINAL_ROT_FACTOR * xfactor;
	rotYMul = NOMINAL_ROT_FACTOR * yfactor;
    }

    /**
     * Sets the rotation x factor.  The factors are used to determine
     * how many radians to rotate the view for each pixel of mouse movement.
     * The view is rotated factor * 0.01 radians for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param xfactor The x movement multiplier
     **/   
    public synchronized void setRotXFactor(double xfactor) {
	rotXFactor = xfactor;
	rotXMul = NOMINAL_ROT_FACTOR * xfactor;
    }

    /**
     * Sets the rotation y factor.  The factors are used to determine
     * how many radians to rotate the view for each pixel of mouse movement.
     * The view is rotated factor * 0.01 radians for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param yfactor The y movement multiplier
     **/   
    public synchronized void setRotYFactor(double yfactor) {
	rotYFactor = yfactor;
	rotYMul = NOMINAL_ROT_FACTOR * yfactor;
    }

    /**
     * Sets the translation x and y factors.  The factors are used to determine
     * how many units to translate the view for each pixel of mouse movement.
     * The view is translated factor * 0.01 units for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param xfactor The x movement multiplier
     * @param yfactor The y movement multiplier
     **/   
    public synchronized void setTransFactors(double xfactor,
						 double yfactor) {
	transXFactor = xfactor;
	transYFactor = yfactor;
	transXMul = NOMINAL_TRANS_FACTOR * xfactor;
	transYMul = NOMINAL_TRANS_FACTOR * yfactor;
    }

    /**
     * Sets the translation x factor.  The factors are used to determine
     * how many units to translate the view for each pixel of mouse movement.
     * The view is translated factor * 0.01 units for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param xfactor The x movement multiplier
     **/   
    public synchronized void setTransXFactor(double xfactor) {
	transXFactor = xfactor;
	transXMul = NOMINAL_TRANS_FACTOR * xfactor;
    }

    /**
     * Sets the translation y factor.  The factors are used to determine
     * how many units to translate the view for each pixel of mouse movement.
     * The view is translated factor * 0.01 units for each pixel of mouse
     * movement.  The default factor is 1.0.
     * @param yfactor The y movement multiplier
     **/   
    public synchronized void setTransYFactor(double yfactor) {
	transYFactor = yfactor;
	transYMul = NOMINAL_TRANS_FACTOR * yfactor;
    }

    /**
     * Sets the zoom factor.  The factor is used to determine how many
     * units to zoom the view for each pixel of mouse movement.
     * The view is zoomed factor * 0.01 units for each pixel of mouse
     * movement.  For proportional zoom, the view is zoomed factor * 1%
     * of the distance from the center of rotation for each pixel of
     * mouse movement.  The default factor is 1.0.
     * @param zfactor The movement multiplier
     */
    public synchronized void setZoomFactor(double zfactor) {
	zoomFactor = zfactor;
	if (proportionalZoom) {
	    zoomMul = NOMINAL_PZOOM_FACTOR * zfactor;
	}
	else {
	    zoomMul = NOMINAL_ZOOM_FACTOR * zfactor;
	}
    }

    /**
     * Returns the x rotation movement multiplier
     * @return The movement multiplier for x rotation
     */
    public double getRotXFactor() {
	return rotXFactor;
    }

    /**
     * Returns the y rotation movement multiplier
     * @return The movement multiplier for y rotation
     */
    public double getRotYFactor() {
	return rotYFactor;
    }

    /**
     * Returns the x translation movement multiplier
     * @return The movement multiplier for x translation
     */
    public double getTransXFactor() {
	return transXFactor;
    }

    /**
     * Returns the y translation movement multiplier
     * @return The movement multiplier for y translation
     */
    public double getTransYFactor() {
	return transYFactor;
    }

    /**
     * Returns the zoom movement multiplier
     * @return The movement multiplier for zoom
     */
    public double getZoomFactor() {
	return zoomFactor;
    }

    /**
     * Enables or disables rotation.  The default is true.
     * @param enabled true or false to enable or disable rotate
     */
    public synchronized void setRotateEnable(boolean enabled) {
	rotateEnabled = enabled;
    }

    /**
     * Enables or disables zoom. The default is true.
     * @param enabled true or false to enable or disable zoom
     */
    public synchronized void setZoomEnable(boolean enabled) {
	zoomEnabled = enabled;
    }

    /**
     * Enables or disables translate. The default is true.
     * @param enabled true or false to enable or disable translate
     */
    public synchronized void setTranslateEnable(boolean enabled) {
	translateEnabled = enabled;
    }

    /**
     * Retrieves the state of rotate enabled
     * @return the rotate enable state
     */
    public boolean getRotateEnable() {
	return rotateEnabled;
    }

    /**
     * Retrieves the state of zoom enabled
     * @return the zoom enable state
     */
    public boolean getZoomEnable() {
	return zoomEnabled;
    }

    /**
     * Retrieves the state of translate enabled
     * @return the translate enable state
     */
    public boolean getTranslateEnable() {
	return translateEnabled;
    }

    boolean rotate(MouseEvent evt) {
	if (rotateEnabled) {
	    if ((leftButton == ROTATE) &&
		(!evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((middleButton == ROTATE) &&
		(evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((rightButton == ROTATE) &&
		(!evt.isAltDown() && evt.isMetaDown())) {
		return true;
	    }
	}
	return false;
    }

    boolean zoom(MouseEvent evt) {
	if (zoomEnabled) {
	    if ((leftButton == ZOOM) &&
		(!evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((middleButton == ZOOM) &&
		(evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((rightButton == ZOOM) &&
		(!evt.isAltDown() && evt.isMetaDown())) {
		return true;
	    }
	}
	return false;
    }

    boolean translate(MouseEvent evt) {
	if (translateEnabled) {
	    if ((leftButton == TRANSLATE) &&
		(!evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((middleButton == TRANSLATE) &&
		(evt.isAltDown() && !evt.isMetaDown())) {
		return true;
	    }
	    if ((rightButton == TRANSLATE) &&
		(!evt.isAltDown() && evt.isMetaDown())) {
		return true;
	    }
	}
	return false;
    }

    /**
     * Sets the minimum radius for the myOrbitBehavior.  The zoom will
     * stop at this distance from the center of rotation.  The default
     * is 0.0.  The minimum will have no affect if the STOP_ZOOM constructor
     * flag is not set.
     * @param r the minimum radius
     * @exception IllegalArgumentException if the radius is less than 0.0
     */
    public synchronized void setMinRadius(double r) {
	if (r < 0.0) {
	    throw new IllegalArgumentException(J3dUtilsI18N.getString("myOrbitBehavior1"));
	}
	minRadius = r;
    }

    /**
     * Returns the minimum orbit radius.  The zoom will stop at this distance
     * from the center of rotation if the STOP_ZOOM constructor flag is set.
     * @return the minimum radius
     */
    public double getMinRadius() {
	return minRadius;
    }


}

