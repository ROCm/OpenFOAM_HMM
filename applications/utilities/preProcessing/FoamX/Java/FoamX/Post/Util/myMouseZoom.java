/*
 *	@(#)myMouseZoom.java 1.31 01/01/11 07:22:24
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

import com.sun.j3d.utils.behaviors.mouse.*;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;


/**
 * myMouseZoom is a Java3D behavior object that lets users control the 
 * Z axis translation of an object via a mouse drag motion with the second
 * mouse button. See MouseRotate for similar usage info.
 */
 
public class myMouseZoom extends MouseBehavior {

    double z_factor = .04;
    Vector3d translation = new Vector3d();
  
    private MouseBehaviorCallback callback = null;

    /**
     * Creates a zoom behavior given the transform group.
     * @param transformGroup The transformGroup to operate on.
     */
    public myMouseZoom(TransformGroup transformGroup) {
	super(transformGroup);
    }

    /**
     * Creates a default mouse zoom behavior.
     **/
    public myMouseZoom(){
	super(0);
    }

    /**
     * Creates a zoom behavior.
     * Note that this behavior still needs a transform
     * group to work on (use setTransformGroup(tg)) and
     * the transform group must add this behavior.
     * @param flags
     */
    public myMouseZoom(int flags) {
	super(flags);
    }

    /**
     * Creates a zoom behavior that uses AWT listeners and behavior
     * posts rather than WakeupOnAWTEvent.  The behavior is added to the
     * specified Component.  A null component can be passed to specify
     * the behavior should use listeners.  Components can then be added
     * to the behavior with the addListener(Component c) method.
     * @param c The Component to add the MouseListener
     * and MouseMotionListener to.
     * @since Java 3D 1.2.1
     */
    public myMouseZoom(Component c) {
	super(c, 0);
    }

    /**
     * Creates a zoom behavior that uses AWT listeners and behavior
     * posts rather than WakeupOnAWTEvent.  The behaviors is added to
     * the specified Component and works on the given TransformGroup.
     * @param c The Component to add the MouseListener and
     * MouseMotionListener to.  A null component can be passed to specify
     * the behavior should use listeners.  Components can then be added
     * to the behavior with the addListener(Component c) method.
     * @param transformGroup The TransformGroup to operate on.
     * @since Java 3D 1.2.1
     */
    public myMouseZoom(Component c, TransformGroup transformGroup) {
	super(c, transformGroup);
    }

    /**
     * Creates a zoom behavior that uses AWT listeners and behavior
     * posts rather than WakeupOnAWTEvent.  The behavior is added to the
     * specified Component.  A null component can be passed to specify
     * the behavior should use listeners.  Components can then be added
     * to the behavior with the addListener(Component c) method.
     * Note that this behavior still needs a transform
     * group to work on (use setTransformGroup(tg)) and the transform
     * group must add this behavior.
     * @param flags interesting flags (wakeup conditions).
     * @since Java 3D 1.2.1
     */
    public myMouseZoom(Component c, int flags) {
	super(c, flags);
    }

    public void initialize() {
	super.initialize();
	if ((flags & INVERT_INPUT) == INVERT_INPUT) {
	    z_factor *= -1;
	    invert = true;
	}
    }
    
    /**
     * Return the y-axis movement multipler.
     **/
    public double getFactor() {
	return z_factor;
    }
  
    /**
     * Set the y-axis movement multipler with factor.
     **/
    public void setFactor( double factor) {
	z_factor = factor;
    }
  

    public void processStimulus (Enumeration criteria) {
	WakeupCriterion wakeup;
	AWTEvent[] events;
 	MouseEvent evt;
// 	int id;
// 	int dx, dy;

        System.out.println("processStimulus");

	while (criteria.hasMoreElements()) {
	    wakeup = (WakeupCriterion) criteria.nextElement();
	    if (wakeup instanceof WakeupOnAWTEvent) {
		events = ((WakeupOnAWTEvent)wakeup).getAWTEvent();
		if (events.length > 0) {
		    evt = (MouseEvent) events[events.length-1];
		    doProcess(evt);
		}
	    }

	    else if (wakeup instanceof WakeupOnBehaviorPost) {
		while (true) {
		    synchronized (mouseq) {
			if (mouseq.isEmpty()) break;
			evt = (MouseEvent)mouseq.remove(0);
			// consolodate MOUSE_DRAG events
			while((evt.getID() == MouseEvent.MOUSE_DRAGGED) &&
			      !mouseq.isEmpty() &&
			      (((MouseEvent)mouseq.get(0)).getID() ==
			       MouseEvent.MOUSE_DRAGGED)) {
			    evt = (MouseEvent)mouseq.remove(0);
			}
		    }
		    doProcess(evt);
		}
	    }
	    
	}
	wakeupOn (mouseCriterion);
    }

    void doProcess(MouseEvent evt) {
	int id;
	int dx, dy;

	processMouseEvent(evt);
	
        System.out.println("doProcess");

	if (((buttonPress)&&((flags & MANUAL_WAKEUP) == 0)) ||
	    ((wakeUp)&&((flags & MANUAL_WAKEUP) != 0))){
	    id = evt.getID();
	    if ((id == MouseEvent.MOUSE_DRAGGED) &&
		evt.isAltDown() && !evt.isMetaDown()){
		
		x = evt.getX();
		y = evt.getY();
		
		dx = x - x_last;
		dy = y - y_last;
		
		if (!reset){
		    transformGroup.getTransform(currXform);
		    
		    translation.z  = dy*z_factor;
		    
		    transformX.set(translation);
		    
		    if (invert) {
			currXform.mul(currXform, transformX);
		    } else {
			currXform.mul(transformX, currXform);
		    }
		    
		    transformGroup.setTransform(currXform);
		    
		    transformChanged( currXform );
		    
		    if (callback!=null)
			callback.transformChanged( MouseBehaviorCallback.ZOOM,
						   currXform );
		    
		}
		else {
		    reset = false;
		}
		
		x_last = x;
		y_last = y;
	    }
	    else if (id == MouseEvent.MOUSE_PRESSED) {
		x_last = evt.getX();
		y_last = evt.getY();
	    }
	}
    }

    
  /**
    * Users can overload this method  which is called every time
    * the Behavior updates the transform
    *
    * Default implementation does nothing
    */
  public void transformChanged( Transform3D transform ) {
  }
 
  /**
    * The transformChanged method in the callback class will
    * be called every time the transform is updated
    */
  public void setupCallback( MouseBehaviorCallback callback ) {
      this.callback = callback;
  }
}

