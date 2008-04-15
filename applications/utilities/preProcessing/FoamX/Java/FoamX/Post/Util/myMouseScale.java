package FoamX.Post.Util;

//import com.sun.j3d.utils.behaviors.vp.*;
import com.sun.j3d.utils.behaviors.mouse.*;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;


/**
 * myMouseScale is a Java3D behavior object that lets users control the 
 * Z axis translation of an object via a mouse drag motion with the second
 * mouse button. See MouseRotate for similar usage info.
 */
 
public class myMouseScale extends MouseZoom
{
    private MouseBehaviorCallback callback = null;

    public myMouseScale(TransformGroup transformGroup) {
	super(transformGroup);
    }

    /**
     * Creates a default mouse zoom behavior.
     **/
    public myMouseScale(){
	super(0);
    }

    /**
     * Creates a zoom behavior.
     * Note that this behavior still needs a transform
     * group to work on (use setTransformGroup(tg)) and
     * the transform group must add this behavior.
     * @param flags
     */
    public myMouseScale(int flags) {
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
    public myMouseScale(Component c) {
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
    public myMouseScale(Component c, TransformGroup transformGroup) {
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
    public myMouseScale(Component c, int flags) {
	super(c, flags);
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

	//processMouseEvent(evt);
	
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
		    
		    double scale = currXform.getScale();

                    scale += 0.01*dy;

System.out.println("setting scale " + scale);

                    currXform.setScale(scale);

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
    * The transformChanged method in the callback class will
    * be called every time the transform is updated
    */
  public void setupCallback( MouseBehaviorCallback callback ) {
      this.callback = callback;
  }
}

