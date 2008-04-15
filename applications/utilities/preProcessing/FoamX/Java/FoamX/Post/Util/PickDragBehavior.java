/*
 *	@(#)PickDragBehavior.java 1.8 01/01/11 07:32:11
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

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.geometry.Sphere;

/**
 * Class:       PickDragBehavior
 * 
 * Description: Used to respond to mouse pick and drag events
 *              in the 3D window.
 *
 * Version:     1.0
 *
 */
public class PickDragBehavior extends Behavior
{
    WakeupCriterion[] mouseEvents;
    WakeupOr mouseCriterion;
    int x, y;
    int x_last, y_last;
    double x_angle, y_angle;
    double x_factor, y_factor;
    Transform3D modelTrans;
    Transform3D transformX;
    Transform3D transformY;
    Vector3d transVector;
    TransformGroup transformGroup;

    public PickDragBehavior(TransformGroup transformGroup)
    {

        this.transformGroup = transformGroup;

        modelTrans = new Transform3D();
        transformX = new Transform3D();
        transformY = new Transform3D();
        transVector = new Vector3d();
    }

    public void initialize()
    {
        x = 0;
        y = 0;
        x_last = 0;
        y_last = 0;
        x_angle = 0;
        y_angle = 0;
        x_factor = .02;
        y_factor = .02;

        mouseEvents = new WakeupCriterion[2];
        mouseEvents[0] = new WakeupOnAWTEvent(MouseEvent.MOUSE_DRAGGED);
        mouseEvents[1] = new WakeupOnAWTEvent(MouseEvent.MOUSE_PRESSED);
        mouseCriterion = new WakeupOr(mouseEvents);
        wakeupOn (mouseCriterion);
    }

    public void processStimulus (Enumeration criteria)
    {
        WakeupCriterion wakeup;
        AWTEvent[] event;
        int id;
        int dx, dy;

        while (criteria.hasMoreElements())
        {
            wakeup = (WakeupCriterion) criteria.nextElement();
            if (wakeup instanceof WakeupOnAWTEvent)
            {
                event = ((WakeupOnAWTEvent)wakeup).getAWTEvent();
                for (int i=0; i<event.length; i++)
                {
                    id = event[i].getID();
                    if (id == MouseEvent.MOUSE_DRAGGED)
                    {
                        x = ((MouseEvent)event[i]).getX();
                        y = ((MouseEvent)event[i]).getY();

                        dx = x - x_last;
                        dy = y - y_last;

                        // rotate
	                if (rotate((MouseEvent)event[i]))
                        {

                            x_angle = dy * y_factor;
                            y_angle = dx * x_factor;

                            transformX.rotX(x_angle);
                            transformY.rotY(y_angle);

                            modelTrans.mul(transformX, modelTrans);
                            modelTrans.mul(transformY, modelTrans);

                        }
                        else if (translate((MouseEvent)event[i]))
                        {
		            transformGroup.getTransform(modelTrans);
                            modelTrans.get(transVector);

                            
                            transVector.x += dx*0.1;
                            transVector.y -= dy*0.1;

                            modelTrans.setTranslation(transVector);
                        }
                        else if (zoom((MouseEvent)event[i]))
                        {
                            double scale = modelTrans.getScale();
                            scale += dy*0.01;

                            scale = Math.max(1E-6, scale);
                            //System.out.println("scale:" + scale);
                            modelTrans.setScale(scale);
                        }

                        transformGroup.setTransform(modelTrans);

                        x_last = x;
                        y_last = y;
                    }
                    else if (id == MouseEvent.MOUSE_PRESSED)
                    {

                        x = x_last = ((MouseEvent)event[i]).getX();
                        y = y_last = ((MouseEvent)event[i]).getY();
                    }
                }
            }
        }
        wakeupOn (mouseCriterion);
    }
    

    boolean rotate(MouseEvent evt)
    {
        if (evt.getModifiers() == evt.BUTTON1_MASK)
        {
            return true;
        }
        else
        {
	    return false;
        }
    }

    boolean zoom(MouseEvent evt)
    {
        if (evt.getModifiers() == evt.BUTTON2_MASK)
        {
            return true;
        }
        else
        {
	    return false;
        }
    }

    boolean translate(MouseEvent evt)
    {
        if (evt.getModifiers() == evt.BUTTON3_MASK)
        {
            return true;
        }
        else
        {
	    return false;
        }
    }
}
