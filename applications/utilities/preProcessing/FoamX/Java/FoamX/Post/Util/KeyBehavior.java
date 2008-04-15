package FoamX.Post.Util;

import java.awt.AWTEvent;
import java.awt.event.*;
import java.util.Enumeration;

import javax.media.j3d.*;
import javax.vecmath.*;

import com.sun.j3d.utils.universe.*;

import FoamX.Post.PostApp;

/**
 *  KeyBehavior is a generic behavior class to take key presses and move a
 *  TransformGroup through a Java3D scene. The actions resulting from the key strokes
 *  are modified by using the Ctrl, Alt and Shift keys.
 *
 *  (version 1.0) reconstructed class to make more generic.
 *
 * MODIFIED:
 * @version 1.0, 25 September 1998 aajc
 * @author  Andrew AJ Cain, Swinburne University, Australia
 *          <acain@it.swin.edu.au>
 *
 * edited from code by:
 *   Gary S. Moss <moss@arl.mil>
 *   U. S. Army Research Laboratory
 *
 * CLASS NAME:
 *   KeyBehavior
 *
 * PUBLIC FEATURES:
 *   // Data
 *
 *   // Constructors
 *
 *   // Methods:
 *
 * COLLABORATORS:
 *
 */
public class KeyBehavior extends Behavior
{
    protected static final double FAST_SPEED = 2.0;
    protected static final double NORMAL_SPEED = 1.0;
    protected static final double SLOW_SPEED = 0.5;

    private TransformGroup transformGroup;
    private Transform3D transform3D;
    private WakeupCondition keyCriterion;
    private PostApp pp_;

    private final static double TWO_PI = (2.0 * Math.PI);
    private double rotateXAmount = Math.PI / 16.0;
    private double rotateYAmount = Math.PI / 16.0;
    private double rotateZAmount = Math.PI / 16.0;

    private double moveRate = 0.3;
    private double speed = NORMAL_SPEED;

    private int forwardKey = KeyEvent.VK_UP;
    private int backKey = KeyEvent.VK_DOWN;
    private int leftKey = KeyEvent.VK_LEFT;
    private int rightKey = KeyEvent.VK_RIGHT;
    private int resetKey = KeyEvent.VK_F;

    public KeyBehavior(PostApp pp, TransformGroup tg )
    {
        transformGroup = tg;
        transform3D = new Transform3D();
        pp_ = pp;
    }

    public void initialize()
    {
        WakeupCriterion[] keyEvents = new WakeupCriterion[2];

        keyEvents[0] = new WakeupOnAWTEvent( KeyEvent.KEY_PRESSED );
        keyEvents[1] = new WakeupOnAWTEvent( KeyEvent.KEY_RELEASED );

        keyCriterion = new WakeupOr( keyEvents );
        wakeupOn( keyCriterion );
    }





    public void processStimulus( Enumeration criteria )
    {
        WakeupCriterion wakeup;
        AWTEvent[] event;
    
        while( criteria.hasMoreElements() )
        {
            wakeup = (WakeupCriterion) criteria.nextElement();
    
            if( !(wakeup instanceof WakeupOnAWTEvent) )
                continue;

            event = ((WakeupOnAWTEvent)wakeup).getAWTEvent();

            for( int i = 0; i < event.length; i++ )
            {
                if( event[i].getID() == KeyEvent.KEY_PRESSED )
                {
                    processKeyEvent((KeyEvent)event[i]);
                }
            }
        }
        wakeupOn( keyCriterion );
    }

    protected void processKeyEvent(KeyEvent event)
    {
        int keycode = event.getKeyCode();

        if(event.isShiftDown()) speed = FAST_SPEED;
        else speed = NORMAL_SPEED;

        if( event.isAltDown() )
            altMove(keycode);
        else if( event.isControlDown())
            controlMove(keycode);
        else
            standardMove(keycode);
    }





    //moves forward backward or rotates left right
    private void standardMove(int keycode)
    {
        if(keycode == forwardKey)
            moveForward();
        else if(keycode == backKey)
            moveBackward();
        else if(keycode == leftKey)
            rotLeft();
        else if(keycode == rightKey)
            rotRight();
        else if(keycode == resetKey)
            resetView();
    }

  //moves left right, rotate up down
  protected void altMove(int keycode)
  {
    if(keycode == forwardKey)
      rotUp();
    else if(keycode == backKey)
      rotDown();
    else if(keycode == leftKey)
      moveLeft();
    else if(keycode == rightKey)
      moveRight();
  }

  //move up down, rot left right
  protected void controlMove(int keycode)
  {
    if(keycode == forwardKey)
      moveUp();
    else if(keycode == backKey)
      moveDown();
    else if(keycode == leftKey)
      rollLeft();
    else if(keycode == rightKey)
      rollRight();
  }

    private void resetView()
    {
        System.out.println("ResetView");
        View view = pp_.getUniverse().getViewer().getView();
        BoundingSphere bounds = (BoundingSphere)pp_.getScene().getBounds();
        PostApp.viewAll(view, bounds);
    }

  private void moveForward()
  {
    doMove(new Vector3d(0.0,0.0, -getMovementRate()));
  }

  private void moveBackward()
  {
    doMove(new Vector3d(0.0,0.0, getMovementRate()));
  }

  private void moveLeft()
  {
    doMove(new Vector3d( -getMovementRate() ,0.0,0.0));
  }

  private void moveRight()
  {
    doMove(new Vector3d(getMovementRate(),0.0,0.0));
  }

  private void moveUp()
  {
    doMove(new Vector3d(0.0, getMovementRate() ,0.0));
  }

  private void moveDown()
  {
    doMove(new Vector3d(0.0, -getMovementRate() ,0.0));
  }

  protected void rotRight()
  {
    doRotateY(getRotateRightAmount());
  }

  protected void rotUp()
  {
    doRotateX(getRotateUpAmount());
  }

  protected void rotLeft()
  {
    doRotateY(getRotateLeftAmount());
  }

  protected void rotDown()
  {
    doRotateX(getRotateDownAmount());
  }

  protected void rollLeft()
  {
     doRotateZ(getRollLeftAmount());
  }

  protected void rollRight()
  {
     doRotateZ(getRollRightAmount());
  }

  protected void doRotateY(double radians)
  {
    transformGroup.getTransform(transform3D);
    Transform3D toMove = new Transform3D();
    toMove.rotY(radians);
    transform3D.mul(toMove);
    transformGroup.setTransform(transform3D);
  }

  protected void doRotateX(double radians)
  {
    transformGroup.getTransform(transform3D);
    Transform3D toMove = new Transform3D();
    toMove.rotX(radians);
    transform3D.mul(toMove);
    transformGroup.setTransform(transform3D);
  }

  protected void doRotateZ(double radians)
  {
    transformGroup.getTransform(transform3D);
    Transform3D toMove = new Transform3D();
    toMove.rotZ(radians);
    transform3D.mul(toMove);
    transformGroup.setTransform(transform3D);
  }

  protected void doMove(Vector3d theMove)
  {
    transformGroup.getTransform(transform3D);
    Transform3D toMove = new Transform3D();
    toMove.setTranslation(theMove);
    transform3D.mul(toMove);
    transformGroup.setTransform(transform3D);
  }

  protected double getMovementRate()
  {
    return moveRate * speed;
  }

  protected double getRollLeftAmount()
  {
    return rotateZAmount * speed;
  }

  protected double getRollRightAmount()
  {
    return -rotateZAmount * speed;
  }

  protected double getRotateUpAmount()
  {
    return rotateYAmount * speed;
  }

  protected double getRotateDownAmount()
  {
    return -rotateYAmount * speed;
  }

  protected double getRotateLeftAmount()
  {
    return rotateYAmount * speed;
  }

  protected double getRotateRightAmount()
  {
    return -rotateYAmount * speed;
  }

  public void setRotateXAmount(double radians)
  {
    rotateXAmount = radians;
  }

  public void setRotateYAmount(double radians)
  {
    rotateYAmount = radians;
  }

  public void setRotateZAmount(double radians)
  {
    rotateZAmount = radians;
  }

  public void setMovementRate(double meters)
  {
    moveRate = meters; // Travel rate in meters/frame
  }

  public void setForwardKey(int key)
  {
    forwardKey = key;
  }

  public void setBackKey(int key)
  {
    backKey = key;
  }

  public void setLeftKey(int key)
  {
    leftKey = key;
  }

  public void setRightKey(int key)
  {
    rightKey = key;
  }
}
