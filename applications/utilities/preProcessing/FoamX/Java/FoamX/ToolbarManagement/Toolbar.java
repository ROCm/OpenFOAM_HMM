/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/
package FoamX.ToolbarManagement;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import FoamX.App;
import FoamX.Util.ReferenceCount;

public class Toolbar
    extends JToolBar
    implements ReferenceCount
{
    //--------------------------------------------------------------------------
    /** Constant of grip layout constraint. */
    private static final String GRIP = "Grip";
    /** Constant of action layout constraint. */
    private static final String ACTION = "Action";
    /** Basic toolbar height. */
    public static final int BASIC_HEIGHT = 34;
    /** 5 pixels is tolerance of toolbar height so toolbar can be high (BASIC_HEIGHT + HEIGHT_TOLERANCE) but it will be set to BASIC_HEIGHT high. */
    static int HEIGHT_TOLERANCE = 5;
    /** Horizontal separator gap. */
    private static final int SEPARARTOR_GAP = 3;
    /** TOP of toolbar empty border. */
    static int TOP    = 2;
    /** LEFT of toolbar empty border. */
    static int LEFT   = 4;
    /** BOTTOM of toolbar empty border. */
    static int BOTTOM = 2;
    /** RIGHT of toolbar empty border. */
    static int RIGHT  = 4;

    /** Is toolbar floatable */
    private boolean floatable_;
    /** Toolbar DnDListener */
    private DnDListener listener_;
    /** Toolbar mouse listener */
    private ToolbarMouseListener mouseListener_;
    /** Display name of the toolbar */
    private String displayName_;
    /** Map of action objects managed by this toolbar. */
    private HashMap actionMap_;
    /** Map of toolbar button objects managed by this toolbar. */
    private HashMap buttonMap_;
    /** Reference count for multi-toolbars. */
    private int refCount_;

    static final long serialVersionUID = 5011742660516204764L;

    //--------------------------------------------------------------------------
    /**
     * Create a new not floatable Toolbar with programmatic name.
     * Display name is set to be the same as name.
     */
    public Toolbar(String name)
    {
        this(name, name, false);
    }

    /**
     * Create a new not floatable Toolbar with specified programmatic name
     * and display name.
     */
    public Toolbar(String name, String displayName)
    {
        this(name, displayName, false);
    }

    /**
     * Create a new <code>Toolbar</code>.
     * @param name a <code>String</code> containing the associated name
     * @param f specified if Toolbar is floatable
     * Display name of the toolbar is set equal to the name.
     */
    public Toolbar(String name, boolean f)
    {
        this(name, name, f);
    }

    /**
     * Create a new <code>Toolbar</code>.
     * @param name a <code>String</code> containing the associated name
     * @param f specifies if Toolbar is floatable
     */
    public Toolbar(String name, String displayName, boolean f)
    {
        super();
        setDisplayName(displayName);
        initAll(name, f);
    }

    //--------------------------------------------------------------------------

    private void initAll(String name, boolean f)
    {
        floatable_     = f;
        mouseListener_ = null;
        actionMap_     = new HashMap();
        buttonMap_     = new HashMap();
        refCount_      = 0;

        setName(name);

        setFloatable(false);
        setBorder(new CompoundBorder(new EtchedBorder(),
                                       new EmptyBorder(TOP, LEFT, BOTTOM, RIGHT)));
        setLayout(new InnerLayout(2, 2));
        addGrip();
    }

    //--------------------------------------------------------------------------
    /**
     * When Toolbar is floatable, ToolbarGrip is added as Grip as first toolbar component
     */
    private void addGrip()
    {
        if (floatable_)
        {
            ToolbarGrip grip = new ToolbarGrip();

            if (mouseListener_ == null)
            {
                mouseListener_ = new ToolbarMouseListener();
            }

            grip.addMouseListener(mouseListener_);
            grip.addMouseMotionListener(mouseListener_);

            add(GRIP, grip);
            //grip.updateUI();
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Compute with HEIGHT_TOLERANCE number of rows for specific toolbar height.
     * @param height of some toolbar
     * @return number of rows
     */
    static public int rowCount(int height)
    {
        int rows       = 1;
        int max_height = (BASIC_HEIGHT + HEIGHT_TOLERANCE);
        while (height> max_height)
        {
            rows++;
            height -= max_height;
        }
        return rows;
    }

    //--------------------------------------------------------------------------
    /**
     * @return Display name of this toolbar. Display name is localizable,
     * on the contrary to the programmatic name.
     */
    public String getDisplayName()
    {
        return displayName_;
    }

    //--------------------------------------------------------------------------
    /**
     * Sets new display name of this toolbar. Display name is localizable,
     * on the contrary to the programmatic name.
     */
    public void setDisplayName(String displayName)
    {
        displayName_ = displayName;
    }

    //--------------------------------------------------------------------------

    public JButton addMultiAction(String actionName, String clientKey, javax.swing.AbstractAction action)
    {
        JButton button = null;

        try
        {
            if (actionMap_.containsKey(actionName))
            {
                ActionDelegate delegate = (ActionDelegate)actionMap_.get(actionName);

                // Add sub-action to delegate object.
                delegate.addSubAction(clientKey, action);

                // Return the appropriate button for this action.
                button = (JButton)buttonMap_.get(actionName);
            }
            else
            {
                // Create a new action delegate object.
                ActionDelegate delegate = new ActionDelegate(clientKey, action);

                // Add action to map.
                actionMap_.put(actionName, delegate);

                // Add button to toolbar connected to delegate.
                button = add(delegate);

                // Add button to map.
                buttonMap_.put(actionName, button);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return button;
    }

    //--------------------------------------------------------------------------

    public void selectClient(String clientKey)
    {
        try
        {
            // Loop over all actions and set the action key,
            Iterator iter = actionMap_.values().iterator();
            while (iter.hasNext())
            {
                ActionDelegate delegate = (ActionDelegate)iter.next();
                delegate.selectAction(clientKey);
            }
            //invalidate();
            //updateUI();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void removeMultiActions(String clientKey)
    {
        try
        {
            // Loop over all actions and set the action key,
            Iterator iter = actionMap_.values().iterator();
            while (iter.hasNext())
            {
                ActionDelegate delegate = (ActionDelegate)iter.next();
                delegate.removeSubAction(clientKey);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- Container Overrides
    //--------------------------------------------------------------------------
    /** Add new Component as ACTION. */
    public Component add(Component comp)
    {
        add(ACTION, comp);
        return comp;
    }

    //--------------------------------------------------------------------------
    /** Removes all ACTION components. */
    public void removeAll()
    {
        super.removeAll();
        addGrip();
    }

    //--------------------------------------------------------------------------
    //---- JToolBar Overrides
    //--------------------------------------------------------------------------

    public JButton add(Action action)
    {
        return super.add(action);
    }

    public void addSeparator()
    {
        addSeparator(new Dimension(SEPARARTOR_GAP, SEPARARTOR_GAP));
    }

    //--------------------------------------------------------------------------
    //---- Drag & Drop Methods
    //--------------------------------------------------------------------------
    /**
     * Set DnDListener to Toolbar.
     * @param DndListener for toolbar
     */
    public void setDnDListener(DnDListener l)
    {
        listener_ = l;
    }

    //--------------------------------------------------------------------------
    /**
     * Fire drag of Toolbar
     * @param dx distance of horizontal dragging
     * @param dy distance of vertical dragging
     * @param type type of toolbar dragging
     */
    protected void fireDragToolbar(int dx, int dy, int type)
    {
        if (listener_ != null)
        {
            listener_.dragToolbar(new DnDEvent(this, getName(), dx, dy, type));
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Fire drop of Toolbar
     * @param dx distance of horizontal dropping
     * @param dy distance of vertical dropping
     * @param type type of toolbar dropping
     */
    protected void fireDropToolbar(int dx, int dy, int type)
    {
        if (listener_ != null)
        {
            listener_.dropToolbar(new DnDEvent(this, getName(), dx, dy, type));
        }
    }

    //--------------------------------------------------------------------------
    /** DnDListener is Drag and Drop listener for Toolbar motion events. */
    public interface DnDListener
    extends java.util.EventListener
    {
        /** Invoked when toolbar is dragged. */
        public void dragToolbar(DnDEvent e);

        /** Invoked when toolbar is dropped. */
        public void dropToolbar(DnDEvent e);
    }

    //--------------------------------------------------------------------------
    //---- ReferenceCount Interface
    //--------------------------------------------------------------------------
    public synchronized int addRef()
    {
        // Return new reference count.
        return ++refCount_;
    }

    public synchronized int release()
    {
        // Return new reference count.
        return --refCount_;
    }

    public synchronized int refCount()
    {
        // Return current reference count.
        return refCount_;
    }

    //--------------------------------------------------------------------------
    //---- DnDEvent Class
    //--------------------------------------------------------------------------

    public static class DnDEvent
    extends EventObject
    {
        /** Type of DnDEvent. Dragging with only one Toolbar. */
        public static final int DND_ONE  = 1;
        /** Type of DnDEvent. Only horizontal dragging with Toolbar and it's followers. */
        public static final int DND_END  = 2;
        /** Type of DnDEvent. Only vertical dragging with whole lines. */
        public static final int DND_LINE = 3;

        /** Name of toolbar where event occured. */
        private String name;
        /** distance of horizontal dragging */
        private int dx;
        /** distance of vertical dragging */
        private int dy;
        /** Type of event. */
        private int type;

        static final long serialVersionUID = 4389530973297716699L;

        public DnDEvent(Toolbar toolbar, String n, int x, int y, int t)
        {
            super(toolbar);

            name = n;
            dx   = x;
            dy   = y;
            type = t;
        }

        /** @return name of toolbar where event occured. */
        public String getName()
        {
            return name;
        }

        /** @return distance of horizontal dragging */
        public int getDX()
        {
            return dx;
        }

        /** @return distance of vertical dragging */
        public int getDY()
        {
            return dy;
        }

        /** @return type of event. */
        public int getType()
        {
            return type;
        }
    }

    //--------------------------------------------------------------------------
    //---- ToolbarMouseListener Class
    //--------------------------------------------------------------------------

    class ToolbarMouseListener
    implements MouseInputListener
    {
        /** Is toolbar dragging now. */
        private boolean dragging_;
        /** Start point of dragging. */
        private Point startPoint_;

        public ToolbarMouseListener()
        {
            dragging_   = false;
            startPoint_ = null;
        }

        /** Invoked when the mouse has been clicked on a component. */
        public void mouseClicked(MouseEvent e)
        {
        }

        /** Invoked when the mouse enters a component. */
        public void mouseEntered(MouseEvent e)
        {
        }

        /** Invoked when the mouse exits a component. */
        public void mouseExited(MouseEvent e)
        {
        }

        /** Invoked when a mouse button has been pressed on a component. */
        public void mousePressed(MouseEvent e)
        {
            startPoint_ = new Point(e.getX(), e.getY());
        }

        /** Invoked when a mouse button has been released on a component. */
        public void mouseReleased(MouseEvent e)
        {
            if (dragging_)
            {
                fireDropToolbar(e.getX() - startPoint_.x,
                                 e.getY() - startPoint_.y,
                                 DnDEvent.DND_ONE);
                dragging_ = false;
            }
        }

        /** Invoked when a mouse button is pressed on a component and then dragged. */
        public void mouseDragged(MouseEvent e)
        {
            int m    = e.getModifiers();
            int type = DnDEvent.DND_ONE;
            if (e.isControlDown())
            {
                type = DnDEvent.DND_LINE;
            }
            else if (((m & InputEvent.BUTTON2_MASK)   != 0) ||
                      ((m & InputEvent.BUTTON3_MASK) != 0))
            {
                type = DnDEvent.DND_END;
            }
            if (startPoint_ == null)
            {
                startPoint_ = new Point(e.getX(), e.getY());
            }
            fireDragToolbar(e.getX() - startPoint_.x,
                             e.getY() - startPoint_.y,
                             type);
            dragging_ = true;
        }

        /** Invoked when the mouse button has been moved on a component (with no buttons no down). */
        public void mouseMoved(MouseEvent e)
        {
        }
    }

    //--------------------------------------------------------------------------
    //---- ToolbarGrip Class
    //--------------------------------------------------------------------------

    private class ToolbarGrip
    extends JPanel
    {
        /** Horizontal gaps. */
        static final int HGAP  = 1;
        /** Vertical gaps. */
        static final int VGAP  = 1;
        /** Step between two grip elements. */
        static final int STEP  = 1;
        /** Width of grip element. */
        static final int WIDTH = 2;

        /** Number of grip elements. */
        int columns_;
        /** Minimum size. */
        Dimension dim_;

        static final long serialVersionUID = -8819972936203315276L;

        /** Create new ToolbarGrip for default number of grip elements. */
        public ToolbarGrip()
        {
            this(2);
        }

        /** Create new ToolbarGrip for specific number of grip elements.
         * @param col number of grip elements
         */
        public ToolbarGrip(int col)
        {
            super();
            columns_ = col;
            int width = (col - 1) * STEP + col * WIDTH + 2 * HGAP;
            dim_ = new Dimension(width, width);
            setBorder(new EmptyBorder(VGAP, HGAP, VGAP, HGAP));
            setToolTipText(Toolbar.this.getDisplayName());
        }

        /** Paint grip to specific Graphics. */
        public void paint(Graphics g)
        {
            Dimension size = getSize();
            int top = VGAP;
            int bottom = size.height - 1 - VGAP;
            int height = bottom - top;
            g.setColor(getBackground());

            for (int i = 0, x = HGAP; i <columns_; i++, x += WIDTH + STEP)
            {
                g.draw3DRect(x, top, WIDTH, height, true); // Grip element is 3D rectangle now
            }
        }

        /** @return minimum size */
        public Dimension getMinimumSize()
        {
            return dim_;
        }

        /** @return preferred size */
        public Dimension getPreferredSize()
        {
            return getMinimumSize();
        }
    }

    //--------------------------------------------------------------------------
    //---- InnerLayout Class
    //--------------------------------------------------------------------------

    private class InnerLayout
    implements LayoutManager
    {
        /** Grip component */
        private Component grip;

        /** Vector of Components */
        private Vector actions;

        /** horizontal gap */
        private int hgap;

        /** vertical gap */
        private int vgap;

        /**
         * Constructs a new InnerLayout.
         */
        public InnerLayout()
        {
            this(5, 5);
        }

        /**
         * Constructs a new InnerLayout with the specified gap values.
         * @param hgap the horizontal gap variable
         * @param vgap the vertical gap variable
         */
        public InnerLayout(int h, int v)
        {
            hgap = h;
            vgap = v;
            actions = new Vector();
        }

        /**
         * Adds the specified component with the specified name to the layout.
         */
        public void addLayoutComponent(String name, Component comp)
        {
            synchronized (comp.getTreeLock())
            {
                if (GRIP.equals(name))
                {
                    grip = comp;
                }
                else if (ACTION.equals(name))
                {
                    actions.addElement(comp);
                }
                else
                {
                    throw new IllegalArgumentException("Cannot add to layout : Unknown constraint : " + name);
                }
            }
        }

        /**
         * Adds the specified component to the layout, using the specified
         * constraint object.
         */
        public void addLayoutComponent(Component comp, Object constraints)
        {
            throw new IllegalArgumentException();
        }

        /**
         * Removes the specified component from the layout.
         */
        public void removeLayoutComponent(Component comp)
        {
            synchronized (comp.getTreeLock())
            {
                if (grip == comp)
                    grip = null;
                else
                    actions.removeElement(comp);
            }
        }

        /**
         * Calculates the preferred size dimensions for the specified panal given
         * the components in the specified parent container.
         */
        public Dimension preferredLayoutSize(Container target)
        {
            synchronized (target.getTreeLock())
            {
                Dimension dim = new Dimension(0, 0);
                int size = actions.size();
                if ((grip != null) && grip.isVisible())
                {
                    Dimension d = grip.getPreferredSize();
                    dim.width += d.width;
                    dim.width += hgap;
                }
                for (int i = 0; i <size; i++)
                {
                    Component comp = (Component)actions.elementAt(i);
                    if (comp.isVisible())
                    {
                        Dimension d = comp.getPreferredSize();
                        dim.width += d.width;
                        dim.height = Math.max(dim.height, d.height);
                        dim.width += hgap;
                    }
                }
                Insets insets = target.getInsets();
                dim.width  += insets.left + insets.right;
                dim.height += insets.top  + insets.bottom;

                return dim;
            }
        }

        /**
         * Calculates the minimum size dimensions for the specified panal given
         * the components in the specified parent container.
         */
        public Dimension minimumLayoutSize(Container target)
        {
            return preferredLayoutSize(target);
        }

        /**
         * Calculates the maximum size dimensions for the specified panal given
         * the components in the specified parent container.
         */
        public Dimension maximumLayoutSize(Container target)
        {
            synchronized (target.getTreeLock())
            {
                Dimension dim = new Dimension(0, 0);
                int size = actions.size();
                if ((grip != null) && grip.isVisible())
                {
                    Dimension d = grip.getPreferredSize();
                    dim.width += d.width;
                    dim.width += hgap;
                }
                Component last = null;
                for (int i = 0; i <size; i++)
                {
                    Component comp = (Component)actions.elementAt (i);
                    if (comp.isVisible())
                    {
                        Dimension d = comp.getPreferredSize();
                        dim.width += d.width;
                        dim.height = Math.max (dim.height, d.height);
                        dim.width += hgap;
                        last = comp;
                    }
                }
                if (last != null)
                {
                    Dimension prefSize = last.getPreferredSize();
                    Dimension maxSize = last.getMaximumSize();
                    if (prefSize != maxSize)
                    {
                        dim.width  = dim.width - prefSize.width + maxSize.width;
                        dim.height = Math.max(dim.height, maxSize.height);
                    }
                }
                Insets insets = target.getInsets();
                dim.width += insets.left + insets.right;
                dim.height += insets.top + insets.bottom;

                return dim;
            }
        }

        /**
         * Lays out the container in the specified panel.
         */
        public void layoutContainer(Container target)
        {
            synchronized (target.getTreeLock())
            {
                Insets insets = target.getInsets();
                Dimension dim = target.getSize();

                int fullHeight = dim.height;
                int bottom = dim.height - insets.bottom - 1;
                int top = 0 + insets.top;
                int left = insets.left;
                int maxPosition = dim.width - (insets.left + insets.right) - hgap;
                int right = dim.width - insets.right;
                maxPosition = right;
                int height = bottom - top;
                int size = actions.size();
                int w, h;
                if ((grip != null) && grip.isVisible())
                {
                    Dimension d = grip.getPreferredSize();
                    grip.setBounds (left, top + vgap, d.width, bottom - top - 2*vgap);
                    left += d.width;
                    left += hgap;
                }
                for (int i = 0; i <size; i++)
                {
                    left += hgap;
                    Component comp = (Component)actions.elementAt (i);
                    Dimension d;
                    Dimension minSize = comp.getMinimumSize();

                    if (i == size - 1)
                        d = comp.getMaximumSize();
                    else
                        d = comp.getPreferredSize();
                    w = d.width;
                    h = Math.min (height, d.height);
                    if ((left <maxPosition) &&
                            (left + w> maxPosition))
    {
                        if (maxPosition - left>= minSize.width)
                            w = maxPosition - left;
                    }

                    comp.setBounds (left, (fullHeight - d.height)/2, w, h);

                    left += d.width;
                }
            }
        }

        /**
         * Returns the alignment along the x axis.
         */
        public float getLayoutAlignmentX(Container target)
        {
            return 0;
        }

        /**
         * Returns the alignment along the y axis.
         */
        public float getLayoutAlignmentY(Container target)
        {
            return (float)0.5;
        }

        /**
         * Invalidates the layout, indicating that if the layout manager has cached
         * information it should ne discarded.
         */
        public void invalidateLayout(Container target)
        {
            // check target components with local vars (grip, actions)
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
