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

import java.beans.*;
import java.awt.event.*;
import java.util.HashMap;
import javax.swing.*;
import javax.swing.event.*;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

public class ActionDelegate
    extends javax.swing.AbstractAction
    implements PropertyChangeListener
{
    //--------------------------------------------------------------------------

    public static final String DEFAULT_ACTION = "Default";

    protected HashMap subActions_;
    protected Action currentAction_;

    //--------------------------------------------------------------------------
    /** ActionDelegate constructor. */
    public ActionDelegate(javax.swing.AbstractAction defaultAction)
    {
        try
        {
            subActions_ = new HashMap();
            subActions_.put(DEFAULT_ACTION, defaultAction);

            // Select the default action.
            selectAction(DEFAULT_ACTION);

            // Use the default action's client properties.
            if (defaultAction.getValue(Action.SMALL_ICON) != null)
            {
                putValue(Action.SMALL_ICON, defaultAction.getValue(Action.SMALL_ICON));
            }
            if (defaultAction.getValue(Action.NAME) != null)
            {
                putValue(Action.NAME, defaultAction.getValue(Action.NAME));
            }
            if (defaultAction.getValue(Action.SHORT_DESCRIPTION) != null)
            {
                putValue(Action.SHORT_DESCRIPTION, defaultAction.getValue(Action.SHORT_DESCRIPTION));
            }
            if (defaultAction.getValue(Action.LONG_DESCRIPTION) != null)
            {
                putValue(Action.LONG_DESCRIPTION, defaultAction.getValue(Action.LONG_DESCRIPTION));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** ActionDelegate constructor. */
    public ActionDelegate(String actionName, javax.swing.AbstractAction action)
    {
        try
        {
            subActions_ = new HashMap();
            subActions_.put(actionName, action);

            // Select the default action.
            selectAction(actionName);

            // Use the default action's client properties.
            if (action.getValue(Action.SMALL_ICON) != null)
            {
                putValue(Action.SMALL_ICON, action.getValue(Action.SMALL_ICON));
            }
            if (action.getValue(Action.NAME) != null)
            {
                putValue(Action.NAME, action.getValue(Action.NAME));
            }
            if (action.getValue(Action.SHORT_DESCRIPTION) != null)
            {
                putValue(Action.SHORT_DESCRIPTION, action.getValue(Action.SHORT_DESCRIPTION));
            }
            if (action.getValue(Action.LONG_DESCRIPTION) != null)
            {
                putValue(Action.LONG_DESCRIPTION, action.getValue(Action.LONG_DESCRIPTION));
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Overridded isEnabled method. Return enabled status of current sub-action. */
    public boolean isEnabled()
    {
        // Return the enabled status of the current sub-action.
        return currentAction_ == null ? false : currentAction_.isEnabled();
    }

    //--------------------------------------------------------------------------

    public void addSubAction(String actionName, Action action) throws FoamXException
    {
        // Check for duplicate name.
        if (subActions_.containsKey(actionName))
        {
            throw new FoamXException("Duplicate sub-action name : " + actionName);
        }

        // Add to map.
        subActions_.put(actionName, action);
    }

    //--------------------------------------------------------------------------

    public void removeSubAction(String actionName) throws FoamXException
    {
        // Check for valid name.
        if (!subActions_.containsKey(actionName))
        {
            throw new FoamXException("Invalid sub-action name : " + actionName);
        }

        // Reset current action if we have invalidated it.
        Action subAction = (Action)subActions_.get(actionName);
        if (currentAction_ == subAction)     // Object reference equality?
        {
            // Unregister for property change events.
            currentAction_.removePropertyChangeListener(this);
            currentAction_ = null;
        }

        // Remove from map.
        subActions_.remove(actionName);
    }

    //--------------------------------------------------------------------------

    public void selectAction(String actionName) throws FoamXException
    {
        // Check for valid name.
        if (!subActions_.containsKey(actionName))
        {
            throw new FoamXException("Invalid sub-action name : " + actionName);
        }

        if (currentAction_ != null)
        {
            // Unregister for property change events.
            currentAction_.removePropertyChangeListener(this);
        }

        // Select specified action.
        currentAction_ = (Action)subActions_.get(actionName);

        // Register for property change events.
        currentAction_.addPropertyChangeListener(this);

        // Set the enabled state of this object to reflect the
        // enabled state of the current sub-action object.
        // Causes property change event from this object.
        setEnabled(currentAction_.isEnabled());
    }

    //--------------------------------------------------------------------------

    public synchronized void actionPerformed(ActionEvent evt)
    {
        try
        {
            // Belt and braces!
            if (currentAction_ != null && currentAction_.isEnabled())
            {
                // Invoke sub-action.
                currentAction_.actionPerformed(evt);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void propertyChange(java.beans.PropertyChangeEvent evt)
    {
        // Set the enabled state of this object to reflect the
        // enabled state of the current sub-action object.
        // Causes property change event from this object.
        Boolean enabled = (Boolean)evt.getNewValue();
        setEnabled(enabled.booleanValue());
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
