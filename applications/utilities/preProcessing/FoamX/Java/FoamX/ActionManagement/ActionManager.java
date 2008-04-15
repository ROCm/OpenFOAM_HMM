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
package FoamX.ActionManagement;

import java.lang.ref.*;
import java.util.*;
import javax.swing.*;

import FoamX.App;

public class ActionManager
    extends Object
{
    //--------------------------------------------------------------------------

    protected Hashtable actionMap_;

    //--------------------------------------------------------------------------
    /** ActionManager constructor. */
    public ActionManager()
    {
        actionMap_ = new Hashtable();
    }

    //--------------------------------------------------------------------------

    public boolean isRegistered(String actionName)
    {
        return actionMap_.containsKey(actionName);
    }

    //--------------------------------------------------------------------------

    public void registerAction(String actionName, Action action)
        throws IllegalArgumentException
    {
        if (actionMap_.containsKey(actionName))
        {
            throw new IllegalArgumentException
            (
                "Duplication action name in ActionManager.registerAction"
                + " method call."
            );
        }

        actionMap_.put(actionName, action);
    }

    //--------------------------------------------------------------------------

    public void unregisterAction(String actionName)
        throws IllegalArgumentException
    {
        if (!actionMap_.containsKey(actionName))
        {
            throw new IllegalArgumentException
            (
                "Duplication action name in ActionManager.unregisterAction"
                + " method call."
            );
        }

        actionMap_.remove(actionName);
    }

    //--------------------------------------------------------------------------

    public Action getAction(String actionName)
        throws IllegalArgumentException
    {
        if (!actionMap_.containsKey(actionName))
        {
            throw new IllegalArgumentException
            (
                "Invalid action name in ActionManager.getAction method call."
            );
        }

        return (Action)actionMap_.get(actionName);
    }

    //--------------------------------------------------------------------------

    public void registerActions(IActionProvider provider)
    {
        try
        {
            // Get provider's prefix string.
            String prefix = provider.getPrefix();
            if (prefix != null)
            {
                prefix += ".";
            } else {
                prefix = "";
            }

            // Loop over all actions and register each.
            javax.swing.Action[] actions = provider.getActions();
            for (int i = 0; i <actions.length; i++)
            {
                Action action = actions[i];
                String actionName =
                    prefix + (String)action.getValue(Action.NAME);

                registerAction(actionName, action);
            }
        }
        catch (IllegalArgumentException ex)
        {
            App.printMessage
            (
                "Illegal action name in ActionManager.registerProvider method"
                + " call."
            );
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}








