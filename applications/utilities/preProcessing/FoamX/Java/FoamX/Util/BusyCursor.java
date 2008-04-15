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
package FoamX.Util;

import java.awt.*;
import javax.swing.*;

public class BusyCursor
    extends Object
{
    //--------------------------------------------------------------------------

    private java.awt.Component component_;
    private java.awt.Cursor prevCursor_;

    //--------------------------------------------------------------------------
    /** BusyCursor constructor. */
    public BusyCursor()
    {
        component_  = null;
        prevCursor_ = null;
    }

    //--------------------------------------------------------------------------
    /** BusyCursor constructor. */
    public BusyCursor(java.awt.Component component)
    {
        component_  = component;
        prevCursor_ = component_.getCursor();
        component_.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

        // Put a restore cursor event into the event dispatch thread event queue.
        // This will restore the previous cursor when the current activity
        // has finished.
        SwingUtilities.invokeLater(new RestoreCursorEvent());
    }

    //--------------------------------------------------------------------------

    protected void restoreCursor()
    {
        // Restore previous cursor.
        if (component_ != null && prevCursor_ != null)
        {
            component_.setCursor(prevCursor_);
        }
        component_  = null;
        prevCursor_ = null;
    }

    //--------------------------------------------------------------------------

    protected void finalize()
    {
        // Restore previous cursor.
        //restoreCursor();
    }

    //--------------------------------------------------------------------------

    class RestoreCursorEvent
    implements Runnable
    {
        public void run()
        {
            //? Container parent = component_.getParent();

            // Restore the original cursor.
            restoreCursor();

            //? parent.repaint();
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
