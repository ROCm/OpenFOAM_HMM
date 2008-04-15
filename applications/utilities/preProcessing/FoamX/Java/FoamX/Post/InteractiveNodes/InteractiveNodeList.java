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

package FoamX.Post.InteractiveNodes;

import javax.swing.table.*;
import javax.swing.event.TableModelListener; 
import javax.swing.event.TableModelEvent; 

import java.util.*;

public class InteractiveNodeList
    extends AbstractTableModel
    implements TableModelListener
{
    private boolean DEBUG = false;

    final String[] columnNames_ = {"Name", "Show"};

    Vector shapes_;

    //--------------------------------------------------------------------------

    public InteractiveNodeList()
    {
        shapes_ = new Vector();
    }


    //--------------------------------------------------------------------------
    // Vector-like interface
    //--------------------------------------------------------------------------

    public void addElement(InteractiveNode shape)
    {
        shapes_.addElement(shape);

        int maxRow = this.getRowCount() - 1;
        int nCols = this.getColumnCount();
        for(int i = 0; i < nCols; i++)
        {
            fireTableCellUpdated(maxRow, i);
        }
    }

    //--------------------------------------------------------------------------

    public void removeElement(InteractiveNode shape)
    {
        int index = shapes_.indexOf(shape);

        this.removeElementAt(index);
    }

    //--------------------------------------------------------------------------

    public void removeElementAt(int index)
    {
        InteractiveNode shape = (InteractiveNode)shapes_.elementAt(index);
        shape.hide();
        shape.remove();
        shapes_.removeElementAt(index);
    }


    //--------------------------------------------------------------------------
    // Shape specifics
    //--------------------------------------------------------------------------

    public InteractiveNode getByName(String name)
    {
        Enumeration enum = shapes_.elements();
        while (enum.hasMoreElements())
        {
            InteractiveNode shape = (InteractiveNode)enum.nextElement();
            if (shape.name().equals(name))
            {
                return shape;
            }
        }
        return null;
    }

    //--------------------------------------------------------------------------

    public void show(String name)
    {
        InteractiveNode shape = this.getByName(name);
        if (shape != null)
        {
            shape.show();
        }
    }

    //--------------------------------------------------------------------------

    public void hide(String name)
    {
        InteractiveNode shape = this.getByName(name);
        if (shape != null)
        {
            shape.hide();
        }
    }

    //--------------------------------------------------------------------------

    public void showAll()
    {
        Enumeration enum = shapes_.elements();
        while (enum.hasMoreElements())
        {
            InteractiveNode shape = (InteractiveNode)enum.nextElement();
            shape.show();
        }
    }

    //--------------------------------------------------------------------------

    public void hideAll()
    {
        Enumeration enum = shapes_.elements();
        while (enum.hasMoreElements())
        {
            InteractiveNode shape = (InteractiveNode)enum.nextElement();
            shape.hide();
        }
    }

    //--------------------------------------------------------------------------

    public void println()
    {
        System.out.println("-- Contents --");

        int i = 0;
        Enumeration enum = shapes_.elements();
        while (enum.hasMoreElements())
        {
            InteractiveNode shape = (InteractiveNode)enum.nextElement();
            System.out.println(i + " " + shape.name() + " " + shape.isShowing());
            i++;
        }

        System.out.println();
    }


    //--------------------------------------------------------------------------
    // AbstractTableModel interface
    //--------------------------------------------------------------------------

    public int getColumnCount()
    {
        return columnNames_.length;
    }

    public int getRowCount()
    {
        return shapes_.size();
    }

    public String getColumnName(int col)
    {
        return columnNames_[col];
    }

    public Object getValueAt(int row, int col)
    {
        InteractiveNode shape = (InteractiveNode)shapes_.elementAt(row);

        if (col == 0)
        {
            return shape.name();
        }
        else
        {
            return new Boolean(shape.isShowing());
        }
    }


    public Class getColumnClass(int c) {
        return getValueAt(0, c).getClass();
    }


    public boolean isCellEditable(int row, int col) {
        if (col == 0) { 
            return false;
        } else {
            return true;
        }
    }


    public void setValueAt(Object value, int row, int col) {
        if (DEBUG)
        {
            System.out.println
            (
                "Setting value at " + row + "," + col
              + " to " + value
              + " (an instance of " 
              + value.getClass() + ")"
            );
        }

        InteractiveNode shape = (InteractiveNode)shapes_.elementAt(row);

        if (col != 1)
        {
            System.out.println("**ERROR: col: " + col);
        }

        Boolean bool = (Boolean)value;
        if (bool.booleanValue() == true)
        {
            shape.show();
        }
        else
        {
            shape.hide();
        }

        if (DEBUG)
        {
            println();
        }

        fireTableCellUpdated(row, col);
    }


    //
    // Implementation of the TableModelListener interface 
    //

    // By default forward all events to all the listeners. 
    public void tableChanged(TableModelEvent e) {
        fireTableChanged(e);
    }
}
