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
package FoamX.Modules.RawDictionaryEditor;

import java.awt.*;
import java.beans.*;
import javax.swing.*;

public class DictionaryTextEditor
    extends javax.swing.JPanel
{
    //--------------------------------------------------------------------------

    private PropertyChangeSupport propertySupport_;
    private JScrollPane scrollPane_;
    private JTextArea editorPane_;
    private String sDictionaryName_;
    private boolean bModified_;

    //--------------------------------------------------------------------------

    /** DictionaryEditor constructor. */
    public DictionaryTextEditor(String sDictionaryName)
    {
        sDictionaryName_ = sDictionaryName;
        propertySupport_ = new PropertyChangeSupport(this);
        bModified_       = false;

        // Initialise GUI.
        initComponents();
    }

    //--------------------------------------------------------------------------

    public String toString()
    {
        // Used to set the title of the parent internal frame window.
        return sDictionaryName_;
    }

    //--------------------------------------------------------------------------

    /**
     * @param listener
     */
    public void addPropertyChangeListener(PropertyChangeListener listener)
    {
        propertySupport_.addPropertyChangeListener(listener);
    }

    //--------------------------------------------------------------------------

    /**
     * @param listener
     */
    public void removePropertyChangeListener(PropertyChangeListener listener)
    {
        propertySupport_.removePropertyChangeListener(listener);
    }

    //--------------------------------------------------------------------------

    private void initComponents()
    {
        scrollPane_ = new JScrollPane();
        editorPane_ = new JTextArea();
        scrollPane_.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane_.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane_.setViewportView(editorPane_);
        editorPane_.setFont(new Font("monospaced", Font.PLAIN, 10));

        // Set the preferred size of the editor window.
        scrollPane_.setPreferredSize(new Dimension(400, 400));

        setBackground(new java.awt.Color(255, 0, 0));
        setLayout(new GridLayout(1, 1));
        add(scrollPane_);
        repaint();
    }

    //--------------------------------------------------------------------------

    public void setEditorText(String sText)
    {
        editorPane_.setText(sText);
        editorPane_.setCaretPosition(0);
        return;
    }

    //--------------------------------------------------------------------------

    public String getEditorText()
    {
        return editorPane_.getText();
    }

    //--------------------------------------------------------------------------

    public boolean isModified()
    {
        return bModified_;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}

