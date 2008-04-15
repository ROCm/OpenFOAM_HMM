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
package FoamX.Editors.TypeEditor;

import java.awt.Toolkit;
import javax.swing.*;
import javax.swing.table.*;

import FoamX.App;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.StringListEditor;

class TypeParameter
{
    //--------------------------------------------------------------------------

    // Common properties.
    public static final int PARAM_TYPE           = 0;
    public static final int PARAM_PATH           = 1;
    public static final int PARAM_NAME           = 2;
    public static final int PARAM_DISPLAYNAME    = 3;
    public static final int PARAM_DESCRIPTION    = 4;
    public static final int PARAM_COMMENT        = 5;
    public static final int PARAM_CATEGORY       = 6;
    public static final int PARAM_HELPURL        = 7;
    public static final int PARAM_ICONURL        = 8;
    public static final int PARAM_OPTIONAL       = 9;
    public static final int PARAM_VISIBLE        = 10;

    // Properties for non-compound types.
    public static final int PARAM_EDITABLE       = 11;
    public static final int PARAM_MINVALUE       = 12;
    public static final int PARAM_MAXVALUE       = 13;
    public static final int PARAM_DEFAULTVALUE   = 14;
    public static final int PARAM_LOOKUPDICT     = 15;
    public static final int PARAM_VALUELIST      = 16;

    // Properties for compound types.
    public static final int PARAM_DICTIONARYPATH = 17;
    public static final int PARAM_NUMELEMENTS    = 18;
    public static final int PARAM_ELEMENTLABELS  = 19;

    private static final String[] parameterNames_ =
    {
        "Type",
        "Path",
        "Foam Name",
        "Display Name",
        "Description",
        "Comment",
        "Category",
        "Help Url",
        "Icon Url",
        "Optional",
        "Visible",
        "Editable",
        "Minimum Value",
        "Maximum Value",
        "Default Value",
        "Value Lookup Dictionary",
        "Constrained Values",
        "Dictionary Path",
        "Vector Size",
        "Vector Element Labels"
    };

    //--------------------------------------------------------------------------

    private TypeDescriptorCache typeDescriptor_;
    private int parameterID_;
    private TableCellEditor customEditor_;
    private TableCellEditor entryEditor_;
    private boolean editable_;

    //--------------------------------------------------------------------------
    /** TypeParameter constructor. */
    TypeParameter(TypeDescriptorCache typeDescriptor, int parameterID)
    {
        try
        {
            typeDescriptor_ = typeDescriptor;
            parameterID_    = parameterID;
            customEditor_   = null;
            entryEditor_   = null;
            editable_       = true;

            // See if this parameter needs a custom editor.
            if (parameterID_ == PARAM_DEFAULTVALUE)
            {
                entryEditor_ = typeDescriptor_.getDefaultValue().getEditor();
            }
            else if
            (
                parameterID_ == PARAM_VALUELIST
             || parameterID_ == PARAM_ELEMENTLABELS
            )
            {
                TypeParameterCustomEditor editor =
                    new TypeParameterCustomEditor();

                // Subscribe to the edit button's action event.
                editor.addActionListener
                (
                    new java.awt.event.ActionListener()
                    {
                        public void actionPerformed
                        (
                            java.awt.event.ActionEvent evt
                        )
                        {
                            showEditDialog();
                        }
                    }
                );
                customEditor_ = editor;
            }
            else if
            (
                parameterID_ == PARAM_OPTIONAL
             || parameterID_ == PARAM_VISIBLE
             || parameterID_ == PARAM_EDITABLE
            )
            {
                JComboBox combo = new JComboBox();
                combo.addItem("True");
                combo.addItem("False");
                combo.setFont(new java.awt.Font("Dialog", 0, 10));
                customEditor_ = new DefaultCellEditor(combo)
                {
                    // Hook into getCellEditorValue to intercept the
                    // string value.
                    public Object getCellEditorValue()
                    {
                        String objValue = (String)super.getCellEditorValue();
                        updateValue(objValue);

                        // Model expects a TypeParameter object.
                        return TypeParameter.this;
                    }
                };
            }

            // See if this parameter is read only.
            if
            (
                parameterID_ == PARAM_TYPE ||
                parameterID_ == PARAM_PATH
            )
            {
                editable_ = false;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    String getName()
    {
        return parameterNames_[parameterID_];
    }

    //--------------------------------------------------------------------------

    String getValueString()
    {
        String strValue = "";

        switch (parameterID_)
        {
            case PARAM_TYPE:
                // Read-only.
                strValue = 
                    FoamXTypeEx.toDisplayString(typeDescriptor_.getType());
                break;
            case PARAM_PATH:
                // Read-only.
                strValue = typeDescriptor_.getPath();
                break;
            case PARAM_NAME:
                strValue = typeDescriptor_.getName();
                break;
            case PARAM_DISPLAYNAME:
                strValue = typeDescriptor_.getDisplayName();
                break;
            case PARAM_DESCRIPTION:
                strValue = typeDescriptor_.getDescription();
                break;
            case PARAM_COMMENT:
                strValue = typeDescriptor_.getComment();
                break;
            case PARAM_CATEGORY:
                strValue = typeDescriptor_.getCategory();
                break;
            case PARAM_HELPURL:
                strValue = typeDescriptor_.getHelpURL();
                break;
            case PARAM_ICONURL:
                strValue = typeDescriptor_.getIconURL();
                break;
            case PARAM_OPTIONAL:
                strValue = typeDescriptor_.isOptional() ? "True" : "False";
                break;
            case PARAM_VISIBLE:
                strValue = typeDescriptor_.isVisible() ? "True" : "False";
                break;
            case PARAM_EDITABLE:
                strValue = typeDescriptor_.isEditable() ? "True" : "False";
                break;
            case PARAM_MINVALUE:
                strValue = typeDescriptor_.getMinValueStr();
                break;
            case PARAM_MAXVALUE:
                strValue = typeDescriptor_.getMaxValueStr();
                break;
            case PARAM_DEFAULTVALUE:
                strValue = typeDescriptor_.getDefaultValueStr();
                break;
            case PARAM_LOOKUPDICT:
                strValue = typeDescriptor_.getLookupDict();
                break;
            case PARAM_VALUELIST:
                // Return a string containing an indication of the value list.
                strValue = listDisplayString(typeDescriptor_.getValueList(), 3);
                break;
            case PARAM_DICTIONARYPATH:
                strValue = typeDescriptor_.getDictionaryPath();
                break;
            case PARAM_NUMELEMENTS:
                strValue = Integer.toString(typeDescriptor_.getNumElements());
                break;
            case PARAM_ELEMENTLABELS:
                // Return a string containing an indication of the
                // element label list.
                strValue = listDisplayString
                (
                    typeDescriptor_.getElementLabels(), 3
                );
                break;
            default:
                break;
        }

        return strValue;
    }

    //--------------------------------------------------------------------------

    boolean updateValue(String newValue)
    {
        boolean bRet = true;

        try
        {
            switch (parameterID_)
            {
                case PARAM_NAME:
                    typeDescriptor_.setName(newValue);
                    break;

                case PARAM_DISPLAYNAME:
                    typeDescriptor_.setDisplayName(newValue);
                    break;

                case PARAM_DESCRIPTION:
                    typeDescriptor_.setDescription(newValue);
                    break;

                case PARAM_COMMENT:
                    typeDescriptor_.setComment(newValue);
                    break;

                case PARAM_CATEGORY:
                    typeDescriptor_.setCategory(newValue);
                    break;

                case PARAM_HELPURL:
                    typeDescriptor_.setHelpURL(newValue);
                    break;

                case PARAM_ICONURL:
                    typeDescriptor_.setIconURL(newValue);
                    break;

                case PARAM_OPTIONAL:
                    boolean opt = newValue.equalsIgnoreCase("True");
                    typeDescriptor_.setOptional(opt);
                    break;

                case PARAM_VISIBLE:
                    boolean vis = newValue.equalsIgnoreCase("True");
                    typeDescriptor_.setVisible(vis);
                    break;

                case PARAM_EDITABLE:
                    boolean edit = newValue.equalsIgnoreCase("True");
                    typeDescriptor_.setEditable(edit);
                    break;

                case PARAM_MINVALUE:
                    typeDescriptor_.setMinValue(newValue);
                    break;

                case PARAM_MAXVALUE:
                    typeDescriptor_.setMaxValue(newValue);
                    break;

                case PARAM_DEFAULTVALUE:
                    // This only does something if defaultValue is a primitive
                    // Otherwise the update is handles by the custom editor
                    typeDescriptor_.getDefaultValue().updateValue(newValue);
                    break;

                case PARAM_LOOKUPDICT:
                    typeDescriptor_.setLookupDict(newValue);
                    break;

                case PARAM_DICTIONARYPATH:
                    typeDescriptor_.setDictionaryPath(newValue);
                    break;

                case PARAM_NUMELEMENTS:
                    typeDescriptor_.setNumElements(Integer.parseInt(newValue));
                    break;

                case PARAM_TYPE:            // Read only parameters.
                case PARAM_VALUELIST:       // Custom edit parameter.
                case PARAM_ELEMENTLABELS:   // Custom edit parameter.
                default:
                    break;
            }
        }
        catch (IllegalArgumentException ex)
        {
            // Aaaarrrgghhh!
            Toolkit.getDefaultToolkit().beep();
            //App.handleAllExceptions(ex);
            bRet = false;
        }

        return bRet;
    }

    //--------------------------------------------------------------------------

    public Object getValue()
    {
        return typeDescriptor_.getDefaultValue();
    }

    TableCellEditor entryEditor()
    {
        return entryEditor_;
    }

    TableCellEditor customEditor()
    {
        return customEditor_;
    }

    //--------------------------------------------------------------------------

    boolean isEditable()
    {
        // Type and key parameters not editable.
        return editable_;
    }

    //--------------------------------------------------------------------------

    void showEditDialog()
    {
        try
        {
            switch (parameterID_)
            {
                case PARAM_DEFAULTVALUE:
                {
                    
                }
                break;

                case PARAM_VALUELIST:
                {
                    // Show string list editor.
                    java.util.Vector values = new java.util.Vector
                    (
                        typeDescriptor_.getValueList().length
                    );
                    for
                    (
                        int i=0;
                        i < typeDescriptor_.getValueList().length;
                        i++
                    )
                    {
                        values.addElement
                        (
                            typeDescriptor_.getValueList()[i]
                        );
                    }
                    StringListEditor dlg = new StringListEditor
                    (
                        App.getRootFrame(),
                        values
                    );
                    dlg.show();
                    if (!dlg.wasCancelled())
                    {
                        String[] valArray = new String[values.size()];
                        for (int i = 0; i <values.size(); i++)
                        {
                            valArray[i] = (String)values.elementAt(i);
                        }
                        typeDescriptor_.setValueList(valArray);
                    }
                }
                break;

                case PARAM_ELEMENTLABELS:
                {
                    // Show string list editor.
                    java.util.Vector values = new java.util.Vector
                    (
                        typeDescriptor_.getElementLabels().length
                    );
                    for
                    (
                        int i=0;
                        i < typeDescriptor_.getElementLabels().length;
                        i++
                    )
                    {
                        values.addElement
                        (
                            typeDescriptor_.getElementLabels()[i]
                        );
                    }
                    StringListEditor dlg = new StringListEditor
                    (
                        App.getRootFrame(),
                        values
                    );
                    dlg.show();
                    if (!dlg.wasCancelled())
                    {
                        String[] valArray = new String[values.size()];
                        for (int i = 0; i <values.size(); i++)
                        {
                            valArray[i] = (String)values.elementAt(i);
                        }
                        typeDescriptor_.setElementLabels(valArray);
                    }
                }
                break;

                default:
                break;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        // Editing has finished.
        customEditor_.stopCellEditing();
    }

    //--------------------------------------------------------------------------

    String listDisplayString(String[] list, int maxEntries)
    {
        String strValue = "";

        for (int i = 1; i <= list.length; i++)
        {
            strValue += list[i-1];
            if (i <list.length && i <maxEntries)
            {
                strValue += ", ";
            }
            else if (i == maxEntries)
            {
                strValue += " ...";
                break;
            }
        }

        return strValue;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



