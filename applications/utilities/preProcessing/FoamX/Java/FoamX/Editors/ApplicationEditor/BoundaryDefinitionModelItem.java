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
package FoamX.Editors.ApplicationEditor;

import javax.swing.Icon;
import java.util.Enumeration;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamXServer.CaseServer.IPatchPhysicalTypeDescriptor;
import FoamXServer.CaseServer.IPatchDescriptor;
import FoamXServer.StringPair;

public class BoundaryDefinitionModelItem
    implements FoamX.Util.FoamXTreeRenderer.FoamXTreeItem
{
    //--------------------------------------------------------------------------

    private static final int NODE_TYPE_ROOT        = 0;
    private static final int NODE_TYPE_PATCHTYPE   = 1;
    private static final int NODE_TYPE_BOUNDARYDEF = 2;
    private static Icon[] icons_;

    private IPatchPhysicalTypeDescriptor patchPhysicalTypeDesc_;
    private int type_;                 // Type of this node - Root, Patch or BoundaryDefinition.
    private String name_;
    private String displayName_;
    private String description_;
    private String parentType_;
    private String patchType_;
    private java.util.Hashtable patchFieldTypeMap_;    // Map between field names and patch field types.

    //--------------------------------------------------------------------------

    static
    {
        try
        {
            icons_ = new Icon[3];
            icons_[NODE_TYPE_ROOT]        = App.getResources().getIcon("PatchRootImage");
            icons_[NODE_TYPE_PATCHTYPE]   = App.getResources().getIcon("PatchTypeImage");
            icons_[NODE_TYPE_BOUNDARYDEF] = App.getResources().getIcon("BoundaryDefinitionImage");
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Constructor for Root BoundaryDefinitionModelItem object. */
    public BoundaryDefinitionModelItem()
    {
        type_        = NODE_TYPE_ROOT;
        name_        = "Patch Types";
        displayName_ = "Patch Types";
    }

    //--------------------------------------------------------------------------
    /** Constructor for Patch Type BoundaryDefinitionModelItem object. */
    public BoundaryDefinitionModelItem(IPatchDescriptor patchType)
    {
        try
        {
            type_        = NODE_TYPE_PATCHTYPE;
            name_        = patchType.name();
            displayName_ = patchType.displayName();
            description_ = patchType.description();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /** Constructor for Boundary Definition Type BoundaryDefinitionModelItem object. */
    public BoundaryDefinitionModelItem(IPatchPhysicalTypeDescriptor patchPhysicalType)
    {
        try
        {
            patchPhysicalTypeDesc_  = patchPhysicalType;
            type_              = NODE_TYPE_BOUNDARYDEF;
            name_              = patchPhysicalType.name();
            displayName_       = patchPhysicalType.displayName();
            description_       = patchPhysicalType.description();
            parentType_         = patchPhysicalType.parentType();
            patchType_         = patchPhysicalType.patchType();
            patchFieldTypeMap_ = new java.util.Hashtable();

            // Get the patch field types.
            StringPair[] patchFieldTypes = patchPhysicalType.patchFieldTypes();
            for (int j=0; j <patchFieldTypes.length; j++)
            {
                addField(patchFieldTypes[j].name, patchFieldTypes[j].value);
                //setPatchFieldType(patchFieldTypes[j].name, patchFieldTypes[j].value);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    // FoamXTreeItem interface methods.

    public Icon getIcon()
    {
        return icons_[type_];
    }
    public String getText()
    {
        return toString();
    }

    //--------------------------------------------------------------------------

    public String toString()
    {
        return displayName_;
    }

    public int getType()
    {
        return type_;
    }
    public String getName()
    {
        return name_;
    }
    public String getDisplayName()
    {
        return displayName_;
    }
    public String getDescription()
    {
        return description_;
    }
    public String getParentType()
    {
        return parentType_;
    }
    public String getPatchType()
    {
        return patchType_;
    }

    public void    setName(String name)
    {
        name_ = name.trim();
    }
    public void    setDisplayName(String displayName)
    {
        displayName_ = displayName.trim();
    }
    public void    setDescription(String description)
    {
        description_ = description.trim();
    }
    public void    setParentType(String parentType)
    {
        parentType_ = parentType.trim();
    }
    public void    setPatchType(String patchType)
    {
        patchType_ = patchType.trim();
    }

    public boolean isRoot()
    {
        return type_ == NODE_TYPE_ROOT;
    }
    public boolean isPatchType()
    {
        return type_ == NODE_TYPE_PATCHTYPE;
    }
    public boolean isBoundaryDefinition()
    {
        return type_ == NODE_TYPE_BOUNDARYDEF;
    }

    //--------------------------------------------------------------------------

    void addField(String fieldName, String patchFieldType)
    {
        try
        {
            if (patchFieldTypeMap_.containsKey(fieldName))
            {
                throw new FoamXException("Duplicate field name.");
            }
            patchFieldTypeMap_.put(fieldName, patchFieldType);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    void deleteField(String fieldName)
    {
        try
        {
            if (!patchFieldTypeMap_.containsKey(fieldName))
            {
                throw new FoamXException("Invalid field name.");
            }
            patchFieldTypeMap_.remove(fieldName);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    void renameField(String oldName, String newName)
    {
        try
        {
            if (!patchFieldTypeMap_.containsKey(oldName))
            {
                throw new FoamXException("Invalid field name.");
            }
            String patchFieldType = (String)patchFieldTypeMap_.get(oldName);
            patchFieldTypeMap_.remove(oldName);
            patchFieldTypeMap_.put(newName, patchFieldType);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public String getPatchFieldType(String fieldName)
    {
        try
        {
            if (!patchFieldTypeMap_.containsKey(fieldName))
            {
                throw new FoamXException("Invalid field name.");
            }
            return (String)patchFieldTypeMap_.get(fieldName);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
        return "";
    }

    //--------------------------------------------------------------------------

    public void setPatchFieldType(String fieldName, String patchFieldType)
    {
        try
        {
            if (!patchFieldTypeMap_.containsKey(fieldName))
            {
                throw new FoamXException("Invalid field name.");
            }
            patchFieldTypeMap_.put(fieldName, patchFieldType);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Update the cached IPatchPhysicalTypeDescriptor object.
     */
    public void updatePatchPhysicalTypeDescriptor()
    {
        try
        {
            if (patchPhysicalTypeDesc_ != null)
            {
                // Update the IPatchPhysicalTypeDescriptor object.
                patchPhysicalTypeDesc_.name(name_);
                patchPhysicalTypeDesc_.displayName(displayName_);
                patchPhysicalTypeDesc_.description(description_);
                patchPhysicalTypeDesc_.parentType(parentType_);
                patchPhysicalTypeDesc_.patchType(patchType_);

                // Set the patch field types.
                Enumeration enum = patchFieldTypeMap_.keys();
                StringPair[] patchFieldTypes = new StringPair[patchFieldTypeMap_.size()];
                int nField = 0;
                while (enum.hasMoreElements())
                {
                    String fieldName            = (String)enum.nextElement();
                    String patchFieldType       = (String)patchFieldTypeMap_.get(fieldName);
                    patchFieldTypes[nField++] = new StringPair(fieldName, patchFieldType);
                }
                patchPhysicalTypeDesc_.patchFieldTypes(patchFieldTypes);
            }
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
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



