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

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import javax.swing.tree.*;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Enumeration;
import java.lang.ref.WeakReference;

import org.omg.CORBA.StringHolder;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamXServer.CaseServer.IApplication;
import FoamXServer.CaseServer.IPatchPhysicalTypeDescriptor;
import FoamXServer.CaseServer.IPatchPhysicalTypeDescriptorHolder;
import FoamXServer.CaseServer.IPatchDescriptor;
//import FoamXServer.CaseServer.IPatchFieldDescriptor;

class BoundaryDefinitionModel
    extends AbstractTableModel
    implements FieldModelListener
{
    //--------------------------------------------------------------------------

    // Table constants.
    public    static final int      FIELD_NAME_INDEX        = 0;
    public    static final int      PATCH_FIELD_TYPE_INDEX  = 1;
    public    static final String   FIELD_NAME_COLUMN       = "Field";
    public    static final String   PATCH_FIELD_TYPE_COLUMN = "Patch Field Type";
    protected static final String[] ColumnNames             = {
        FIELD_NAME_COLUMN, PATCH_FIELD_TYPE_COLUMN };

    private WeakReference appModel_;         // Weak reference to the application class model object.
    private java.util.Hashtable patchNodeMap_;          // Map of tree node objects for each defined patch type.
    private java.util.Hashtable boundaryDefMap_;        // Map of BoundaryDefinitionModelItem objects.
    private BoundaryDefinitionModelItem currentboundaryDef_;    // Current BoundaryDefinitionModelItem object.
    private DefaultTreeModel treeModel_;             // Tree model to support rendering of this model via a tree control.
    private DefaultMutableTreeNode rootNode_;              // Tree root node.
    private java.util.Hashtable patchFieldTypeNameMap_; // Map between patch true and display type names.

    //--------------------------------------------------------------------------

    BoundaryDefinitionModel(ApplicationModel appModel)
    {
        try
        {
            appModel_         = new java.lang.ref.WeakReference(appModel);
            patchNodeMap_          = new Hashtable();
            boundaryDefMap_        = new Hashtable();
            currentboundaryDef_    = null;
            patchFieldTypeNameMap_ = new Hashtable();

            // Initialise the tree model.
            rootNode_  = new DefaultMutableTreeNode(new BoundaryDefinitionModelItem(), true);
            treeModel_ = new DefaultTreeModel(rootNode_);

            // Initialise the patch field type name map.
            String[] patchFieldTypes = getAppClassModel().getPatchFieldTypes();
            for (int i = 0; i <patchFieldTypes.length; i++)
            {
                FoamXServer.ITypeDescriptor patchFieldType = getAppClassModel().getPatchFieldType(patchFieldTypes[i]);
                patchFieldTypeNameMap_.put(patchFieldTypes[i], patchFieldType.displayName());
                patchFieldTypeNameMap_.put(patchFieldType.displayName(), patchFieldTypes[i]);
            }

            // Initialise patch type nodes.
            initialisePatchTypes();

            // Initialise rest of boundary types.
            initialisePatchPhysicalTypes();

            // Listen out for changes to the field model.
            getAppClassModel().getFieldDefinitionModel().addFieldModelListener(this);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    DefaultTreeModel getTreeModel()
    {
        return treeModel_;
    }

    BoundaryDefinitionModelItem getCurrentBoundaryModelItem()
    {
        return currentboundaryDef_;
    }

    //--------------------------------------------------------------------------

    void fireBoundaryDefinitionDataChanged()
    {
        // Update the boundary type desriptor objects.
        Enumeration enum = boundaryDefMap_.elements();
        while (enum.hasMoreElements())
        {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode)enum.nextElement();
            treeModel_.nodeChanged(node);
        }

        fireTableDataChanged();
    }

    //--------------------------------------------------------------------------

    void select(DefaultMutableTreeNode node)
    {
        BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();
        currentboundaryDef_ = nodeInfo;

        // Since the selected boundary type has changed, the table's data has changed.
        fireTableDataChanged();
    }

    //--------------------------------------------------------------------------

    DefaultMutableTreeNode addBoundaryDefinitionNode(String patchPhysicalTypeName, DefaultMutableTreeNode parentNode)
    {
        DefaultMutableTreeNode newNode = null;

        try
        {
            if (!boundaryDefMap_.containsKey(patchPhysicalTypeName))
            {
                // Create a new Boundary Type Descriptor object.
                IPatchPhysicalTypeDescriptorHolder bndTypeDescHolder = new IPatchPhysicalTypeDescriptorHolder();
                getAppClassModel().getAppClassDescriptor().addPatchPhysicalType(patchPhysicalTypeName, bndTypeDescHolder);

                // Wrap it in a BoundaryDefinitionModelItem object and create node object.
                BoundaryDefinitionModelItem nodeInfo = new BoundaryDefinitionModelItem(bndTypeDescHolder.value);
                newNode = new DefaultMutableTreeNode(nodeInfo, true);

                // Add to map.
                boundaryDefMap_.put(patchPhysicalTypeName, newNode);

                // Get parent details.
                BoundaryDefinitionModelItem parentNodeInfo = (BoundaryDefinitionModelItem)parentNode.getUserObject();

                // Initialise the new object.
                nodeInfo.setDisplayName(patchPhysicalTypeName);
                nodeInfo.setDescription(patchPhysicalTypeName);
                nodeInfo.setPatchType(parentNodeInfo.getPatchType());

                // See if this is a new sub-type.
                if (parentNodeInfo.isBoundaryDefinition())
                {
                    nodeInfo.setParentType(parentNodeInfo.getName());
                }

                // Add fields.
                for (int i = 0; i <getAppClassModel().getFieldDefinitionModel().getRowCount(); i++)
                {
                    FieldDescriptorCache fieldDesc = getAppClassModel().getFieldDefinitionModel().getFieldDescriptor(i);

                    // Set default patch field type.
                    String[] patchFieldTypes = getAppClassModel().getPatchFieldTypes();
                    nodeInfo.addField(fieldDesc.getFieldName(), patchFieldTypes[0]);
                }

                // Add node to tree.
                treeModel_.insertNodeInto(newNode, parentNode, 0);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return newNode;
    }

    //--------------------------------------------------------------------------

    boolean deleteBoundaryDefinitionNode(DefaultMutableTreeNode nodeItem)
    {
        boolean bRet = false;

        try
        {
            // See if the given node has any children or is a patch type.
            BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)nodeItem.getUserObject();
            if (nodeItem.getChildCount() == 0 && nodeInfo.isBoundaryDefinition())
            {
                // Delete boundary type in the application class.
                getAppClassModel().getAppClassDescriptor().deletePatchPhysicalType(nodeInfo.getName());

                // Remove from map and tree.
                boundaryDefMap_.remove(nodeInfo.getName());
                treeModel_.removeNodeFromParent(nodeItem);
                bRet = true;
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return bRet;
    }

    //--------------------------------------------------------------------------

    void updateApplication(IApplication app)
    {
        try
        {
            // Update the boundary type desriptor objects.
            Enumeration enum = boundaryDefMap_.elements();
            while (enum.hasMoreElements())
            {
                DefaultMutableTreeNode      node     = (DefaultMutableTreeNode)enum.nextElement();
                BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();
                nodeInfo.updatePatchPhysicalTypeDescriptor();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    void validate() throws ApplicationExeption
    {
        // TODO : Validate boundary definitions.
    }

    //--------------------------------------------------------------------------

    private ApplicationModel getAppClassModel()
    {
        return (ApplicationModel)appModel_.get();
    }

    //--------------------------------------------------------------------------

    private void initialisePatchTypes()
    {
        try
        {
            String[] patchTypes = getAppClassModel().getPatchTypes();
            for (int i = 0; i <patchTypes.length; i++)
            {
                IPatchDescriptor       patchType = getAppClassModel().getPatchType(patchTypes[i]);
                BoundaryDefinitionModelItem nodeInfo  = new BoundaryDefinitionModelItem(patchType);
                DefaultMutableTreeNode      newNode   = new DefaultMutableTreeNode(nodeInfo, true);

                treeModel_.insertNodeInto(newNode, rootNode_, i);

                patchNodeMap_.put(patchTypes[i], newNode);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void initialisePatchPhysicalTypes()
    {
        try
        {
            // Get defined boundary type names from the application class object.
            String[] patchPhysicalTypes = getAppClassModel().getAppClassDescriptor().patchPhysicalTypes();
            for (int i=0; i <patchPhysicalTypes.length; i++)
            {
                // Get Boundary Type Descriptor object.
                IPatchPhysicalTypeDescriptorHolder bndTypeDescHolder = new IPatchPhysicalTypeDescriptorHolder();
                getAppClassModel().getAppClassDescriptor().getPatchPhysicalType(patchPhysicalTypes[i], bndTypeDescHolder);
                IPatchPhysicalTypeDescriptor bndTypeDesc = bndTypeDescHolder.value;

                // Make sure this boundary definition has not been encountered before.
                if (boundaryDefMap_.containsKey(bndTypeDesc.name()))
                {
                    throw new FoamXException("Duplicate boundary definition in BoundaryDefinitionModel::initialisePatchPhysicalTypes.");
                }

                // Wrap it in a BoundaryDefinitionModelItem object and create node object.
                BoundaryDefinitionModelItem nodeInfo = new BoundaryDefinitionModelItem(bndTypeDesc);
                DefaultMutableTreeNode      newNode  = new DefaultMutableTreeNode(nodeInfo, true);

                // Add to map.
                boundaryDefMap_.put(patchPhysicalTypes[i], newNode);
            }

            // Now that the map is initialised, build the tree.
            Enumeration enum = boundaryDefMap_.elements();
            while (enum.hasMoreElements())
            {
                DefaultMutableTreeNode      node     = (DefaultMutableTreeNode)enum.nextElement();
                BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();

                // Make sure this patch type is defined.
                if (!patchNodeMap_.containsKey(nodeInfo.getPatchType()))
                {
                    throw new FoamXException("Unknown patch type in BoundaryDefinitionModel::initialisePatchPhysicalTypes.");
                }

                // If we have a super type, make sure we have it in the map.
                if (nodeInfo.getParentType().length()> 0 && !boundaryDefMap_.containsKey(nodeInfo.getParentType()))
                {
                    throw new FoamXException("Invalid super type in boundary definition in BoundaryDefinitionModel::initialisePatchPhysicalTypes.");
                }

                // Get the parent node.
                DefaultMutableTreeNode parentNode = null;
                if (nodeInfo.getParentType().length()> 0)
                {
                    // Sub-type boundary definition.
                    parentNode = (DefaultMutableTreeNode)boundaryDefMap_.get(nodeInfo.getParentType());
                }
                else
                {
                    // Top-level boundary definition.
                    parentNode = (DefaultMutableTreeNode)patchNodeMap_.get(nodeInfo.getPatchType());
                }

                // Add new node to tree.
                treeModel_.insertNodeInto(node, parentNode, 0);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- FieldModelListener Interface
    //--------------------------------------------------------------------------

    public void fieldAdded(FieldModelEvent evt)
    {
        try
        {
            Enumeration enum = boundaryDefMap_.elements();
            while (enum.hasMoreElements())
            {
                DefaultMutableTreeNode      node     = (DefaultMutableTreeNode)enum.nextElement();
                BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();

                // Set default patch field type.
                String[] patchFieldTypes = getAppClassModel().getPatchFieldTypes();
                nodeInfo.addField(evt.getFieldName(), patchFieldTypes[0]);
            }
            // When the field table model changes, so does this one.
            fireTableDataChanged();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void fieldRemoved(FieldModelEvent evt)
    {
        try
        {
            Enumeration enum = boundaryDefMap_.elements();
            while (enum.hasMoreElements())
            {
                DefaultMutableTreeNode      node     = (DefaultMutableTreeNode)enum.nextElement();
                BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();
                nodeInfo.deleteField(evt.getFieldName());
            }
            // When the field table model changes, so does this one.
            fireTableDataChanged();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void fieldRenamed(FieldModelEvent evt)
    {
        try
        {
            Enumeration enum = boundaryDefMap_.elements();
            while (enum.hasMoreElements())
            {
                DefaultMutableTreeNode      node     = (DefaultMutableTreeNode)enum.nextElement();
                BoundaryDefinitionModelItem nodeInfo = (BoundaryDefinitionModelItem)node.getUserObject();
                //nodeInfo.renameField(evt.getOldFieldName(), evt.getFieldName());
            }
            // When the field table model changes, so does this one.
            fireTableDataChanged();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- AbstractTableModel methods
    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        // If we are not displaying a boundary definition, show no rows.
        if (currentboundaryDef_ != null && currentboundaryDef_.isBoundaryDefinition())
        {
            return getAppClassModel().getFieldDefinitionModel().getRowCount();
        }
        else
        {
            return 0;
        }
    }

    //--------------------------------------------------------------------------

    public int getColumnCount()
    {
        return ColumnNames.length;
    }

    //--------------------------------------------------------------------------

    public String getColumnName(int columnIndex)
    {
        return ColumnNames[columnIndex];
    }

    //--------------------------------------------------------------------------

    public Class getColumnClass(int columnIndex)
    {
        // Everything is of type String.
        return String.class;
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        // Patch field type column is editable.
        return (columnIndex == FIELD_NAME_INDEX) ? false : true;
    }

    //--------------------------------------------------------------------------

    public Object getValueAt(int rowIndex, int columnIndex)
    {
        // Get field details.
        FieldDescriptorCache fieldDesc = getAppClassModel().getFieldDefinitionModel().getFieldDescriptor(rowIndex);

        switch (columnIndex)
        {
            case FIELD_NAME_INDEX:
                return fieldDesc.getFieldName();
            case PATCH_FIELD_TYPE_INDEX:
                // Return the patch field type for this field.
                if (currentboundaryDef_ != null && currentboundaryDef_.isBoundaryDefinition())
                {
                    // Convert from true patch field type to display name.
                    String patchFieldType = currentboundaryDef_.getPatchFieldType(fieldDesc.getFieldName());
                    return (String)patchFieldTypeNameMap_.get(patchFieldType);
                }
                break;
            default:
                // Buggeration.
                break;
        }

        return null;
    }

    //--------------------------------------------------------------------------

    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
        // Get field details.
        FieldDescriptorCache fieldDesc = getAppClassModel().getFieldDefinitionModel().getFieldDescriptor(rowIndex);

        switch (columnIndex)
        {
            case FIELD_NAME_INDEX:
                break;
            case PATCH_FIELD_TYPE_INDEX:
                if (currentboundaryDef_ != null && currentboundaryDef_.isBoundaryDefinition())
                {
                    // Set the patch field type for the current BoundaryDefinition object.
                    String displayName = (String)aValue;

                    // Convert from display name to true patch field type.
                    String patchFieldType = (String)patchFieldTypeNameMap_.get(displayName);
                    currentboundaryDef_.setPatchFieldType(fieldDesc.getFieldName(), patchFieldType);
                }
                break;
            default:
                // Buggeration.
                break;
        }
    }
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



