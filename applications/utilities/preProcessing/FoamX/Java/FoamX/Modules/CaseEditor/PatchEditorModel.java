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
package FoamX.Modules.CaseEditor;

import javax.swing.tree.*;
import javax.swing.table.*;
import java.util.*;

import org.omg.CORBA.StringHolder;

import FoamX.App;
import FoamX.Editors.ApplicationEditor.BoundaryDefinitionModelItem;
import FoamX.Exceptions.FoamXException;

import FoamXServer.CaseServer.IPatchPhysicalTypeDescriptor;
import FoamXServer.CaseServer.IPatchPhysicalTypeDescriptorHolder;
import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.ITypeDescriptorHolder;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;

public class PatchEditorModel
    extends AbstractTableModel
{
    //--------------------------------------------------------------------------

    public  static final int      FIELDNAME_INDEX   = 0;
    public  static final int      PATCHFIELD_INDEX  = 1;
    public  static final String   FIELDNAME_COLUMN  = "Field";
    public  static final String   PATCHFIELD_COLUMN = "Patch Field Type";
    private static final String[] s_columnNames     =
    {
        FIELDNAME_COLUMN, PATCHFIELD_COLUMN
    };

    private String patchName_;
    private ICaseServer caseServer_;    // Reference to case object.
    private CaseEditorModule module_;   // Reference to parent module.
    private String[] fields_;           // List of field names.
    private int numFields_;             // Number of field in the table.
    // Map of field name to patch field type names.
    private java.util.Hashtable patchFieldMap_;

    // Tree model used to populate the boundary type selector dialog.
    private DefaultTreeModel treeModel_;
  // Tree root node.
    private DefaultMutableTreeNode rootNode_;
    // Map of DefaultMutableTreeNode objects, one for each boundary type.
    private java.util.Hashtable boundaryDefMap_;

    // Currently selected boundary type.
    private BoundaryDefinitionModelItem currentBoundaryDef_;

    //--------------------------------------------------------------------------

    public PatchEditorModel(ICaseServer caseServer, String patchName)
    {
        try
        {
            // Store patch name and reference to the case server.
            caseServer_    = caseServer;
            patchName_     = patchName;
            fields_        = caseServer_.application().fields();
            numFields_     = fields_.length;
            patchFieldMap_ = new Hashtable();

            // Initialise the boundary types tree strcuture.
            initialisePatchPhysicalTypes();

            // Initialise the patch parameter table.
            refreshPatchTable();
        }
        catch (Exception ex)
        {
             App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public DefaultTreeModel getTreeModel()
    {
        return treeModel_;
    }

    //--------------------------------------------------------------------------

    BoundaryDefinitionModelItem getBoundaryDefinition()
    {
        return currentBoundaryDef_;
    }

    //--------------------------------------------------------------------------

    void setBoundaryDefinition(BoundaryDefinitionModelItem boundaryDef)
    {
        try
        {
            // Store reference to the selected boundary type.
            currentBoundaryDef_ = boundaryDef;

            // Update the boundary type in the case server.
            caseServer_.setPatchPhysicalType
            (
                patchName_,
                currentBoundaryDef_.getName()
            );

            // Update the table.
            refreshPatchTable();
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
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
            // Initialise the tree model. We need to built the map of all boundary definitions
            // to make sure that all super types can be resolved before we build the tree.
            rootNode_ = new DefaultMutableTreeNode
            (
                new BoundaryDefinitionModelItem(),
                true
            );
            treeModel_          = new DefaultTreeModel(rootNode_);
            boundaryDefMap_     = new Hashtable();
            currentBoundaryDef_ = null;

            // Get defined boundary type names from the application class object.
            String[] patchPhysicalTypes =
                caseServer_.application().patchPhysicalTypes();
            for (int i=0; i <patchPhysicalTypes.length; i++)
            {
                // Get Boundary Type Descriptor object.
                IPatchPhysicalTypeDescriptorHolder bndTypeDescHolder =
                    new IPatchPhysicalTypeDescriptorHolder();
                caseServer_.application().getPatchPhysicalType
                (
                    patchPhysicalTypes[i],
                    bndTypeDescHolder
                );
                IPatchPhysicalTypeDescriptor bndTypeDesc = bndTypeDescHolder.value;

                // See if this boundary definition has been encountered before.
                if (!boundaryDefMap_.containsKey(bndTypeDesc.name()))
                {
                    // Create a new node.
                    BoundaryDefinitionModelItem nodeInfo =
                        new BoundaryDefinitionModelItem(bndTypeDesc);
                    DefaultMutableTreeNode newNode  =
                        new DefaultMutableTreeNode(nodeInfo, true);

                    // Add to map.
                    boundaryDefMap_.put(bndTypeDesc.name(), newNode);
                }
            }

            // Now that the map is initialised, build the tree.
            Enumeration iter = boundaryDefMap_.elements();
            while (iter.hasMoreElements())
            {
                DefaultMutableTreeNode node =
                    (DefaultMutableTreeNode)iter.nextElement();
                BoundaryDefinitionModelItem nodeInfo =
                    (BoundaryDefinitionModelItem)node.getUserObject();

                // Get the parent node.
                DefaultMutableTreeNode parentNode = null;
                if
                (
                    nodeInfo.getParentType().length() != 0
                 && boundaryDefMap_.containsKey(nodeInfo.getParentType())
                )
                {
                    // Sub-type boundary definition.
                    parentNode = (DefaultMutableTreeNode)boundaryDefMap_.get
                    (
                        nodeInfo.getParentType()
                    );
                }
                else
                {
                    // Top-level boundary definition.
                    parentNode = rootNode_;
                }

                // Add new node to tree.
                treeModel_.insertNodeInto(node, parentNode, 0);
            }

            // Select the specified boundary definition.
            StringHolder holder = new StringHolder();
            caseServer_.getPatchPhysicalType(patchName_, holder);
            String patchPhysicalType = holder.value;

            // Find the node corresponding to the specified boundary type.
            if (!boundaryDefMap_.containsKey(patchPhysicalType))
            {
                throw new FoamXException("Invalid boundary definition name.");
            }
            DefaultMutableTreeNode nodeItem =
                (DefaultMutableTreeNode)boundaryDefMap_.get(patchPhysicalType);
            currentBoundaryDef_ =
                (BoundaryDefinitionModelItem)nodeItem.getUserObject();
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
             App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void refreshPatchTable()
    {
        try
        {
            if (currentBoundaryDef_ == null)
            {
                // Reset each entry.
                for (int i = 0; i <numFields_; i++)
                {
                    patchFieldMap_.put(fields_[i], "-");
                }
            }
            else
            {
                // Loop over all fields.
                for (int i = 0; i <numFields_; i++)
                {
                    // Get the patch field type for this field.
                    String patchFieldType =
                        currentBoundaryDef_.getPatchFieldType(fields_[i]);

                    ITypeDescriptorHolder patchFieldDescHolder =
                        new ITypeDescriptorHolder();
                    caseServer_.foamProperties().getPatchFieldType
                    (
                        patchFieldType,
                        patchFieldDescHolder
                    );

                    // Set the patch field type string.
                    String patchFieldName =
                        patchFieldDescHolder.value.displayName();
                    patchFieldMap_.put(fields_[i], patchFieldName);
                }
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (Exception e)
        {
            App.handleException(e);
        }
    }

    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        return numFields_;
    }

    //--------------------------------------------------------------------------

    public int getColumnCount()
    {
        return s_columnNames.length;
    }

    //--------------------------------------------------------------------------

    public String getColumnName(int columnIndex)
    {
        return s_columnNames[columnIndex];
    }

    //--------------------------------------------------------------------------

    public Class getColumnClass(int columnIndex)
    {
        // All columns are of type String.
        return String.class;
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        return false;
    }

    //--------------------------------------------------------------------------

    public Object getValueAt(int rowIndex, int columnIndex)
    {
        Object obj = null;

        try
        {
            // Check for valid row index.
            if (rowIndex>= numFields_)
            {
                throw new Exception("Invalid row index");
            }

            if (columnIndex == FIELDNAME_INDEX) 
            {
                // Return a new string containing the field name.
                obj = fields_[rowIndex];
            }
            else if (columnIndex == PATCHFIELD_INDEX) 
            {
                obj = (String)patchFieldMap_.get(fields_[rowIndex]);
            }
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
        }
        catch (FoamXException fxEx)
        {
            App.handleException(fxEx);
        }
        catch (Exception ex)
        {
             App.handleAllExceptions(ex);
        }

        return obj;
    }

    //--------------------------------------------------------------------------

    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



