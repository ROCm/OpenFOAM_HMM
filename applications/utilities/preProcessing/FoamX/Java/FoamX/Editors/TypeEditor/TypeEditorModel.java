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

import java.lang.ref.*;
import java.util.Hashtable;
import java.util.Enumeration;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;

import javax.swing.*;
import javax.swing.table.*;
import javax.swing.tree.*;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;

import FoamX.App;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.ApplicationEditor.ApplicationExeption;
import FoamX.Util.FoamXTreeRenderer;

import FoamXServer.CaseServer.IApplication;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.FoamXType;
import FoamXServer.ITypeDescriptor;
import FoamXServer.ITypeDescriptorHolder;

class TypeEditorModel
    extends AbstractTableModel
    implements java.awt.event.ItemListener
{
    //--------------------------------------------------------------------------

    private IApplication app_;
    private java.util.Hashtable dictionaryMap_;         // Map of application classes top-level dictionary objects.
    private TypeNodeInfo currentType_;
    private JComboBox addSubTypeCombo_;

    // Tree related objects.
    private WeakReference treeRef_;               // Weak reference to tree object. Required to prevent circular reference lock.
    private DefaultTreeModel treeModel_;
    private DefaultMutableTreeNode rootNode_;

    // Table related objects.
    private WeakReference tableRef_;               // Weak reference to table object. Required to prevent circular reference lock.
    private java.util.Vector typeParameters_;         // Vector of TypeParameter objects.
    private TypeParameterEditor typeParameterEditor_;
    public  static final int      PARAMETERNAME_INDEX  = 0;
    public  static final int      PARAMETERVALUE_INDEX = 1;
    public  static final String   PARAMETERNAME_TEXT   = "Parameter";
    public  static final String   PARAMETERVALUE_TEXT  = "Value";
    private static final String[] s_columnNames        = {
        PARAMETERNAME_TEXT, PARAMETERVALUE_TEXT };

    // Actions provided by this model.
    private AddDictionaryAction addDictionaryAction_;
    private DeleteDictionaryAction deleteDictionaryAction_;
    private AddSubTypeAction[]     addSubTypeActions_;    // One action for each sub-type.
    private DeleteSubTypeAction deleteSubTypeAction_;

    //--------------------------------------------------------------------------
    /** TypeEditorModel constructor. */
    public TypeEditorModel(IApplication app, JTree tree, JTable table)
    {
        try
        {
            app_      = app;
            dictionaryMap_ = new Hashtable();
            currentType_   = null;

            //-------------------------------------------------------------------------------------
            // Initialise tree objects.
            treeRef_      = new WeakReference(tree);
            rootNode_     = new DefaultMutableTreeNode(new TypeNodeInfo("Dictionaries"));
            treeModel_    = new DefaultTreeModel(rootNode_);

            // Set tree model and renderer.
            getTree().setModel(treeModel_);
            getTree().setCellRenderer(new FoamXTreeRenderer());

            // Set single selection mode for tree.
            getTree().getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);

            // Listen out for selection events.
            getTree().getSelectionModel().addTreeSelectionListener
            (
                new TreeSelectionListener()
                {
                    public void valueChanged(TreeSelectionEvent evt)
                    {
                        OnTreeSelectionChanged(evt);
                    }
                }
            );

            //-------------------------------------------------------------------------------------
            // Initialise table objects.
            tableRef_            = new WeakReference(table);
            typeParameters_      = new java.util.Vector(10);
            typeParameterEditor_ = new TypeParameterEditor();
            getTable().setModel(this);

            // Set cell renderer and editor for value column.
            getTable().getColumn(PARAMETERVALUE_TEXT).setCellRenderer(new TypeParameterRenderer());
            getTable().getColumn(PARAMETERVALUE_TEXT).setCellEditor(typeParameterEditor_);

            //-------------------------------------------------------------------------------------
            // Create action objects.
            addDictionaryAction_ = new AddDictionaryAction();
            addDictionaryAction_.setEnabled(false);

            deleteDictionaryAction_ = new DeleteDictionaryAction();
            deleteDictionaryAction_.setEnabled(false);

            deleteSubTypeAction_ = new DeleteSubTypeAction();
            deleteSubTypeAction_.setEnabled(false);

            addSubTypeActions_ = new AddSubTypeAction[FoamXTypeEx.NUM_TYPES];
            Enumeration enum = FoamXTypeEx.EnumerateTypes();
            while (enum.hasMoreElements())
            {
                FoamXType type  = (FoamXType)enum.nextElement();
                int       index = type.value();     // Use type code as index.
                addSubTypeActions_[index] = new AddSubTypeAction(type);
                addSubTypeActions_[index].setEnabled(false);
            }

            //-------------------------------------------------------------------------------------
            // Initialise the model by reading the application class dictionary definitions.
            initialiseDictionaries();
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
    public JTree getTree()
    {
        return (JTree)treeRef_.get();
    }
    public JTable getTable()
    {
        return (JTable)tableRef_.get();
    }

    //--------------------------------------------------------------------------

    void initialiseToolbar(JToolBar toolbar)
    {
        try
        {
            addToolbarButton(toolbar, addDictionaryAction_);
            addToolbarButton(toolbar, deleteDictionaryAction_);
            addToolbarButton(toolbar, deleteSubTypeAction_);

            // Initialise the add sub-type combo box.
            addSubTypeCombo_ = new JComboBox();
            addSubTypeCombo_.setFont(toolbar.getFont());      // Use same font as toolbar.
            addSubTypeCombo_.setToolTipText("Add sub-type");

            addSubTypeCombo_.addItem("Add Sub-Type...");
            Enumeration enum = FoamXTypeEx.EnumerateTypes();
            while (enum.hasMoreElements())
            {
                FoamXType type  = (FoamXType)enum.nextElement();
                int       index = type.value();     // Use type code as index.
                //addToolbarButton(toolbar, addSubTypeActions_[index]);
                addSubTypeCombo_.addItem(addSubTypeActions_[index].getValue(Action.SHORT_DESCRIPTION));
            }
            toolbar.add(addSubTypeCombo_);

            // Register for combo box change events.
            addSubTypeCombo_.addItemListener(this);
            addSubTypeCombo_.setEnabled(false);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void addToolbarButton(JToolBar toolBar, Action action)
    {
        // Add new button to the toolbar using the action object.
        javax.swing.JButton button = toolBar.add(action);

        // Set the tooltip text.
        button.setToolTipText((String)action.getValue(Action.SHORT_DESCRIPTION));
        button.setText("");
        button.setFont(toolBar.getFont());      // Use same font as toolbar.
    }

    //--------------------------------------------------------------------------

    public void itemStateChanged(final java.awt.event.ItemEvent evt)
    {
        int actionIndex = addSubTypeCombo_.getSelectedIndex();

        // Invoke the required action.
        if (actionIndex> 0)
        {
            addSubTypeActions_[actionIndex].actionPerformed(null);
        }

        // Reselect the title item.
        addSubTypeCombo_.setSelectedIndex(0);
    }

    //--------------------------------------------------------------------------

    void OnTreeSelectionChanged(TreeSelectionEvent evt)
    {
        try
        {
            // Cancel any pending edit.
            typeParameterEditor_.cancelCellEditing();

            // Get the node object for the current selection.
            if (evt.isAddedPath())    // Has the new node been added to the selection?
            {
                TreePath tp = evt.getPath();
                if (tp != null)
                {
                    DefaultMutableTreeNode nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    TypeNodeInfo           nodeInfo = (TypeNodeInfo)nodeItem.getUserObject();

                    // Update the action enabled/disabled status.
                    addDictionaryAction_.setEnabled(nodeInfo.isRootNode());
                    deleteDictionaryAction_.setEnabled(nodeInfo.isTopLevelType());
                    deleteSubTypeAction_.setEnabled(!nodeInfo.isRootNode() && !nodeInfo.isTopLevelType());

                    boolean canAddSubtype = !nodeInfo.isRootNode() && nodeInfo.isCompound();

                    // If this is a list or a fixedList type, then only one sub-type is allowed.
                    if (canAddSubtype &&
                         (nodeInfo.getTypeDescriptorCache().getType().value() == FoamXType._Type_List ||
                           nodeInfo.getTypeDescriptorCache().getType().value() == FoamXType._Type_FixedList))
                    {
                        canAddSubtype = (nodeInfo.getTypeDescriptorCache().getNumElements() == 0);
                    }

                    // Set the add sub type action status.
                    Enumeration enum = FoamXTypeEx.EnumerateTypes();
                    while (enum.hasMoreElements())
                    {
                        FoamXType type  = (FoamXType)enum.nextElement();
                        int       index = type.value();     // Use type code as index.
                        addSubTypeActions_[index].setEnabled(canAddSubtype);
                    }

                    // Set the combo box status.
                    addSubTypeCombo_.setEnabled(canAddSubtype);

                    // Commit the changes to the previous type, if applicable.
                    if (currentType_ != null && !currentType_.isRootNode())
                    {
                        currentType_.updateTypeDescriptor();
                    }

                    // Take a reference to the new type.
                    currentType_ = nodeInfo;

                    // Update the type parameter list.
                    updateParameters();
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    void initialiseContextMenu(DefaultMutableTreeNode typeNode, JPopupMenu contextMenu)
    {
        try
        {
            TypeNodeInfo nodeInfo = (TypeNodeInfo)typeNode.getUserObject();

            // Construct an appropriate context menu,
            contextMenu.removeAll();
            java.awt.Font font = contextMenu.getFont();

            if (nodeInfo.isRootNode())
            {
                contextMenu.add(addDictionaryAction_).setFont(font);
            }
            else if (nodeInfo.isCompound())
            {
                // Construct sub-menu for adding entries of a specific type.
                JMenu subMenu = new JMenu("Add Sub Type");
                contextMenu.add(subMenu);
                subMenu.setFont(contextMenu.getFont());
                Enumeration enum = FoamXTypeEx.EnumerateTypes();
                while (enum.hasMoreElements())
                {
                    FoamXType type  = (FoamXType)enum.nextElement();
                    int       index = type.value();     // Use type code as index.
                    subMenu.add(addSubTypeActions_[index]).setFont(font);
                }
            }

            if (nodeInfo.isTopLevelType())
            {
                contextMenu.add(deleteDictionaryAction_).setFont(font);
            }
            else if (!nodeInfo.isRootNode())
            {
                contextMenu.add(deleteSubTypeAction_).setFont(font);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private String getNewDictionaryKey()
    {
        String dictKey;
        int i = dictionaryMap_.size() + 1;

        do
        {
            dictKey = "Dictionary" + Integer.toString(i++);
        }
        while (dictionaryMap_.containsKey(dictKey));

        return dictKey;
    }

    //--------------------------------------------------------------------------

    private void initialiseDictionaries()
    {
        try
        {
            // Get list of dictionary keys from the application class and build the tree.
            String[] dictKeys = app_.dictionaries();
            for (int i=0; i <dictKeys.length; i++)
            {
                ITypeDescriptorHolder typeDescHolder = new ITypeDescriptorHolder();
                app_.getDictionary(dictKeys[i], typeDescHolder);

                // Add node for this top level dictionary.
                DefaultMutableTreeNode dictNode = addType(typeDescHolder.value);

                // Add to map.
                dictionaryMap_.put(dictKeys[i], dictNode);
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
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

    }

    //--------------------------------------------------------------------------

    private DefaultMutableTreeNode addType(ITypeDescriptor typeDesc)
    {
        DefaultMutableTreeNode typeNode = null;

        try
        {
            // Add a node for this top level type
            TypeNodeInfo nodeInfo = new TypeNodeInfo(typeDesc, true);
            typeNode = new DefaultMutableTreeNode(nodeInfo);
            treeModel_.insertNodeInto(typeNode, rootNode_, rootNode_.getChildCount());

            // Add sub types to this node if it has any.
            if (nodeInfo.isCompound())
            {
                addSubTypes(typeNode, typeDesc);
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return typeNode;
    }

    //--------------------------------------------------------------------------

    private void addSubTypes(DefaultMutableTreeNode parentNode, ITypeDescriptor typeDesc)
    {
        try
        {
            // Loop over all sub types.
            ITypeDescriptor[] subTypes = typeDesc.subTypes();
            for (int i=0; i <subTypes.length; i++)
            {
                // Get TypeDescriptor for this entry.
                ITypeDescriptor subType = subTypes[i];

                // Add node for this entry.
                TypeNodeInfo           nodeInfo = new TypeNodeInfo(subType, false);
                DefaultMutableTreeNode typeNode = new DefaultMutableTreeNode(nodeInfo);
                treeModel_.insertNodeInto(typeNode, parentNode, parentNode.getChildCount());

                // If this entry is another compound type, add its sub types to the tree.
                if (nodeInfo.isCompound())
                {
                    addSubTypes(typeNode, subType);
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void updateParameters()
    {
        try
        {
            typeParameters_.removeAllElements();

            if (currentType_ != null && !currentType_.isRootNode())
            {
                TypeDescriptorCache typeDesc =
                    currentType_.getTypeDescriptorCache();

                // Add non-editable type and key properties.
                typeParameters_.add
                (
                    new TypeParameter(typeDesc, TypeParameter.PARAM_TYPE)
                );

                typeParameters_.add
                (
                    new TypeParameter(typeDesc, TypeParameter.PARAM_PATH)
                );

                // Add top-level type properties.
                if (currentType_.isTopLevelType())
                {
                    typeParameters_.add
                    (
                        new TypeParameter
                        (
                            typeDesc,
                            TypeParameter.PARAM_DICTIONARYPATH
                        )
                    );
                }

                // Add common properties.
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_NAME));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_DISPLAYNAME));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_DESCRIPTION));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_COMMENT));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_CATEGORY));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_HELPURL));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_ICONURL));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_OPTIONAL));
                typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_VISIBLE));

                // Add non-compound properties.
                if (!currentType_.isCompound())
                {
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_EDITABLE));
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_DEFAULTVALUE));
                    if (currentType_.getTypeDescriptorCache().getType() == FoamXType.Type_Scalar ||
                         currentType_.getTypeDescriptorCache().getType() == FoamXType.Type_Label)
                    {
                        typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_MINVALUE));
                        typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_MAXVALUE));
                    }
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_VALUELIST));
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_LOOKUPDICT));
                }
                else if (currentType_.getTypeDescriptorCache().getType() == FoamXType.Type_FixedList)
                {
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_NUMELEMENTS));
                    typeParameters_.add(new TypeParameter(typeDesc, TypeParameter.PARAM_ELEMENTLABELS));
                }
            }

            // Make sure the table is refreshed.
            getTable().updateUI();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public void updateModel()
    {
        // Commit the changes to the previous type, if applicable.
        if (currentType_ != null && !currentType_.isRootNode())
        {
            currentType_.updateTypeDescriptor();
        }
    }

    //--------------------------------------------------------------------------

    void validate() throws ApplicationExeption
    {
        //throw new ApplicationExeption("Dictionaries", "No Fields Defined.");
    }

    //--------------------------------------------------------------------------
    //---- AbstractTableModel methods
    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        return typeParameters_.size();
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
        // Check for value column.
        if (columnIndex == PARAMETERVALUE_INDEX)
        {
            return TypeParameter.class;
        }

        // Everything else is of type String.
        return String.class;
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        // Check for value column.
        if (columnIndex != PARAMETERVALUE_INDEX) return false;

        // Check that the type parameter is editable..
        TypeParameter typeParam = (TypeParameter)typeParameters_.get(rowIndex);

        return typeParam.isEditable();
    }

    //--------------------------------------------------------------------------

    public Object getValueAt(int rowIndex, int columnIndex)
    {
        Object        value     = null;
        TypeParameter typeParam = (TypeParameter)typeParameters_.get(rowIndex);

        switch (columnIndex)
        {
            case PARAMETERNAME_INDEX:
                value = typeParam.getName();
                break;
            case PARAMETERVALUE_INDEX:
                value = typeParam;
                break;
            default:
                // Bugger and damnation.
                break;
        }

        return value;
    }

    //--------------------------------------------------------------------------

    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
        switch (columnIndex)
        {
            case PARAMETERVALUE_INDEX:
                //TypeParameter typeParam = (TypeParameter)typeParameters_.get(rowIndex);
                break;
            case PARAMETERNAME_INDEX:
            default:
                // Bugger.
                break;
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    private class TypeNodeInfo
    implements FoamX.Util.FoamXTreeRenderer.FoamXTreeItem
    {
        private TypeDescriptorCache typeDescCache_ = null;       // Type descriptor cache for this entry.
        private String              name_          = null;       // Name of this node if not a dictionary entry.
        private boolean             rootNode_      = false;      // Root node flag.
        private boolean             topLevelType_  = false;      // Top level type indicator.
        private boolean             compound_      = false;      // Compound type indicator.
        private Icon                icon_          = null;       // Appropriate icon for this node.

        //-----------------------------------------------------------------------------------------
        /** Constructor for a root DictionaryEntryNodeInfo object. */
        public TypeNodeInfo(String name)
        {
            name_     = name;
            rootNode_ = true;
            icon_     = App.getResources().getIcon("DictionaryImage");
        }

        //-----------------------------------------------------------------------------------------
        /** Constructor for DictionaryEntryNodeInfo object. */
        public TypeNodeInfo(TypeDescriptorCache typeDesc, boolean topLevelType)
        {
            typeDescCache_ = typeDesc;
            topLevelType_  = topLevelType;
            compound_      = typeDescCache_.isCompound();
            icon_          = typeDescCache_.getIcon();
        }

        //-----------------------------------------------------------------------------------------
        /** Constructor for DictionaryEntryNodeInfo object. */
        public TypeNodeInfo(ITypeDescriptor typeDesc, boolean topLevelType)
        {
            // Construct a new TypeDescriptorCache object to manage the raw ITypeDescriptor object.
            typeDescCache_ = new TypeDescriptorCache(typeDesc, false);
            topLevelType_  = topLevelType;
            compound_      = typeDescCache_.isCompound();
            icon_          = typeDescCache_.getIcon();
        }

        //-----------------------------------------------------------------------------------------

        public Icon getIcon()
    {
        return icon_;
    }
        public String getText()
    {
        return getName();
    }

        //-----------------------------------------------------------------------------------------

        public String toString()
    {
        return getName();
    }
        public TypeDescriptorCache getTypeDescriptorCache()
    {
        return typeDescCache_;
    }
        public boolean isRootNode()
    {
        return rootNode_;
    }
        public boolean isTopLevelType()
    {
        return topLevelType_;
    }
        public boolean isCompound()
    {
        return compound_;
    }

        //-----------------------------------------------------------------------------------------

        public String getName()
        {
            // If a name has been specified, return it.
            if (name_ != null) return name_;

            // Generate a suitable name to display.
            String name = typeDescCache_.getName();
            if (name.length() != 0) return name;

            // Use the display name.
            name = typeDescCache_.getDisplayName();
            if (name.length() != 0) return name;

            // Otherwise , use the type name.
            name = typeDescCache_.toString();
            return name;
        }

        //--------------------------------------------------------------------------
        /** Update the cached ITypeDescriptor object. */
        public void updateTypeDescriptor()
        {
            typeDescCache_.updateTypeDescriptor(false);      // Do not recurse.
        }
    }

    //--------------------------------------------------------------------------

    private class AddDictionaryAction
    extends AbstractAction
    {
           AddDictionaryAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("NewDictionaryImage"));
            putValue(Action.NAME, "Add Dictionary");
            putValue(Action.SHORT_DESCRIPTION, "Add Dictionary");
            putValue(Action.LONG_DESCRIPTION, "Add Dictionary");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Construct a default dictionary key.
                String dictKey = getNewDictionaryKey();
                if (dictKey != null && dictKey.length()> 0 && !dictionaryMap_.containsKey(dictKey))
                {
                    // Create a new top level dictionary type descriptor object.
                    ITypeDescriptorHolder dictTypeDescHolder = new ITypeDescriptorHolder();
                    app_.addDictionary(dictKey, dictTypeDescHolder);

                    // Add new node to model.
                    DefaultMutableTreeNode nodeItem = addType(dictTypeDescHolder.value);

                    // Add to map.
                    dictionaryMap_.put(dictKey, nodeItem);

                    // If addition successful, select the new node.
                    TreePath childTP = new TreePath(nodeItem.getPath());
                    getTree().expandPath(childTP);
                    getTree().getSelectionModel().clearSelection();
                    getTree().getSelectionModel().setSelectionPath(childTP);
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
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    private class DeleteDictionaryAction
    extends AbstractAction
    {
        DeleteDictionaryAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("DelDictionaryImage"));
            putValue(Action.NAME, "Delete Dictionary");
            putValue(Action.SHORT_DESCRIPTION, "Delete Dictionary");
            putValue(Action.LONG_DESCRIPTION, "Delete Dictionary");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = getTree().getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    DefaultMutableTreeNode nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    TypeNodeInfo           nodeInfo = (TypeNodeInfo)nodeItem.getUserObject();

                    // Make sure this node is a top-level dictionary and is in the map.
                    if (nodeInfo.isTopLevelType() && dictionaryMap_.containsKey(nodeInfo.getName()))
                    {
                        // Delete the top level dictionary type descriptor object.
                        app_.deleteDictionary(nodeInfo.getName());

                        // Remove from map.
                        dictionaryMap_.remove(nodeInfo.getName());

                        // Remove node from tree.
                        treeModel_.removeNodeFromParent(nodeItem);
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
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    private class AddSubTypeAction
    extends AbstractAction
    {
        private FoamXType dictEntryType_;

        AddSubTypeAction(FoamXType dictEntryType)
        {
            dictEntryType_ = dictEntryType;

            putValue(Action.SMALL_ICON, FoamXTypeEx.getIcon(dictEntryType_));
            putValue(Action.NAME, "Add " + FoamXTypeEx.toDisplayString(dictEntryType_));
            putValue(Action.SHORT_DESCRIPTION, "Add " + FoamXTypeEx.toDisplayString(dictEntryType_));
            putValue(Action.LONG_DESCRIPTION, "Add " + FoamXTypeEx.toDisplayString(dictEntryType_));
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = getTree().getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    TypeNodeInfo           nodeInfo = (TypeNodeInfo)nodeItem.getUserObject();

                    // Make sure this node is a compound dictionary entry.
                    if (!nodeInfo.isRootNode() && nodeInfo.isCompound())
                    {
                        // Add a new sub-type through the TypeDescriptorCache object,
                        TypeDescriptorCache parentCache = nodeInfo.getTypeDescriptorCache();
                        TypeDescriptorCache entryCache  = parentCache.addSubType(dictEntryType_);

                        // If successful, create node for this dictionary entry.
                        if (entryCache != null)
                        {
                            TypeNodeInfo           subNodeInfo = new TypeNodeInfo(entryCache, false);
                            DefaultMutableTreeNode subNodeItem = new DefaultMutableTreeNode(subNodeInfo);

                            // Add new node to tree.
                            treeModel_.insertNodeInto(subNodeItem, nodeItem, nodeItem.getChildCount());

                            // If addition successful, select the new node.
                            TreePath childTP = new TreePath(subNodeItem.getPath());
                            getTree().expandPath(childTP);
                            getTree().getSelectionModel().clearSelection();
                            getTree().getSelectionModel().setSelectionPath(childTP);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    private class DeleteSubTypeAction
    extends AbstractAction
    {
        DeleteSubTypeAction()
        {
            putValue(Action.SMALL_ICON, App.getResources().getIcon("DelDictionaryEntryImage"));
            putValue(Action.NAME, "Delete Sub Type");
            putValue(Action.SHORT_DESCRIPTION, "Delete Sub Type");
            putValue(Action.LONG_DESCRIPTION, "Delete Sub Type");
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                // Get currently selected node.
                TreePath tp = getTree().getSelectionModel().getSelectionPath();
                if (tp != null)
                {
                    // Get node info.
                    DefaultMutableTreeNode nodeItem   = (DefaultMutableTreeNode)tp.getLastPathComponent();
                    TypeNodeInfo           nodeInfo   = (TypeNodeInfo)nodeItem.getUserObject();
                    TypeDescriptorCache    entryCache = nodeInfo.getTypeDescriptorCache();

                    // Make sure this node is a dictionary entry.
                    if (!nodeInfo.isRootNode() && !nodeInfo.isTopLevelType())
                    {
                        // Get the parent node info.
                        DefaultMutableTreeNode parentNodeItem   = (DefaultMutableTreeNode)nodeItem.getParent();
                        TypeNodeInfo           parentNodeInfo   = (TypeNodeInfo)parentNodeItem.getUserObject();
                        TypeDescriptorCache    parentEntryCache = parentNodeInfo.getTypeDescriptorCache();

                        // Remove this sub-type from the parent TypeDescriptorCache.
                        if (parentEntryCache.removeSubType(entryCache))
                        {
                            // Remove node from tree.
                            treeModel_.removeNodeFromParent(nodeItem);

                            // If deletion successful, select the parent node.
                            getTree().getSelectionModel().clearSelection();
                            getTree().getSelectionModel().setSelectionPath(tp.getParentPath());
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}




