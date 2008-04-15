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

import java.lang.ref.WeakReference;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
import javax.swing.table.*;
import javax.swing.event.*;

import FoamX.App;
import FoamXServer.CaseServer.IApplication;
import FoamXServer.CaseServer.IGeometricFieldDescriptor;
import FoamXServer.CaseServer.IGeometricFieldDescriptorHolder;
import FoamXServer.CaseServer.IGeometryDescriptor;
import FoamXServer.DimensionSet;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.ITypeDescriptor;

class FieldDefinitionModel
    extends AbstractTableModel
{
    //--------------------------------------------------------------------------

    public  static final int      FIELDNAME_INDEX  = 0;
    public  static final int      FIELDDESC_INDEX  = 1;
    public  static final int      FIELDTYPE_INDEX  = 2;
    public  static final int      GEOMTYPE_INDEX   = 3;
    public  static final int      DIMENSION_INDEX  = 4;
    public  static final String   FIELDNAME_COLUMN = "Name";
    public  static final String   FIELDDESC_COLUMN = "Description";
    public  static final String   FIELDTYPE_COLUMN = "Type";
    public  static final String   GEOMTYPE_COLUMN  = "Geometry Type";
    public  static final String   DIMENSION_COLUMN = "Dimension";
    private static final String[] s_columnNames    = {
        FIELDNAME_COLUMN, FIELDDESC_COLUMN, FIELDTYPE_COLUMN, GEOMTYPE_COLUMN, DIMENSION_COLUMN };

    private WeakReference appModel_;      // Weak reference to the application class model object.
    private EventListenerList listenerList_;
    private java.util.Vector fieldDefList_;       // List of FieldDescriptorCache objects.
    private java.util.Hashtable fieldDefMap_;        // Map of FieldDescriptorCache objects.
    private java.util.Hashtable fieldTypeNameMap_;   // Map between true and display field type names.
    private java.util.Hashtable geomTypeNameMap_;    // Map between true and display geometry type names.

    //--------------------------------------------------------------------------

    FieldDefinitionModel(ApplicationModel appModel)
    {
        try
        {
            appModel_    = new java.lang.ref.WeakReference(appModel);
            fieldDefList_     = new java.util.Vector(5);
            fieldDefMap_      = new java.util.Hashtable();
            fieldTypeNameMap_ = new java.util.Hashtable();
            geomTypeNameMap_  = new java.util.Hashtable();

            // Create listener list.
            listenerList_ = new EventListenerList();

            // Initialise the field model.
            String[] fieldList = getAppClassModel().getAppClassDescriptor().fields();
            for (int i = 0; i <fieldList.length; i++)
            {
                IGeometricFieldDescriptorHolder fieldDescHolder = new IGeometricFieldDescriptorHolder();
                getAppClassModel().getAppClassDescriptor().getField(fieldList[i], fieldDescHolder);
                addField(fieldDescHolder.value);
            }

            // Initialise the field type name map.
            String[] foamTypes = getAppClassModel().getFoamTypes();
            for (int i = 0; i <foamTypes.length; i++)
            {
                ITypeDescriptor fieldType = getAppClassModel().getFoamType(foamTypes[i]);
                fieldTypeNameMap_.put(foamTypes[i], fieldType.displayName());
                fieldTypeNameMap_.put(fieldType.displayName(), foamTypes[i]);
            }

            // Initialise the geometry type name map.
            String[] geometryTypes = getAppClassModel().getGeometryTypes();
            for (int i = 0; i <geometryTypes.length; i++)
            {
                IGeometryDescriptor geomType = getAppClassModel().getGeometryType(geometryTypes[i]);
                geomTypeNameMap_.put(geometryTypes[i], geomType.displayName());
                geomTypeNameMap_.put(geomType.displayName(), geometryTypes[i]);
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

    void addField(IGeometricFieldDescriptor fieldDescriptor)
    {
        // Create a new field descriptor cache object and add to model.
        FieldDescriptorCache fieldDesc = new FieldDescriptorCache(fieldDescriptor);

        // Add to map and listy things.
        fieldDefMap_.put(fieldDesc.getFieldName(), fieldDesc);
        fieldDefList_.add(fieldDesc);

        // Fire an event to update the table.
        fireTableRowsInserted(fieldDefList_.size() - 1, fieldDefList_.size() - 1);
        fireFieldAdded(fieldDesc.getFieldName());
    }

    //--------------------------------------------------------------------------

    void addNewField()
    {
        try
        {
            // Create a name for the new field.
            String fieldName = generateFieldName();

            // Create a new field descriptor object.
            IGeometricFieldDescriptorHolder fieldDescHolder = new IGeometricFieldDescriptorHolder();
            getAppClassModel().getAppClassDescriptor().addField(fieldName, fieldDescHolder);
            addField(fieldDescHolder.value);
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

    void removeField(int rowIndex)
    {
        try
        {
            // Remove field descriptor cache from application class, list and map.
            FieldDescriptorCache fieldDesc = (FieldDescriptorCache)fieldDefList_.get(rowIndex);
            getAppClassModel().getAppClassDescriptor().deleteField(fieldDesc.getFieldName());
            fieldDefMap_.remove(fieldDesc.getFieldName());
            fieldDefList_.remove(rowIndex);

            // Fire an event to update the table.
            fireTableRowsDeleted(rowIndex, rowIndex);
            fireFieldRemoved(fieldDesc.getFieldName());
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

    FieldDescriptorCache getFieldDescriptor(int index)
    {
        FieldDescriptorCache fieldDesc = (FieldDescriptorCache)fieldDefList_.get(index);
        return fieldDesc;
    }

    //--------------------------------------------------------------------------

    void validate() throws ApplicationExeption
    {
        // Make sure at least one field is defined.
        if (fieldDefList_.size() == 0)
        {
            throw new ApplicationExeption("Fields", "No Fields Defined.");
        }

        // Make sure all fields are properly defined.
        Enumeration enum = fieldDefList_.elements();
        while (enum.hasMoreElements())
        {
            FieldDescriptorCache fieldDesc = (FieldDescriptorCache)enum.nextElement();
            if (fieldDesc.getFieldName().length()        == 0) throw new ApplicationExeption("Fields", "Invalid Field Name.");
            if (fieldDesc.getFieldDescription().length() == 0) throw new ApplicationExeption("Fields", "Invalid Field Description.");
            if (fieldDesc.getFoamTypeName().length()    == 0) throw new ApplicationExeption("Fields", "Invalid Field FoamXType.");
            if (fieldDesc.getGeometryTypeName().length() == 0) throw new ApplicationExeption("Fields", "Invalid Field Geometry FoamXType.");
        }
    }

    //--------------------------------------------------------------------------

    void updateApplication(IApplication app)
    {
        try
        {
            // Update the field desriptor objects.
            for (int i = 0; i <fieldDefList_.size(); i++)
            {
                FieldDescriptorCache fieldDesc = (FieldDescriptorCache)fieldDefList_.get(i);
                fieldDesc.updateFieldDescriptor();
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private ApplicationModel getAppClassModel()
    {
        return (ApplicationModel)appModel_.get();
    }

    //--------------------------------------------------------------------------

    private String generateFieldName()
    {
        String fieldNameRoot = "Field";
        String fieldName;
        int i = fieldDefList_.size() + 1;
        do
        {
            fieldName = fieldNameRoot + Integer.toString(i++);
        }
        while (fieldDefMap_.containsKey(fieldName));

        return fieldName;
    }

    //--------------------------------------------------------------------------
    //---- FieldModelListener Methods
    //--------------------------------------------------------------------------

    void addFieldModelListener(FieldModelListener l)
    {
        listenerList_.add(FieldModelListener.class, l);
    }

    //--------------------------------------------------------------------------

    void removeFieldModelListener(FieldModelListener l)
    {
        listenerList_.remove(FieldModelListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireFieldAdded(String fieldName)
    {
        // Create event object.
        FieldModelEvent evt = new FieldModelEvent(this, fieldName);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FieldModelListener.class)
            {
                ((FieldModelListener)listeners[i+1]).fieldAdded(evt);
            }              
        }
    }

    //--------------------------------------------------------------------------

    protected void fireFieldRemoved(String fieldName)
    {
        // Create event object.
        FieldModelEvent evt = new FieldModelEvent(this, fieldName);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FieldModelListener.class)
            {
                ((FieldModelListener)listeners[i+1]).fieldRemoved(evt);
            }              
        }
    }

    //--------------------------------------------------------------------------

    protected void fireFieldRenamed(String newName, String oldName)
    {
        // Create event object.
        FieldModelEvent evt = new FieldModelEvent(this, newName, oldName);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FieldModelListener.class)
            {
                ((FieldModelListener)listeners[i+1]).fieldRenamed(evt);
            }              
        }
    }

    //--------------------------------------------------------------------------
    //---- AbstractTableModel methods
    //--------------------------------------------------------------------------

    public int getRowCount()
    {
        return fieldDefList_.size();
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
        // Check for dimension column.
        if (columnIndex == DIMENSION_INDEX)
        {
            return DimensionSet.class;
        }

        // Everything else is of type String.
        return String.class;
    }

    //--------------------------------------------------------------------------

    public boolean isCellEditable(int rowIndex, int columnIndex)
    {
        // All columns are editable
        return true;
    }

    //--------------------------------------------------------------------------

    public Object getValueAt(int rowIndex, int columnIndex)
    {
        FieldDescriptorCache fieldDesc = (FieldDescriptorCache)fieldDefList_.get(rowIndex);
        switch (columnIndex)
        {
            case FIELDNAME_INDEX:
                return fieldDesc.getFieldName();
            case FIELDDESC_INDEX:
                return fieldDesc.getFieldDescription();
            case FIELDTYPE_INDEX:
                // Convert from the true type name to the display type name.
                String typeName = (String)fieldTypeNameMap_.get(fieldDesc.getFoamTypeName());
                return typeName;
            case GEOMTYPE_INDEX:
                // Convert from the true geometry name to the display geometry name.
                String geomName = (String)geomTypeNameMap_.get(fieldDesc.getGeometryTypeName());
                return geomName;
            case DIMENSION_INDEX:
                return fieldDesc.getFieldDimensions();
            default:
                // Bugger and damnation.
                break;
        }

        return null;
    }

    //--------------------------------------------------------------------------

    public void setValueAt(Object aValue, int rowIndex, int columnIndex)
    {
        FieldDescriptorCache fieldDesc = (FieldDescriptorCache)fieldDefList_.get(rowIndex);
        switch (columnIndex)
        {
            case FIELDNAME_INDEX:
                String oldName = fieldDesc.getFieldName();
                String newName = (String)aValue;
                if (!oldName.equals(newName))
                {
                    fireFieldRenamed(newName, oldName);
                }
                fieldDesc.setFieldName((String)aValue);                    // Change name.
                fireTableDataChanged();
                break;
            case FIELDDESC_INDEX:
                fieldDesc.setFieldDescription((String)aValue);
                fireTableDataChanged();
                break;
            case FIELDTYPE_INDEX:
                // Convert from the display type name to the true type name.
                String typeName = (String)fieldTypeNameMap_.get((String)aValue);
                fieldDesc.setFieldTypeName(typeName);
                fireTableDataChanged();
                break;
            case GEOMTYPE_INDEX:
                // Convert from the display geometry name to the true geometry name.
                String geomName = (String)geomTypeNameMap_.get((String)aValue);
                fieldDesc.setGeometryTypeName(geomName);
                fireTableDataChanged();
                break;
            case DIMENSION_INDEX:
                fieldDesc.setFieldDimensions((DimensionSet)aValue);
                fireTableDataChanged();
                break;
            default:
                // Bugger.
                break;
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



