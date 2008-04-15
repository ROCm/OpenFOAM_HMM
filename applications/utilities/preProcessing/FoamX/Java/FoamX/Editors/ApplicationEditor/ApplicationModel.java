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

import java.util.Vector;
import java.util.Hashtable;

import FoamX.App;
import FoamX.FoamXException.*;
import FoamXServer.CaseServer.IApplication;
import FoamXServer.CaseServer.IApplicationHolder;
import FoamXServer.CaseServer.IFoamProperties;
import FoamXServer.CaseServer.IGeometryDescriptor;
import FoamXServer.CaseServer.IGeometryDescriptorHolder;
import FoamXServer.CaseServer.IPatchDescriptor;
import FoamXServer.CaseServer.IPatchDescriptorHolder;
import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.ITypeDescriptor;
import FoamXServer.ITypeDescriptorHolder;
import FoamXServer.ValidationError;

class ApplicationModel
{
    //--------------------------------------------------------------------------

    private IFoamProperties foamSys_;
    private IApplication app_;
    private boolean modified_;

    // Simple application class properties.
    private String appName_;
    private String appDescription_;
    private String appExeName_;
    private Vector appParameters_;
    private Vector appModules_;

    private FieldDefinitionModel fieldModel_;
    private BoundaryDefinitionModel boundaryModel_;

    // Field, Geometry, Patch and PatchField type information,
    private Hashtable fieldTypeMap_;         // Map of ITypeDescriptor objects.
    private Hashtable geometryTypeMap_;      // Map of IGeometryDescriptor objects.
    private Hashtable patchTypeMap_;         // Map of IPatchDescriptor objects.
    private Hashtable patchFieldTypeMap_;    // Map of ITypeDescriptor objects.

    // Simple lists of global FoamX type names.
    private String[]                availableModules_;     // Available FoamX modules.
    private String[]                foamTypes_;           // Field type names.
    private String[]                geometryTypes_;        // Geometry type names.
    private String[]                patchTypes_;           // Defined patch types.
    private String[]                patchFieldTypes_;      // Defined patch field types.

    //--------------------------------------------------------------------------

    ApplicationModel(IFoamProperties foamSys)
    {
        try
        {
            // Extract FoamX global data.
            foamSys_ = foamSys;

            // Initialise available module list.
            availableModules_ = foamSys_.availableModules();

            // Initialise patch types - wall, patch, symmetry, etc.
            patchTypeMap_ = new java.util.Hashtable();
            patchTypes_   = foamSys_.patchTypes();
            for (int i = 0; i <patchTypes_.length; i++)
            {
                // Get descriptor for this patch type.
                IPatchDescriptorHolder holder = new IPatchDescriptorHolder();
                foamSys_.getPatchType(patchTypes_[i], holder);
                patchTypeMap_.put(patchTypes_[i], holder.value);
            }

            // Initialise patch field types - zeroGradient, fixedValue, etc.
            patchFieldTypeMap_ = new java.util.Hashtable();
            patchFieldTypes_   = foamSys_.patchFieldTypes();
            for (int i = 0; i <patchFieldTypes_.length; i++)
            {
                // Get descriptor for this patch field type.
                ITypeDescriptorHolder holder = new ITypeDescriptorHolder();
                foamSys_.getPatchFieldType(patchFieldTypes_[i], holder);
                patchFieldTypeMap_.put(patchFieldTypes_[i], holder.value);
            }

            // Field Types - Construct a map between the true type name and the display name.
            fieldTypeMap_ = new java.util.Hashtable();
            foamTypes_   = foamSys_.foamTypes();
            for (int i = 0; i <foamTypes_.length; i++)
            {
                // Get descriptor for this field type.
                ITypeDescriptorHolder holder = new ITypeDescriptorHolder();
                foamSys_.getFoamType(foamTypes_[i], holder);
                fieldTypeMap_.put(foamTypes_[i], holder.value);
            }

            // Geometry Types - Construct a map between the true geometry type name and the display name.
            geometryTypeMap_ = new java.util.Hashtable();
            geometryTypes_   = foamSys_.geometryTypes();
            for (int i = 0; i <geometryTypes_.length; i++)
            {
                // Get descriptor for this geometry type.
                IGeometryDescriptorHolder holder = new IGeometryDescriptorHolder();
                foamSys_.getGeometryType(geometryTypes_[i], holder);
                geometryTypeMap_.put(geometryTypes_[i], holder.value);
            }

            // Reset modified flag.
            modified_ = false;
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

    void loadApplication(String appName)
    {
        try
        {
            // Open the specified application class object.
            IApplicationHolder appHolder = new IApplicationHolder();
            foamSys_.getApplication(appName, appHolder);

            // Store reference to app class object.
            app_ = appHolder.value;

            // Extract application class details.
            // Use specified name rather than the one defined (Overrides cloned application class name).
            appName_   = appName;
            appDescription_ = app_.description();
            appExeName_     = app_.executableName();

            // Initialise the parameter list.
            String[] params = app_.parameters();
            appParameters_ = new Vector(params.length);
            for (int i = 0; i <params.length; i++)
            {
                appParameters_.add(params[i]);
            }

            // Initialise the modules list.
            String[] modules = app_.modules();
            appModules_     = new Vector(modules.length);
            for (int i = 0; i <modules.length; i++)
            {
                appModules_.add(modules[i]);
            }

            // Create and initialise the field model.
            fieldModel_ = new FieldDefinitionModel(this);

            // Create and initialise the boundary definition model.
            boundaryModel_ = new BoundaryDefinitionModel(this);

            // Reset modified flag.
            modified_ = false;
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

    void saveApplication()
    {
        try
        {
            // Make sure we have an application class object to save to.
            if (app_ == null)
            {
                throw new FoamXException("Invalid application class in ApplicationModel::saveApplication.");
            }

            // Update the application class details.
            app_.name(appName_);
            app_.description(appDescription_);
            app_.executableName(appExeName_);

            // Update the parameter list.
            String[] params = new String[appParameters_.size()];
            appParameters_.copyInto(params);
            app_.parameters(params);

            // Update the modules list.
            String[] modules = new String[appModules_.size()];
            appModules_.copyInto(modules);
            app_.modules(modules);

            // Update the field information.
            fieldModel_.updateApplication(app_);
            
            // Update the boundary types.
            boundaryModel_.updateApplication(app_);

            // Finally, save the application class details to file.
            app_.save();
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

    public void validate() throws ApplicationExeption
    {
        // Validate the general properties.
        if (appName_.length()   == 0) throw new ApplicationExeption("General", "Invalid application class name.");
        if (appDescription_.length() == 0) throw new ApplicationExeption("General", "Invalid application class description.");
        if (appExeName_.length()     == 0) throw new ApplicationExeption("General", "Invalid application class executable name.");

        // Validate the field definition model.
        fieldModel_.validate();
        
        // Validate the boundary definition model.
        boundaryModel_.validate();

        // Validate the dictionary type descriptors.
        try
        {
            app_.validate();
        }
        catch (ValidationError vErr)
        {
            throw new ApplicationExeption
            (
                "Dictionaries",
                vErr.errorMessage
            );
        }
        catch (FoamXError fxErr)
        {
            throw new ApplicationExeption
            (
                "Dictionaries",
                fxErr.errorMessage
            );
        }
        catch (FoamXIOError ioErr)
        {
            throw new ApplicationExeption
            (
                "Dictionaries",
                ioErr.errorMessage
            );
        }
    }

    //--------------------------------------------------------------------------

    boolean isModified()
    {
        return modified_;
    }

    //--------------------------------------------------------------------------

    IApplication getAppClassDescriptor()
    {
        return app_;
    }

    String getAppClassName()
    {
        return appName_;
    }
    void                    setAppClassName(String name)
    {
        appName_ = name.trim();
    }

    String getAppDescription()
    {
        return appDescription_;
    }
    void                    setAppDescription(String desc)
    {
        appDescription_ = desc.trim();
    }

    String getAppExeName()
    {
        return appExeName_;
    }
    void                    setAppExeName(String name)
    {
        appExeName_ = name.trim();
    }

    Vector getAppParameters()
    {
        return appParameters_;
    }
    Vector getAppModules()
    {
        return appModules_;
    }

    FieldDefinitionModel getFieldDefinitionModel()
    {
        return fieldModel_;
    }

    BoundaryDefinitionModel getBoundaryDefinitionModel()
    {
        return boundaryModel_;
    }

    String[]                getAvailableModules()
    {
        return availableModules_;
    }
    String[]                getFoamTypes()
    {
        return foamTypes_;
    }
    String[]                getGeometryTypes()
    {
        return geometryTypes_;
    }
    String[]                getPatchTypes()
    {
        return patchTypes_;
    }
    String[]                getPatchFieldTypes()
    {
        return patchFieldTypes_;
    }

    //--------------------------------------------------------------------------

    ITypeDescriptor getFoamType(String fieldTypeName)
    {
        return (ITypeDescriptor)fieldTypeMap_.get(fieldTypeName);
    }

    IGeometryDescriptor getGeometryType(String geometryTypeName)
    {
        return (IGeometryDescriptor)geometryTypeMap_.get(geometryTypeName);
    }

    IPatchDescriptor getPatchType(String patchTypeName)
    {
        return (IPatchDescriptor)patchTypeMap_.get(patchTypeName);
    }

    ITypeDescriptor getPatchFieldType(String patchFieldTypeName)
    {
        return (ITypeDescriptor)patchFieldTypeMap_.get(patchFieldTypeName);
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}



