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

import FoamX.App;
import FoamXServer.CaseServer.IGeometricFieldDescriptor;
import FoamXServer.DimensionSet;

public class FieldDescriptorCache
    extends Object
{
    //--------------------------------------------------------------------------

    private IGeometricFieldDescriptor fieldDescriptor_;
    private String fieldName_;
    private String fieldDescription_;
    private String fieldTypeName_;
    private String geometryTypeName_;
    private DimensionSet fieldDimension_;

    //--------------------------------------------------------------------------
    /** FieldDescriptorCache constructor. */
    public FieldDescriptorCache(IGeometricFieldDescriptor fieldDescriptor)
    {
        try
        {
            fieldDescriptor_  = fieldDescriptor;
            fieldName_        = fieldDescriptor_.name();
            fieldDescription_ = fieldDescriptor_.description();
            fieldTypeName_    = fieldDescriptor_.fieldTypeName();        // True field type name.
            geometryTypeName_ = fieldDescriptor_.geometryTypeName();     // True geometry type name.
            fieldDimension_   = fieldDescriptor_.dimensions();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public String getFieldName()
    {
        return fieldName_;
    }
    public String getFieldDescription()
    {
        return fieldDescription_;
    }
    public String getFoamTypeName()
    {
        return fieldTypeName_;
    }
    public String getGeometryTypeName()
    {
        return geometryTypeName_;
    }
    public DimensionSet getFieldDimensions()
    {
        return fieldDimension_;
    }

    public void setFieldName(String fieldName)
    {
        fieldName_        = fieldName.trim();
    }
    public void setFieldDescription(String fieldDescription)
    {
        fieldDescription_ = fieldDescription.trim();
    }
    public void setFieldTypeName(String fieldTypeName)
    {
        fieldTypeName_    = fieldTypeName.trim();
    }
    public void setGeometryTypeName(String geometryTypeName)
    {
        geometryTypeName_ = geometryTypeName.trim();
    }
    public void setFieldDimensions(DimensionSet fieldDimension)
    {
        fieldDimension_   = fieldDimension;
    }

    //--------------------------------------------------------------------------
    /**
     * Update the cached IGeometricFieldDescriptor object.
     */
    public void updateFieldDescriptor()
    {
        try
        {
            // Update the IGeometricFieldDescriptor object.
            fieldDescriptor_.name(fieldName_);
            fieldDescriptor_.description(fieldDescription_);
            fieldDescriptor_.fieldTypeName(fieldTypeName_);
            fieldDescriptor_.geometryTypeName(geometryTypeName_);
            fieldDescriptor_.dimensions(fieldDimension_);
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


