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

import java.text.DecimalFormat;
import java.util.Vector;
import java.net.URL;
import javax.swing.Icon;

import FoamXServer.FoamXError;
import FoamXServer.FoamXIOError;
import FoamXServer.FoamXType;
import FoamXServer.IDictionaryEntryHolder;
import FoamXServer.ITypeDescriptor;
import FoamXServer.ITypeDescriptorHolder;

import FoamX.App;
import FoamX.Exceptions.FoamXException;
import FoamX.Util.FoamXAny;
import FoamX.Util.FoamXTypeEx;
import FoamX.Editors.DictionaryEntryEditor.EntryCache.DictionaryEntryCache;

public class TypeDescriptorCache
{
    //--------------------------------------------------------------------------

    // Reference to the ITypeDescriptor object.
    private ITypeDescriptor typeDesc_;

    // Common properties.
    private FoamXType type_;       // The type of this parameter. Read only.
    private boolean compound_;     // Compound type flag. Read only.
    private String path_;          // Path of this parameter relative to root
                                    // dictionary. Read only.

    private String name_;          // Foam parameter name.
    private String displayName_;   // Short display name.
    private String description_;   // Long descriptive name.
    private String comment_;       // Dictionary comment.
    private String category_;      // GUI category.
    private String helpURL_;       // URL to relevant help file.
    private String iconURL_;       // URL to a nice icon for this entry.
    private boolean optional_;     // Flag indicating whether this parameter
                                    // is optional. Defaults to false.
    private boolean visible_;      // Flag indicating whether this parameter is
                                    // visible within FoamX. Defaults to true.

    // Properties for non-compound types.
    private boolean editable_;     // Flag indicating whether this parameter
                                    // is editable within FoamX. (Default:true)
    private FoamXAny minValue_;    // Minimum value for this parameter.
    private FoamXAny maxValue_;    // Maximum value for this parameter.

    
    private DictionaryEntryCache defaultValue_; // Default value for this
                                                 // parameter.

    private String lookupDict_;    // Name of dictionary to use as a lookup
                                    // table.
    private Vector valueList_;     // List of permissible values for this
                                    // parameter.

    // Properties for compound types.
    private String dictionaryPath_;    // Top-level dictionary location
                                        // relative to case root.
    private int numElements_;          // Number of elements in vector space
                                        // or list.
    private String[] elementLabels_;   // Labels of list or vector space
                                        // elements (eg, x, y, z).
    private Vector subTypes_;          // List of sub-types for compound types.

    //--------------------------------------------------------------------------
    /** Construct a TypeDescriptorCache from a ITypeDescriptor reference. */
    public TypeDescriptorCache(ITypeDescriptor typeDesc, boolean addSubTypes)
    {
        try
        {
            // Store ITypeDescriptor reference.
            typeDesc_ = typeDesc;

            // Get type attributes from ITypeDescriptor object.
            type_        = typeDesc_.type();
            compound_    = typeDesc_.isCompoundType();
            path_        = typeDesc_.path();

            name_        = typeDesc_.name();
            displayName_ = typeDesc_.displayName();
            description_ = typeDesc_.description();
            comment_     = typeDesc_.comment();
            category_    = typeDesc_.category();
            helpURL_     = typeDesc_.helpURL();
            iconURL_     = typeDesc_.iconURL();
            optional_    = typeDesc_.optional();
            visible_     = typeDesc_.visible();
            editable_    = typeDesc_.editable();

            if (!compound_)
            {
                minValue_ = new FoamXAny(typeDesc_.minValue());
                maxValue_ = new FoamXAny(typeDesc_.maxValue());

                IDictionaryEntryHolder dictEntryHolder =
                    new IDictionaryEntryHolder();
                typeDesc_.getDefaultValue(dictEntryHolder);

                lookupDict_ = typeDesc_.lookupDict();

                // Construct constrained value list.
                FoamXServer.FoamXAny[] valueList = typeDesc_.valueList();
                valueList_ = new Vector(valueList.length);
                for (int i = 0; i <valueList.length; i++)
                {
                    // Wrap each Any object in a FoamXAny object.
                    valueList_.add(i, new FoamXAny(valueList[i]));
                }

                defaultValue_ =
                    DictionaryEntryCache.New(dictEntryHolder.value, this);
            }
            else
            {
                dictionaryPath_ = typeDesc_.dictionaryPath();
                numElements_    = typeDesc_.numElements();
                elementLabels_  = typeDesc_.elementLabels();
                subTypes_       = new Vector(10);

                // Recursively add all of the sub-types if required.
                if (addSubTypes)
                {
                    ITypeDescriptor[] subTypes = typeDesc.subTypes();
                    subTypes_.setSize(subTypes.length);
                    for (int i = 0; i <subTypes.length; i++)
                    {
                        TypeDescriptorCache subType =
                            new TypeDescriptorCache(subTypes[i], true);
                        subTypes_.set(i, subType);
                    }
                }
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    // Access to the underlying type descriptor object.
    public ITypeDescriptor TypeDescriptor()
    {
        return typeDesc_;
    }

    // Get methods.
    public FoamXType getType()
    {
        return type_;
    }
    public boolean isCompound()
    {
        return compound_;
    }
    public String getPath()
    {
        return path_;
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
    public String getComment()
    {
        return comment_;
    }
    public String getCategory()
    {
        return category_;
    }
    public String getHelpURL()
    {
        return helpURL_;
    }
    public String getIconURL()
    {
        return iconURL_;
    }
    public boolean isOptional()
    {
        return optional_;
    }
    public boolean isVisible()
    {
        return visible_;
    }

    public boolean isEditable()
    {
        return editable_;
    }
    public FoamXAny getMinValue()
    {
        return minValue_;
    }
    public FoamXAny getMaxValue()
    {
        return maxValue_;
    }
    public String getMinValueStr()
    {
        return minValue_.toString();
    }
    public String getMaxValueStr()
    {
        return maxValue_.toString();
    }
    public DictionaryEntryCache getDefaultValue()
    {
        return defaultValue_;
    }
    public String getDefaultValueStr()
    {
        return defaultValue_.toString();
    }
    public String getLookupDict()
    {
        return lookupDict_;
    }

    public String getDictionaryPath()
    {
        return dictionaryPath_;
    }
    public int getNumElements()
    {
        return numElements_;
    }
    public String[] getElementLabels()
    {
        return elementLabels_;
    }

    // Set methods.
    public void setName(String name)
        throws IllegalArgumentException
    {
        name_ = name;
    }
    public void setDisplayName(String displayName)
        throws IllegalArgumentException
    {
        displayName_ = displayName;
    }
    public void setDescription(String description)
       throws IllegalArgumentException
    {
        description_ = description;
    }
    public void setComment(String comment)
        throws IllegalArgumentException
    {
        comment_ = comment;
    }
    public void setCategory(String category)
        throws IllegalArgumentException
    {
        category_ = category;
    }
    public void setHelpURL(String helpURL)
        throws IllegalArgumentException
    {
        helpURL_ = helpURL;
    }
    public void setIconURL(String iconURL)
        throws IllegalArgumentException
    {
        iconURL_ = iconURL;
    }
    public void setOptional(boolean optional)
        throws IllegalArgumentException
    {
        optional_ = optional;
    }
    public void setVisible(boolean visible)
        throws IllegalArgumentException
    {
        visible_ = visible;
    }

    public void setEditable(boolean editable)
        throws IllegalArgumentException
    {
        editable_ = editable;
    }
    public void setMinValue(String val)
        throws IllegalArgumentException
    {
        minValue_.setValue(val);
    }
    public void setMaxValue(String val)
        throws IllegalArgumentException
    {
        maxValue_.setValue(val);
    }
    public void setDefaultValue(String val) 
        throws IllegalArgumentException 
    {
        defaultValue_.updateValue(val);
    }
    public void setLookupDict(String lookupDict)
        throws IllegalArgumentException
    {
        lookupDict_  = lookupDict;
    }

    public void setDictionaryPath(String path)
        throws IllegalArgumentException
    {
        dictionaryPath_ = path;
    }
    public void setNumElements(int num)
        throws IllegalArgumentException
    {
        numElements_    = num;
    }
    public void setElementLabels(String[] labels)
        throws IllegalArgumentException
    {
        elementLabels_  = labels;
    }

    //--------------------------------------------------------------------------

    public String[] getValueList()
    {
        // Convert any list into string array.
        String[] valueList = null;

        if (valueList_ != null)
        {
            valueList = new String[valueList_.size()];
            for (int i=0; i <valueList_.size(); i++)
            {
                FoamXAny fxAny = (FoamXAny)valueList_.get(i);
                valueList[i]   = new String(fxAny.toString());
            }
        }
        return valueList;
    }

    //--------------------------------------------------------------------------

    public void setValueList(String[] vals)
    {
        // Replace the any list.
        valueList_ = new Vector(vals.length);

        for (int i = 0; i <vals.length; i++)
        {
            // Ignore any invalid values.
            try
            {
                FoamXAny fxAny = new FoamXAny(type_);
                fxAny.setValue(vals[i]);
                valueList_.add(i, fxAny);
            }
            catch (IllegalArgumentException ex)
            {
                App.handleAllExceptions(ex);
            }
        }
    }

    //--------------------------------------------------------------------------

    public TypeDescriptorCache getSubType(int index)
    {
        TypeDescriptorCache subType = null;

        try
        {
            if (!compound_)
            {
                throw new FoamXException
                (
                    "No sub-types on non-compound type descriptors."
                );
            }

            if (index>= subTypes_.size())
            {
                throw new FoamXException("Invalid sub-type index.");
            }

            // Return the sub-type specified by the index.
            subType = (TypeDescriptorCache)subTypes_.get(index);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return subType;
    }

    //--------------------------------------------------------------------------

    public TypeDescriptorCache addSubType(FoamXType type)
    {
        TypeDescriptorCache subType = null;

        try
        {
            if (!compound_)
            {
                throw new FoamXException
                (
                    "No sub-types on non-compound type descriptors."
                );
            }

            // Create a new TypeDescriptor sub object.
            ITypeDescriptorHolder typeDescHolder = new ITypeDescriptorHolder();
            typeDesc_.addSubType(type, typeDescHolder);

            // If successful, create a new TypeDescriptorCache wrapper
            // object and add to collection.
            subType = new TypeDescriptorCache(typeDescHolder.value, false);
            subTypes_.addElement(subType);

            //? numElements_ = typeDesc_.numElements();
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

        return subType;
    }

    //--------------------------------------------------------------------------

    public boolean removeSubType(TypeDescriptorCache subType)
    {
        boolean bRet = false;

        try
        {
            if (!compound_)
            {
                throw new FoamXException
                (
                    "No sub-types on non-compound type descriptors."
                );
            }

            // Remove the specified sub-type.
            typeDesc_.removeSubType(subType.TypeDescriptor());

            // If successful, remove from sub-type list.
            subTypes_.remove(subType);

            //? numElements_ = typeDesc_.numElements();

            // Return success.
            bRet = true;
        }
        catch (FoamXIOError ioErr)
        {
            App.handleException(ioErr);
        }
        catch (FoamXError fxErr)
        {
            App.handleException(fxErr);
            bRet = false;
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
            bRet = false;
        }

        return bRet;
    }

    //--------------------------------------------------------------------------
    /** Return an appropriate string representation of this type. */
    public String toString()
    {
        return FoamXTypeEx.toDisplayString(type_);
    }

    //--------------------------------------------------------------------------
    /** Return an appropriate icon for the given typeDescriptor. */
    public static Icon getIcon(ITypeDescriptor typeDescriptor)
    {
        Icon icon;

        String iconURL = typeDescriptor.iconURL();

        if (iconURL.length() > 0)
        {
            URL url = typeDescriptor.getClass().getResource(iconURL);

            // Check that this Url is valid.
            if (url == null)
            {
                App.printMessage("Cannot resolve URL from " + iconURL);

                // Return default icon for this type.
                icon = FoamXTypeEx.getIcon(typeDescriptor.type());
            }
            else
            {
                // Return icon specified by the URL.
                icon = App.getResources().getIcon(url);
            }
        }
        else
        {
            // Return default icon for this type.
            icon = FoamXTypeEx.getIcon(typeDescriptor.type());
        }

        return icon;
    }

    //--------------------------------------------------------------------------
    /** Return an appropriate icon for this type. */
    public Icon getIcon()
    {
        Icon icon;

        if (iconURL_ != null && iconURL_.length() > 0)
        {
            URL url = getClass().getResource(iconURL_);

            // Check that this Url is valid.
            if (url == null)
            {
                App.printMessage("Cannot resolve URL from " + iconURL_);

                // Return default icon for this type.
                icon = FoamXTypeEx.getIcon(type_);
            }
            else
            {
                // Return icon specified by the URL.
                icon = App.getResources().getIcon(url);
            }
        }
        else
        {
            // Return default icon for this type.
            icon = FoamXTypeEx.getIcon(type_);
        }

        return icon;
    }

    //--------------------------------------------------------------------------
    /**
     * Update the cached ITypeDescriptor object.
     */
    public void updateTypeDescriptor(boolean updateSubTypes)
    {
        try
        {
            // Update the TypeDescriptor object
            typeDesc_.name(name_);
            typeDesc_.displayName(displayName_);
            typeDesc_.description(description_);
            typeDesc_.comment(comment_);
            typeDesc_.category(category_);
            typeDesc_.helpURL(helpURL_);
            typeDesc_.iconURL(iconURL_);
            typeDesc_.optional(optional_);
            typeDesc_.visible(visible_);
            typeDesc_.editable(editable_);

            if (!compound_)
            {
                typeDesc_.minValue(minValue_.getAny());
                typeDesc_.maxValue(maxValue_.getAny());
                //***HGW typeDesc_.setDefaultValue(defaultValue_);
                typeDesc_.lookupDict(lookupDict_);

                // Update the constrained value list.
                FoamXServer.FoamXAny[] valueList =
                    new FoamXServer.FoamXAny[valueList_.size()];
                for (int i = 0; i <valueList_.size(); i++)
                {
                    FoamXAny fxAny = (FoamXAny)valueList_.get(i);
                    valueList[i] = fxAny.getAny();
                }
                typeDesc_.valueList(valueList);
            }
            else
            {
                typeDesc_.dictionaryPath(dictionaryPath_);
                typeDesc_.numElements(numElements_);
                typeDesc_.elementLabels(elementLabels_);

                // Recursively update the sub-types if required.
                if (updateSubTypes)
                {
                    for (int i = 0; i <subTypes_.size(); i++)
                    {
                        TypeDescriptorCache subType =
                            (TypeDescriptorCache)subTypes_.get(i);
                        subType.updateTypeDescriptor(updateSubTypes);
                    }
                }
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
}


