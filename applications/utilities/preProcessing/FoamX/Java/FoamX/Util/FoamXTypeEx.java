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
package FoamX.Util;

import javax.swing.JComboBox;
import javax.swing.Icon;

import FoamXServer.FoamXType;

import FoamX.App;

public class FoamXTypeEx
    extends Object
{
    //--------------------------------------------------------------------------

    public static final int NUM_TYPES = 20;

    // An appropriate icon for each type.
    private static Icon[]           icons_      = new Icon[NUM_TYPES];
    // An appropriate description for each type.
    private static String[]         names_      = new String[NUM_TYPES];
    // Named pair types for combo boxes.
    private static FoamXNamedType[] namedTypes_ = new FoamXNamedType[NUM_TYPES];

    // Map between the integer type code and the index into the icon and name 
    // arrays.
    public static final int UNDEFINED_TYPE    = FoamXType._Type_Undefined;
    public static final int BOOLEAN_TYPE      = FoamXType._Type_Boolean;
    public static final int LABEL_TYPE        = FoamXType._Type_Label;
    public static final int SCALAR_TYPE       = FoamXType._Type_Scalar;
    public static final int CHAR_TYPE         = FoamXType._Type_Char;
    public static final int WORD_TYPE         = FoamXType._Type_Word;
    public static final int STRING_TYPE       = FoamXType._Type_String;

    public static final int ROOTDIR_TYPE      = FoamXType._Type_RootDir;
    public static final int ROOTANDCASE_TYPE  = FoamXType._Type_RootAndCase;
    public static final int CASENAME_TYPE     = FoamXType._Type_CaseName;
    public static final int HOSTNAME_TYPE     = FoamXType._Type_HostName;
    public static final int FILE_TYPE         = FoamXType._Type_File;
    public static final int DIRECTORY_TYPE    = FoamXType._Type_Directory;
    public static final int TIME_TYPE         = FoamXType._Type_Time;

    public static final int DIMENSIONSET_TYPE = FoamXType._Type_DimensionSet;
    public static final int VECTORSPACE_TYPE  = FoamXType._Type_FixedList;
    public static final int LIST_TYPE         = FoamXType._Type_List;
    public static final int DICTIONARY_TYPE   = FoamXType._Type_Dictionary;
    public static final int SELECTION_TYPE    = FoamXType._Type_Selection;
    public static final int COMPOUND_TYPE     = FoamXType._Type_Compound;

    //--------------------------------------------------------------------------

    static
    {
        try
        {
            // Load icons.
            icons_[UNDEFINED_TYPE]    = App.getResources().getIcon("UndefinedImage");
            icons_[BOOLEAN_TYPE]      = App.getResources().getIcon("BooleanImage");
            icons_[LABEL_TYPE]        = App.getResources().getIcon("LabelImage");
            icons_[SCALAR_TYPE]       = App.getResources().getIcon("ScalarImage");
            icons_[CHAR_TYPE]         = App.getResources().getIcon("CharImage");
            icons_[WORD_TYPE]         = App.getResources().getIcon("WordImage");
            icons_[STRING_TYPE]       = App.getResources().getIcon("StringImage");

            icons_[ROOTDIR_TYPE]      = App.getResources().getIcon("StringImage");
            icons_[ROOTANDCASE_TYPE]  = App.getResources().getIcon("StringImage");
            icons_[CASENAME_TYPE]     = App.getResources().getIcon("StringImage");
            icons_[HOSTNAME_TYPE]     = App.getResources().getIcon("StringImage");
            icons_[FILE_TYPE]         = App.getResources().getIcon("StringImage");
            icons_[DIRECTORY_TYPE]    = App.getResources().getIcon("StringImage");
            icons_[TIME_TYPE]         = App.getResources().getIcon("StringImage");

            icons_[DIMENSIONSET_TYPE] = App.getResources().getIcon("DimensionSetImage");
            icons_[VECTORSPACE_TYPE]  = App.getResources().getIcon("FixedListImage");
            icons_[LIST_TYPE]         = App.getResources().getIcon("ListImage");
            icons_[DICTIONARY_TYPE]   = App.getResources().getIcon("DictionaryImage");
            icons_[SELECTION_TYPE]    = App.getResources().getIcon("SelectionImage");
            icons_[COMPOUND_TYPE]     = App.getResources().getIcon("CompoundImage");

            // Load display names.
            names_[UNDEFINED_TYPE]    = "Undefined";
            names_[BOOLEAN_TYPE]      = "Boolean";
            names_[LABEL_TYPE]        = "Label";
            names_[SCALAR_TYPE]       = "Scalar";
            names_[CHAR_TYPE]         = "Char";
            names_[WORD_TYPE]         = "Word";
            names_[STRING_TYPE]       = "String";

            names_[ROOTDIR_TYPE]      = "RootDir";
            names_[ROOTANDCASE_TYPE]  = "RootAndCase";
            names_[CASENAME_TYPE]     = "CaseName";
            names_[HOSTNAME_TYPE]     = "HostName";
            names_[FILE_TYPE]         = "File";
            names_[DIRECTORY_TYPE]    = "Directory";
            names_[TIME_TYPE]         = "Time";

            names_[DIMENSIONSET_TYPE] = "DimensionSet";
            names_[VECTORSPACE_TYPE]  = "FixedList";
            names_[LIST_TYPE]         = "List";
            names_[DICTIONARY_TYPE]   = "Dictionary";
            names_[SELECTION_TYPE]    = "Selection";
            names_[COMPOUND_TYPE]     = "Compound";

            // Load name-type pair objects.
            namedTypes_[UNDEFINED_TYPE]    = new FoamXNamedType(FoamXType._Type_Undefined);
            namedTypes_[BOOLEAN_TYPE]      = new FoamXNamedType(FoamXType._Type_Boolean);
            namedTypes_[LABEL_TYPE]        = new FoamXNamedType(FoamXType._Type_Label);
            namedTypes_[SCALAR_TYPE]       = new FoamXNamedType(FoamXType._Type_Scalar);
            namedTypes_[CHAR_TYPE]         = new FoamXNamedType(FoamXType._Type_Char);
            namedTypes_[WORD_TYPE]         = new FoamXNamedType(FoamXType._Type_Word);
            namedTypes_[STRING_TYPE]       = new FoamXNamedType(FoamXType._Type_String);

            namedTypes_[ROOTDIR_TYPE]      = new FoamXNamedType(FoamXType._Type_RootDir);
            namedTypes_[ROOTANDCASE_TYPE]  = new FoamXNamedType(FoamXType._Type_RootAndCase);
            namedTypes_[CASENAME_TYPE]     = new FoamXNamedType(FoamXType._Type_CaseName);
            namedTypes_[HOSTNAME_TYPE]     = new FoamXNamedType(FoamXType._Type_HostName);
            namedTypes_[FILE_TYPE]         = new FoamXNamedType(FoamXType._Type_File);
            namedTypes_[DIRECTORY_TYPE]    = new FoamXNamedType(FoamXType._Type_Directory);
            namedTypes_[TIME_TYPE]         = new FoamXNamedType(FoamXType._Type_Time);

            namedTypes_[DIMENSIONSET_TYPE] = new FoamXNamedType(FoamXType._Type_DimensionSet);
            namedTypes_[VECTORSPACE_TYPE]  = new FoamXNamedType(FoamXType._Type_FixedList);
            namedTypes_[LIST_TYPE]         = new FoamXNamedType(FoamXType._Type_List);
            namedTypes_[DICTIONARY_TYPE]   = new FoamXNamedType(FoamXType._Type_Dictionary);
            namedTypes_[SELECTION_TYPE]    = new FoamXNamedType(FoamXType._Type_Selection);
            namedTypes_[COMPOUND_TYPE]     = new FoamXNamedType(FoamXType._Type_Compound);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public static void initialiseCombo(JComboBox combo)
    {
        for (int i=0; i <namedTypes_.length; i++)
        {
            combo.addItem(namedTypes_[i]);
        }
    }

    //--------------------------------------------------------------------------
    /** Convert the type to a nice displayable string. */
    public static String toDisplayString(FoamXType fxType)
    {
        return names_[fxType.value()];
    }

    //--------------------------------------------------------------------------
    /** Return an appropriate icon for this type. */
    public static Icon getIcon(FoamXType fxType)
    {
        return icons_[fxType.value()];
    }

    //--------------------------------------------------------------------------

    public static java.util.Enumeration EnumerateTypes()
    {
        return new TypeEnumerator();
    }

    //--------------------------------------------------------------------------
    /** Static member class to define a Name/Type-code pair. */
    public static class FoamXNamedType
    {
        private String name_;
        private int typeCode_;
        public FoamXNamedType(int typeCode)
        {
            // Use the FoamXTypeEx's toDisplayString method to convert type to
            // string.
            name_     = toDisplayString(FoamXType.from_int(typeCode));
            typeCode_ = typeCode;
        }
        public String toString()
    {
        return name_;
    }
        public String getName()
    {
        return name_;
    }
        public int getTypeCode()
    {
        return typeCode_;
    }
    }

    //--------------------------------------------------------------------------
    /** Static member class to define an enumerator for looping over all known
     *  types.
     *  Returns a Type object for each iteration.
     */
    public static class TypeEnumerator
    implements java.util.Enumeration
    {
        private int currType_ = FoamXTypeEx.UNDEFINED_TYPE;

        public boolean hasMoreElements()
        {
            return currType_ != FoamXTypeEx.COMPOUND_TYPE;
        }

        public Object nextElement()
        {
            return FoamXType.from_int(++currType_);
        }
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


