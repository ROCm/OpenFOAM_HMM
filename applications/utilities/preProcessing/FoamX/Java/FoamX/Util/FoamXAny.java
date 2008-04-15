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

import java.text.DecimalFormat;
import java.text.MessageFormat;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

import FoamXServer.DimensionSet;
import FoamXServer.DimensionSetHelper;
import FoamXServer.FoamXType;

public class FoamXAny
{
    //-------------------------------------------------------------------------

    // Reference to the managed Any object.
    private FoamXServer.FoamXAny fxAny_;

    //-------------------------------------------------------------------------
    /** FoamXAny constructor. */
    public FoamXAny(FoamXType type)
    {
        // Create an uninitialised any object to manage.
        fxAny_ = new FoamXServer.FoamXAny
        (
            type,
            org.omg.CORBA.ORB.init().create_any()
        );
    }

    //-------------------------------------------------------------------------
    /** FoamXAny constructor. */
    public FoamXAny(FoamXServer.FoamXAny any)
    {
        fxAny_ = new FoamXServer.FoamXAny(any.type, any.value);
    }

    //-------------------------------------------------------------------------
    /** Returns a string representation of the value with the Any object. */
    public String toString()
    {
        return anyToString();
    }

    //-------------------------------------------------------------------------
    /** Returns the raw Any object. */
    public FoamXServer.FoamXAny getAny()
    {
        return fxAny_;
    }

    //-------------------------------------------------------------------------
    /** Sets the value within the Any object by parsing the given string. */
    public void setValue(String val) throws IllegalArgumentException
    {
        try
        {
            // Don't bother if the string is empty.
            if (val == null || val.length() == 0)
            {
                // Reset the Any object to null state.
                fxAny_.value = org.omg.CORBA.ORB.init().create_any();
            }
            else
            {
                switch (fxAny_.type.value())
                {
                    case FoamXType._Type_Boolean:
                        try
                        {
                            Boolean bVal = Boolean.valueOf(val);
                            fxAny_.value.insert_boolean(bVal.booleanValue());
                        }
                        catch (NumberFormatException ex)
                        {
                            throw new IllegalArgumentException
                            (
                                val + " is not a boolean"
                            );
                        }
                        break;

                    case FoamXType._Type_Label:
                        try
                        {
                            int lVal = Integer.parseInt(val);
                            fxAny_.value.insert_long(lVal);
                        }
                        catch (NumberFormatException ex)
                        {
                            throw new IllegalArgumentException
                            (
                                val + " is not an integer"
                            );
                        }
                        break;

                    case FoamXType._Type_Scalar:
                        try
                        {
                            double dVal = Double.parseDouble(val);
                            fxAny_.value.insert_double(dVal);
                        }
                        catch (NumberFormatException ex)
                        {
                            throw new IllegalArgumentException
                            (
                                val + " is not a number"
                            );
                        }
                        break;

                    case FoamXType._Type_Char:
                        if (val.length() != 1)
                        {
                            throw new IllegalArgumentException
                            (
                                val + " is not a single character"
                            );
                        }
                        else
                        {
                            fxAny_.value.insert_string(val);
                        }
                        break;

                    case FoamXType._Type_Word:
                    case FoamXType._Type_String:

                    case FoamXType._Type_RootDir:
                    case FoamXType._Type_RootAndCase:
                    case FoamXType._Type_CaseName:
                    case FoamXType._Type_HostName:
                    case FoamXType._Type_File:
                    case FoamXType._Type_Directory:
                    case FoamXType._Type_Time:

                        fxAny_.value.insert_string(val);
                        break;

                    case FoamXType._Type_DimensionSet:
                        throw new IllegalArgumentException
                        (
                            "Cannot create a dimensionSet in FoamXAny.setValue."
                        );

                    default:
                        throw new IllegalArgumentException
                        (
                            "Invalid type " + fxAny_.type.value()
                          + " in FoamXAny.setValue."
                        );
                }
            }
        }
        catch (org.omg.CORBA.BAD_OPERATION ex)
        {
            App.handleAllExceptions(ex);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //-------------------------------------------------------------------------
    /** Static formatters */
    static private java.text.DecimalFormat intFormat_;
    static private java.text.DecimalFormat fixedFormat_;
    static private java.text.DecimalFormat expFormat_;

    //-------------------------------------------------------------------------
    /** Initialisation for the static formatters */
    static
    {
        // Get and set decimal notation format object for labels and integers.
        intFormat_ = new java.text.DecimalFormat();
        intFormat_.setGroupingUsed(false);
        intFormat_.setMaximumFractionDigits(6);

        // Get and set decimal notation format object for fixed-format doubles
        fixedFormat_ = new java.text.DecimalFormat("####.#####");
        fixedFormat_.setGroupingUsed(false);
        fixedFormat_.setMinimumFractionDigits(1);
        fixedFormat_.setMaximumFractionDigits(5);
        fixedFormat_.setMinimumIntegerDigits(1);
        fixedFormat_.setMaximumIntegerDigits(5);
        fixedFormat_.setDecimalSeparatorAlwaysShown(true);

        // Get and set decimal notation format object for
        // exponential-format doubles
        expFormat_ = new java.text.DecimalFormat("0.#####E0");
        expFormat_.setGroupingUsed(false);
        expFormat_.setMinimumFractionDigits(1);
        expFormat_.setMaximumFractionDigits(5);
        expFormat_.setMinimumIntegerDigits(1);
        expFormat_.setMaximumIntegerDigits(1);
        expFormat_.setDecimalSeparatorAlwaysShown(true);
    }

    //-------------------------------------------------------------------------
    /** Format a long as a String */
    static public String format(long iVal)
    {
        return intFormat_.format(iVal);
    }

    //-------------------------------------------------------------------------
    /** Format a double as a String */
    static public String format(double sVal)
    {
        if (sVal == 0)
        {
            return "0.0";
        }
        else if (Math.abs(sVal) >= 1e3 || Math.abs(sVal) < 1e-3)
        {
            return expFormat_.format(sVal).replace('E', 'e');
        }
        else
        {
            return fixedFormat_.format(sVal);
        }
    }

    //-------------------------------------------------------------------------
    /** Format a DimensionSet as a String */
    static public String format(DimensionSet dimSet)
    {
        java.lang.Object[] args = new java.lang.Object[]
        {
            new Double(dimSet.mass),
            new Double(dimSet.length),
            new Double(dimSet.time),
            new Double(dimSet.temperature),
            new Double(dimSet.moles),
            new Double(dimSet.current),
            new Double(dimSet.luminousIntensity)
        };
        return MessageFormat.format("[{0} {1} {2} {3} {4} {5} {6}]", args);
    }


    //-------------------------------------------------------------------------
    /** Convert a Corba Any value into a nice string. */
    protected String anyToString()
    {
        String text = "";   // Return a blank string by default.

        try
        {
            // See if we have a null Any object.
            if
            (
                fxAny_.value.type().kind().value()
             != org.omg.CORBA.TCKind._tk_null
            )
            {
                switch (fxAny_.type.value())
                {
                    case FoamXType._Type_Boolean:
                        text = new Boolean
                        (
                            fxAny_.value.extract_boolean()
                        ).toString();
                        break;

                    case FoamXType._Type_Label:
                        text = format(fxAny_.value.extract_long());
                        break;
                    case FoamXType._Type_Scalar:
                        text = format(fxAny_.value.extract_double());
                        break;

                    case FoamXType._Type_Char:
                    case FoamXType._Type_Word:
                    case FoamXType._Type_String:

                    case FoamXType._Type_RootDir:
                    case FoamXType._Type_RootAndCase:
                    case FoamXType._Type_CaseName:
                    case FoamXType._Type_HostName:
                    case FoamXType._Type_File:
                    case FoamXType._Type_Directory:
                    case FoamXType._Type_Time:
                        text = fxAny_.value.extract_string();
                        break;

                    case FoamXType._Type_DimensionSet:
                        text = format(DimensionSetHelper.extract(fxAny_.value));
                        break;

                    default:
                        // Buggeration and damnation.
                        throw new Exception
                        (
                            "Invalid type in FoamXAny.anyToString."
                        );
                }
            }
        }
        catch (org.omg.CORBA.BAD_OPERATION ex)
        {
            App.handleAllExceptions(ex);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }

        return text;
    }
}


//-----------------------------------------------------------------------------
