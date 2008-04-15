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
package FoamX.Options;

import java.io.*;
import java.util.*;

import FoamX.App;

import FoamX.Exceptions.FoamXException;

public class OptionsManager
    extends java.util.Properties
{
    //--------------------------------------------------------------------------

    private String foamXSystemPath_;
    private String foamXSystemConfigPath_;
    private String foamXUserConfigPath_;

    //--------------------------------------------------------------------------
    /** OptionsManager constructor. */
    public OptionsManager()
    {
        try
        {
            // Get the FoamX config path.
            // To override default, use -DFoamX.ConfigPath=path on command line.
            // Would really like to just pick up the FOAMX_USER_CONFIG
            // environment variable. But, Java thinks that environment
            // variables are not portable and does not allow access to them!
            // Pile of crap that it is....

            // See if it has been overridden on the command line.
            foamXSystemPath_ = java.lang.System.getProperty
            (
                "FoamX.SystemPath", foamXSystemPath_
            );


            // Set the default path.
            foamXSystemConfigPath_ =
                foamXSystemPath_
                + "applications/utilities/preProcessing/FoamX/config";

            // See if it has been overridden on the command line.
            foamXSystemConfigPath_ = java.lang.System.getProperty
            (
                "FoamX.SystemConfigPath", foamXSystemConfigPath_
            );


            // Get user config path
            foamXUserConfigPath_ = java.lang.System.getProperty
            (
                "FoamX.UserConfigPath"
            );
            if (foamXUserConfigPath_ == null)
            {
                throw new FoamXException
                (
                    "No FoamX.UserConfigPath property specified in arguments."
                );
            }

            // Initialise the default properties.
            Properties defProps = new Properties();

            // TODO : Load default values from resource file.
            //defProps.load();
            defaults = defProps;

            // See if we have a custom properties file and load it if we have.
            loadFromFile();
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    public String foamXSystemPath()
    {
        return foamXSystemPath_;
    }

    //--------------------------------------------------------------------------

    public String foamXSystemConfigPath()
    {
        return foamXSystemConfigPath_;
    }

    //--------------------------------------------------------------------------

    public String foamXUserConfigPath()
    {
        return foamXUserConfigPath_;
    }

    //--------------------------------------------------------------------------

    public void loadFromFile()
    {
        try
        {
            // Load the properties object from file.
            String propFileName = foamXUserConfigPath_ + "/FoamXClient.cfg";
            File propFile = new File(propFileName);
            if (propFile.exists())
            {
                load(new FileInputStream(propFileName));

                // Add the config paths to the properties bundle.
                put("FoamX.SystemConfigPath", foamXSystemConfigPath_);
                put("FoamX.UserConfigPath", foamXUserConfigPath_);
            }
            else
            {
                throw new Exception
                (
                    "FoamX client properties file not found at '"
                  + propFileName + "'."
                );
            }
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

//    //--------------------------------------------------------------------------
//
//    public void saveToFile()
//    {
//        try
//        {
//            // Persist the properties object to file.
//            String propFileName = foamXUserConfigPath_ + "/FoamXClient.cfg";
//            store
//            (
//                new FileOutputStream(propFileName),
//                " FoamX Client Configuration."
//            );
//        }
//        catch (Exception ex)
//        {
//            App.handleAllExceptions(ex);
//        }
//    }
//
//    //--------------------------------------------------------------------------
//
//    public void showEditDialog(java.awt.Frame parent)
//    {
//        OptionsDlg dlg = new OptionsDlg(parent, this);
//        dlg.show();
//        if (!dlg.wasCancelled())
//        {
//            saveToFile();
//        }
//    }
//
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}
