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

import javax.swing.*;
import java.io.*;
import java.util.*;

import FoamX.App;

import FoamX.Exceptions.FoamXException;

public class ArgsManager
    extends java.util.Properties
{
    //--------------------------------------------------------------------------

    private String foamXSystemPath_;
    private String foamXSystemConfigPath_;
    private String foamXUserConfigPath_;

    //--------------------------------------------------------------------------
    /** ArgsManager constructor. */
    public ArgsManager(String[] args)
    {
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
}
