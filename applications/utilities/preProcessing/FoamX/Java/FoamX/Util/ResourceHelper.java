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

import java.util.*;
import java.net.URL;
import javax.swing.ImageIcon;

import FoamX.App;

public class ResourceHelper
    extends Object
{
    //--------------------------------------------------------------------------

    private ResourceBundle resources_;

    //--------------------------------------------------------------------------
    /** ResourceHelper constructor. */
    public ResourceHelper(String bundleName) throws MissingResourceException
    {
        try
        {
            // Get the application resource bundle object.
            resources_ = ResourceBundle.getBundle(bundleName, Locale.getDefault());
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }
    }

    //--------------------------------------------------------------------------
    /** Extract a simple string from the resource bundle. */
    public String getResourceString(String key) throws MissingResourceException
    {
        String str = null;

        try
        {
            str = resources_.getString(key);
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }

        return str;
    }

    //--------------------------------------------------------------------------

    public String[] getResourceStringArray(String key) throws MissingResourceException
    {
        String[] strArray = null;

        try
        {
            String str = resources_.getString(key);
            strArray = tokenizeString(str);
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }

        return strArray;
    }

    //--------------------------------------------------------------------------

    public URL getResourceURL(String key) throws MissingResourceException
    {
        URL url = null;

        try
        {
            // Get the Url string from the resource file.
            String name = resources_.getString(key);
            url = getClass().getResource(name);

            // Check that this Url is valid.
            if (url == null)
            {
                throw new MissingResourceException
                (
                    "Invalid Url",
                    "URL",
                    key
                );
            }
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }

        return url;
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the Icon associated with the name from the resources.
     * The resouce should be in the path.
     * @param  key Name of the resource specifying the icon.
     * @return The ImageIcon object or null if the icon is not found.
     */
    public ImageIcon getIcon(String key) throws MissingResourceException
    {
        ImageIcon icon = null;

        try
        {
            // Get the Url string from the resource file.
            String name = resources_.getString(key);
            URL url = getClass().getResource(name);

            // Check that this Url is valid.
            if (url == null)
            {
                throw new MissingResourceException
                (
                    "Invalid Url",
                    "URL",
                    key
                );
            }

            // Load the icon at the Url.
            icon = new ImageIcon(url);

            if (icon == null)
            {
                throw new MissingResourceException
                (
                    "Invalid Icon Url",
                    "ImageIcon",
                    key
                );
            }
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }

        return icon;
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the Icon from the specified URL.
     */
    public ImageIcon getIcon(URL url) throws MissingResourceException
    {
        // Check that this Url is valid.
        if (url == null)
        {
            throw new MissingResourceException
            (
                "Invalid Url",
                "URL",
                "null"
            );
        }

        ImageIcon icon = null;

        try
        {
            // Load the icon at the Url.
            icon = new ImageIcon(url);

            if (icon == null)
            {
                throw new MissingResourceException
                (
                    "Invalid Url",
                    "ImageIcon",
                    url.toString()
                );
            }
        }
        catch (MissingResourceException mre)
        {
            App.handleException(mre);
        }

        return icon;
    }

    //--------------------------------------------------------------------------
    /**
     * Take the given string and chop it up into a series of
     * strings on whitespace boundries. This is useful for
     * getting an array of strings out of the resource file.
     */
    protected String[] tokenizeString(String input)
    {
        Vector          v = new Vector();
        StringTokenizer t = new StringTokenizer(input);
        String cmd[];

        while (t.hasMoreTokens())
        {
            v.addElement(t.nextToken());
        }
        cmd = new String[v.size()];
        for (int i=0; i <cmd.length; i++)
        {
            cmd[i] = (String)v.elementAt(i);
        }

        return cmd;
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
}


