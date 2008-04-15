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

public class VarString
{
    private String expandedStr_ = "";

    //--------------------------------------------------------------------------
    /** VarString constructor. */
    public VarString(String str, Hashtable vars)
        throws Exception
    {
        int i = 0;

        while (i < str.length())
        {
            char c = str.charAt(i);

            if (c == '$')
            {
                int varEnd = parseVar(str, i+1);
                String varName = str.substring(i+1, varEnd);
                if (!vars.containsKey(varName))
                {
                    throw new Exception
                    (
                       "Unknown variable name " + varName + " in " + str
                    );
                }
                expandedStr_ += vars.get(varName);
                i = varEnd;
            }
            else
            {
                expandedStr_ += c;
                i++;
            }
        }
    }


    //--------------------------------------------------------------------------
    /** Get expanded string */
    public String string()
    {
         return expandedStr_;
    }

    //--------------------------------------------------------------------------
    /** Find end of variable name. Returns index of first non identifier 
      * character
      */
    private int parseVar(String str, int i)
    {
       // First character after $ is always excepted.
       i++;
       while ((i < str.length()) && Character.isLetterOrDigit(str.charAt(i)))
       {
            i++;
       }
       return i;
    }
}
