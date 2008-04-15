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
package FoamX.Editors.DictionaryEntryEditor;

import FoamX.App;

public class CompoundEntryEvent
    extends java.util.EventObject
{
    public static final int TYPE_PATCH      = 0;
    public static final int TYPE_DICTIONARY = 1;

    protected int type_;
    protected String name_;

    /** CompoundEntryEvent constructor. */
    public CompoundEntryEvent(Object source, int type, String name)
    {
        
        super(source);
        type_ = type;
        name_ = name;
    }

    public int getType()
    {
        return type_;
    }

    public String getName()
    {
        return name_;
    }
}

