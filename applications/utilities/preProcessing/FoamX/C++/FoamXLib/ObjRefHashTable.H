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

Class
    Foam::ObjRefHashTable

Description
    Hashtable to hold CORBA implementation (ie, RefCountServantBase derived)
    object references. Overrides the clear and erase method to properly release
    the object references. All references are released in Hashtable destructor
    via call to clear method.

\*---------------------------------------------------------------------------*/

#ifndef ObjRefHashTable_H
#define ObjRefHashTable_H

// Foam header files.
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                    Class ObjRefHashTable Declaration
\*---------------------------------------------------------------------------*/

template<class ObjImplClass>
class ObjRefHashTable
:
    public Foam::HashTable<ObjImplClass>
{

public:

    // Destructor

        ~ObjRefHashTable()
        {
            clear();
        }


    // Member Functions
 
        //- Erase an entry specified by given iterator.
        bool erase(typename Foam::HashTable<ObjImplClass>::iterator& iter)
        {
            // Release Corba implementation object.
            if (iter != this->end())
            {
                iter()->_remove_ref();
            }

            // Invoke base class erase method.
            return Foam::HashTable<ObjImplClass>::erase(iter);
        }

        //- Clear all entries from table.
        void clear()
        {
            // Release all Corba implementation objects.
            for
            (
                typename Foam::HashTable<ObjImplClass>::iterator iter = 
                    this->begin();
                iter != this->end();
                ++iter
            )
            {
                iter()->_remove_ref();
            }

            // Invoke base class clear method.
            Foam::HashTable<ObjImplClass>::clear();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
