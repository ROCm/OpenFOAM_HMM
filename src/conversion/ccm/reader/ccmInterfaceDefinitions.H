/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Containers for holding STARCCM interface definitions

\*---------------------------------------------------------------------------*/
#ifndef ccmInterfaceDefinitions_H
#define ccmInterfaceDefinitions_H

#include "Map.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ccm
{

/*---------------------------------------------------------------------------*\
                     Class ccm::interfaceEntry Declaration
\*---------------------------------------------------------------------------*/

//- A STARCCM interface definition is a pair of boundary ids
class interfaceEntry
{
public:
    // Public Data

        //- The internal interface id
        label id;

        //- The first boundary
        label bnd0;

        //- The second boundary
        label bnd1;

    // Constructors

        //- Construct null
        interfaceEntry()
        :
            id(-1),
            bnd0(-1),
            bnd1(-1)
        {}


        //- Construct empty interface definition
        interfaceEntry(const label index)
        :
            id(index),
            bnd0(-1),
            bnd1(-1)
        {}


        //- Construct from components
        interfaceEntry
        (
            const label index,
            const label boundary0,
            const label boundary1
        )
        :
            id(index),
            bnd0(boundary0),
            bnd1(boundary1)
        {}


    // Access

        //- Check for in-place interfaces
        static bool isInPlace(const std::string& configurationType)
        {
            return configurationType == "IN_PLACE";
        }


        //- True if the boundary id is in this interface
        bool inInterface(label bndId) const
        {
            return bndId == bnd0 || bndId == bnd1;
        }


        //- Canonical name for boundary 0
        word canonicalName0() const
        {
            return "Interface" + ::Foam::name(id) + "_0";
        }

        //- Canonical name for boundary 1
        word canonicalName1() const
        {
            return "Interface" + ::Foam::name(id) + "_1";
        }

        //- Canonical name for boundary
        word canonicalName(label bndId) const
        {
            if (bndId == bnd0)
            {
                return canonicalName0();
            }
            else if (bndId == bnd1)
            {
                return canonicalName1();
            }
            else
            {
                return word::null;
            }
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const interfaceEntry& entry
        )
        {
            os  << "(" << entry.bnd0 << " " << entry.bnd1 << ")";
            return os;
        }

};


/*---------------------------------------------------------------------------*\
                  Class ccm::interfaceDefinitions Declaration
\*---------------------------------------------------------------------------*/

//- A list of available interface definitions
class interfaceDefinitions
:
    public Map<interfaceEntry>
{
public:
    // Constructor

        //- Null construct
        interfaceDefinitions()
        :
            Map<interfaceEntry>()
        {}


        //- Scan available interface entries for one matching this boundary id
        bool add(interfaceEntry entry)
        {
            if
            (
                entry.id >= 0
             && entry.bnd0 >= 0 && entry.bnd1 >= 0
             && entry.bnd0 != entry.bnd1
            )
            {
                this->set(entry.id, entry);
                return true;
            }
            else
            {
                return false;
            }
        }


        //- Scan available interface entries for one matching this boundary id
        bool isInterface(label bndId)
        {
            forAllConstIter(Map<interfaceEntry>, *this, iter)
            {
                if (iter().inInterface(bndId))
                {
                    return true;
                }
            }
            return false;
        }


        //- Scan interface entries for one matching this boundary id
        //  return the canonical name
        word interfaceName(label bndId)
        {
            word ifname;
            forAllConstIter(Map<interfaceEntry>, *this, iter)
            {
                ifname = iter().canonicalName(bndId);
                if (!ifname.empty())
                {
                    break;
                }
            }

            return ifname;
        }


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const interfaceDefinitions& defs
        )
        {
            os  << static_cast<const Map<interfaceEntry>& >(defs)
                << nl;

            return os;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
