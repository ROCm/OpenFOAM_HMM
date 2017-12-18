/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::functionObjects::fieldSelection

Description
    Helper class to manage field selections

SourceFiles
    fieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldSelection_H
#define functionObjects_fieldSelection_H

#include "HashSet.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fieldSelection Declaration
\*---------------------------------------------------------------------------*/

class fieldSelection
:
    public HashSet<wordRe>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fieldSelection(const fieldSelection&) = delete;


protected:

    // Protected member data

        //- Reference to the database
        const objectRegistry& obr_;

        //- Current field selection
        wordHashSet selection_;


    // Protected Member Functions

        //- Add registered objects of a given type
        template<class Type>
        void addRegistered(wordHashSet& set) const;


public:

    //- Construct from object registry
    fieldSelection(const objectRegistry& obr);


    //- Destructor
    virtual ~fieldSelection();


    // Member Functions

        //- Return the current field selection
        const wordHashSet& selection() const
        {
            return selection_;
        }

        //- Read the fieldSelection data from dictionary
        virtual bool read(const dictionary& dict);

        //- Return whether the field names contain a pattern
        virtual bool containsPattern() const;

        //- Clear the current selection
        virtual void clearSelection();

        //- Update the selection using current contents of obr_
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldSelectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
