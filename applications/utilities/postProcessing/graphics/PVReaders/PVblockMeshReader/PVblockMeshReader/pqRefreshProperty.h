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
    pqRefreshProperty

Description
    Custom refresh button (ParaView blockMesh reader)

SourceFiles
    pqRefreshProperty.cxx

\*---------------------------------------------------------------------------*/
#ifndef pqRefreshProperty_h
#define pqRefreshProperty_h

#include "pqPropertyWidget.h"

// Forward declarations (ParaView)
class vtkSMIntVectorProperty;


/*---------------------------------------------------------------------------*\
                      Class pqRefreshProperty Declaration
\*---------------------------------------------------------------------------*/

class pqRefreshProperty
:
    public pqPropertyWidget
{
    Q_OBJECT;
    typedef pqPropertyWidget Superclass;

    // Private data

        //- Refresh (bool property - as push button)
        vtkSMIntVectorProperty* refresh_;


protected slots:

    // Protected Member Functions

    //- Trigger refresh
    void refreshPressed();


public:

    //- Construct from components
    pqRefreshProperty
    (
        vtkSMProxy* proxy,
        vtkSMProperty* prop,
        QWidget* parent = nullptr
    );


    //- Destructor
    virtual ~pqRefreshProperty();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
