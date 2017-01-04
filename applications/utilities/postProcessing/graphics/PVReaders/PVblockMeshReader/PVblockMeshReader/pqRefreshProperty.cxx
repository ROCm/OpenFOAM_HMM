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

\*---------------------------------------------------------------------------*/

#include "pqRefreshProperty.h"

#include <QPushButton>
#include <QGridLayout>

#include "pqApplicationCore.h"
#include "pqView.h"
#include "vtkSMDocumentation.h"
#include "vtkSMIntVectorProperty.h"
#include "vtkSMSourceProxy.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope
static void setButtonProperties
(
    QAbstractButton* b,
    vtkSMIntVectorProperty* prop,
    bool initChecked = true
)
{
    QString tip;

    vtkSMDocumentation* doc = prop->GetDocumentation();
    if (doc)
    {
        const char* txt = doc->GetDescription();
        if (txt)
        {
            tip = QString(txt).simplified();
        }
    }

    b->setText(prop->GetXMLLabel());
    if (tip.size())
    {
        b->setToolTip(tip);
    }
    b->setFocusPolicy(Qt::NoFocus); // avoid dotted border

    // initial checked state
    if (initChecked)
    {
        b->setChecked(prop->GetElement(0));
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pqRefreshProperty::refreshPressed()
{
    // Update everything
    refresh_->Modified();

    vtkSMSourceProxy::SafeDownCast(this->proxy())->UpdatePipeline();

    // Render all views
    pqApplicationCore::instance()->render();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pqRefreshProperty::pqRefreshProperty
(
    vtkSMProxy* proxy,
    vtkSMProperty* prop,
    QWidget* parent
)
:
    Superclass(proxy, parent),
    refresh_(vtkSMIntVectorProperty::SafeDownCast(prop))
{
    // Replace with our UI content
    this->setShowLabel(false);

    QGridLayout* form = new QGridLayout(this);

    QPushButton* b = new QPushButton(this);
    setButtonProperties(b, refresh_, false);
    form->addWidget(b, 0, 0, Qt::AlignLeft);

    connect(b, SIGNAL(clicked()), this, SLOT(refreshPressed()));
    refresh_->SetImmediateUpdate(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pqRefreshProperty::~pqRefreshProperty()
{}


// ************************************************************************* //
