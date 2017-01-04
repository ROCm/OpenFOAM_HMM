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

#include "pqShowPointNumbersProperty.h"

#include <QCheckBox>
#include <QGridLayout>

#include "pqApplicationCore.h"
#include "pqView.h"
#include "vtkSMDocumentation.h"
#include "vtkSMIntVectorProperty.h"

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

void pqShowPointNumbersProperty::showPointNumbers(bool checked)
{
    showPointNumbers_->SetElement(0, checked);

    // Update the active view
    if (this->view())
    {
        this->view()->render();
    }

    // OR: update all views
    // pqApplicationCore::instance()->render();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pqShowPointNumbersProperty::pqShowPointNumbersProperty
(
    vtkSMProxy* proxy,
    vtkSMProperty* prop,
    QWidget* parent
)
:
    Superclass(proxy, parent),
    showPointNumbers_(vtkSMIntVectorProperty::SafeDownCast(prop))
{
    // Replace with our UI content
    this->setShowLabel(false);

    QGridLayout* form = new QGridLayout(this);

    QCheckBox* b = new QCheckBox(this);
    setButtonProperties(b, showPointNumbers_);
    form->addWidget(b, 0, 0, Qt::AlignLeft);

    connect(b, SIGNAL(toggled(bool)), this, SLOT(showPointNumbers(bool)));
    showPointNumbers_->SetImmediateUpdate(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pqShowPointNumbersProperty::~pqShowPointNumbersProperty()
{}


// ************************************************************************* //
