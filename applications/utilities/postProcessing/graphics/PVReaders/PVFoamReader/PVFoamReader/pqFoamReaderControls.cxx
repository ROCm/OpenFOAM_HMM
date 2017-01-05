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

#include "pqFoamReaderControls.h"

#include <QCheckBox>
#include <QFrame>
#include <QGridLayout>
#include <QPushButton>

#include "pqApplicationCore.h"
#include "pqPipelineRepresentation.h"
#include "pqView.h"
#include "vtkSMDocumentation.h"
#include "vtkSMIntVectorProperty.h"
#include "vtkSMPropertyGroup.h"
#include "vtkSMSourceProxy.h"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope
static QAbstractButton* setButtonProperties
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

    return b;
}


static vtkSMIntVectorProperty* lookupIntProp
(
    vtkSMPropertyGroup* group,
    const char* name
)
{
    vtkSMProperty* prop = group->GetProperty(name);

    if (prop)
    {
        return vtkSMIntVectorProperty::SafeDownCast(prop);
    }

    return nullptr;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pqFoamReaderControls::updatePartsStatus()
{
    vtkSMProperty* prop = this->proxy()->GetProperty("PartArrayStatus");
    if (prop)
    {
        this->proxy()->UpdatePropertyInformation(prop);
    }
}


void pqFoamReaderControls::updatePartsStatus(bool)
{
    updatePartsStatus();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void pqFoamReaderControls::refreshPressed()
{
    // Update everything
    refresh_->Modified();

    vtkSMSourceProxy::SafeDownCast(this->proxy())->UpdatePipeline();

    // Update all views
    pqApplicationCore::instance()->render();
}


void pqFoamReaderControls::cacheMesh(bool checked)
{
    cacheMesh_->SetElement(0, checked);
}


void pqFoamReaderControls::showPatchNames(bool checked)
{
    showPatchNames_->SetElement(0, checked);

    // update the active view
    if (this->view())
    {
        this->view()->render();
    }
    // OR: update all views
    // pqApplicationCore::instance()->render();
}


void pqFoamReaderControls::showGroupsOnly(bool checked)
{
    showGroupsOnly_->SetElement(0, checked);
    updatePartsStatus();
}


void pqFoamReaderControls::includeSets(bool checked)
{
    includeSets_->SetElement(0, checked);
    updatePartsStatus();
}


void pqFoamReaderControls::includeZones(bool checked)
{
    includeZones_->SetElement(0, checked);
    updatePartsStatus();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pqFoamReaderControls::pqFoamReaderControls
(
    vtkSMProxy* proxy,
    vtkSMPropertyGroup* group,
    QWidget* parent
)
:
    Superclass(proxy, parent),
    refresh_(lookupIntProp(group, "Refresh")),
    showPatchNames_(lookupIntProp(group, "ShowPatchNames")),
    showGroupsOnly_(lookupIntProp(group, "ShowGroupsOnly")),
    includeSets_(lookupIntProp(group, "IncludeSets")),
    includeZones_(lookupIntProp(group, "IncludeZones")),
    cacheMesh_(lookupIntProp(group, "CacheMesh"))
{
    typedef vtkSMIntVectorProperty intProp;

    QGridLayout* form = new QGridLayout(this);

    // ROW
    // ~~~
    int row = 0;

    if (refresh_)
    {
        QPushButton* b = new QPushButton(this);
        setButtonProperties(b, refresh_, false);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        connect(b, SIGNAL(clicked()), this, SLOT(refreshPressed()));
        refresh_->SetImmediateUpdate(true);
    }

    intProp* zeroTime = lookupIntProp(group, "ZeroTime");
    if (zeroTime)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, zeroTime);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), zeroTime);
    }

    // LINE
    // ~~~~
    ++row;
    {
        QFrame* hline = new QFrame(this);
        hline->setFrameStyle(QFrame::HLine | QFrame::Sunken);
        form->addWidget(hline, row, 0, 1, 4);
    }

    // ROW
    // ~~~
    ++row;

    if (includeSets_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, includeSets_);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), includeSets_);

        connect(b, SIGNAL(toggled(bool)), this, SLOT(includeSets(bool)));
        includeSets_->SetImmediateUpdate(true);
    }

    if (showGroupsOnly_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showGroupsOnly_);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), showGroupsOnly_);

        connect(b, SIGNAL(toggled(bool)), this, SLOT(showGroupsOnly(bool)));
        showGroupsOnly_->SetImmediateUpdate(true);
    }


    // ROW
    // ~~~
    ++row;

    if (includeZones_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, includeZones_);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), includeZones_);

        connect(b, SIGNAL(toggled(bool)), this, SLOT(includeZones(bool)));
        includeZones_->SetImmediateUpdate(true);
    }

    if (showPatchNames_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, showPatchNames_);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        connect(b, SIGNAL(toggled(bool)), this, SLOT(showPatchNames(bool)));
        showPatchNames_->SetImmediateUpdate(true);
    }

    // LINE
    // ~~~~
    ++row;
    {
        QFrame* hline = new QFrame(this);
        hline->setFrameStyle(QFrame::HLine | QFrame::Sunken);
        form->addWidget(hline, row, 0, 1, 4);
    }

    // ROW
    // ~~~
    ++row;

    intProp* interpolate = lookupIntProp(group, "InterpolateFields");
    if (interpolate)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, interpolate);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), interpolate);
    }

    intProp* extrapolate = lookupIntProp(group, "ExtrapolatePatches");
    if (extrapolate)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, extrapolate);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), extrapolate);
    }

    // LINE
    // ~~~~
    ++row;
    {
        QFrame* hline = new QFrame(this);
        hline->setFrameStyle(QFrame::HLine | QFrame::Sunken);
        form->addWidget(hline, row, 0, 1, 4);
    }

    // ROW
    // ~~~
    ++row;

    intProp* updateGui = lookupIntProp(group, "UpdateGUI");
    if (updateGui)
    {
        QPushButton* b = new QPushButton(this);
        setButtonProperties(b, updateGui, false);
        form->addWidget(b, row, 0, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(clicked()), updateGui);
    }

    intProp* usePolyhedron = lookupIntProp(group, "UseVTKPolyhedron");
    if (usePolyhedron)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, usePolyhedron);
        form->addWidget(b, row, 1, Qt::AlignLeft);

        addPropertyLink(b, "checked", SIGNAL(toggled(bool)), usePolyhedron);
    }

    if (cacheMesh_)
    {
        QCheckBox* b = new QCheckBox(this);
        setButtonProperties(b, cacheMesh_);
        form->addWidget(b, row, 2, Qt::AlignLeft);

        connect(b, SIGNAL(toggled(bool)), this, SLOT(cacheMesh(bool)));
        cacheMesh_->SetImmediateUpdate(true);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pqFoamReaderControls::~pqFoamReaderControls()
{}


// ************************************************************************* //
