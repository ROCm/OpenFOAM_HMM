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
package FoamX.CaseManagement;

import javax.swing.Action;
import javax.swing.Icon;
import javax.swing.JFrame;

import FoamXServer.CaseServer.ICaseServer;
import FoamXServer.CaseBrowser.ICaseBrowser;

import FoamX.WindowManagement.FoamXInternalFrame;

public interface IModuleHost
{
    // Miscellaneous methods.
    public JFrame getRootFrame();
    public ICaseServer getCaseServer();
    public ICaseBrowser getCaseBrowser();

    // Reporting methods.
    public void printMessage(String message);

    // Tree management methods.
    public void setRootLabel(String label);
    public String getRootLabel();
    public int addNode(int parentNodeID, String nodeLabel, Icon icon);
    public void removeNode(int nodeID);

    // Add action as popup to treeNode
    public void addNodeAction(int nodeID, Action action, String context);
    // Retrieve stored action
    public Action getAction(String actionName);

    // Window management methods.
    public int createWindow(javax.swing.JComponent component);
    public FoamXInternalFrame getWindowFrame(int windowID);
    public void closeWindow(int windowID);
    public void showWindow(int windowID);
    public void hideWindow(int windowID);
}


