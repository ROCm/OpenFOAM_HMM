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
import java.io.*;
import java.awt.Point;

import org.omg.CORBA.StringHolder;
import FoamXServer.CaseBrowser.ICaseBrowser;
import FoamXServer.FoamXIOError;

import FoamX.App;
import FoamX.Exceptions.FoamXException;

public class FileEditor
    implements FileListener
{
    String fileName_ = "";
    File workFile_;
    ICaseBrowser caseBrowser_;
    SimpleEditor simpleEditor_;

    public FileEditor
    (
        boolean modal,
        boolean followMode,
        ICaseBrowser caseBrowser,
        String fileName,
        int startLineNo,
        int endLineNo,
        int columnNo
    )
    {
        try
        {
            fileName_ = fileName;
            caseBrowser_ = caseBrowser;
            workFile_= new File(fileName_);

            // Get file. Use case browser if nessecary.

            boolean nfsFile = workFile_.canRead();

            if (!nfsFile && (fileName_.length() > 0))
            {
                if (caseBrowser_ == null)
                {
                    throw new FoamXException
                    (
                        "Cannot open file " + fileName_ + " for editing."
                    );
                }

                // Generate temporary file for the editors to edit.
                workFile_ = java.io.File.createTempFile
                (
                    workFile_.getName(),
                    ""
                );

                App.printMessage
                (
                    "Working on temporary file:" + workFile_.getName()
                );

                // Quick and dirty test for availability of file
                if
                (
                    !copyFromCaseBrowser
                    (
                        caseBrowser_,
                        fileName_,
                        new File("/dev/null")
                    )
                )
                {
                    throw new FoamXException
                    (
                        "Cannot open file " + fileName_ + " for editing."
                    );
                }
            }


            String editorMethod = App.getOptions().getProperty
            (
                "FoamX.Editor",
                "internal"
            );
            if (editorMethod.equals("internal"))
            {
                internalEditor(modal, followMode, nfsFile);
            }
            else
            {
                externalEditor
                (
                    modal,
                    followMode,
                    nfsFile, 
                    editorMethod, 
                    startLineNo, 
                    endLineNo, 
                    columnNo
                );
            }
        }
        catch (java.io.IOException ex)
        {
            App.handleAllExceptions(ex);
        }
        catch (FoamXException ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------

    private void internalEditor
    (
        boolean modal,
        boolean followMode,
        boolean nfsFile
    )
    {
        // Start editor.
        App.printMessage("Starting internal editor on " + fileName_);

        simpleEditor_ = new SimpleEditor();

        // Center on screen
        Point topLeft = App.getRootFrame().getLocationOnScreen();
//        System.out.println("Size:" + App.getRootFrame().getSize());
//        System.out.println("Position:" + topLeft.x + "  " + topLeft.y);
//        System.out.println("Here");

//        java.awt.Dimension screenSize =
//            java.awt.Toolkit.getDefaultToolkit().getScreenSize();
        simpleEditor_.setLocation(topLeft);

        if (!nfsFile)
        {
            // Make sure we are warned of any changes (FileListener events)
            simpleEditor_.addFileListener(this);
        }

        if (followMode)
        {
            simpleEditor_.followFile(caseBrowser_, workFile_);
        }
        else
        {
            simpleEditor_.openFile(workFile_);
        }

        simpleEditor_.show();
    }


    //--------------------------------------------------------------------------
    /** Starts External editor.
      * Since no events blocks until editing finished
      */
    private void externalEditor
    (
        boolean modal,
        boolean followMode,
        boolean nfsFile,
        String editorMethod,
        int startLineNo,
        int endLineNo,
        int columnNo
    )   throws FoamXException
    {
        if (!nfsFile)
        {
            copyFromCaseBrowser(caseBrowser_, fileName_, workFile_);
        }

        // Start editor on workFile_

        Process editProcess = null;

        // External: replace $FILE etc. in editorMethod with actual values
        try
        {
            // Execute the desired browser passing it the required URL.
            int hasFile = editorMethod.indexOf("$FILE");
            if (hasFile == -1)
            {
                throw new Exception
                (
                   "FoamX.Editor property malformed " + editorMethod +
                   ".\nIt should at least contain $FILE " +
                   "(expands to file name)"
                );
            }

            Hashtable vars = new Hashtable(3);
            vars.put("FILE", workFile_);

            if (startLineNo <= 0)
            {
                startLineNo = 1;
            }
            vars.put("START", new Integer(startLineNo));
            if (endLineNo <= 0)
            {
                endLineNo = startLineNo;
            }
            vars.put("END", new Integer(endLineNo));
            vars.put("COL", new Integer(columnNo));
            VarString editCmd = new VarString(editorMethod, vars);

            // Start editor.
            App.printMessage("Starting external editor : " + editCmd.string());
            editProcess = Runtime.getRuntime().exec(editCmd.string());
        }
        catch (java.io.IOException ex)
        {
            App.handleAllExceptions(ex);
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }



        // If nessecary copy working file back to actual one.

        try
        {
            if (!nfsFile && (editProcess != null))
            {
                // Wait for editing to finish before copying back
                editProcess.waitFor();

                copyToCaseBrowser(workFile_, caseBrowser_, fileName_);

                workFile_.delete();
                workFile_ = null;
            }
            else if (modal)
            {
                // Wait for editing to finish
                editProcess.waitFor();
            }
        }
        catch (java.lang.InterruptedException ex)
        {
            App.handleAllExceptions(ex);
        }
    }


    //--------------------------------------------------------------------------
    /** Copies file inFileName, located at caseBrowser, into local outFile */
    static public boolean copyFromCaseBrowser
    (
        ICaseBrowser caseBrowser,
        String inFileName,
        File outFile
    )
    {
        try
        {
            StringHolder holder = new StringHolder();
            caseBrowser.readFile(inFileName, holder);
            String contents = holder.value;

            DataOutputStream out =
            new DataOutputStream
            (
                new BufferedOutputStream
                (
                    new FileOutputStream(outFile)
                )
            );
            out.writeBytes(contents);
            out.close();
            return true;
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(new FoamXException(ioErr.errorMessage));
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
        return false;
    }


    //--------------------------------------------------------------------------
    /** Copies local file inFile to caseBrowser as outFileName */
    static public void copyToCaseBrowser
    (
        File inFile,
        ICaseBrowser caseBrowser,
        String outFileName
    )
    {
        try
        {
            // Limited to editing files < 2Gb only
            int size = (int)inFile.length();

            FileInputStream inStream = new FileInputStream(inFile);

            int nRead = 0;
            byte[] data = new byte[size];
            while(nRead < size)
            {
                nRead += inStream.read(data, nRead, size-nRead);
            }
            inStream.close();

            caseBrowser.writeFile(outFileName, new String(data));
        }
        catch (FoamXIOError ioErr)
        {
            App.handleAllExceptions(new FoamXException(ioErr.errorMessage));
        }
        catch (Exception ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    //--------------------------------------------------------------------------
    //---- FileListener Interface
    //--------------------------------------------------------------------------

    // SimpleEditor wants to read file
    public void fileOpened(FileEvent evt)
    {
        File file = evt.file();

        try
        {
            if(file.getCanonicalPath().equals(workFile_.getCanonicalPath()))
            {
                copyFromCaseBrowser(caseBrowser_, fileName_, workFile_);
            }
        }
        catch (java.io.IOException ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    // SimpleEditor has closed a file.
    public void fileClosed(FileEvent evt)
    {
        File file = evt.file();

        try
        {
            if(file.getCanonicalPath().equals(workFile_.getCanonicalPath()))
            {
                // event can only happen once per FileEditor so cleanup
                simpleEditor_.removeFileListener(this);
                workFile_.delete();
                caseBrowser_ = null;
            }
        }
        catch (java.io.IOException ex)
        {
            App.handleAllExceptions(ex);
        }
    }

    // SimpleEditor has saved a file.
    public void fileChanged(FileEvent evt)
    {
        File file = evt.file();

        try
        {
            if(file.getCanonicalPath().equals(workFile_.getCanonicalPath()))
            {
                copyToCaseBrowser(workFile_, caseBrowser_, fileName_);
            }
        }
        catch (java.io.IOException ex)
        {
            App.handleAllExceptions(ex);
        }
    }
}
