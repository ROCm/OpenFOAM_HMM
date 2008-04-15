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

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.Hashtable;

import javax.swing.*;
import javax.swing.text.*;
import javax.swing.event.*;
import javax.swing.undo.*;

import FoamX.App;
import FoamXServer.CaseBrowser.ICaseBrowser;


public class SimpleEditor extends JFrame implements FileListener
{
    JTextPane textPane_;
    DefaultStyledDocument lsd_;
    Hashtable actions_;
    JFileChooser fileChooser_;
    File file_;

    // Listen for changes to followed file
    FileMonitor monitor_;

    // Casebrowser reference
    ICaseBrowser caseBrowser_;

    boolean followMode_;

    // Objects interested in file status
    protected EventListenerList listenerList_;

    //undo helpers
    protected UndoAction undoAction_;
    protected RedoAction redoAction_;
    protected UndoManager undo_ = new UndoManager();

    public SimpleEditor()
    {
        //some initial setup
        super("Simple Editor");

        monitor_ = null;

        caseBrowser_ = null;

        followMode_ = false;

        //Create the document for the text area.
        lsd_ = new DefaultStyledDocument();

        initComponents();

        // Create listener list.
        listenerList_ = new EventListenerList();

        //Start watching for undoable edits and caret changes.
        lsd_.addUndoableEditListener(new MyUndoableEditListener());
    }

    //--------------------------------------------------------------------------
    /** Access textPane */
    public JTextPane getTextPane()
    {
        return textPane_;
    }


    //--------------------------------------------------------------------------
    /** Initialize gui
      */
    private void initComponents()
    {
        //Create the text pane and configure it.
        textPane_ = new JTextPane(lsd_);
        textPane_.setCaretPosition(0);
        textPane_.setMargin(new Insets(5,5,5,5));
        textPane_.setFont(new java.awt.Font("Monospaced", 0, 12));

        JScrollPane scrollPane = new JScrollPane(textPane_);
        scrollPane.setPreferredSize(new Dimension(660, 400));

        setContentPane(scrollPane);

        //Set up the menu bar.
        createActionTable(textPane_);

        fileChooser_ =  new JFileChooser();

        JMenu fileMenu = createfileMenu();
        JMenu editMenu = createEditMenu();
        JMenuBar mb = new JMenuBar();
        mb.add(fileMenu);
        mb.add(editMenu);
        setJMenuBar(mb);

        pack();

        setDefaultCloseOperation
        (
            javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE
        );
        addWindowListener
        (
            new WindowAdapter()
            {
                public void windowClosing(WindowEvent evt)
                {
                    exitEditor();
                }
            }
        );
    }


    //--------------------------------------------------------------------------
    /** Open file */
    public void openFile(File inFile)
    {
        try
        {
            // Inform file status listeners
            fireOpenFile(inFile);

            // Read data from file
            int size = (int)inFile.length();

            FileInputStream inStream = new FileInputStream(inFile);

            int nRead = 0;
            byte[] data = new byte[size];
            while(nRead < size)
            {
                nRead += inStream.read(data, nRead, size-nRead);
            }
            inStream.close();

            textPane_.setText(new String(data));
            file_ = inFile;

            // Clear undoes
            undo_.discardAllEdits();

            // Does not work. Todo
            SimpleEditor.this.setTitle(file_.getName());
        }
        catch (Exception ex)
        {
            System.out.println("Cannot open file " + inFile);
        }
    }

    //--------------------------------------------------------------------------
    /** Follow file */
    public void followFile(ICaseBrowser caseBrowser, File inFile)
    {
        try
        {
            caseBrowser_ = caseBrowser;

            openFile(inFile);

            followMode_ = true;

            // Disable editing
            getTextPane().setEditable(false);

            System.out.println("Positioning cursor");
            JTextPane textPane = getTextPane();
            textPane.setCaretPosition(textPane.getText().length());

            // Install file change listener

            File[] files = new File[1];
            files[0] = inFile;

            int sleep = Integer.parseInt
            (
                App.getOptions().getProperty("FoamX.FileSleep", "1000")
            );

            monitor_ = new FileMonitor(caseBrowser_, files, sleep);

            monitor_.addFileListener(this);

            Thread monitorThread = new Thread(monitor_);

            monitorThread.start();
        }
        catch (Exception ex)
        {
            System.out.println("Cannot open file " + inFile);
        }
    }

    //--------------------------------------------------------------------------
    /** Save file */
    protected void saveFile(File outFile)
    {
        try
        {
            FileOutputStream outStream = new FileOutputStream(outFile);

            outStream.write(textPane_.getText().getBytes());
            outStream.close();

            // Clear undoes
            undo_.discardAllEdits();

            // Inform file status listeners
            fireSaveFile(outFile);
        }
        catch (Exception ex)
        {
            System.out.println("Cannot save file " + outFile);
        }
    }

    //--------------------------------------------------------------------------
    /** Exit panel */
    protected void exitEditor()
    {
        if (undo_.canUndo())
        {
            if
            (
                JOptionPane.showConfirmDialog
                (
                    SimpleEditor.this,
                    "File has changed. Are you sure you want to quit?",
                    "Confirm Quit Operation",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE
                )
                !=
                JOptionPane.OK_OPTION
            )
            {
                return;
            }
        }

        if (monitor_ != null)
        {
            //System.out.println("Uninstalling fileMonitor on " + file_);

            monitor_.removeFileListener(this);
            monitor_.stop();
            monitor_ = null;
        }

        // Inform file status listeners
        fireCloseFile(file_);

        // Delete my frame
        setVisible(false);
        dispose();
    }


    //This one listens for edits that can be undone.
    protected class MyUndoableEditListener
                    implements UndoableEditListener
    {
        public void undoableEditHappened(UndoableEditEvent e)
        {
            //Remember the edit and update the menus.
            undo_.addEdit(e.getEdit());
            undoAction_.updateUndoState();
            redoAction_.updateRedoState();
        }
    }


    //Create the file menu.
    protected JMenu createfileMenu()
    {
        JMenu menu = new JMenu("File");

        ReOpenFileAction reOpenFileAction = new ReOpenFileAction();
        menu.add(reOpenFileAction);

        OpenFileAction openFileAction = new OpenFileAction();
        menu.add(openFileAction);

        SaveFileAction saveFileAction = new SaveFileAction();
        menu.add(saveFileAction);

        menu.addSeparator();

        ExitAction exitAction = new ExitAction();
        menu.add(exitAction);

        return menu;
    }

    //--------------------------------------------------------------------------
    //Create the edit menu.
    protected JMenu createEditMenu()
    {
        JMenu menu = new JMenu("Edit");

        //Undo and redo are actions of our own creation.
        undoAction_ = new UndoAction();
        menu.add(undoAction_);

        redoAction_ = new RedoAction();
        menu.add(redoAction_);

        menu.addSeparator();

        //These actions come from the default editor kit.
        //Get the ones we want and stick them in the menu.
        menu.add(getActionByName(DefaultEditorKit.cutAction));
        menu.add(getActionByName(DefaultEditorKit.copyAction));
        menu.add(getActionByName(DefaultEditorKit.pasteAction));

        menu.addSeparator();

        menu.add(getActionByName(DefaultEditorKit.selectAllAction));
        return menu;
    }

    //--------------------------------------------------------------------------

    //The following two methods allow us to find an
    //action provided by the editor kit by its name.
    private void createActionTable(JTextComponent textComponent)
    {
        actions_ = new Hashtable();
        Action[] actionsArray = textComponent.getActions();
        for (int i = 0; i < actionsArray.length; i++) {
            Action a = actionsArray[i];
            actions_.put(a.getValue(Action.NAME), a);
        }
    }


    //--------------------------------------------------------------------------

    private Action getActionByName(String name)
    {
        return (Action)(actions_.get(name));
    }

    //--------------------------------------------------------------------------
    
    class OpenFileAction extends AbstractAction
    {
        public OpenFileAction()
        {
            super("Open");
            setEnabled(true);
        }
          
        public void actionPerformed(ActionEvent e)
        {
            if (file_ != null)
            {
                File parentFile = file_.getParentFile();
                if (parentFile == null)
                {
                    parentFile = new File(".");
                }

                fileChooser_.setCurrentDirectory(parentFile);
                fileChooser_.setSelectedFile(file_);
            }
            if
            (
                fileChooser_.showOpenDialog(SimpleEditor.this)
             != JFileChooser.APPROVE_OPTION
            )
            {
                return;
            }
            openFile(fileChooser_.getSelectedFile());
        }
    }

    //--------------------------------------------------------------------------
    
    class ReOpenFileAction extends AbstractAction
    {
        public ReOpenFileAction()
        {
            super("ReOpen");
            setEnabled(true);
        }
          
        public void actionPerformed(ActionEvent e)
        {
            openFile(file_);
        }
    }

    //--------------------------------------------------------------------------
          
    class SaveFileAction extends AbstractAction
    {
        public SaveFileAction()
        {
            super("Save");
            setEnabled(true);
        }
          
        public void actionPerformed(ActionEvent e)
        {
            File parentFile = file_.getParentFile();
            if (parentFile == null)
            {
                parentFile = new File(".");
            }

            fileChooser_.setCurrentDirectory(parentFile);
            fileChooser_.setSelectedFile(file_);
 
            if
            (
                fileChooser_.showSaveDialog(SimpleEditor.this)
             != JFileChooser.APPROVE_OPTION
            )
            {
                return;
            }
            saveFile(fileChooser_.getSelectedFile());
        }
    }

    //--------------------------------------------------------------------------
          
    class ExitAction extends AbstractAction
    {
        public ExitAction()
        {
            super("Exit");
            setEnabled(true);
        }
          
        public void actionPerformed(ActionEvent e)
        {
            exitEditor();
        }
    }

    //--------------------------------------------------------------------------
          
    class UndoAction extends AbstractAction
    {
        public UndoAction()
        {
            super("Undo");
            setEnabled(false);
        }
          
        public void actionPerformed(ActionEvent e)
        {
            try {
                undo_.undo();
            }
            catch (CannotUndoException ex)
            {
                System.out.println("Unable to undo: " + ex);
                ex.printStackTrace();
            }
            updateUndoState();
            redoAction_.updateRedoState();
        }
          
        protected void updateUndoState()
        {
            if (undo_.canUndo())
            {
                setEnabled(true);
                putValue(Action.NAME, undo_.getUndoPresentationName());
            }
            else
            {
                setEnabled(false);
                putValue(Action.NAME, "Undo");
            }
        }      
    }    

    //--------------------------------------------------------------------------

    class RedoAction extends AbstractAction
    {
        public RedoAction()
        {
            super("Redo");
            setEnabled(false);
        }

        public void actionPerformed(ActionEvent e)
        {
            try
            {
                undo_.redo();
            }
            catch (CannotRedoException ex)
            {
                System.out.println("Unable to redo: " + ex);
                ex.printStackTrace();
            }
            updateRedoState();
            undoAction_.updateUndoState();
        }

        protected void updateRedoState()
        {
            if (undo_.canRedo())
            {
                setEnabled(true);
                putValue(Action.NAME, undo_.getRedoPresentationName());
            }
             else
            {
                setEnabled(false);
                putValue(Action.NAME, "Redo");
            }
        }
    }    

    //--------------------------------------------------------------------------
    //---- FileListener Interface
    //--------------------------------------------------------------------------

    public void addFileListener(FileListener l)
    {
        listenerList_.add(FileListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeFileListener(FileListener l)
    {
        listenerList_.remove(FileListener.class, l);
    }

    //--------------------------------------------------------------------------


    protected void fireOpenFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileOpened(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireCloseFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileClosed(evt);
            }
        }
    }

    //--------------------------------------------------------------------------

    protected void fireSaveFile(File file)
    {
        // Create event object.
        FileEvent evt = new FileEvent(this, file);

        // Process the listeners last to first, notifying those that are interested in this event.
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == FileListener.class)
            {
                ((FileListener)listeners[i+1]).fileChanged(evt);
            }
        }
    }

//    //--------------------------------------------------------------------------
//    /** Starts SimpleEditor. All the reading to/from casebrowser is handled
//      * by FileListener events
//      */
//    public class KeepOpeningThread implements Runnable
//    {
//        SimpleEditor editor_;
//        File file_;
//
//        public KeepOpeningThread(SimpleEditor editor, File file)
//        {
//            super();
//            editor_ = editor;
//            file_ = file;
//        }
//
//        public void run()
//        {
//            try
//            {
//                while (true)
//                {
//                    System.out.println("Trying to load " + file_);
//                    editor_.openFile(file_);
//
//                    System.out.println("Positioning cursor");
//                    JTextPane textPane = editor_.getTextPane();
//                    String text = textPane.getText();
//                    textPane.setCaretPosition(text.length());
//
//                    Thread.currentThread().sleep(1000);
//                }
//            }
//            catch(InterruptedException iex) {}
//        }
//    }


    //--------------------------------------------------------------------------
    //---- FileListener Interface.
    //     Callback from FileMonitor (only if in follow mode)
    //--------------------------------------------------------------------------

    public void fileOpened(FileEvent evt)
    {}

    public void fileClosed(FileEvent evt)
    {}

    public void fileChanged(FileEvent evt)
    {
        File file = evt.file();

        //System.out.println("Followed file changed " + file);

        openFile(file);

        JTextPane textPane = getTextPane();
        textPane.setCaretPosition(textPane.getText().length());
    }

    //--------------------------------------------------------------------------

    //The standard main method.
    public static void main(String[] args)
    {
        if (args.length != 1)
        {
            System.out.println("Usage: java SimpleEditor file");
            System.exit(1);
        }
        final SimpleEditor frame =
            new SimpleEditor();

        File file = new File(args[0]);

        frame.openFile(file);

        frame.show();
    }
}
