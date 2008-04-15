package FoamX.Post.Table;


import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.*;

import FoamX.Post.InteractiveNodes.InteractiveNodeList;
import FoamX.Post.InteractiveNodes.InteractiveNode;

import javax.swing.event.TableModelListener; 
import javax.swing.event.TableModelEvent; 

public class NodeTable extends JFrame implements TableModelListener
{
    private boolean DEBUG = true;

    InteractiveNodeList shapes_;
    JTable table_;

    //--------------------------------------------------------------------------
    

    public NodeTable(InteractiveNodeList shapes)
    {
        super("NodeTable");

        shapes_ = shapes;

        table_ = new JTable(shapes_);

        //AbstractTableModel model = shapes_.getModel();
        shapes_.addTableModelListener(this);

        table_.setPreferredScrollableViewportSize(new Dimension(200, 150));

        //Create the scroll pane and add the table to it. 
        JScrollPane scrollPane = new JScrollPane(table_);

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new GridLayout(2, 1));
        mainPanel.add(scrollPane);

        JPanel buttonPanel = new JPanel();
        JButton removeButton = new JButton("Remove");
        JButton editButton = new JButton("Edit");

        buttonPanel.setLayout(new GridLayout(1, 2));
        buttonPanel.add(removeButton);
        buttonPanel.add(editButton);
        mainPanel.add(buttonPanel);

        //Add the main panel to this window.
        getContentPane().add(mainPanel, BorderLayout.CENTER);

        removeButton.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    removeButton_ActionPerformed(evt);
                }
            }
        );
        editButton.addActionListener
        (
            new ActionListener()
            {
                public void actionPerformed(ActionEvent evt)
                {
                    editButton_ActionPerformed(evt);
                }
            }
        );

        addWindowListener
        (
            new WindowAdapter()
            {
                public void windowClosing(WindowEvent e)
                {
                    //System.exit(0);
                }
            }
        );
    }

    //--------------------------------------------------------------------------

    /** remove button handler */
    private void removeButton_ActionPerformed(ActionEvent evt)
    {
        InteractiveNodeList model = (InteractiveNodeList)table_.getModel();
        int startRow = table_.getSelectionModel().getMinSelectionIndex();
        int endRow = table_.getSelectionModel().getMaxSelectionIndex();
        if (startRow < 0)
        {
            return;
        }
        for (int row = startRow; row <= endRow; row++)
        {
            model.removeElementAt(row);
        }

        table_.getSelectionModel().clearSelection();
        initTableSize(table_, model);
        table_.repaint();
    }

    //--------------------------------------------------------------------------

    /** edit button handler */
    private void editButton_ActionPerformed(ActionEvent evt)
    {
        InteractiveNodeList model = (InteractiveNodeList)table_.getModel();
        int startRow = table_.getSelectionModel().getMinSelectionIndex();
        int endRow = table_.getSelectionModel().getMaxSelectionIndex();
        if (startRow < 0)
        {
            return;
        }
        for (int row = startRow; row <= endRow; row++)
        {
            System.out.println("Edit of row " + row);
        }
    }

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------

    /** Size table itself */
    static private void initTableSize(JTable table, InteractiveNodeList model)
    {
        int height = table.getRowHeight() * model.getRowCount();

        table.setPreferredSize
        (
            new java.awt.Dimension(table.getWidth(), height)
        );
        table.revalidate();
    } 

    //--------------------------------------------------------------------------

    public void tableChanged(TableModelEvent e) {
        //int row = e.getFirstRow();
        //int column = e.getColumn();

        initTableSize(table_, shapes_);
        table_.repaint();
    }


    //--------------------------------------------------------------------------
    
//    public static void main(String[] args) {
//        NodeTable frame = new TableDemo();
//        frame.pack();
//        frame.setVisible(true);
//    }
}
