package FoamX.Post;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.util.*;

import FoamX.Post.InteractiveNodes.InteractiveNode;
import FoamX.Post.InteractiveNodes.InteractiveNodeList;
import FoamX.Post.InteractiveNodes.FeatureLinesNode;
import FoamX.Post.Shapes.FeatureLines;


public class MouseModePanel extends JPanel
{
    static JFrame frame;

    static String zeroString = "rotate&pick";
    static String oneString = "create loop from patch0";
    static String twoString = "load arguments";
    static String threeString = "print shapes";

    protected EventListenerList listenerList_;
    InteractiveNodeList shapes_;


    public MouseModePanel(InteractiveNodeList shapes)
    {
        shapes_ = shapes;

        // Create listener list.
        listenerList_ = new EventListenerList();


        // Create the radio buttons.
        JRadioButton zeroButton = new JRadioButton(zeroString);
        zeroButton.setMnemonic(KeyEvent.VK_B);
        zeroButton.setActionCommand(zeroString);
        zeroButton.setSelected(true);

        JRadioButton oneButton = new JRadioButton(oneString);
        oneButton.setMnemonic(KeyEvent.VK_R);
        oneButton.setActionCommand(oneString);

        JRadioButton twoButton = new JRadioButton(twoString);
        twoButton.setMnemonic(KeyEvent.VK_R);
        twoButton.setActionCommand(twoString);

        JRadioButton threeButton = new JRadioButton(threeString);
        threeButton.setMnemonic(KeyEvent.VK_R);
        threeButton.setActionCommand(threeString);

        // Group the radio buttons.
        ButtonGroup group = new ButtonGroup();
        group.add(zeroButton);
        group.add(oneButton);
        group.add(twoButton);
        group.add(threeButton);

        // Register a listener for the radio buttons.
        RadioListener zeroListener = new RadioListener(0);
        zeroButton.addActionListener(zeroListener);
        RadioListener oneListener = new RadioListener(1);
        oneButton.addActionListener(oneListener);
        RadioListener twoListener = new RadioListener(2);
        twoButton.addActionListener(twoListener);
        RadioListener threeListener = new RadioListener(3);
        threeButton.addActionListener(threeListener);


        // Put the radio buttons in a column in a panel
        JPanel radioPanel = new JPanel();
        radioPanel.setLayout(new GridLayout(0, 1));
        radioPanel.add(zeroButton);
        radioPanel.add(oneButton);
        radioPanel.add(twoButton);
        radioPanel.add(threeButton);

        setLayout(new BorderLayout());
        add(radioPanel, BorderLayout.WEST);
        setBorder(BorderFactory.createEmptyBorder(20,20,20,20));
    }

    /** Listens to the radio buttons. */
    class RadioListener implements ActionListener
    {
        int choice_;

        public RadioListener(int choice)
        {
            super();
            choice_ = choice;
        }

        public void actionPerformed(ActionEvent e) {
//            System.out.println(e.getActionCommand());

//            InteractiveNode shape = shapes_.getByName("featureLines_0");
//            System.out.println("Found shape=" + shape);
//
//            if (e.getActionCommand() == oneString)
//            {
//                FeatureLinesNode node = (FeatureLinesNode)shape;
//                FeatureLines fLines = node.getFeatureLines();
//
//                int[] edges = fLines.getSelectedEdges();
//                for(int i = 0; i < edges.length; i++)
//                {
//                    System.out.println("i=" + i + "  edges=" + edges[i]);
//                }
//
//                //shape.hide();
//            }
//            else if (e.getActionCommand() == zeroString)
//            {
//                //shape.show();
//            }

//            fireChoice(new ChoiceEvent(this, choice_));
            fireChoice(e);
        }
    }


    //--------------------------------------------------------------------------
    //---- ChoiceListener Methods
    //--------------------------------------------------------------------------
    public void addChoiceListener(ChoiceListener l)
    {
        listenerList_.add(ChoiceListener.class, l);
    }

    //--------------------------------------------------------------------------

    public void removeChoiceListener(ChoiceListener l)
    {
        listenerList_.remove(ChoiceListener.class, l);
    }

    //--------------------------------------------------------------------------

    protected void fireChoice(ActionEvent e)
    {
        Object[] listeners = listenerList_.getListenerList();
        for (int i=listeners.length-2; i>= 0; i-=2)
        {
            if (listeners[i] == ChoiceListener.class)
            {
                ((ChoiceListener)listeners[i+1]).actionPerformed(e);
            }
        }
    }


    //--------------------------------------------------------------------------
    //---- Main
    //--------------------------------------------------------------------------
    public static void main(String s[])
    {
        frame = new JFrame("Mouse Action");
        frame.addWindowListener
        (
            new WindowAdapter()
            {
                public void windowClosing(WindowEvent e) {System.exit(0);}
            }
        );
 
        frame.getContentPane().add
        (
            new MouseModePanel(new InteractiveNodeList()),
            BorderLayout.CENTER
        );
        frame.pack();
        frame.setVisible(true);
    }
}
