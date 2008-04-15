package FoamX.Util;

import javax.swing.JOptionPane;

import FoamX.App;

public class Assert
{
    private static void perr(String msg)
    {
        JOptionPane.showMessageDialog
        (
            null, 
            "Assertion Failed",
            msg,
            JOptionPane.ERROR_MESSAGE
        );

        App.printMessage(msg);
    }

    public final static void True(boolean exp)
    {
        if (!exp) perr("Assertion failed");
    }

    public final static void False(boolean exp)
    {
        if (exp) perr("Assertion failed");
    }

    public final static void True(boolean exp, String msg)
    {
        if (!exp) perr("Assertion failed: " + msg);
    }

    public final static void False(boolean exp, String msg)
    {
        if (exp) perr("Assertion failed: " + msg);
    }
}
