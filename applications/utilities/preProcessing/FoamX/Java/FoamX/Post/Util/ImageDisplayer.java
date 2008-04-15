/*
 *	@(#)ImageDisplayer.java 1.3 01/01/11 07:38:16
 *
 * Copyright (c) 1996-2001 Sun Microsystems, Inc. All Rights Reserved.
 *
 * Sun grants you ("Licensee") a non-exclusive, royalty free, license to use,
 * modify and redistribute this software in source and binary code form,
 * provided that i) this copyright notice and license appear on all copies of
 * the software; and ii) Licensee does not utilize the software in a manner
 * which is disparaging to Sun.
 *
 * This software is provided "AS IS," without a warranty of any kind. ALL
 * EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY
 * IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR
 * NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL NOT BE
 * LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A RESULT OF USING, MODIFYING
 * OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL SUN OR ITS
 * LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR DIRECT,
 * INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES, HOWEVER
 * CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE USE OF
 * OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES.
 *
 * This software is not designed or intended for use in on-line control of
 * aircraft, air traffic, aircraft navigation or aircraft communications; or in
 * the design, construction, operation or maintenance of any nuclear
 * facility. Licensee represents and warrants that it will not use or
 * redistribute the Software for such purposes.
 */
package FoamX.Post.Util;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.awt.event.*;

public class ImageDisplayer extends JFrame implements ActionListener {
    BufferedImage bImage;

    private class ImagePanel extends JPanel {
	public void paint(Graphics g) {
	    g.setColor(Color.black);
	    g.fillRect(0, 0, getSize().width, getSize().height);
	    g.drawImage(bImage, 0, 0, this);
	}

	private ImagePanel() {
	    setPreferredSize(new Dimension(bImage.getWidth(),
					   bImage.getHeight()));
	}
    }

    private JMenuItem printItem;
    private JMenuItem closeItem;

    public void actionPerformed (ActionEvent event) {
	Object target = event.getSource();

	if (target == printItem) {
	    new ImagePrinter(bImage).print();
	}
	else if (target == closeItem) {
	    this.removeAll();
	    this.setVisible(false);
	    bImage = null;
	}
    }

    private JMenuBar createMenuBar() {
	JMenuBar menuBar = new JMenuBar();
	JMenu fileMenu = new JMenu("File");
	printItem = new JMenuItem("Print...");
	printItem.addActionListener(this);
	closeItem = new JMenuItem("Close");
	closeItem.addActionListener(this);
	fileMenu.add(printItem);
	fileMenu.add(new JSeparator());
	fileMenu.add(closeItem);
	menuBar.add(fileMenu);
	return menuBar;
    }

    public ImageDisplayer(BufferedImage bImage) {
	this.bImage = bImage;
	this.setTitle("Off-screen Canvas3D Snapshot");

	// Create and initialize menu bar
	this.setJMenuBar(createMenuBar());

	// Create scroll pane, and embedded image panel
	ImagePanel imagePanel = new ImagePanel();
	JScrollPane scrollPane = new JScrollPane(imagePanel);

        int width = Math.min(700, bImage.getWidth());
        int height = Math.min(700, bImage.getHeight());

//	scrollPane.getViewport().setPreferredSize(new Dimension(700, 700));
	scrollPane.getViewport().setPreferredSize(new Dimension(width, height));

	// Add scroll pane to the frame and make it visible
	this.getContentPane().add(scrollPane);
	this.pack();
	this.setVisible(true);
    }
}
