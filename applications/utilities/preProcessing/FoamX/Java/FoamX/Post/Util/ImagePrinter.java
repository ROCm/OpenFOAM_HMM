/*
 *	@(#)ImagePrinter.java 1.3 01/01/11 07:38:16
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

import java.awt.print.*;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.awt.geom.AffineTransform;

public class ImagePrinter implements Printable, ImageObserver {
    BufferedImage bImage;

    public int print(Graphics g, PageFormat pf, int pi)
	throws PrinterException {


	if (pi >= 1) {
	    return Printable.NO_SUCH_PAGE;
	}

	Graphics2D g2d = (Graphics2D)g;
	//g2d.translate(pf.getImageableX(), pf.getImageableY());
	AffineTransform t2d = new AffineTransform();
	t2d.translate(pf.getImageableX(), pf.getImageableY());
	double xscale  = pf.getImageableWidth() / (double)bImage.getWidth();
	double yscale  = pf.getImageableHeight() / (double)bImage.getHeight();
	double scale = Math.min(xscale, yscale);
	t2d.scale(scale, scale);
	try {
	    g2d.drawImage(bImage,t2d, this);
	}
	catch (Exception ex) {
	    ex.printStackTrace();
	    return Printable.NO_SUCH_PAGE;
	}
	return Printable.PAGE_EXISTS;
    }

    void print() {
	PrinterJob printJob = PrinterJob.getPrinterJob();
	PageFormat pageFormat = printJob.defaultPage();
	pageFormat.setOrientation(PageFormat.LANDSCAPE);
	pageFormat = printJob.validatePage(pageFormat);
	printJob.setPrintable(this, pageFormat);
	if (printJob.printDialog()) {
	    try {
		printJob.print();
	    }
	    catch (PrinterException ex) {
		ex.printStackTrace();
	    }
	}
    }

    public boolean imageUpdate(Image img,
			       int infoflags,
			       int x,
			       int y,
			       int width,
			       int height) {
	return false;
    }

    ImagePrinter(BufferedImage bImage) {
	this.bImage = bImage;
    }
}
