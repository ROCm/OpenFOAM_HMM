/*
 *      Axis.java 1.0 98/11/25
 *
 * Copyright (c) 1998 Sun Microsystems, Inc. All Rights Reserved.
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

/*
 * Getting Started with the Java 3D API
 * written in Java 3D
 *
 * This program demonstrates:
 *   1. writing a visual object class
 *      In this program, Axis class defines a visual object
 *      This particular class extends Shape3D
 *      See the text for a discussion.
 *   2. Using LineArray to draw 3D lines.
 */

//package geometry;


import javax.media.j3d.*;
import javax.vecmath.*;

    public class Axis extends Shape3D{

        ////////////////////////////////////////////
        //
        // create axis visual object
        //
        public Axis() {
        
            this.setGeometry(createGeometry());
            this.setCapability(Shape3D.ALLOW_BOUNDS_READ);
        }

        private Geometry createGeometry(){
            // create line for X axis
            IndexedLineArray axisLines = new IndexedLineArray
            (
                18,
                GeometryArray.COORDINATES,
                30
            );

            axisLines.setCoordinate( 0, new Point3f(-1.0f, 0.0f, 0.0f));
            axisLines.setCoordinate( 1, new Point3f( 1.0f, 0.0f, 0.0f));
            axisLines.setCoordinate( 2, new Point3f( 0.9f, 0.1f, 0.1f));
            axisLines.setCoordinate( 3, new Point3f( 0.9f,-0.1f, 0.1f));
            axisLines.setCoordinate( 4, new Point3f( 0.9f, 0.1f,-0.1f));
            axisLines.setCoordinate( 5, new Point3f( 0.9f,-0.1f,-0.1f));
            axisLines.setCoordinate( 6, new Point3f( 0.0f,-1.0f, 0.0f));
            axisLines.setCoordinate( 7, new Point3f( 0.0f, 1.0f, 0.0f));
            axisLines.setCoordinate( 8, new Point3f( 0.1f, 0.9f, 0.1f));
            axisLines.setCoordinate( 9, new Point3f(-0.1f, 0.9f, 0.1f));
            axisLines.setCoordinate(10, new Point3f( 0.1f, 0.9f,-0.1f));
            axisLines.setCoordinate(11, new Point3f(-0.1f, 0.9f,-0.1f));
            axisLines.setCoordinate(12, new Point3f( 0.0f, 0.0f,-1.0f));
            axisLines.setCoordinate(13, new Point3f( 0.0f, 0.0f, 1.0f));
            axisLines.setCoordinate(14, new Point3f( 0.1f, 0.1f, 0.9f));
            axisLines.setCoordinate(15, new Point3f(-0.1f, 0.1f, 0.9f));
            axisLines.setCoordinate(16, new Point3f( 0.1f,-0.1f, 0.9f));
            axisLines.setCoordinate(17, new Point3f(-0.1f,-0.1f, 0.9f));

            axisLines.setCoordinateIndex( 0, 0);
            axisLines.setCoordinateIndex( 1, 1);
            axisLines.setCoordinateIndex( 2, 2);
            axisLines.setCoordinateIndex( 3, 1);
            axisLines.setCoordinateIndex( 4, 3);
            axisLines.setCoordinateIndex( 5, 1);
            axisLines.setCoordinateIndex( 6, 4);
            axisLines.setCoordinateIndex( 7, 1);
            axisLines.setCoordinateIndex( 8, 5);
            axisLines.setCoordinateIndex( 9, 1);
            axisLines.setCoordinateIndex(10, 6);
            axisLines.setCoordinateIndex(11, 7);
            axisLines.setCoordinateIndex(12, 8);
            axisLines.setCoordinateIndex(13, 7);
            axisLines.setCoordinateIndex(14, 9);
            axisLines.setCoordinateIndex(15, 7);
            axisLines.setCoordinateIndex(16,10);
            axisLines.setCoordinateIndex(17, 7);
            axisLines.setCoordinateIndex(18,11);
            axisLines.setCoordinateIndex(19, 7);
            axisLines.setCoordinateIndex(20,12);
            axisLines.setCoordinateIndex(21,13);
            axisLines.setCoordinateIndex(22,14);
            axisLines.setCoordinateIndex(23,13);
            axisLines.setCoordinateIndex(24,15);
            axisLines.setCoordinateIndex(25,13);
            axisLines.setCoordinateIndex(26,16);
            axisLines.setCoordinateIndex(27,13);
            axisLines.setCoordinateIndex(28,17);
            axisLines.setCoordinateIndex(29,13);


            axisLines.setCapability(IndexedLineArray.ALLOW_COORDINATE_READ);
            axisLines.setCapability(IndexedLineArray.ALLOW_FORMAT_READ);
            axisLines.setCapability(IndexedLineArray.ALLOW_COUNT_READ);
            axisLines.setCapability(IndexedLineArray.ALLOW_COORDINATE_INDEX_READ);

            return axisLines;

        } // end of Axis createGeometry()


    } // end of class Axis
