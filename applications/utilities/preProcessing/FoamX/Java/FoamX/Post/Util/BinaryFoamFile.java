package FoamX.Post.Util;

//import java.net.URL;
//import java.net.URLConnection;
import java.awt.Component;
import java.io.*;
import java.util.*;

import javax.media.j3d.*;
//import javax.vecmath.Color3f;

import javax.vecmath.Point3f;

import org.j3d.loaders.stl.*;

public class BinaryFoamFile
{

    public BinaryFoamFile()
    {
    }

    public static void readAll(String fName)
    {
        System.out.println("fName=" + fName);
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            while(true)
            {
                int a = LittleEndianConverter.read4ByteBlock(bufferedStream);
                //byte a = dataStream.readByte();

                System.out.println("a=" + a);
            }
        }
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
    }


    // read integer array
    public static int[] readIntArray(String fName)
    {
        System.out.println("fName=" + fName);

        int[] values = null;
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            int nElems =
                LittleEndianConverter.read4ByteBlock(bufferedStream);
            System.out.println("nInt=" + nElems);

            values = new int[nElems];

            for(int i = 0; i < nElems; i++)
            {
                values[i] = 
                    LittleEndianConverter.read4ByteBlock(bufferedStream);
            }
            stream.close();
        }            
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
        return values;
    }

    // read array of integer arrays
    public static int[][] readIntIntArray(String fName)
    {
        System.out.println("fName=" + fName);

        int[][] values = null;
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            int nArrays =
                LittleEndianConverter.read4ByteBlock(bufferedStream);
            System.out.println("nArrays=" + nArrays);

            values = new int[nArrays][];

            for(int arrayI = 0; arrayI < nArrays; arrayI++)
            {
                int nElems =
                    LittleEndianConverter.read4ByteBlock(bufferedStream);
                //System.out.println("nInt=" + nElems);

                values[arrayI] = new int[nElems];

                for(int i = 0; i < nElems; i++)
                {
                    values[arrayI][i] = 
                        LittleEndianConverter.read4ByteBlock(bufferedStream);
                }
            }
            stream.close();
        }            
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
        return values;
    }


    // read float array
    public static float[] readFloatArray(String fName)
    {
        System.out.println("fName=" + fName);

        float[] values = null;
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            int nElems =
                LittleEndianConverter.read4ByteBlock(bufferedStream);
            System.out.println("nFloat=" + nElems);

            values = new float[nElems];
            byte[] readBuffer = new byte[nElems*4];
            int[] dataBuffer = new int[nElems];

            LittleEndianConverter.read
            (
                readBuffer, dataBuffer, 0, nElems, bufferedStream
            );

            for(int i = 0; i < nElems; i++)
            {
                values[i] = Float.intBitsToFloat(dataBuffer[i]);
            }
            stream.close();
        }            
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
        return values;
    }


    // read point3f array
    public static Point3f[] readPoint3fArray(String fName)
    {
        System.out.println("fName=" + fName);

        Point3f[] values = null;
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            int nElems =
                LittleEndianConverter.read4ByteBlock(bufferedStream);
            System.out.println("nPoint3f=" + nElems);

            values = new Point3f[nElems];

            int nFloats = 3*nElems;
            byte[] readBuffer = new byte[nFloats*4];
            int[] dataBuffer = new int[nFloats];

            LittleEndianConverter.read
            (
                readBuffer, dataBuffer, 0, nFloats, bufferedStream
            );

            int floati = 0;
            for(int i = 0; i < nElems; i++)
            {
                float x = Float.intBitsToFloat(dataBuffer[floati++]);
                float y = Float.intBitsToFloat(dataBuffer[floati++]);
                float z = Float.intBitsToFloat(dataBuffer[floati++]);

                values[i] = new Point3f(x, y, z);
            }
        }            
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
        return values;
    }

    // read point array as floats
    public static float[] readFloatPointArray(String fName)
    {
        System.out.println("fName=" + fName);

        float[] values = null;
        try
        {
            FileInputStream stream =
                new FileInputStream(fName);
            BufferedInputStream bufferedStream =
                new BufferedInputStream(stream);
            DataInputStream dataStream =
                new DataInputStream(bufferedStream);

            int nElems =
                LittleEndianConverter.read4ByteBlock(bufferedStream);
            System.out.println("nFloatPoint=" + nElems);


            int nFloats = 3*nElems;
            values = new float[3*nElems];
            byte[] readBuffer = new byte[nFloats*4];
            int[] dataBuffer = new int[nFloats];

            LittleEndianConverter.read
            (
                readBuffer, dataBuffer, 0, nFloats, bufferedStream
            );

            for(int i = 0; i < nFloats; i++)
            {
                values[i] = Float.intBitsToFloat(dataBuffer[i]);
            }
        }            
        catch(FileNotFoundException ex)
        {
            System.out.println("file " + fName + " not found.");
        }
        catch( IOException e )
        {
            System.out.println("IOException caught");
        }
        return values;
    }

    public static void main(String[] args)
    {
        System.out.println("args[0] = " + args[0]);

//        float[] vals = foamFile.readFloatArray();
//        for(int i = 0; i < vals.length; i++)
//        {
//            System.out.println("i=" + i + "  val=" + vals[i]);
//        }
        Point3f[] vals = BinaryFoamFile.readPoint3fArray(args[0]);
        for(int i = 0; i < vals.length; i++)
        {
            System.out.println("i=" + i + "  val=" + vals[i]);
        }
    }
}

