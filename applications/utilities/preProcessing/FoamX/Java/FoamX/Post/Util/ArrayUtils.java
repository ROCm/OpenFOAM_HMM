/*
 * ArrayUtils.java
 */


package FoamX.Post.Util;

import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;


public class ArrayUtils
{
    public static float[] mag(float[] a)
    {
        int nVecs = a.length / 3;

        float[] result = new float[nVecs];

        int coordI = 0;
        for(int i = 0; i < nVecs; i++)
        {
            float x = a[coordI++];
            float y = a[coordI++];
            float z = a[coordI++];

            result[i] = (float)Math.sqrt(x*x + y*y + z*z);
        }
        return result;
    }

    public static float[] add(float[] a, float b)
    {
        float[] result = new float[a.length];

        for(int i = 0; i < a.length; i++)
        {
            result[i] = a[i] + b;
        }
        return result;
    }

    public static float[] add(float[] a, float[] b)
    {
        float[] result = new float[a.length];

        for(int i = 0; i < a.length; i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    public static float[] sub(float[] a, float b)
    {
        float[] result = new float[a.length];

        for(int i = 0; i < a.length; i++)
        {
            result[i] = a[i] - b;
        }
        return result;
    }

    public static float[] sub(float[] a, float[] b)
    {
        float[] result = new float[a.length];

        for(int i = 0; i < a.length; i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    public static float[] div(float[] a, float b)
    {
        float[] result = new float[a.length];

        for(int i = 0; i < a.length; i++)
        {
            result[i] = a[i] / b;
        }
        return result;
    }

    public static float[] div(float[] a, float[] b)
    {
        if (a.length != 3*b.length)
        {
            System.out.println("Error: div: length not ok:" + a.length + " " + b.length);
            return null;
        }

        int nVecs = a.length / 3;

        float[] result = new float[nVecs];

        int coordI = 0;
        for(int i = 0; i < nVecs; i++)
        {
            result[coordI] = a[coordI] / b[i];
            coordI++;
            result[coordI] = a[coordI] / b[i];
            coordI++;
            result[coordI] = a[coordI] / b[i];
            coordI++;
        }

        return result;
    }

    // debug print
    public static void print(int[] array)
    {
        for(int i = 0; i < array.length; i++)
        {
            System.out.println(i + ":" + array[i]);
        }
    }

    // debug print
    public static void print(int[][] array)
    {
        for(int i = 0; i < array.length; i++)
        {
            System.out.print(i + ":");
            for(int j = 0; j < array[i].length; j++)
            {
                System.out.print
                (
                    " " + array[i][j]
                );
            }
            System.out.println();
        }
    }

    // debug print
    public static void print(float[] array)
    {
        for(int i = 0; i < array.length; i++)
        {
            System.out.println(i + ":" + array[i]);
        }
    }

    // debug print
    public static void print(Object[] array)
    {
        for(int i = 0; i < array.length; i++)
        {
            System.out.println(i + ":" + array[i]);
        }
    }
} // end of class ArrayUtils
