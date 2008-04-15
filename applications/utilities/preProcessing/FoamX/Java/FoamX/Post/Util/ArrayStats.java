/*
 * ArrayStats.java
 */


package FoamX.Post.Util;

import javax.media.j3d.*;
import javax.vecmath.*;


public class ArrayStats
{
    float min_;
    float max_;
    float avg_;

    public ArrayStats(float[] values)
    {
        if (values.length == 0)
        {
            System.out.println("**ERROR: ArrayStats: length == 0");
            return;
        }
        min_ = values[0];
        max_ = values[0];
        float sum = 0.0f;
        for(int i = 1; i < values.length; i++)
        {
            float val = values[i];
            min_ = Math.min(min_, val);
            max_ = Math.max(max_, val);
            sum += val;
        }
        avg_ = sum / values.length;
    }

    public float getMin()
    {
        return min_;
    }

    public float getMax()
    {
        return max_;
    }

    public float getAvg()
    {
        return avg_;
    }

} // end of class ArrayStats
