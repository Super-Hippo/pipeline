package edu.stevens.dhutchis.accumuloiter;



import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author dhutchis
 */


public final class myGlobal{
    public static List<String> data ; // x-axis=#sequences, y-axis=time-to-get-answer

    private myGlobal()
    {
        data = new ArrayList<>();
    }
}
