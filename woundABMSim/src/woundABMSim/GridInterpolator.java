package woundABMSim;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.util.List;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import java.lang.Double;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;

/**
 * GridInterpolator
 * This class includes methods (and the two associated sub classes, BicubicInterpolator and CubicInterpolator) to read 
 * a CSV file to a 2D double variable, perform a bicubic interpolation of a 2D double variable, and write 2D double 
 * variable to a CSV file (output file name: originialName_gridNum.csv).
 * 
 * Adapted from code at...
 * https://stackoverflow.com/questions/9668821/array-interpolation
 * http://commons.apache.org/proper/commons-csv/user-guide.html)
 * https://stackoverflow.com/questions/4397907/updating-specific-cell-csv-file-using-java
 * https://stackoverflow.com/questions/23403552/writing-2d-int-array-to-file-java
 * 
 * Author: Arlynn C Baker
 * Created and validated: 1/11/2019
 * Modified: 
 */

public class GridInterpolator {
    
    public GridInterpolator() {
        // Constructor
    }
    
    public void writeToCSV(String csvOutPath, double[][] data) throws Exception{// csvOutPath example: "path/name.csv";
        File file = new File(csvOutPath);									// Create a new File
        
        if (!file.exists()) {												// If file does not exists create it and
            file.createNewFile();											// fill it with zeros
            CSVReader cr = new CSVReader(new FileReader(file), ',');
            List<String[]> csv = cr.readAll();
            String[] placeHolder = new String[data.length];
            for (int row = 0; row < data.length; row++) { 
                for (int col = 0; col < data.length; col++) {
                    placeHolder[col] = "0";
                }
                csv.add(row, placeHolder);
            }
            cr.close();
            
            CSVWriter writer = new CSVWriter(new FileWriter(file), ',');    // Write CSV information to the new file
            writer.writeAll(csv);
            writer.flush();
            writer.close();
        }
        
        CSVReader reader = new CSVReader(new FileReader(file), ',');        // Read file by each integer row and column
        List<String[]> csvBody = reader.readAll();							// (0 - by 1 - to data length)
        for (int row = 0; row < data.length; row++) {
            for (int col = 0; col < data.length; col++) {
                csvBody.get(row)[col] = String.valueOf(data[row][col]);     // Get CSV cell and replace with data value
                reader.close();												// by row and column
                
                CSVWriter writer = new CSVWriter(new FileWriter(file), ',');// Write to CSV file
                writer.writeAll(csvBody);
                writer.flush();
                writer.close();
            }
        }
    }
    
    public double[][] csvToDouble(String csvPath) throws Exception{
        try {
            Reader in = new FileReader(csvPath);                            // Get the size of the CSV (total # of rows)
            Iterable<CSVRecord> recordSize = CSVFormat.EXCEL.parse(in);
            int gridSize = 0;
            for (CSVRecord record : recordSize) {
                    gridSize = record.size();                               
            }
            
            Reader fr = new FileReader(csvPath);							// Read each value in the CSV and store it
            Iterable<CSVRecord> records = CSVFormat.EXCEL.parse(fr);		// by row and column in a 2D double
            int row = 0;													// variable (source) which is returned
            double[][] source = new double[gridSize][gridSize];
            for (CSVRecord record : records) {
                for (int col = 0; col < gridSize; col++){
                    source[row][col] = Double.parseDouble(record.get(col));
                }
                row = row+1;
            }
            return source;
        }
        
        catch (Exception ex) {                                              // If an exception is thrown, print the
            System.out.println("Exception thrown:"+ex);						// error and return null
        }
        return null;
    }
    
    public void interpolateCSV(String csvPath, int gridNum) throws Exception{
        try {
            double[][] source = csvToDouble(csvPath);                       // Get input CSV as a 2D double
            
            double stepSize = (source.length-1)/1.0/(gridNum-1);            // Perform a bicubic interpolation on source
            BicubicInterpolator bi = new BicubicInterpolator();				// by calculating a row and column stepSize
            double[][] sourceInt = new double[gridNum][gridNum];			// based on the source length and desired
            for (int i = 0; i < gridNum; i++) {								// grid number
                double idx = i*stepSize;
                for(int j = 0; j < gridNum; j++) {
                    double idy = j*stepSize;
                    sourceInt[i][j] = bi.getValue(source, idx, idy);
                }
            }
            
            String csvOutPath = csvPath.replace(csvPath.substring(csvPath.lastIndexOf("_")+1, csvPath.indexOf(".")),
            		Integer.toString(gridNum));								// Output file name reads: name_gridNum.csv
            
            writeToCSV(csvOutPath, sourceInt);                              // Write to CSV
        }
        
        catch (IOException ex) {                                            // If an exception is thrown, print it
            System.out.println("Exception thrown:"+ex);
        }
    }

    public static class CubicInterpolator {
        public static double getValue(double[] p, double x) {
            int xi = (int) x;												// Perform cubic interpolation
            x -= xi;
            double p0 = p[Math.max(0, xi - 1)];
            double p1 = p[xi];
            double p2 = p[Math.min(p.length - 1,xi + 1)];
            double p3 = p[Math.min(p.length - 1, xi + 2)];
            return p1 + 0.5 * x * (p2 - p0 + x * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + x * (3.0 * (p1 - p2) + p3 - p0)));
        }
    }

    public static class BicubicInterpolator extends CubicInterpolator {
        private double[] arr = new double[4];								// Perform bicubic interpolation

        public double getValue(double[][] p, double x, double y) {
            int xi = (int) x;
            x -= xi;
            arr[0] = getValue(p[Math.max(0, xi - 1)], y);
            arr[1] = getValue(p[xi], y);
            arr[2] = getValue(p[Math.min(p.length - 1,xi + 1)], y);
            arr[3] = getValue(p[Math.min(p.length - 1, xi + 2)], y);
            return getValue(arr, x+ 1);
        }
    }
}