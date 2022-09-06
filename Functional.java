package toShare;
import java.util.*;
import java.io.*;

public class Functional{

	public static void main(String[] args) throws IOException {
		/**double[][] spline = cubicSpline(Functional.importData("glassAlpha1V12Fitted.pot"));
		System.out.println(evaluateSpline(spline, 1.15));*/
		/**double rho = 0.68;
		double[][] g2 = standardize(importData("rsa.g2"), 0.05, 0, 20);
		double[][] h = new double[g2.length][2];
		for (int i = 0; i < g2.length; i++) {
			h[i][0] = g2[i][0];
			h[i][1] = g2[i][1] - 1;
		}
		for (double k = 0; k < 20; k = k + 0.05) {
			double sk = 1 + fourier2D(h, k)*rho;
			System.out.println(k + " " + sk);
		}*/
                double[][] f = importData("Quintic/Targ.g2");
		for (int i = 0; i < f.length; i++) {
			f[i][1] -= 1;
		}
                //for (double k = 200; k < 250; k += 0.05) {
                for (double k = 50; k < 60; k += 0.05) {        
			//double ftilde = inverseFourier2D(f, k);
			double ftilde = 1 + 0.8*fourier2D(f, k);
                        System.out.println(k + " " + ftilde);
                }
	}

	/**imports single column data
 *      @param fileName String, the file name of the data source
 *      @return double[][] a n*2 array of the xy-data */
        public static double[] importData1D(String fileName) throws IOException {
                double[] tempResult = new double[20000];
                Scanner s = null;
                int line = 0; /** this will be the number of lines in the file */
                try {
                        String path = "../InverseSk/" + fileName;
                        s = new Scanner(new BufferedReader(new FileReader(path)));
                        while (s.hasNextDouble() && line < 20000) {
                                double x = s.nextDouble();
                                tempResult[line] = x;
                                line ++;
                        }
                }
                finally {
                        s.close();
                }
                double[] result = new double[line];
                for (int i = 0; i < line; i++) {
                        result[i] = tempResult[i];
                }
                return result;
        }

	/**imports xy-valued data
 *   	@param fileName String, the file name of the data source
 *    	@return double[][] a n*2 array of the xy-data */
	public static double[][] importData(String fileName) throws IOException {
		double[][] tempResult = new double[20000][2];
		Scanner s = null;
		int line = 0; /** this will be the number of lines in the file */
		try {
                        String path = "../InverseSk/" + fileName;
                        s = new Scanner(new BufferedReader(new FileReader(path)));
                        while (s.hasNextDouble() && line < 20000) {
                                double x = s.nextDouble();
                                double y = s.nextDouble();
                                tempResult[line][0] = x;
				tempResult[line][1] = y;
				line ++;  
                        }
                }
                finally {
                        s.close();
                }
		double[][] result = new double[line][2];
		for (int i = 0; i < line; i++) {
			result[i][0] = tempResult[i][0];
			result[i][1] = tempResult[i][1];	
		}
		return result;		
	}
	
	/**imports xy-valued data from a certain line to a certain line
 *      @param fileName String, the file name of the data source
 *      @param startLine int, the starting line to import. The first line is line 0
 *      @param endLine int
 *      @return double[][] a n*2 array of the xy-data */
        public static double[][] importData(String fileName, int startLine, int endLine) throws IOException {
                double[][] tempResult = new double[30000][2];
                Scanner s = null;
                int line = 0; /** this will be the number of lines in the file */
                try {
                        String path = "../InverseSk/" + fileName;
                        s = new Scanner(new BufferedReader(new FileReader(path)));
                        while (s.hasNextLine() && line < startLine) {
                                s.nextLine();
                                line ++;
                        }
                        while (s.hasNextDouble() && line < endLine) {
                                int dataRow = line - startLine;//row in the exported data
                                double x = s.nextDouble();
                                double y = s.nextDouble();
				s.nextLine();
                                tempResult[dataRow][0] = x;
                                tempResult[dataRow][1] = y;
                                line ++;
                        }
                }
                finally {
                        s.close();
                }
                double[][] result = new double[endLine - startLine][2];
                for (int i = 0; i < result.length; i++) {
                        result[i][0] = tempResult[i][0];
                        result[i][1] = tempResult[i][1];
                }
                return result;
        }

        /**imports xy-valued data from a certain line to a certain line
 *  	@param fileName String, the file name of the data source
 *      @param startLine int, the starting line to import. The first line is line 0
 *      @param endLine int
 *      @return double[][] a n*3 array of the xyz-data */
        public static double[][] importData3D(String fileName, int startLine, int endLine) throws IOException {
        	double[][] tempResult = new double[30000][3];
        	Scanner s = null;
        	int line = 0; /** this will be the number of lines in the file */
        	try {
        		String path = "../InverseSk/" + fileName;
        		s = new Scanner(new BufferedReader(new FileReader(path)));
        		while (s.hasNextLine() && line < startLine) {
        			s.nextLine();
        			line ++;
        		}
        		while (s.hasNextDouble() && line < endLine) {
        			int dataRow = line - startLine;//row in the exported data
        			double x = s.nextDouble();
        			double y = s.nextDouble();
        			double z = s.nextDouble();
        			tempResult[dataRow][0] = x;
        			tempResult[dataRow][1] = y;
        			tempResult[dataRow][2] = z;
        			line ++;
        		}
        	}
        	finally {
        		s.close();
        	}
        	double[][] result = new double[endLine - startLine][3];
        	for (int i = 0; i < result.length; i++) {
        		result[i][0] = tempResult[i][0];
        		result[i][1] = tempResult[i][1];
        		result[i][2] = tempResult[i][2];
        	}
        	return result;
        }



	/**imports xyz-valued data
 *      @param fileName String, the file name of the data source
 *      @return double[][] a n*3 array of the xy-data */
        public static double[][] importData3D(String fileName) throws IOException {
                double[][] tempResult = new double[5000][3];
                Scanner s = null;
                int line = 0; /** this will be the number of lines in the file */
                try {
                        String path = "../InverseSk/" + fileName;
                        s = new Scanner(new BufferedReader(new FileReader(path)));
                        while (s.hasNextDouble() && line < 5000) {
                                double x = s.nextDouble();
                                double y = s.nextDouble();
				double z = s.nextDouble();
				s.nextLine();
                                tempResult[line][0] = x;
                                tempResult[line][1] = y;
				tempResult[line][2] = z;
                                line ++;
                        }
                }
                finally {
                        s.close();
                }
                double[][] result = new double[line][3];
                for (int i = 0; i < line; i++) {
                        result[i][0] = tempResult[i][0];
                        result[i][1] = tempResult[i][1];
			result[i][2] = tempResult[i][2];
                }
                return result;
        }


	/** computes the L2 distance squared between two functions
 * 	@param f1 double[][], a function
 * 	@param f2 double[][], another function */
	public static double l2sq(double[][] f1, double[][] f2) {
		double binSize = 0.05;
		double start = 0.005;
		double end = 10;
		f1 = standardize(f1, binSize, start, end);
		f2 = standardize(f2, binSize, start, end);
		double result = 0;
		for (int i = 0; i < f2.length; i++) {
			result = result + Math.pow(f1[i][1]-f2[i][1],2);
		}
		return result;
	}

	/** computes the L2 distance squared between two functions, with custom start and end
 *      @param f1 double[][], a function
 *      @param f2 double[][], another function
 *      @param start double, the left point of the range
 *      @param end double, the right point of the range */
        public static double l2sq(double[][] f1, double[][] f2, double start, double end) {
                double binSize = 0.05;
                f1 = standardize(f1, binSize, start, end);
                f2 = standardize(f2, binSize, start, end);
                double result = 0;
                for (int i = 0; i < f2.length; i++) {
                        result = result + Math.pow(f1[i][1] - f2[i][1], 2)*2*Math.PI*f1[i][0]*binSize;
                }
                return result;
        }

	
	/**standardizes the bin positions of a function: valued at multiples of a given bin size
 * 	@param f double[][], a function
 * 	@param binSize double, the bin size needed
 * 	@param start double, the left point of the range
 * 	@param end double, the right point of the range
 * 	@return double[][], the standardized function */
	public static double[][] standardize(double[][] f, double binSize, double start, double end) {
		int resultSize = (int)((end - start)/(double)binSize) + 1;
		double[][] result = new double[resultSize][2];
		double fStart = f[0][0];
		double fEnd = f[f.length-1][0];
		for (int i = 0; i < resultSize; i++) {
			double x =  start + i*binSize;
			result[i][0] = x;
			if (x <= fStart) {
				/**double x1 = f[0][0];
                                double y1 = f[0][1];
                                double x2 = f[1][0];
                                double y2 = f[1][1];
                                double derivative = (y2 - y1)/(x2 - x1);
				result[i][1] = f[i][1] - derivative*(x1 - x);*/
				//result[i][1] = x/(fStart)*f[0][1];
				//result[i][1] = x*x/(fStart*fStart)*f[0][1];
				result[i][1] = f[0][1];
			}
			else if (x < fEnd) {
				int j = 0;
				while (f[j][0] <= x) {
					j++;
				}
				double x2 = f[j][0];
				double y2 = f[j][1];
				double x1 = f[j-1][0];
				double y1 = f[j-1][1];
				double y = y1 + (x - x1)*(y2 - y1)/(x2 - x1);
				result[i][1] = y;
			}
			else {
				result[i][1] = f[f.length-1][1];
			}
		}
		return result;
	}

	/**evaluate a function at a given point, assume even grid!!!
 * 	@param f double[][], a function
 * 	@param x double, the input
 * 	@param binSize double
 * 	@return double, the output function value*/
	public static double evaluate(double[][] f, double x) {
		double fStart = f[0][0];
		double fEnd = f[f.length-1][0];
		double binSize = f[1][0] - f[0][0];
		if (x <= fStart) {
			return f[0][1];
		}
		else if (x < fEnd) {
			int j = (int)((x-f[0][0])/binSize);
			double x2 = f[j+1][0];
			double y2 = f[j+1][1];
			double x1 = f[j][0];
			double y1 = f[j][1];
			double y = y1 + (x - x1)*(y2 - y1)/(x2 - x1);
			return y;
		}
		else {
			return f[f.length-1][1];
		}
	}

	/**create a cubic spline a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
	 * REF: https://en.wikipedia.org/wiki/Spline_(mathematics)
	 * @param f double[][], a function
	 * @return double[][], the spline parameters. The columns of output are (x, a, b, c, d)
	 */
	public static double[][] cubicSpline(double[][] f) {
		int n = f.length - 1; //number of bins
		double[][] result = new double[n + 1][5];
		for (int i = 0; i < n + 1; i++) {
			result[i][0] = f[i][0];
			result[i][1] = f[i][1];
		}
		double h = f[1][0] - f[0][0];//uniform bins
		double[] alpha = new double[n + 1];
		for (int i = 1; i < n; i++) {
			alpha[i] = 3.0/h*(result[i + 1][1] + result[i - 1][1] - 2*result[i][1]);
		}
		double[] l = new double[n + 1];
		l[0] = 1;
		l[n] = 0;
		double[] mu = new double[n + 1];
		double[] z = new double[n + 1];
		for (int i = 1; i < n; i++) {
			l[i] = 4*h - h*(mu[i - 1]);
			mu[i] = h/l[i];
			z[i] = (alpha[i] - h*z[i - 1])/l[i];
		}
		for (int i = n - 1; i >= 0; i--) {
			result[i][3] = z[i] - mu[i]*result[i + 1][3];
			result[i][2] = (result[i + 1][1] - result[i][1])/h - h*(result[i + 1][3] + 2*result[i][3])/3.0;
			result[i][4] = (result[i + 1][3] - result[i][3])/(3.0*h);
		}
		return result;
	}
	
	/**
	 * evaluates a cubic spline at point x
	 * @param f, double[][5], the spline
	 * @param x, the evaluation point
	 * @return double
	 */
	public static double evaluateSpline(double[][] f, double x) {
		double fStart = f[0][0];
		double fEnd = f[f.length-1][0];
		double binSize = f[1][0] - f[0][0];
		int i;//bin index where x falls!
		if (x <= fStart) {
			i = 0;
		}
		else if (x < fEnd) {
			i = (int)((x-f[0][0])/binSize);
		}
		else {
			i = f.length - 1;
		}
		double xi = f[i][0];
		double ai = f[i][1];
		double bi = f[i][2];
		double ci = f[i][3];
		double di = f[i][4];
		double dx = x - xi;
		double dxsq = dx*dx;
		return ai + bi*dx + ci*dx*dx + di*dx*dxsq;
	}
	
	/** Fourier transform for 3D radial functions Phys Rep eqn. 10
 *      @param fr double[][], the orginal function
 *      @param k double, the Fourier space input
 *      @return double, the inverse fourier transform ftilda(k) */
        public static double fourier3D(double[][] fr, double k) throws IOException {
                double sum = 0;
                for (int i = 0; i < fr.length - 1; i++) {
                        double r = fr[i][0];
                        double rpdr = fr[i+1][0]; /**k+dk*/
			double besselHalf = Math.sqrt(2.0/(Math.PI*k*r))*Math.sin(k*r);
                        sum = sum + (rpdr - r) * r * r * fr[i][1] * besselHalf/Math.sqrt(k*r);
                }
                return sum*Math.pow(2*Math.PI, 1.5);
        }	

	/**inverse Fourier transform for 3D radial functions Phys Rep eqn. 10
 *	@param fk double[][], the fourier transformed function
 *	@param r double, the direct space input
 *	@return double, the inverse fourier transform f(r) */
	public static double inverseFourier3D(double[][] fk, double r) throws IOException {
		double sum = 0;
		for (int i = 0; i < fk.length - 1; i++) {
			double k = fk[i][0];
			double kpdk = fk[i+1][0]; /**k+dk*/
			if (k!=0 && r!=0) {
				double besselHalf = Math.sqrt(2.0/(Math.PI*k*r))*Math.sin(k*r);	
				sum = sum + (kpdk - k) * k * fk[i][1] * besselHalf/Math.sqrt(k*r);
			}
			else {
				sum = sum + (kpdk - k) * k * fk[i][1] * 0.797885;
			}
		}
		return sum/Math.pow(2*Math.PI, 1.5);				
	}

		/** Fourier transform for 2D radial functions Phys Rep eqn. 10
 *      @param fr double[][], the orginal function
 *      @param k double, the Fourier space input
 *      @return double, the inverse fourier transform ftilda(k) */
        public static double fourier2D(double[][] fr, double k) throws IOException {
                double sum = 0;
                double[][] besselData = importData("BesselJ0.dat");
                for (int i = 0; i < fr.length - 1; i++) {
                        double r = fr[i][0];
                        double rpdr = fr[i+1][0]; /**k+dk*/
                        sum = sum + (rpdr - r) * r * fr[i][1] * evaluate(besselData, k*r);
                }
                return sum*2*Math.PI;
	}	

	public static double[][] fourier2D(double[][] fr, double binSize, double range) throws IOException {
		double[][] result = new double[(int)(range/binSize)][2];
		for (int i = 0; i < result.length; i++) {
			double k = (i + 1)*binSize;
			result[i][0] = k;
			result[i][1] = fourier2D(fr, k);
		}
		return result;
	}

	/**inverse Fourier transform for 2D radial functions Phys Rep eqn. 10
 *	@param fk double[][], the fourier transformed function
 *	@param r double, the direct space input
 *	@return double, the inverse fourier transform f(r) */
	public static double inverseFourier2D(double[][] fk, double r) throws IOException {
		double sum = 0;
		double[][] besselData = importData("BesselJ0.dat");
		for (int i = 0; i < fk.length - 1; i++) {
			double k = fk[i][0];
			double kpdk = fk[i+1][0]; /**k+dk*/			
			sum = sum + (kpdk - k) * k * fk[i][1] * evaluate(besselData, k*r);
		}
		return sum/(2*Math.PI);				
	}

	public static double fourierOfUniform2D(double a, Vector3D k) {
		double sum = 0;
		double height = 1/(a*a);
		int bins = 20;
		double binSize = a/(double)bins;
		double dxdy = binSize*binSize;
		for (int i = 0; i < bins; i++) {
			double rx = -0.5*a + (0.5 + i)*binSize;
			for (int j = 0; j < bins; j++) {
				double rj = rx = -0.5*a + (0.5 + j)*binSize;
				Vector3D rVec = new Vector3D(rx, rj, 0);
				double kr = rVec.innerProd(k);
				sum += height*Math.cos(kr)*dxdy;
			}
		}
		return sum*sum;
	}	


        /**a quicker sine function than Math.sin
 *      @param angle double, angle in radians
 *      @return double*/
        public static double sin(double angle) {
                if (angle >= 0) {
                        angle = angle%(2.0*Math.PI);
                        if (0 <= angle && angle <= Math.PI/2.0) {
                                double x = angle/(Math.PI/2.0);
                                double x2 = x*x;
                                double xsin = ((((.00015148419 * x2 - .00467376557) * x2 + .07968967928) * x2 - .64596371106) * x2 + 1.57079631847) * x;
                                return xsin;
                        }
                        else if (angle <= Math.PI) {
                                return sin(Math.PI - angle);
                        }
                        else {
                                return -sin(angle - Math.PI);
                        }
                }
                else {
                        return -sin(-angle);
                }
        }

        /**a quicker consine
 *      @param angle double, angle in radians
 *      @return double*/
        public static double cos(double angle) {
                return sin(angle + Math.PI/2.0);
        }

	/**computes the median of an array of doubles
 * 	@param numArr double[]
 * 	@return double*/
	public static double median(double[] numArray) {
                Arrays.sort(numArray);
                if (numArray.length % 2 == 0) {
                        return ((double)numArray[numArray.length/2] + (double)numArray[numArray.length/2 - 1])/2;
                }
                else {
                        return (double) numArray[numArray.length/2];
                }
        }

	/**prints an array like I want
 * 	@param f double[][]
 * 	@return String */
	public static String funcToString(double[][] f) {
		String result = "";
		for (int i = 0; i < f.length; i++) {
			for (int j = 0; j < f[0].length; j++) {
				result = result + f[i][j] + " ";
			}
			result = result + "\n";
		}
		return result;
	}

	/**prints an array like I want
 *      @param f double[]
 *       @return String */
        public static String funcToString(double[] f) {
                String result = "";
                for (int i = 0; i < f.length; i++) {
                        result = result + f[i];
                        result = result + "\n";
                }
                return result;
        }	
}
