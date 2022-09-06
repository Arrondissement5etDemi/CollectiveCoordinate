package toShare;
import java.util.*;
import java.io.*;
import org.apache.commons.math3.special.BesselJ;
import org.apache.commons.math3.special.Gamma;


/** A class with functions to sample radial distribution function */
public class RadialStat {
	
	private static double pi = Math.PI;

	/**compute the average g2 for an array of boxes
 *      @param boxArr Box[], the array of boxes
 *      @param binSize double
 *      @return double[][], the g2 average*/
        public static double[][] g2Equil(Box[] boxArr, double binSize) {
                Box b = boxArr[0];
		double range = b.getD()/2.0;
                double[][] g2Temp = g22d(b, binSize, range);
		double[][] allData = new double[g2Temp.length][boxArr.length];
                double[][] g2 = new double[g2Temp.length][2];
		for (int i = 0; i < g2.length; i++) {
                        g2[i][0] = g2Temp[i][0];
                        g2[i][1] = g2Temp[i][1];
                }
		/**for (int i = 0; i < boxArr.length; i++) {
                        g2Temp = g22d(b, 0.05, range);
                        for (int j = 0; j < g2.length; j++) {
                                allData[j][i] = g2Temp[j][1];
                        }
                }
                for (int j = 0; j < g2.length; j++) {
                        g2[j][1] = Functional.median(allData[j]);
                }*/
                for (int i = 1; i < boxArr.length; i++) {
                        b = boxArr[i];
                        g2Temp = g22d(b, binSize, range);
                        for (int j = 0; j < g2.length; j++) {
                                g2[j][1] = g2[j][1] + g2Temp[j][1];
                        }
                }
                for (int i = 0; i < g2.length; i++) {
                        g2[i][1] = g2[i][1]/boxArr.length;
                }
                return g2;
        }

	/**computes g2 for 2D structures using binning
 *      @param b Box, the Particle configuration that we consider
 *      @param numBinsPerDiam int, how many bins are there in the diameter of the particles
 *      @param upperRad double, the upper radius in multiples of diam of the range we consider in binning
 *      @return double[][], the g2 function */
        public static double[][] g22d(Box b, double binSize, double upperRad) {
                Particle[] partiArr = b.toArray();
                int n = b.getN();
                double d = b.getD();
                double rho = //(double)4.0*n/(Math.PI*d*d);
				(double)n/(d*d);//number density

                int numBins = (int) (upperRad/binSize);
                int[] binningResult = new int[numBins];
                for (int i = 0; i < n; i++) {
                        Particle center = partiArr[i];
			for (int j = i + 1; j < n; j++) {
				Particle thatParti = partiArr[j];
				double distScaled = b.minDist(center, thatParti);
				//double distScaled = center.distanceto(thatParti);
				int binToEnter = (int)(distScaled/binSize);
				if (binToEnter < binningResult.length) {
                         		binningResult[binToEnter]++;
                        	}
			}
                }
                double[][] result = new double[numBins][2];
                for (int i = 0; i < numBins; i++) {
                        result[i][0] = (i + 0.5) * binSize;
                        double shellVolume = 2.0*Math.PI*result[i][0]*binSize;
                        result[i][1] = (double)binningResult[i]*2/((double)(n)*rho*shellVolume);
                }
                return result;
        }

	/**computes the autocovariance function chi for a box of decorated spheres
 * 	@param Box Box, Particle configuration
 * 	@param a double, sphere radius that decorate the particles
 * 	@param binSize double
 * 	@param range double
 * 	@return chi, the autocovariance function*/
	public static double[][] s2(Box b, double a, double binSize, double range) {
		int testN = 30000;
		int nBins = (int)(range/binSize);
		double[][] result = new double[nBins][2];
		double d = b.getD();
		double phi2 = b.getN()/(d*d)*Math.PI*a*a;
		Box testB1 = new Box(testN, d);
		Particle[] testB1Array = testB1.toArray();
		double[] randomAngles = new double[testN];
		for (int i = 0; i < testN; i++) {
			randomAngles[i] = Box.getRandomNumberInRange(0, 2*Math.PI); 
		}
		for (int i = 0; i < nBins; i++) {
			double r = i*binSize;
			result[i][0] = r;
			double[][] testB2Coord = new double[testN][2];
			for (int j = 0; j < testN; j++) {
				testB2Coord[j][0] = (testB1Array[j].getx() + r*Math.cos(randomAngles[j]) + d)%d;
				testB2Coord[j][1] = (testB1Array[j].gety() + r*Math.sin(randomAngles[j]) + d)%d;
			}
			Box testB2 = new Box(testN, d, testB2Coord);
			double[] neighborList1 = Exclusion.nearestNeighborList(testB1, b, range);
			double[] neighborList2 = Exclusion.nearestNeighborList(testB2, b, range);
			double bothInPhase2Count = 0;
			for (int j = 0; j < testN; j++) {
				if (neighborList1[j] < a && neighborList2[j] < a) {
					bothInPhase2Count++;
				}
			}
			result[i][1] = bothInPhase2Count/(double)testN - phi2*phi2;
		}
		return result;
	}

	/**computes the spreadability S(t) from the autocovariance function
 * 
 *
 * */
	public static double spreadability(double t, double[][] chi, double phi2, double diffus, double dim) {
		if (t == 0) {
                        return 0;
                }
		double phi1 = 1.0 - phi2;
		double wd =  Math.pow(Math.PI, dim/2)/Gamma.gamma(1 + dim/2);
		double prefactor = dim*wd/(Math.pow(4.0*Math.PI*diffus*t, dim/2.0)*phi2);
		double phi1MinusSt = 0;
		for (int i = 0; i < chi.length - 1; i++) {
			double r = chi[i][0];
			double dr = chi[i + 1][0] - chi[i][0];
			phi1MinusSt += Math.pow(r, dim - 1.0)*chi[i][1]*Math.exp(-r*r/(4*diffus*t))*dr;
		}
		return phi1 - prefactor*phi1MinusSt;
	}
	

	/**computes the directo correlation function
 *	@param g2 double[][]
 *	@param sk double[][]
 *	@param rho double, density
 *	@return double[][], the direct correlation function c(r)*/
	public static double[][] directCorr(double[][] g2, double[][] sk, double rho) throws IOException {
		double binSize = 0.05;
		double range = g2[g2.length - 1][0];
		g2 = Functional.standardize(g2, binSize, 0.5*binSize, range + 0.5*binSize);
		sk = Functional.standardize(sk, binSize, 0, 100);
		double[][] bracket = new double[sk.length][2];
		for (int i = 0; i < sk.length; i++) {
			bracket[i][0] = sk[i][0];
			double s = sk[i][1];
			bracket[i][1] = //(s - 1.0)*(s - 1.0)/s;//for Heinen's method only ... if not, 
				(s-1)/(rho*s);
		}
		double[][] result = new double[g2.length][2];
		for (int i = 0; i < g2.length; i++) {
			result[i][0] = g2[i][0];
			result[i][1] = //g2[i][1] - 1 - 
				Functional.inverseFourier2D(bracket, g2[i][0]);
			System.out.println(result[i][0] + " " + result[i][1]);
		}
		return result;
	}

}

	
