package toShare;
import java.util.*;
import java.io.*;

public class Exclusion {
	public static void main(String[] args) throws IOException {
		/**double rho = 0.390625;
		int n = 2500;
		double d = Math.sqrt((double)n/rho);
		double[][] coord = InitialGuess.importData("RSA2DIS_BFGS.snap");
		Box b = new Box(n, d, coord);
		for (double r = 0.01; r <= 3; r += 0.01) {
			System.out.println(r + " " + lz(b, 1, r));
		}*/
		int n = 2500;
		double rho = 0.390625;
		double d = Math.sqrt((double)n/rho);
		Box[] boxArr = new Box[433];
		for (int i = 0; i < boxArr.length; i++) {
			double[][] coord;
			if (i < 133) {
				coord = Functional.importData("pgAlpha1/targ/collectiveCoordN2500*_" + (i + 1) +".out", 1201, 3701);
			}
			else if (i < 269) {
				coord = Functional.importData("pgAlpha1/targ/collectiveCoordN2500_" + (i - 133 + 1) +".out", 2401, 4901);
			}
			else {
				coord = Functional.importData("pgAlpha1/targ/collectiveCoordN2500**_" + (i - 269 + 1) +".out", 2501, 5001);
			}
			//double[][] coord = InitialGuess.importData("pgAlpha1/ewald2_" + (i + 33) +".out", 4, 2504);
			boxArr[i] = new Box(n, d, coord);
		}
		double[] rGrid = new double[100];
                for (int i = 0; i < rGrid.length; i++) {
                        rGrid[i] = i*0.05;
                }
                double[][] hv = Exclusion.hvDirect(boxArr, rGrid, 3000000, 0.02);
		System.out.println(Functional.funcToString(hv));
		double[][] ev = Exclusion.ev(boxArr, 0.05, 5);
		System.out.println(Functional.funcToString(ev));
		/**double[][] hpp = hp(boxArr, 0.05, 3);
                for (int i = 0; i < hpp.length; i++) {
                        System.out.println(hpp[i][0]+" "+hpp[i][1]);
                }*/
	}

	/**particle exclusion function, monodisperse
 * 	@param b Box
 * 	@param r double, the radius of the test window
 * 	@return double E_p(r)*/
	public static double ep(Box b, double r) {
		double[] ndl = nearestDistanceList(b);	
		int voidCount = 0;
		int n = b.getN();
		for (int i = 0; i < n; i++) {
			if (r < ndl[i]) {
				voidCount ++;
			}
		}
		return (double)voidCount/(double)n;
	}

	/**generate a list of the nearest neighbor distances wrt every particle
 * 	@param b Box
 * 	@return double[], the neareast distance list*/
	public static double[] nearestDistanceList(Box b) {
		Particle p, q;
		Particle[] partiArr = b.toArray();
		int n = b.getN();
		double[] result = new double[n];
		for (int i = 0; i < n; i++) {
			p = partiArr[i];
			result[i] = b.getD();
			for (int j = 0; j < n; j++) {
				if (j != i){
					q = partiArr[j];
					double distPQ = b.minDist(p, q);
					if (distPQ < result[i]) {
						result[i] = distPQ;
					}
				}
			}
		}
		return result;	
	}

	/**generate a list of nearest neighbor distances from every test particle to every real particle
	 * @param textB Box
	 * @param b Box, box of particles
	 * @return double[] the nearest distance list
	 */
        public static double[] nearestNeighborList(Box testB, Box b, double range) {
                int nProbe = testB.getN();
                int n = b.getN();
                double d = b.getD();
                double[] result = new double[nProbe];
                int nCells = (int)(d*2/range);
                double cellLength = d/nCells;
                Particle[][][] cellList = new Particle[nCells][nCells][100];
                int[][] cellCount = new int[nCells][nCells];
                int[][] cellOfProbeParticle = new int[nProbe][2];
                for (int i = 0; i < n; i++) {
                        Particle p = b.particleAt(i);
                        int ix = (int)(p.getx()/cellLength);
                        int iy = (int)(p.gety()/cellLength);
                        cellList[ix][iy][cellCount[ix][iy]] = p;
                        cellCount[ix][iy]++;
                }
                for (int i = 0; i < nProbe; i++) {
                        Particle p = testB.particleAt(i);
                        int ix = (int)(p.getx()/cellLength);
                        int iy = (int)(p.gety()/cellLength);
                        cellOfProbeParticle[i][0] = ix;
                        cellOfProbeParticle[i][1] = iy;
                }
                for (int i = 0; i < nProbe; i++) {
                        Particle p = testB.particleAt(i);
                        double nearestDist = d;
                        int ix = cellOfProbeParticle[i][0];
                        int iy = cellOfProbeParticle[i][1];
                        result[i] = d;
                        for (int jx = ix - 1; jx < ix + 2; jx++) {
                        for (int jy = iy - 1; jy < iy + 2; jy++) {
                                int jxMod = (jx + nCells)%nCells;
                                int jyMod = (jy + nCells)%nCells;
                                for (int j = 0; j < cellCount[jxMod][jyMod]; j++) {
                                        Particle q = cellList[jxMod][jyMod][j];
                                        double distPQ = b.minDist(p, q);
                                        if (distPQ < result[i]) {
                                                result[i] = distPQ;
                                        }
                                }
                        }
                        }
                }
                return result;
        }

	/**directly computes hv for a box
	 * 
	 * @param b
	 * @param rGrid
	 * @param testN
	 * @param binSize
	 * @return
	 */
	public static double[][] hvDirect(Box b, double[] rGrid, int testN, double binSize) {
                double d = b.getD();
                Box testB = new Box(testN, d);
                double[] nl = nearestNeighborList(testB, b, rGrid[rGrid.length - 1]);
                double halfBin = 0.5*binSize;
                double[][] result = new double[rGrid.length][2];
                for (int j = 0; j < rGrid.length; j++) {
                        double r = rGrid[j];
                        result[j][0] = r;
                        int count = 0;
                        for (int i = 0; i < testN; i++) {
                                if (nl[i] > r - halfBin && nl[i] <= r + halfBin) {
                                        count++;
                                }
                        }
                        result[j][1] = (double)count/((double)testN*binSize);
                }
                return result;
        }

	/**directly compute hv for an array of boxes*/
	public static double[][] hvDirect(Box[] boxArr, double[] rGrid, int testN, double binSize) {
		Box b = boxArr[0];
		double[][] hvtemp = hvDirect(b, rGrid, testN, binSize);
		double[][] hv = new double[hvtemp.length][2];
		for (int i = 0; i < hv.length; i++) {
			hv[i][0] = hvtemp[i][0];
			hv[i][1] = hvtemp[i][1];
		}
		for (int i = 1; i < boxArr.length; i++) {
			System.out.println(i);
			b = boxArr[i];
			hvtemp = hvDirect(b, rGrid, testN, binSize);
			for (int j = 0; j < hv.length; j++) {
				hv[j][1] = hv[j][1] + hvtemp[j][1];
			}
		}
		for (int i = 0; i < hv.length; i++) {
			hv[i][1] = hv[i][1]/boxArr.length;
		}
		return hv;
	}


	/**directly computes hv for a box
 *      *
        * @param b
 *      * @param rGrid
 *      * @param binSize
 *      * @return
 *      */
        public static double[][] hpDirect(Box b, double[] rGrid, double binSize) {
                int n = b.getN();
		double d = b.getD();
                double[] nl = nearestDistanceList(b);
                double halfBin = 0.5*binSize;
                double[][] result = new double[rGrid.length][2];
                for (int j = 0; j < rGrid.length; j++) {
                        double r = rGrid[j];
                        result[j][0] = r;
                        int count = 0;
                        for (int i = 0; i < n; i++) {
                                if (nl[i] > r - halfBin && nl[i] <= r + halfBin) {
                                        count++;
                                }
                        }
                        result[j][1] = (double)count/((double)n*binSize);
                }
                return result;
        }

	public static double[][] hpDirect(Box[] boxArr, double[] rGrid, double binSize) {
                Box b = boxArr[0];
                double[][] hptemp = hpDirect(b, rGrid, binSize);
                double[][] hp = new double[hptemp.length][2];
                for (int i = 0; i < hp.length; i++) {
                        hp[i][0] = hptemp[i][0];
                        hp[i][1] = hptemp[i][1];
                }
                for (int i = 1; i < boxArr.length; i++) {
                        System.out.println(i);
                        b = boxArr[i];
                        hptemp = hpDirect(b, rGrid, binSize);
                        for (int j = 0; j < hp.length; j++) {
                                hp[j][1] = hp[j][1] + hptemp[j][1];
                        }
                }
                for (int i = 0; i < hp.length; i++) {
                        hp[i][1] = hp[i][1]/boxArr.length;
                }
                return hp;
        }


	/**generate a histogram of the nearest neighbor distances wrt every particle
 *      @param b Box
 *      @param binsize double
 *      @param range double
 *      @return double[], the neareast distance histogram*/
        public static double[][] hp(Box b, double binSize, double range) {
		int nBin = (int)(range/binSize);
		double[][] result = new double[nBin][2];
		for (int i = 0; i < nBin; i++) {
			result[i][0] = (0.5 + (double)i)*binSize;
			result[i][1] = 0;
		}
                Particle p, q;
                Particle[] partiArr = b.toArray();
                int n = b.getN();
                for (int i = 0; i < n; i++) {
                        p = partiArr[i];
                        double nearestDist = b.getD();
                        for (int j = 0; j < n; j++) {
                                if (j != i){
                                        q = partiArr[j];
                                        double distPQ = b.minDist(p, q);
                                        if (distPQ < nearestDist) {
                                                nearestDist = distPQ;
                                        }
                                }
                        }
			int binToEnter = (int)(nearestDist/binSize);
			if (binToEnter < result.length) {
				result[binToEnter][1] = result[binToEnter][1] + 1.0;
			}
                }
		for (int i = 0; i < nBin; i++) {
                        result[i][1] = result[i][1]/((double)n*binSize);
                }
                return result;
        }

	/**computes the pore size distribution function, E_V(r)
 * 	@param b Box
 * 	@param binsize double
 * 	@param range double
 * 	@return double[], E_V(r)*/
	public static double[][] ev(Box b, double binSize, double range) {
		int nBin = (int)(range/binSize);
                double[][] h = new double[nBin][2];
                for (int i = 0; i < nBin; i++) {
                        h[i][0] = (0.5 + (double)i)*binSize;
                        h[i][1] = 0;
                }
                Particle[] partiArr = b.toArray();
                int n = b.getN();
		double d = b.getD();
		int probeN = 3000000;
		Box probe = new Box(probeN, d);/**a random box of probing particles*/
		double[] nl = nearestNeighborList(probe, b, range);
		for (int i = 0; i < probeN; i++) {
                        int binToEnter = (int)(nl[i]/binSize);
                        if (binToEnter < nBin) {
                                h[binToEnter][1] = h[binToEnter][1] + 1.0;
                        }
                }
	 	for (int i = 0; i < nBin; i++) {
                        h[i][1] = h[i][1]/((double)probe.getN()*binSize);
                }
		double[][] result = new double[nBin][2];
		for (int i = 0; i < nBin; i++) {
			result[i][0] = (double)i*binSize;
			double accu = 0;
			for (int j = 0; j < i; j++) {
				accu += h[j][1]*binSize;
			}
                        result[i][1] = 1.0 - accu;
                }
                return result;
	}

	/**generate the particle nearest neighbor distribution function for an array of Boxes
 * 	@param boxArr Box[]
 * 	@param binsize double
 * 	@param range double 
 * 	@return double[] */
	public static double[][] hp(Box[] boxArr, double binSize, double range) {
		Box b = boxArr[0];
                double[][] e2temp = hp(b, binSize, range);
                double[][] e2 = new double[e2temp.length][2];
                for (int i = 0; i < e2.length; i++) {
                        e2[i][0] = e2temp[i][0];
                        e2[i][1] = e2temp[i][1];
                }
                for (int i = 1; i < boxArr.length; i++) {
                        b = boxArr[i];
                        e2temp = hp(b, binSize, range);
                        for (int j = 0; j < e2.length; j++) {
                                e2[j][1] = e2[j][1] + e2temp[j][1];
                        }
                }
                for (int i = 0; i < e2.length; i++) {
                        e2[i][1] = e2[i][1]/boxArr.length;
                }
                return e2;
	}

        /**computes the pore size distribution function e_V(r) for an array of Boxes
 *      @param boxArr Box[]
 *      @param binsize double
 *      @param range double
 *      @return double[]*/
	public static double[][] ev(Box[] boxArr, double binSize, double range) {
                Box b = boxArr[0];
                double[][] e2temp = ev(b, binSize, range);
                double[][] e2 = new double[e2temp.length][2];
                for (int i = 0; i < e2.length; i++) {
                        e2[i][0] = e2temp[i][0];
                        e2[i][1] = e2temp[i][1];
                }
                for (int i = 1; i < boxArr.length; i++) {
			System.out.println(i);
                        b = boxArr[i];
                        e2temp = ev(b, binSize, range);
                        for (int j = 0; j < e2.length; j++) {
                                e2[j][1] = e2[j][1] + e2temp[j][1];
                        }
                }
                for (int i = 0; i < e2.length; i++) {
                        e2[i][1] = e2[i][1]/boxArr.length;
                }
                return e2;
        }
	
	/**computes the lineal path function L_z(r) at a fixed variable r
 *	@param b Box
 *	@param rad double, radius of the spheres
 *	@param r the rod length variable
 *	@return double*/
	public static double lz(Box b, double rad, double r) {
		int trials = 20000;
		int n = b.getN();
		double d = b.getD();
		Particle[] partiArr = b.toArray();
		int success = trials; // number of trial rods wholly contained in the matrix phase, i.e. don't overlap with any particle
		for (int i = 0; i < trials; i++) {
			/**generate a rod with length r*/
			double x1 = Box.getRandomNumberInRange(0, d);
			double y1 = Box.getRandomNumberInRange(0, d);
			double theta = Box.getRandomNumberInRange(0, 2*Math.PI);
			double x2 = x1 + r*Math.cos(theta);
			double y2 = y1 + r*Math.sin(theta);
			/**check if any particle overlaps with the rod*/
			for (int j = 0; j < n; j++) {
				Particle p = partiArr[j];
				if (p.overlapWithLine(x1, x2, y1, y2, rad, d)) {
					success--;
					break;
				}
			}
		}
		return (double)success/(double)trials;
	}
}


