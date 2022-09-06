package toShare;
import java.util.*;
import java.io.*;
import org.apache.commons.math3.special.BesselJ;
import org.apache.commons.math3.special.Gamma;

public class StrucFac {

	public static double pi = Math.PI;

	public static void main(String[] args) throws IOException, InterruptedException {
		int n = 2500;
		double rho = 0.390625;
		double d = Math.sqrt(n/rho);
		Box[] boxArr = new Box[20];
		for (int i = 0; i < boxArr.length; i++) {
			double[][] coord 
			= Functional.importData("gesConfig/alpha_2_config_" + i + ".txt", 2, 2502);
			for (int j = 0; j < coord.length; j++) {
				for (int k = 0; k < 2; k++) {
					coord[j][k] /= 10.0;
				}
			}
                        boxArr[i] = new Box(n, d, coord);
                }
		double[][] sk = skEquilParallel(boxArr);
                System.out.println(Functional.funcToString(sk));
	}

	/**compute the average structure factor for an array of boxes, parallelized version
 * 	@param boxArr Box[], the array of boxes
 * 	@param nThread int, the number of threads
 *      @return double[][], the structure factor average*/
	public static double[][] skEquilParallel(Box[] boxArr) throws InterruptedException {
		Thread[] threadArr = new Thread[boxArr.length];
		StrucFacRunnable[] runnableArr = new StrucFacRunnable[boxArr.length];
		for (int i = 0; i < boxArr.length; i++) {
			runnableArr[i] = new StrucFacRunnable(boxArr[i]);
			threadArr[i] = new Thread(runnableArr[i]);
			threadArr[i].start();
		}
		for (int i = 0; i < boxArr.length; i++) {
			threadArr[i].join();
		}
		double[][] result = runnableArr[0].getS();
		for (int i = 1; i < boxArr.length; i ++) {
			for (int j = 0; j < result.length; j++) {
                        	result[j][1] += runnableArr[i].getS()[j][1];
			}
                }
		for (int j = 0; j < result.length; j++) {
                        result[j][1] /= (double)boxArr.length;
                }
		return result;
	}

	/**compute the average structure factor for an array of boxes, serial
 * 	@param boxArr Box[], the array of boxes
 * 	@return double[][], the structure factor average*/
	public static double[][] skEquil(Box[] boxArr) {
		Box b = boxArr[0];
		double[][] skTemp = strucFacReport2D(b);
                double[][] sk = new double[skTemp.length][2];
		for (int i = 0; i < sk.length; i++) {
                        sk[i][0] = skTemp[i][0];
                        sk[i][1] = skTemp[i][1];
                }
		for (int i = 1; i < boxArr.length; i++) {
                        b = boxArr[i];
                        skTemp = strucFacReport2D(b);
                        for (int j = 0; j < sk.length; j++) {
                                sk[j][1] = sk[j][1] + skTemp[j][1];
                        }
                }
		for (int i = 0; i < sk.length; i++) {
                        sk[i][1] = sk[i][1]/boxArr.length;
		}
		return sk;
	}

	/**generates a less bumpy angular averaged s(k) with binning...
 * 	@param raw double[][], the crude list of sk data without binning
 *	@param binSize the bin size of the smoothed sk
 *	@param plotRange the max k value for the sk plot
 * 	@return double[][], the table of angular averaged s(k) with binning */
	private static double[][] skAverage(double[][] raw, double binSize, double plotRange) {
		int outputSize = (int) (plotRange/binSize);
		double[][] result = new double[outputSize][2];
		double kMini = raw[0][0];
		int currentBin = 0;
		double currentLength = kMini;
		result[0][0] = currentLength;
		int countInBin = 0;
		for (int i = 0; i < raw.length; i++) {
			double entryLength = raw[i][0];
			if (entryLength - currentLength < binSize) {
				countInBin++;
				result[currentBin][1] += raw[i][1];
			}
			else {
				result[currentBin][1] = result[currentBin][1]/countInBin;
				if (currentBin < outputSize - 1) {//preventing index overflow
					currentBin ++;
					currentLength = entryLength;
					result[currentBin][0] = currentLength;
					result[currentBin][1] = raw[i][1];
					countInBin = 1;
				}
				else {
					break;
				}
			}
		}
		result[currentBin][1] = result[currentBin][1]/countInBin;
		for (int i = currentBin; i < result.length; i++) {
			result[i][0] = result[currentBin][0] + binSize*(i - currentBin);
		}
		return result;
	}
	
	/**generates the table of ANGULAR AVERAGED s(k) for various k in 2D
*       @param b Box, which specifies the configuration of the particles
*       @return double[][], the table of angular averaged s(k) */
        public static double[][] strucFacReport2D(Box b) {
                double[][] raw = new double[2850][2];//25000
                Vector3D k = new Vector3D(0.0,0.0,0.0);
                double length = 0;
                double kMini = 2*Math.PI/b.getD();/**the length of the smallest meaningful k*/

                double n = 0;/**the square of distance*/
                for (int i = 0; i < raw.length; i++) {
                        boolean meaningfulN = false;
                        int multiplicityN = 0;
                        while (! meaningfulN) {
                                n++;
                                length = Math.sqrt(n)*kMini;
                                int searchRange = (int)Math.floor(Math.sqrt(n));
                                for (int jx = -searchRange; jx <= searchRange; jx++) {
					int jy = (int)Math.round(Math.sqrt((double)(n - jx*jx)));
                                        if (jx*jx + jy*jy == n) {
                                                meaningfulN = true;
                                                k.setX(kMini*jx);
                                                k.setY(kMini*jy);
                                                raw[i][1] += sk(b,k);
						k.setY(-kMini*jy);
						raw[i][1] += sk(b,k);
                                                multiplicityN = multiplicityN + 2;
                                        }
                                }
                        }
                        raw[i][0] = length;
                        raw[i][1] = raw[i][1]/multiplicityN;
                }
                return skAverage(raw, 0.05, 15);
        }

	 /**generates the table of ANGULAR AVERAGED s(k) for various k in 2D
 *      @param b Box, which specifies the configuration of the particles
 *      @param rawLength int, the range to be computed for S(k)
 *      @return double[][], the table of angular averaged s(k) */
        public static double[][] strucFacReport2D(Box b, double bigK) {
                double[][] raw = new double[25000][2];//25000
                Vector3D k = new Vector3D(0.0,0.0,0.0);
                double length = 0;
                double kMini = 4*Math.PI/b.getD();/**the length of the smallest meaningful k*/
                double n = 0;/**the square of distance*/
		int vecCount = 0;
                while (length < bigK + 0.1) {
                        boolean meaningfulN = false;
                        int multiplicityN = 0;
                        while (! meaningfulN) {
                                n++;
                                length = Math.sqrt(n)*kMini;
                                int searchRange = (int)Math.floor(Math.sqrt(n));
                                for (int jx = -searchRange; jx <= searchRange; jx++) {
                                        int jy = (int)Math.round(Math.sqrt((double)(n - jx*jx)));
                                        if (jx*jx + jy*jy == n) {
                                                meaningfulN = true;
                                                k.setX(kMini*jx);
                                                k.setY(kMini*jy);
                                                raw[vecCount][1] = raw[vecCount][1] + sk(b,k);
                                                k.setY(-kMini*jy);
                                                raw[vecCount][1] = raw[vecCount][1] + sk(b,k);
                                                multiplicityN = multiplicityN + 2;
                                        }
                                }
                        }
                        raw[vecCount][0] = length;
                        raw[vecCount][1] = raw[vecCount][1]/multiplicityN;
                	vecCount++;
		}
		double[][] result = new double[vecCount][2];
		for (int i = 0; i < vecCount; i++) {
			for (int j = 0; j < 2; j++) {
				result[i][j] = raw[i][j];
			}
		}
                return Functional.standardize(skAverage(raw, 0.05, bigK), 0.05, kMini, bigK);
        }

        /**generates the table of SPHERE AVERAGED s(k) for various k in 3D
 *      @param b Box, which specifies the configuration of the particles
 *      @return double[][], the table of sphere averaged s(k), where k is reported as multiples of pi */
        public static double[][] strucFacReport3D(Box b) {
                double[][] result = new double[601][2];
                Vector3D k = new Vector3D(0.0,0.0,0.0);
                double length = 0;
                double kMini = 2*pi/b.getD();//the length of the smallest meaningful k

                double n = 0;
                for (int i = 0; i < result.length; i++) {
                        boolean meaningfulN = false;
                        int multiplicityN = 0;
                        while (! meaningfulN) {
                                n++;
                                length = Math.sqrt(n)*kMini;
                                int searchRange = (int)Math.floor(Math.sqrt(n));
                                for (int jx = -searchRange; jx <= searchRange; jx++) {
                                for (int jy = -searchRange; jy <= searchRange; jy++) {
                                for (int jz = -searchRange; jz <= searchRange; jz++) {
                                        if (jx*jx + jy*jy + jz*jz == n) {
                                                meaningfulN = true;
                                                k.setX(kMini*jx);
                                                k.setY(kMini*jy);
                                                k.setZ(kMini*jz);
                                                result[i][1]=result[i][1] + sk(b,k);
                                                multiplicityN++;
                                        }
                                }
                                }
                                }
                        }
                        result[i][0] = length;
                        result[i][1] = result[i][1]/multiplicityN;
                }
                return result;
        }


	/**computes the structure factor s(k) for a Cell config
 * 	@param b Box, which specifies the configuration of the particles
 * 	@param k Vector3D, the argument wavevector
 * 	@return double, the structure factor */
	public static double sk(Box b, Vector3D k) {
		int n = b.getN(); //the numberr of particles
		Particle[] pa = b.toArray();
		Complex sum = new Complex(0.0);
		double sumRe = 0;
		double sumIm = 0;
		for (int i = 0; i < n; i++) {
			Particle r = pa[i];
			Vector3D rVec = new Vector3D(r.getx(),r.gety(),0);
			double arg = -k.innerProd(rVec);
			sumRe += Functional.cos(arg);
			sumIm += Functional.sin(arg);
		}
		double result = (sumRe*sumRe + sumIm*sumIm)/(double)n;
		return result;
	}

	public static double[][] dS(Box b, int ind, Vector3D k, double delta) {
		Vector3D rVec = b.particleAt(ind).toVector3D();
		double arg = -k.innerProd(rVec);
		double prefac = 2.0/b.getN()*(Functional.cos(arg) - Functional.sin(arg));
		double[][] result = new double[1][2];
		result[0][0] = prefac*k.getx()*delta;
		result[0][1] = prefac*k.gety()*delta;
		return result;
	}
	
	/**computes a term of the structure factor
 * 	@param k Vector3D, the argument wavevector
 * 	@param r Particle, the particle in the term
 * 	@return Complex, the term */
	private static Complex term(Vector3D k, Particle r) {
		Vector3D rVec = new Vector3D(r.getx(),r.gety(),0);
		Complex arg = new Complex(-k.innerProd(rVec));
		Complex exponent = Complex.i().times(arg);
		return exponent.expThis();	
	}

	/**converts S(k) to the spectral density \tilde{chi}(k)
 *	See Phys Rep Torquato2018 eq. (142) 
 *  	 * @param s double[][] structure factor
 *  	 * @param rho double number density
 *  	 * @param a double, radius of the spheres that decorate the points
 *  	 * @param dim double, dimension
 *  	 * @return double[][]
 * 	 */
	public static double[][] sToChi(double[][] s, double rho, double a, double dim) {
		double[][] result = new double[s.length][2];
		double v1 = Math.pow(Math.PI, dim/2)*Math.pow(a, dim)/Gamma.gamma(1 + dim/2);
		double phi2 = v1*rho;
		for (int i = 0; i < s.length; i++) {
			double k = s[i][0];
			double sk = s[i][1];
			double alpha2Tilde = 
				Math.pow(2.0, dim)*Math.pow(Math.PI, dim/2)*Gamma.gamma(1 + dim/2)*Math.pow(BesselJ.value(dim/2.0, k*a), 2)/Math.pow(k, dim);
			if (k == 0) {
				alpha2Tilde = v1;
			}
			result[i][0] = k;
			result[i][1] = phi2*alpha2Tilde*sk;
		}
		return result;
	}
}

class StrucFacRunnable implements Runnable {
        Box b;
        double[][] s;

        /**constructor, accepting all the relevant input parameters
         * @param boxArr box[], array of box
         */
        public StrucFacRunnable(Box bb) {
                b = bb;
        }

        /**fetch out result
         * @return boxArr
         */
        public double[][] getS() {
                return s;
        }

        public void run() {
                s = StrucFac.strucFacReport2D(b);
        }

}

