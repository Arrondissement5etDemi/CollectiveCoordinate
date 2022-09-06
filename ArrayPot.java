package toShare;
import java.util.function.*;
import java.io.*;

public class ArrayPot implements PairPot {
	public static double binSize = 0.05;
	double range = 16;
	double cutoff = 16;
	public static double skRange = 30;
	static double hardCore = 1;
	double[][] arr;
	double[][] fk;
	double[][] spline;//the cubic spline of the potential
	PairPot p;

	public ArrayPot(double[][] f, double r) {
		range = r;
		cutoff = r;
		double[][] arrTemp = Functional.standardize(f, binSize, 0, range);
		double cutoffPot = Functional.evaluate(arrTemp, cutoff);
		arr = new double[arrTemp.length][2];
	        for (int i = 0; i < arr.length; i++) {
			arr[i][0] = arrTemp[i][0];
			if (arr[i][0] < cutoff) {
                        	arr[i][1] = arrTemp[i][1] - cutoffPot;
			}
			else {
				arr[i][1] = 0;
			}
        	}
        	spline = Functional.cubicSpline(arr);
	}

	/**build an ArrayPot given the fourier transform of the potential*/
	public ArrayPot(double[][] f, double rho, int n) throws IOException {
		double kmin = f[0][0];
                fk = Functional.standardize(f, 0.05, kmin, skRange + kmin);
		int size = (int)Math.round(range/binSize);
		double[][] a = new double[size][2];
                for (int i = 0; i < size; i++) {
                        a[i][0] = (i + 1)*binSize;
                        a[i][1] = Functional.inverseFourier3D(fk, a[i][0]);
                }
                arr = Functional.standardize(a, binSize, 0.05, range);
		/**protect the repulsive wall*/
                for (int i = 0; i < arr.length; i++) {
                        if (arr[i][0] <= 1.05) {
                                arr[i][1] = arr[21][1] + (21 - (double)i)*(arr[21][1] - arr[22][1]);
                        }
                }
                /**modify by cutoff*/
                double cutoffPot = arr[arr.length - 1][1];
                for (int i = 0; i < arr.length; i++) {
                        arr[i][1] = arr[i][1] - cutoffPot;
                }
	}

	public ArrayPot(ToDoubleFunction<Double> f, double r) {
		range = r;
		cutoff = r;
		int size = (int)(range/binSize) + 1;
                arr = new double[size][2];
                for (int i = 0; i < size; i++) {
                        arr[i][0] = (i + 0.5)*binSize;
                        arr[i][1] = f.applyAsDouble(arr[i][0]);
                }
                arr = Functional.standardize(arr, binSize, 0.5*binSize, range);
                spline = Functional.cubicSpline(arr);
	}

	public ArrayPot(PairPot pot, double r) {
		range = r;
                cutoff = r;
		int size = (int)(range/binSize) + 1;
		arr = new double[size][2];
		for (int i = 0; i < size; i++) {
			arr[i][0] = (i + 0.5)*binSize;
			arr[i][1] = pot.pairE(arr[i][0]);
		}
		arr = Functional.standardize(arr, binSize, 0.5*binSize, range);
		spline = Functional.cubicSpline(arr);
		p = pot;
	}


	public ArrayPot(double[][] sk, double rho, double t, int n) throws IOException { /**by Fourier transform of an S(k) target*/
		int size = (int)(range/binSize);
		double kmin = 2*Math.PI/Math.sqrt((double)n/rho);
		sk = Functional.standardize(sk, 0.05, kmin, skRange + kmin);
		fk = new double[sk.length][2];
		for (int i = 0; i < sk.length; i ++) {
			fk[i][0] = sk[i][0];
			fk[i][1] = 1.0/sk[i][1] - 1.0;
		}
                double[][] a = new double[size][2];
                for (int i = 0; i < size; i++) {
                        a[i][0] = (i + 1)*binSize;
                        a[i][1] = t*Functional.inverseFourier3D(fk, a[i][0])/rho;
                }
		arr = Functional.standardize(a, binSize, 0.05, range);
		double cutoffPot = arr[arr.length - 1][1];
		for (int i = 0; i < arr.length; i++) {
                        arr[i][1] = arr[i][1] - cutoffPot;
			/**if (arr[i][0] <= 0.2) {
				arr[i][1] = arr[6][1] + (6 - (double)i);
			}*/
                }	
	}

	/**constructor using HNC from target cT*/
	public ArrayPot(double[][] g2, double[][] cT, double t, double r) throws IOException {
		range = r;
		g2 = Functional.standardize(g2, binSize, 0.5*binSize, range + 0.5*binSize);
		cT = Functional.standardize(cT, binSize, 0.5*binSize, range + 0.5*binSize);
		double[][] f = new double[g2.length][2];
		for (int i = 0; i < g2.length; i++) {
			f[i][0] = g2[i][0];
			if (g2[i][1] < 1e-6) {
				g2[i][1] = 1e-6;
			}	
			f[i][1] = (g2[i][1] - 1 - cT[i][1] - Math.log(g2[i][1]))*t;
		}
		range = f[f.length - 1][0];
                cutoff = range;
                double[][] arrTemp = f;
                double cutoffPot = Functional.evaluate(f, cutoff);
                arr = new double[arrTemp.length][2];
                for (int i = 0; i < arr.length; i++) {
                        arr[i][0] = arrTemp[i][0];
                        if (arr[i][0] < cutoff) {
                                arr[i][1] = arrTemp[i][1] - cutoffPot;
                        }
                        else {
                                arr[i][1] = 0;
                        }
                }
		/**protect the repulsive wall*/
		int indEndCore = (int)(hardCore/binSize);
		System.out.println(arr[indEndCore][0] + " " + arr[indEndCore][1]);
                for (int i = 0; i < arr.length; i++) {
                        if (i < indEndCore) {
                                arr[i][1] = arr[indEndCore][1] + (indEndCore - (double)i)*(arr[indEndCore][1] - arr[indEndCore + 1][1]);
                        }
                }
		spline = Functional.cubicSpline(arr);
	}

	public static void setBinSize (double b) {
		binSize = b;
	}

	/**gets range
 * 	@return double*/
	public double getRange() {
		return range;
	}

	/**computes the pair energy
	 * @param r double, the pair distance
	 * @return double, the pair energy */
	public double pairE(double r) {
		/**if (1.0 < r && r < 1.05 && p != null) {
			return p.pairE(r);
		}
		if (r < cutoff) {
			return Functional.evaluateSpline(spline, r);
		}
		return 0; */	
		if (r < arr[arr.length - 1][0]) {
			if (r > arr[0][0]) {
				int j = (int)((r - arr[0][0])/binSize);
				double y2 = arr[j + 1][1];
				double x1 = arr[j][0];
				double y1 = arr[j][1];
				double y = y1 + (r - x1)*(y2 - y1)/binSize;
				return y;
			}
			return arr[0][1];
		}
		return 0;
	}

	/**compute the derivative
 * 	@param r double, the pair distance
 * 	@return double, the derivative of the potential */
	public double derivative(double r) {
		double delta = 0.005;
		return (pairE(r+delta)-pairE(r-delta))/(2*delta);
	}

	public void setArr(double[][] f) {
		arr = Functional.standardize(f, binSize, 0, range);
	}

	public static void setHardCore(double d) {
		hardCore = d;
	}

	public double[][] getArr() {
		return arr;
	}

	/**generate a new guess of potential given an equilibrated g2 with the last guess
 *      @param boxArr Box[]
 *      @param t double, temperature
 *      @param target double[][], the target g2
 *      @return void, the next guess of the potential parameters */
        public void nextGuess(Box[] boxArr, double t, double[][] target) {
                /**generate equilibrated g2*/
                double[][] g2 = Functional.standardize(RadialStat.g2Equil(boxArr, binSize), binSize, 0.5*binSize, range);
		target = Functional.standardize(target, binSize, 0.5*binSize, range);
		for (int i = 0; i < g2.length; i++) {
			System.out.println(g2[i][0]+" "+g2[i][1]);
		}
                double g2Dist = Functional.l2sq(g2, target, 0, range);
                System.out.println("g2 distance = "+g2Dist);
                /**generate the new potential guess*/
                double[][] potPlot = new double[g2.length][2];
                for (int i = 0; i < g2.length; i++) {
                        double r = g2[i][0];
                        potPlot[i][0] = r;
                        if (target[i][1] == 0 || g2[i][1] == 0) {
                                potPlot[i][1] = Double.NaN;
                        }
                        else {
                                potPlot[i][1] = this.pairE(r) + 0.5*t*Math.log(g2[i][1]/target[i][1]);
                        }
                }
		setArr(potPlot);
		/**fix the discontinuity at cutoff*/
                double cutoffPot = arr[arr.length - 1][1];
                double[][] arrTemp = new double[arr.length][2];
                for (int i = 0; i < arr.length; i++) {
                        arr[i][0] = arr[i][0];
                        arr[i][1] = arr[i][1] - cutoffPot;
                }
                /**fix the NaN*/
		int indEndCore = (int)(hardCore/binSize) + 1;
                for (int i = 0; i < arr.length; i++) {
                        if (arr[i][0] <= hardCore) {
                                arr[i][1] = arr[indEndCore][1] + (indEndCore - (double)i)*(arr[indEndCore][1] - arr[indEndCore + 1][1]);
                        }
                }
                spline = Functional.cubicSpline(arr);	
        }

	/**amartya's algorithm about inversion of S(k)
 * 	@param boxArr Box[], array of box
 * 	@param t double, temperature
 * 	@param target double[][] S(k) target
 * 	@param rho number density*/
	public void nextGuess(Box[] boxArr, double t, double[][] target, double rho) throws IOException{
                /**generate equilibrated S(k)*/
		double kmin = fk[0][0];
                double[][] sk = Functional.standardize(StrucFac.skEquil(boxArr), 0.05, kmin, skRange + kmin);
                double[][] s0 = Functional.standardize(target, 0.05, kmin, skRange + kmin);
                for (int i = 0; i < sk.length; i++) {
                        System.out.println(sk[i][0]+" "+sk[i][1]);
                }
                double skDist = Functional.l2sq(sk, s0, 0, skRange);
                System.out.println("sk distance = " + skDist);
                /**generate the new guess of the Fourier transform*/
		System.out.println("new vtilde:");
		for (int i = 0; i < fk.length; i++) {
			if (sk[i][0] < 13) {
				fk[i][1] = fk[i][1] + 10*t*(sk[i][1] - s0[i][1])/(sk[i][1]*s0[i][1]*rho);
			}
			System.out.println(fk[i][0]+" "+fk[i][1]);
		}
		/**inverse Fourier of the potential back to direct space*/
		for (int i = 0; i < arr.length; i++) {
                        arr[i][1] = Functional.inverseFourier3D(fk, arr[i][0]);
                }
		/**protect the repulsive wall*/
		for (int i = 0; i < arr.length; i++) {
                        if (arr[i][0] <= 0.55) {
                                arr[i][1] = arr[11][1] + (11 - (double)i)*(arr[11][1] - arr[12][1]);
                        }
                }
		/**modify by cutoff*/
		double cutoffPot = arr[arr.length - 1][1];
		for (int i = 0; i < arr.length; i++) {
                        arr[i][1] = arr[i][1] - cutoffPot;
                } 
       }

	
	/**Heinen's iterative HNC
 * 	@param boxArr Box[]
 * 	@param t double, temperature
 * 	@param g2Targ double[][]
 *	@param skTarg double[][]
 *	@param cT double[][], the target direct correlation function */
	public void nextGuess(Box[] boxArr, double t, double[][] g2Targ, double[][] skTarg, double[][] cT, double rho) throws IOException, InterruptedException {
		/**compute g2*/
		double[][] g2 = Functional.standardize(RadialStat.g2Equil(boxArr, binSize), binSize, 0.5*binSize, range);
                g2Targ = Functional.standardize(g2Targ, binSize, 0.5*binSize, range);
                for (int i = 0; i < g2.length; i++) {
                        System.out.println(g2[i][0]+" "+g2[i][1]);
                }
                double g2Dist = Functional.l2sq(g2, g2Targ, 0, range)*rho;
                System.out.println("g2 distance = "+g2Dist);
		/**compute sk*/
		double kmin = 2*Math.PI/boxArr[0].getD();
		double[][] sk = Functional.standardize(StrucFac.skEquilParallel(boxArr), 0.05, kmin, skRange + kmin);
                double[][] s0 = Functional.standardize(skTarg, 0.05, kmin, skRange + kmin);
                for (int i = 0; i < sk.length; i++) {
                        System.out.println(sk[i][0]+" "+sk[i][1]);
                }
                double skDist = Functional.l2sq(sk, s0, 0, skRange)/(rho*4.0*Math.PI*Math.PI);
		System.out.println("sk distance = " + skDist);
		/**compute c_i*/
		double[][] ci = RadialStat.directCorr(g2, sk, rho);
		/**update v(r)*/
		for (int i = 0; i < arr.length; i++) {
			double r = arr[i][0];
			double gTr = Functional.evaluate(g2Targ, r);
			double gir = Functional.evaluate(g2, r);
			if (gTr <= 1e-6) {
                                gTr = 1e-6;
                        }
			if (gir <= 1e-6) {
				gir = 1e-6;
			}
			arr[i][1] = arr[i][1] + 0.01*t*(gTr - gir - Functional.evaluate(cT, r) + Functional.evaluate(ci, r) + Math.log(gir/gTr));
		}
		/**fix the discontinuity at cutoff*/
                double cutoffPot = arr[arr.length - 1][1];
                for (int i = 0; i < arr.length; i++) {
                        arr[i][1] = arr[i][1] - cutoffPot;
                }
		/**protect the repulsive wall*/
		int indEndCore = (int)(hardCore/binSize) + 1;
                for (int i = 0; i < arr.length; i++) {
                        if (i < indEndCore) {
                                arr[i][1] = arr[indEndCore][1] + (indEndCore - (double)i)*(arr[indEndCore][1] - arr[indEndCore + 1][1]);
                        }
                }
		spline = Functional.cubicSpline(arr);
	}
	

	public String toString() {
		String result = "";
		for (int i = 0; i < arr.length; i++) {
			result = result + arr[i][0] + " " + arr[i][1] + "\n";
		}
		return result;
	}
}
