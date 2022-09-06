package toShare;
import java.io.*;
import java.util.*;
import java.util.function.*;
import static java.lang.Math.*;

public class CollectiveCoord {
	static double bigK = 3;//the maximum constrained K
	static double bigB = 0.0;//the target S(k) for small k is ~bigB*k^\alphaForS
	static double alphaForS = 1;//the target S(k) for small k is ~bigB*k^\alphaForS
	static ToDoubleFunction<Double> v = k -> Math.pow(10.0/k - 1.0, 3);//weighing function that weighs small k heavier than large k
	static ArrayList<Vector3D> allVecs = new ArrayList<Vector3D>();//a list that stores all constrained vectors

	public static void main(String[] args) throws IOException {
		int n = 200;//number of particles per configuration
		double rho = 0.390625;//number density
		double d = Math.pow(n/rho, 1.0/2.0);//box side length
		double chi = bigK*bigK/(16.0*PI*rho);//compute the fraction of contrained degrees of freedom. \rho\chi=v_1(K)/(2*d*(2*Pi)^d).
		System.out.println("chi = " + chi);//print chi
		Box b = new Box(n, d);//the initial configuration
		Box[] boxArr = {b};
		double kMini = 2.0*PI/d;//magninude of the minimum k vector
		getAllK(kMini, bigK);//get all the constrained vectors
		int round = 300;//number of iterations
		evolve(b, 1e-7, round, "BFGS");//start optimization
		System.out.println(b);//print the configuration
		double[][] g2 = RadialStat.g2Equil(boxArr, 0.05);//compute g2
                System.out.println(Functional.funcToString(g2));//print g2
		double[][] sk = StrucFac.skEquil(boxArr);//compute S
		System.out.println(Functional.funcToString(sk));//print S	
	}

	/**generates the list of all vectors that have length smaller than bigK
*       @param rho double, number denstiy
*       @param bigK double
*       @return double[][], the table of angular averaged s(k) */
        public static void getAllKHorizontal(Vector3D k, double kMini, double bigK){
		if (k.length() <= bigK && !(k.gety() == 0 && k.getx() < 0)) {
			allVecs.add(k);
			for (int i = -1; i <= 1; i += 2) {
				Vector3D neighborVec = new Vector3D(k.getx() + i*kMini, k.gety(), 0);
				if (!allVecs.contains(neighborVec)) {
					getAllKHorizontal(neighborVec, kMini, bigK);
				}
			}
		}
		else {
			return;
		}
	}

	public static void getAllK(double kMini, double bigK){
                double ky = 0;
		while (ky < bigK) {
			getAllKHorizontal(new Vector3D(0, ky, 0), kMini, bigK);
                	ky += kMini;
		}
        }

	/**the fictitious energy*/
        public static double phi(Box b) {
                double[][] sk = s(b);
                double result = 0;
                for (int i = 0; i < sk.length; i++) {
                        double k = allVecs.get(i).length();
                        if (k != 0) {
				double term = v.applyAsDouble(k)*Math.pow(sk[i][2] - bigB*Math.pow(k, alphaForS), 2);//for targets like S(k)~B*k^\alpha
				//double term = v.applyAsDouble(k)*Math.pow(sk[i][2], 2);//for stealthy targets
                        	result += term;
			}
                }
                return result;
        }

	/**computes \Sum_i exp^(ikr_i)*/
	public static Complex cc(Box b, Vector3D k) {
		int n = b.getN();
		double re = 0;
		double im = 0;
		for (int i = 0; i < n; i++) {
			double arg = b.particleAt(i).toVector3D().innerProd(k);
			re += Functional.cos(arg);
			im -= Functional.sin(arg);
		}
		return new Complex(re, im);
	}

	/**computes \Sum_i exp^(ikr_i) for all k<K*/
	public static Complex[] cc(Box b) {
		Complex[] result = new Complex[allVecs.size()];
                for (int i = 0; i < allVecs.size(); i++) {
			result[i] = cc(b, allVecs.get(i));
		}
		return result;
        }

	public static double[][] s(Box b) {
		Complex[] ccResult = cc(b);
		double n = (double)b.getN();
		double[][] result = new double[allVecs.size()][3];
		for (int i = 0; i < allVecs.size(); i++) {
			result[i][0] = allVecs.get(i).getx();
			result[i][1] = allVecs.get(i).gety();
			result[i][2] = ccResult[i].normSq()/n;
		}
		return result;
	}

	/**computes the gradient of a potential wrt particle coordinates
 *      @param b Box,
 *      @param delta double, a small number used to compute the gradient
 *  @return double[][] the gradient*/
        public static double[][] gradient(Box b) {
                int n = b.getN();
                double[][] result = new double[n][2];
		Complex[] ccOriginal = cc(b);
                for (int i = 0; i < n; i++) {
                        Vector3D pI = b.particleAt(i).toVector3D();
                        for (int j = 0; j < allVecs.size(); j++) {
                                Vector3D k = allVecs.get(j);
                                if (k.length() != 0) {
                                        double targ = bigB*Math.pow(k.length(), alphaForS);
                                        double weight = v.applyAsDouble(k.length());
					Complex ccK = ccOriginal[j];
					double re = ccK.getRe();
					double im = ccK.getIm();
					double kr = pI.innerProd(k);
					double temp = re*Functional.sin(kr) + im*Functional.cos(kr);
					double s = (re*re + im*im)/n;
					temp *= (weight*(s - targ));
					result[i][0] += temp*k.getx();
					result[i][1] += temp*k.gety();
                                }
                        }
			result[i][0] *= -4.0/(double)n;
			result[i][1] *= -4.0/(double)n;
                }
                return result;
        }

	/**convert a 2D array into a n*dim rank-1 column vector
	 * @param b Box
	 * @param delta, double
	 * @param dim, dimension
	 * @return
	 */
	private static double[][] arrayToVec(double[][] arr) {
		int n = arr.length;
		int dim = arr[0].length;
		double[][] result = new double[n*dim][1];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < dim; j++) {
				result[dim*i + j][0] = arr[i][j];
			}
		}
		return result;
	}
	
	/**
	 * find the inherent structure with Gradient descent, SR1 or BFGS!
	 * @param b Box
	 * @param alpha double, the step size parameter
	 * @param round int, number of iterations
	 * @param method String, BFGS or other quasi-Newton
	 */
	public static void evolve(Box b, double alpha, int round, String method) {
		int n = b.getN(); //number of particles
		double d = b.getD();
		double[][] coord = b.getCoord(); //particle coords
		int dim = coord[0].length; //dimension
		double[][] h = identityMat(n*dim);//approximation of inverse of Hessian
		double[][] gradVecNew = arrayToVec(gradient(b));
		double[][] gradVec = new double[n*dim][1];
		double e = phi(b); 
		double enew = e - 1.0;
		double alphaK;
		double dxVecNorm;
		for (int i = 0; i < round; i++) {
			if (alpha < 1E-13) {
                                System.out.println("");
                                continue;
                        }
			alphaK = 0;
			//compute the gradient at the guess config
			for (int j = 0; j < gradVec.length; j++) {
				gradVec[j][0] = gradVecNew[j][0];
			}
			//compute the amount to move, dx, and move the particles (line search)
			double[][] direction = multiplyMatrices(h, gradVec);
			double[][] dxVec = new double[n*dim][1];
			int lineSearchCount = 0;
			double eOld = phi(b);
			double[][] coordOld = b.getCoord();
			while (lineSearchCount <= 100) {
				/**perform line search and move the particles in the computed direction*/
				direction = multiplyMatrices(h, gradVec);
				for (int j = 0; j < n*dim; j++) {
					dxVec[j][0] = -alpha*direction[j][0];
				}
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < dim; k++) {
						coord[j][k] = (coord[j][k] + dxVec[dim*j + k][0] + d)%d;
					}
				}
				b.setCoord(coord);
				enew = phi(b);
				//System.out.println(enew);
				/**compute the new gradient*/
				if (enew < eOld) {
					alphaK += alpha;
					eOld = enew;
					coordOld = b.getCoord();
					lineSearchCount++;
					if (lineSearchCount >= 3) {
                         		       alpha *= 2;
					       lineSearchCount = 1;
                        		}
				}
				else if (lineSearchCount > 0) {
					for (int j = 0; j < n; j++) {
                                               for (int k = 0; k < dim; k++) {
                                                        coord[j][k] = coordOld[j][k];
                                                }
                                        }
                                	b.setCoord(coord);
					alpha /= 2;
					break;
				}
				else {
					for (int j = 0; j < n; j++) {
                                                for (int k = 0; k < dim; k++) {
                                                        coord[j][k] = coordOld[j][k];
                                                }
                                        }
                                        b.setCoord(coord);
					if (alpha >= 1E-13) {
						alpha /= 2;
						h = identityMat(n*dim);
					}
					else {
						h = identityMat(n*dim);
						break;
					}
				}
			}
			for (int j = 0; j < n*dim; j++) {
                                dxVec[j][0] = -alphaK*direction[j][0];
                        }
			dxVecNorm = multiplyMatrices(transpose(dxVec), dxVec)[0][0];	
			gradVecNew = arrayToVec(gradient(b));
			//compute the change in gradient, yVec
			double[][] yVec = new double[n*dim][1];
			for (int j = 0; j < n*dim; j++) {
				yVec[j][0] = gradVecNew[j][0] - gradVec[j][0]; 
			}
			double gradNorm = multiplyMatrices(transpose(gradVecNew), gradVecNew)[0][0];
			System.out.println(i + " " + eOld + " " + gradNorm + " " + alpha);
			//if dxVec == 0, skip updating H
			double yVecNorm = multiplyMatrices(transpose(yVec), yVec)[0][0];
			if (yVecNorm == 0) {
				h = identityMat(n*dim);
				continue;
			}
			//update h with Sherman-Morrison, awww too much work
			if (method == "GD") {
				continue;
			}
			if (method == "BFGS") {
				double dxTTimesY = multiplyMatrices(transpose(dxVec), yVec)[0][0];
				double yTTimesHTimesY = multiplyMatrices(transpose(yVec), multiplyMatrices(h, yVec))[0][0]; 
				double[][] dxTimesDxT = multiplyMatrices(dxVec, transpose(dxVec));
				double[][] secondTerm = dxTimesDxT;
				for (int j = 0; j < n*dim; j++) {
					for (int k = 0; k < n*dim; k++) {
						secondTerm[j][k] *= (dxTTimesY + yTTimesHTimesY)/(dxTTimesY*dxTTimesY);
					}
				}
				double[][] hTimesYTimesDxT = multiplyMatrices(multiplyMatrices(h, yVec), transpose(dxVec));
				double[][] dxTimesYTTimesH = transpose(hTimesYTimesDxT);
				double[][] thirdTerm = new double[n*dim][n*dim];
				for (int j = 0; j < n*dim; j++) {
					for (int k = 0; k < n*dim; k++) {
						thirdTerm[j][k] = (hTimesYTimesDxT[j][k] + dxTimesYTTimesH[j][k])/dxTTimesY;
					}
				}
				for (int j = 0; j < n*dim; j++) {
					for (int k = 0; k < n*dim; k++) {
						h[j][k] = h[j][k] + secondTerm[j][k] - thirdTerm[j][k];
					}
				}
			}
			if (method == "SR1") {
				double[][] hTimesY = multiplyMatrices(h, yVec);
				double[][] dxMinusHTimesY = new double[n*dim][1];
				for (int j = 0; j < n*dim; j++) {
					dxMinusHTimesY[j][0] = dxVec[j][0] - hTimesY[j][0];
                                }
				double[][] numerator = multiplyMatrices(dxMinusHTimesY, transpose(dxMinusHTimesY));
				double denominator = multiplyMatrices(transpose(dxMinusHTimesY), yVec)[0][0];
				double[][] secondTerm = new double[n*dim][n*dim];
				for (int j = 0; j < n*dim; j++) {
                                        for (int k = 0; k < n*dim; k++) {
                                                secondTerm[j][k] = numerator[j][k]/denominator;
                                        }
                                }
				for (int j = 0; j < n*dim; j++) {
                                        for (int k = 0; k < n*dim; k++) {
                                                h[j][k] = h[j][k] + secondTerm[j][k];
                                        }
                                }
			}
		}
	}
	
	
	/**
	 * multiplies 2 matrices
	 * @param firstMatrix
	 * @param secondMatrix
	 * @return the product
	 */
	private static double[][] multiplyMatrices(double[][] firstMatrix, double[][] secondMatrix) {
        int r1 = firstMatrix.length;
        int c1 = firstMatrix[0].length;
        int c2 = secondMatrix[0].length;
		double[][] product = new double[r1][c2];
        for(int i = 0; i < r1; i++) {
            for (int j = 0; j < c2; j++) {
                for (int k = 0; k < c1; k++) {
                    product[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                }
            }
        }
        return product;
    }

	/**
	 * trasposes a matrix
	 * @param mat
	 * @return
	 */
	private static double[][] transpose(double[][] mat) {
		int r = mat.length;
		int c = mat[0].length;
		double[][] result = new double[c][r];
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				result[i][j] = mat[j][i];
			}
		}
		return result;
	}
	
	/**
	 * A size*size identity matrix
	 * @param size
	 * @return
	 */
	private static double[][] identityMat(int size) {
		double[][] result = new double[size][size];
		for (int i = 0; i < size; i++) {
			result[i][i] = 1.0;
		}
		return result;
	}
}

