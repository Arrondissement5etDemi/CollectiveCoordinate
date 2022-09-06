package toShare;
import java.util.*;
import java.util.Random.*;
import java.util.concurrent.ThreadLocalRandom;

//this is a d*d square box with n Particles in the periodic boundary condition.
public class Box {
	int n; //the number of particles
	double d; //the dimension of the box
	Particle[] partiArr;
	double cutoff; //the cutoff distance for calculating the energy at a particle
	static ThreadLocalRandom r = ThreadLocalRandom.current();

	/** constructs a box with randomized particle distribution 
 *   	@param nGiven int, the given numer of particles
 *    	@param dGiven double, the given size of the box
 *  	@return Box, the constructed box */
	public Box(int nGiven, double dGiven) {
		n = nGiven;
		d = dGiven;
		partiArr = new Particle[n];
		for (int i = 0; i < n; i++) {
			/**generate a candidate particle p*/
			double x = getRandomNumberInRange(0, d);
			double y = getRandomNumberInRange(0, d);
			Particle p = new Particle(x, y);
			partiArr[i] = p;	
		}
	}

	/** constructs a box with randomized particle distribution, but we impose a minimum neighbor distance 
 * 	@param nGiven int, the given numer of particles
 * 	@param dGiven double, the given size of the box
 * 	@param co     double, the given cutoff distance
 * 	@return Box, the constructed box */
	public Box(int nGiven, double dGiven, double multi) {
		n = nGiven;
		d = dGiven;
		cutoff = multi;
		partiArr = new Particle[n];
		double[] xArr = new double[n];
		double[] yArr = new double[n];
		for (int i = 0; i < n; i++) {
			double x = getRandomNumberInRange(0, d);
			xArr[i] = x;
		}
		for (int i = 0; i < n; i++) {
                        double y = getRandomNumberInRange(0, d);
                        yArr[i] = y;
                }
		shuffleArray(yArr);
		for (int i = 0; i < n; i++) {
			/**generate a candidate particle p*/
			Particle p = new Particle(xArr[i], yArr[i]);
			partiArr[i] = p;	
		}
	}

	/**generates a Box of 2D square lattice or 2D RSA
 *      @param nGiven int, the given numer of particles
 *	@param dGiven double, the given size of the box
 *	@param co     double, the given diameter for RSA particles
 *	@param sq     String, must be "sq"
 *	@return Box, a Box with n particles for the 2D lattice */
        public Box(int nGiven, double dGiven, double co, String sq) {
		n = nGiven;
		d = dGiven;
		partiArr = new Particle[n];
		if (sq == "sq") {
                	int partiPerSide = (int)Math.ceil(Math.sqrt((double)nGiven)); 
                	double neighborDist = dGiven/(double)partiPerSide;
                	partiArr = new Particle[n];
                	int index = 0;
               		for (int i = 0; i < partiPerSide && index < n; i++) {
                        for (int j = 0; j < partiPerSide && index < n; j++) {
                               	partiArr[index] = new Particle(i*neighborDist,j*neighborDist);
                               	index++;
                       	}
               		}
		}
		else if (sq=="rsa") {
			int i = 0;
			while (i < n) {
				double x = getRandomNumberInRange(0,d);
                        	double y = getRandomNumberInRange(0,d);
                        	Particle p = new Particle(x,y);
				boolean clash = false;
				for (int j = 0; j < i; j++) {
					Particle q = partiArr[j];
					double distPQ = p.distanceto(q, d);
					if (distPQ <= co) {
						clash = true;
						break;
					}
				}
				if (!clash) {
					partiArr[i] = p;
					i++;
				}
			}
		}
		else {		
                        throw new IllegalArgumentException("check parameter sq in the constructor");
                }
        }

	/** a box whose particle positions are predetermined
 *      @param nGiven int, the given numer of particles
 *      @param dGiven double, the given size of the box
 *      @return Box, the constructed box */
	public Box(int nGiven, double dGiven, double[][] coordinates) {
		n = nGiven;
                d = dGiven;
		partiArr = new Particle[n];
		for (int i = 0; i < n; i++) {
			double x = coordinates[i][0];
			double y = coordinates[i][1];
			Particle p = new Particle(x, y);
			partiArr[i] = p;
		}
	}

	/**clone a box*/
	public Box(Box b) {
		n = b.getN();
		d = b.getD();
		cutoff = b.getCo();
		partiArr = b.toArray();
	}

	/**superimpose two boxes, assuming that they have the same d*/
	public static Box superimpose(Box b1, Box b2) {
		int n1 = b1.getN();
		int n2 = b2.getN();
		int newN = n1 + n2;
		double[][] coord = new double[newN][2];
		for (int i = 0; i < n1; i++) {
			Particle p = b1.particleAt(i);
			coord[i][0] = p.getx();
			coord[i][1] = p.gety();
		}
		for (int i = n1; i < newN; i++) {
			Particle p = b2.particleAt(i - n1);
			coord[i][0] = p.getx();
                        coord[i][1] = p.gety();
		}
		return new Box(newN, b1.getD(), coord);
	}

	/**quadraple a box*/
	public Box quadraple() {
		int nn = n*4;
		double dd = d*2.0;
		double[][] coord = new double[nn][2];
		double[][] coordOriginal = this.getCoord();
		for (int i = 0; i < nn; i++) {
			switch (i/n) {
				case 0 : 
					coord[i][0] = coordOriginal[i%n][0];
					coord[i][1] = coordOriginal[i%n][1];
					break;
				case 1 :
					coord[i][0] = coordOriginal[i%n][0] + d;
					coord[i][1] = coordOriginal[i%n][1];
					break;
				case 2 :
					coord[i][0] = coordOriginal[i%n][0];
					coord[i][1] = coordOriginal[i%n][1] + d;
					break;
				case 3 :
					coord[i][0] = coordOriginal[i%n][0] + d;
					coord[i][1] = coordOriginal[i%n][1] + d;
					break;
			}
		}
		return new Box(nn, dd, coord);
	}

	/**resets the box to a square lattice
 * 	@return void*/
	public void resetToSq() {
		Box sqTemp = new Box(n, d, cutoff, "sq");
		Particle[] thatArr = sqTemp.toArray();
		for (int i = 0; i < n; i++) {
			partiArr[i] = thatArr[i];
		}		
	}

	/**reset the box to an ideal gas
 * 	@return void*/
	public void reset() {
                Box sqTemp = new Box(n, d, cutoff);
                Particle[] thatArr = sqTemp.toArray();
                for (int i = 0; i < n; i++) {
                        partiArr[i] = thatArr[i];
                }
        }

	/**reset partiArr
 * 	@return void*/
	public void setArr(Particle[] arr) {
		if (arr.length == n) {
			partiArr = arr;
		}
		else {
			System.out.println("Input partiArr must have the same number of particles as n");
		}
	}

	/**modify the particle at index i of partiArr 
 * 	@param i int
 * 	@param p Particle
 * 	@return void*/
	public void modifyParticle(int i, Particle p) {
		partiArr[i] = p;
	}

	/**modify the x coord of particle at index i of partiArr
 * 	@param i int
 * 	@param newX double
 * 	@return void*/
	public void modifyParticleX(int i, double newX) {
                partiArr[i].setX(newX);
        }

	
        /**modify the y coord of particle at index i of partiArr
 *      @param i int
 *      @param newY double
 *      @return void */
	public void modifyParticleY(int i, double newY) {
                partiArr[i].setY(newY);
        }

	public void setCoord(double[][] coord) {
		for (int i = 0; i < n; i++) {
			double x = coord[i][0];
			double y = coord[i][1];
			Particle p = new Particle(x, y);
			partiArr[i] = p;
		}
	}

	/** gets d, the dimension of the box
 * 	@return double */
	public double getD() {
		return d;
	}

	/** gets n, the number of particles
 *      @return int */
        public int getN() {
                return n;
        }

	/** gets co, the cutoff
 *      @return int */
        public double getCo() {
                return cutoff;
        }

	/**get the ith particle
 * 	@param i int
 * 	@return Particle */
	public Particle particleAt(int i) {
		return partiArr[i];
	}

        /**equilibrates a box under this potential, a neighborlist has been implemented in case the range of the potential is much smaller than the box size
 *      @param pot PairPot, a pair potential
 *      @param t double, the scaled temperature
 *      @param sweeps int, the number of sweeps
 *      @param displace double, displacement in each MC step
 *      @param hardcore double, imposed hard core
 *      @return void */	
	public void equilibrate(PairPot pot, double t, int sweeps, double displace, double hardcore) {
		double criterion = pot.getRange() + displace * 30;
		Particle[][] nl = new Particle[n][n];/**neighbor list*/
                int[] indList = new int[n];/**list of the number of neighbors for each particle*/
                for (int i = 0; i < sweeps; i++) {
			/**update neighbor list every 5 sweeps*/
			if (i%5 == 0) {
				nl = new Particle[n][n];
				indList = new int[n];
				for (int j = 0; j < n; j++) {
					Particle pj = partiArr[j];
					for (int k = 0; k < j; k++) {	
						Particle pk = partiArr[k];
						double distJK = minDist(pj, pk);
						if (distJK <= criterion) {
                                                        nl[j][indList[j]] = pk;
                                                        indList[j]++;
							nl[k][indList[k]] = pj;
							indList[k]++;
                                                }
					}
				}
			}
			/**compute energies according to neighborlist*/
                        for (int j = 0; j < n; j++) {
				Particle pj = partiArr[j];
				double oldx = pj.getx();
				double oldy = pj.gety();
				double oldE = partiE(j, nl, pot, hardcore);
                                move(j, displace);
                                double newE = partiE(j, nl, pot, hardcore);
                                if (!Boltzmann.accept(newE, oldE, t)) {
					pj.setX(oldx);
					pj.setY(oldy);
                                }
                        }
                }
        }

	/**equilibrates a box under this potential in the microcanonical ensemble
 *      @param pot PairPot, a pair potential
 *      @param e double, the energy
 *      @param sweeps int, the number of sweeps
 *      @param displace double, displacement in each MC step
 *      @param hardcore double, imposed hard core
 *      @return double, the energy corresponding to the microcanonical energy */
        public double equilibrateMicro(ArrayPot pot, double e, int sweeps, double displace, double hardcore) {
                double criterion = pot.getRange() + displace * 30;
                Particle[][] nl = new Particle[n][n];/**neighbor list*/
                int[] indList = new int[n];/**list of the number of neighbors for each particle*/
		double boxE = this.getEnergy(pot, hardcore);
		double demonEnergy = e - boxE;
		double demonAverage = demonEnergy;
                for (int i = 0; i < sweeps; i++) {
			if (i%10 == 0) {
                                System.out.println(demonEnergy);
                        }
                        /**update neighbor list every 5 sweeps*/
                        if (i%5 == 0) {
                                nl = new Particle[n][n];
                                indList = new int[n];
                                for (int j = 0; j < n; j++) {
                                        Particle pj = partiArr[j];
                                        for (int k = 0; k < j; k++) {
                                                Particle pk = partiArr[k];
                                                double distJK = minDist(pj, pk);
                                                if (distJK <= criterion) {
                                                        nl[j][indList[j]] = pk;
                                                        indList[j]++;
                                                        nl[k][indList[k]] = pj;
                                                        indList[k]++;
                                                }
                                        }
                                }
                        }
                        /**compute energies according to neighborlist*/
                        for (int j = 0; j < n; j++) {
                                Particle pj = partiArr[j];
                                double oldx = pj.getx();
                                double oldy = pj.gety();
                                double oldE = partiE(j, nl, pot, hardcore);
                                move(j, displace);
                                double newE = partiE(j, nl, pot, hardcore);
				double dE = newE - oldE;
				if (dE < demonEnergy || dE < 0) {//accept and give the excess energy to the demon
					demonEnergy = demonEnergy - dE;
				}
				else {//reject
					pj.setX(oldx);
                                        pj.setY(oldy);	
				}
                        }
			demonAverage += demonEnergy;
                }
		double t = demonAverage/sweeps;
		return t;
        }

	/**modifies a box such that its energy is equal to a target value
 * 	@param pot PairPot, the pair potential
 * 	@param hardcore double, imposed hard core
 * 	@param e double, the target energy
 * 	@param tol double, tolerance of the energy difference from the target
 * 	@return void */
	public void setEnergy(PairPot pot, double hardcore, double eTarg, double tol) {
		double totalE = this.getEnergy(pot, hardcore);
		double de = Math.abs(totalE - eTarg);//difference of the actual energy and the target energy
		while (de > tol) {
			System.out.println(de);
			for (int j = 0; j < n; j++) {
				Particle pj = partiArr[j];
                                double oldx = pj.getx();
                                double oldy = pj.gety();
                                double oldE = partiE(j, pot, hardcore);
                                move(j, 0.01);
                                double newE = partiE(j, pot, hardcore);
				double newTotalE = totalE + (newE - oldE)/n;
				double newDe = Math.abs(newTotalE - eTarg);
				if (newDe < de) {//accept the move and update totalE and de
					totalE = newTotalE;
					de = newDe;
				}
				else {//reject the move
					pj.setX(oldx);
                                        pj.setY(oldy);
				}
			}	
		}
		System.out.println("energy set to E = " + this.getEnergy(pot, hardcore));
	}

	/**computes the sum energy of all particles in the box
 * 	@param pot PairPot, the pair potential
 * 	@param hardcore double, imposed hard core
 * 	@return double, the ensemble energy */
	public double getEnergy(PairPot pot, double hardcore) {
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum = sum + partiE(i, pot, hardcore);
		}
		return sum/(2.0*n);
	}

	/**gets the energy felt by a single particle in the box
 * 	@param ind integer, the index of the particle in partiArr 
 * 	@param pot PairPot, a pair potential
 * 	@param hardcore double, imposed hard core
 * 	@return double, the energy felt by the particle
 * 	@assume Box side length > 2 * cutcoff */
	public double partiE(int ind, PairPot pot, double hardcore) {
		/**get the particle in question*/
		Particle thisParti = partiArr[ind];	
		/**For every particle in the box, compute their countribution to the particle*/
		double result = 0;
		for (int i = 0; i < n ; i++) {
			if (i != ind) {
				Particle thatParti = partiArr[i];
				double distToDupli = minDist(thisParti, thatParti);
				if (distToDupli < hardcore) {
					return Double.POSITIVE_INFINITY;
				}
				/**get the pair energy*/
				double pairEnergy = pot.pairE(distToDupli);
				result = result + pairEnergy;
			}
		}
		return result;
	}

	public double partiE(int ind, PairPot pot) {//without hardcore
                /**get the particle in question*/
                Particle thisParti = partiArr[ind];
                /**For every particle in the box, compute their countribution to the particle*/
                double result = 0;
                for (int i = 0; i < ind ; i++) {
                       double distToDupli = minDist(thisParti, partiArr[i]);
                       result = result + pot.pairE(distToDupli);;
                }
		for (int i = ind + 1; i < n; i++) {
                       double distToDupli = minDist(thisParti, partiArr[i]);
                       result = result + pot.pairE(distToDupli);;
                }
                return result;
        }

        /**gets the energy felt by a single particle in the box, given its neighborList
 *      @param j int, the index of the particle in partiArr
 *	@param nl Particle[n][n], the neighbor list of the particle
 *	@param pot PairPot, the pair potential
 *	@param hardcore double, imposed hard core
 *      @return double, the energy felt by the particle
 *      @assume Box side length > 2 * cutcoff */
        public double partiE(int j, Particle[][] nl, PairPot pot, double hardcore) {
                Particle pj = partiArr[j];
                /**For every particle in the neighbor list, compute their countribution to the particle*/
                double result = 0;
		Particle pk = nl[j][0];
                for (int k = 0; pk != null; k++) {
			double distJK = minDist(pj, pk);
			if (distJK < hardcore) {
                                return Double.POSITIVE_INFINITY;
                        }
                        Double pairEnergy = pot.pairE(distJK);
			result = result + pairEnergy;
			pk = nl[j][k + 1];
                }
                return result;
        }


	/**gets the neighbor list of a particle
 * 	@param ind  int, the index of the particle in partiArr
 * 	@param criterion double, how close should 2 particles be to be considered neighbors?
 * 	@return List, a list of the neighbors of the particle indexed by ind */
	public List<Particle> neighborList(int ind, double criterion) {
		List<Particle> list = new ArrayList<Particle>();
                Particle thisParti = partiArr[ind];
		/**For every particle in the box, see if they are near the particle*/
		for (int i = 0; i < n ; i++) {
                        if (i != ind) {
                                Particle thatParti = partiArr[i];
                                double distToDupli = minDist(thisParti, thatParti);
                                /**are you a neighbor?*/
                                if (distToDupli <= criterion) {
                                        list.add(thatParti);
                                }
                        }
                }
		return list;
	}

	 /**the minimum distance between two particles in PBC
 *      @param p Particle
 *      @param q Particle, another particle
 *      @return Particle, the minimum image wrt p of q among all its replicates in 2D */
        public MiniResult minimage(Particle p, Particle q) {
                double xp = p.getx();
                double yp = p.gety();
                double xq = q.getx();
                double yq = q.gety();
                int[] positions = new int[2];
                double[] partialDist = new double[2];
                double inf = Double.POSITIVE_INFINITY;
                for (int i = 0; i < 2; i++) {
                        partialDist[i] = inf;
                }
                for (int i = -1; i < 2; i++) {
                        double dx = Math.abs(xp - (xq + i*d));
                        if (dx < partialDist[0]) {
                                partialDist[0] = dx;
                                positions[0] = i;
                        }
                }
                for (int i = -1; i < 2; i++) {
                        double dy = Math.abs(yp - (yq + i*d));
                        if (dy < partialDist[1]) {
                                partialDist[1] = dy;
                                positions[1] = i;
                        }
                }
                Particle imageq = new Particle(xq + positions[0]*d, yq + positions[1]*d);
                double dist = Math.sqrt(partialDist[0]*partialDist[0] + partialDist[1]*partialDist[1]);
                MiniResult result = new MiniResult(imageq, dist);
                return result;
        }

	/**the minimum distance between two particles in PBC
 *     @param thisP Particle
 *     @param another Particle, another particle 
 *     @return the minimum distance to that particle among all its replicates in 2D */
	public double minDist(Particle thisP, Particle another) {
                double thatX = another.getx();
                double thatY = another.gety();
                double miniDx = minIn3(thatX,thatX + d,thatX - d, thisP.getx());
                double miniDy = minIn3(thatY,thatY + d,thatY - d, thisP.gety());
                return Math.sqrt(miniDx*miniDx + miniDy*miniDy);
        }

	public static double minIn3(double a, double b, double c, double center) {
		return Math.min(Math.abs(a-center),Math.min(Math.abs(b-center),Math.abs(c-center)));
	}

	/**moves a random particle in a random direction by a random distance between 0 and maxDist
 * 	@param ind int, the index of the particle moved
 *	@param maxDist double, the maximum distance to move for a movement
 * 	@return void, moves the particle with index ind */
	public void move(int ind, double maxDist) {
		Particle p =  partiArr[ind];
		double x = p.getx();
                double y = p.gety();
		double dx = getRandomNumberInRange(-maxDist,maxDist);
		double dy = getRandomNumberInRange(-maxDist,maxDist);
		partiArr[ind].setX(positiveModulo(x + dx, d));
                partiArr[ind].setY(positiveModulo(y + dy, d));
	}
	
	/**shuffle the particles in a box
 * 	@return void**/
	public void shuffle() {
		for (int i = 0; i < n; i++) {
                        int randomIndexToSwap = r.nextInt(n);
                        Particle temp = partiArr[randomIndexToSwap];
                        partiArr[randomIndexToSwap] = partiArr[i];
                        partiArr[i] = temp;
                }
	}

	/** returns the array of particles
 * 	@return a newly constructed array of particles that is a deep copy of partiArr */
	public Particle[] toArray() {
		Particle[] partiArrCopy = new Particle[n];
		for (int i = 0; i < n; i++) {
			double x = partiArr[i].getx();
			double y = partiArr[i].gety();
			Particle pCopy = new Particle(x,y);
			partiArrCopy[i] = pCopy;
		}
		return partiArrCopy;
	}

	/**returns the particle coordinates*/
	public double[][] getCoord() {
		double[][] coord = new double[n][2];
		for (int i = 0; i < n; i++) {
                        coord[i][0] = partiArr[i].getx();
			coord[i][1] = partiArr[i].gety();
                }
		return coord;
	}

	/** returns the contant of partiArray as a string
 * 	@return the content of partiArr */
	public String toString() {
		String result = "";
		for (int i = 0; i < n; i++) {
			result = result + partiArr[i].toString() + "\n";
		}
		return result;
	}

	//auxillary function
	public static double getRandomNumberInRange(double min, double max) {
		if (min >= max) {
			throw new IllegalArgumentException("max must be greater than min");
	        }	
       		return r.nextDouble()*(max - min) + min;
	}
	
		
	private static double positiveModulo(double x, double d) {
		double result = x%d;
		if (result < 0) {
			result = result + d;
		}
		return result;
	}

	private static double expoRandomNumber(double lambda) {
		//get a uniform random variable 
		double u = r.nextDouble();
		//create the exponentially distributed random number
		double result = Math.log(1 - u)/(-lambda);
		return result;
	}

	public static void shuffleArray(double[] a) {
        	int n = a.length;
        	for (int i = 0; i < n; i++) {
            		int change = i + r.nextInt(n - i);
            		swap(a, i, change);
        	}
    	}

    	private static void swap(double[] a, int i, int change) {
        	double helper = a[i];
        	a[i] = a[change];
        	a[change] = helper;
    	}	
}

class MiniResult {
        Particle image;
        double dist;
        public MiniResult(Particle p, double d) {
                image = p;
                dist = d;
        }
        public void setImage(Particle p) {
                image = p;
        }
        public void setDist(double d) {
                dist = d;
        }
        public Particle getImage() {
                return image;
        }
        public double getDist() {
                return dist;
        }
}

