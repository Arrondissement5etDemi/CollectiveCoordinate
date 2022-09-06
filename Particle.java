package toShare;
import java.util.*;

public class Particle {
        //data
	private double x,y;

	//constructor
        public Particle(double xx, double yy) {
		x = xx;
                y = yy;
        }

        //gets the coordinates
        public double getx() {
		return x;
	}
    
	public double gety() {
		return y;
	}

	//modifiers
	public void setX(double xx) {
		x = xx;
	}

	public void setY(double yy) {
		y = yy;
	}


	//distance to another city
	public double distanceto(Particle another) {
		double anotherX = another.getx();
                double anotherY = another.gety();
                return Math.sqrt(Math.pow((x - anotherX),2) + Math.pow((y - anotherY),2));
	}

	/**the minimum distance between two particles in PBC
 *      @param another Particle, another particle 
 *      @param sideLength double, the side length of the unit cell
 *      @return the minimum distance to that particle among all its replicates */
	public double distanceto(Particle another, double sideLength) {
                double thatX = another.getx();
                double thatY = another.gety();
                double miniDx = minIn3(thatX,thatX+sideLength,thatX-sideLength,x);
                double miniDy = minIn3(thatY,thatY+sideLength,thatY-sideLength,y);
                return Math.sqrt(miniDx*miniDx + miniDy*miniDy);
        }

	/**parallelogram box
	@param another Particle, another particle
 *      @param a double, the base length of the unit cell
 *      @param ma double, the slope of the unit cell
 *      @param b double, the height of the unit cell
 *      @return the minimum distance to that particle among all its replicates */
	public double distanceto(double thisa, double thisma, double thisb, Particle another, double a, double ma, double b) {
		double thisRealX = thisa*x + thisma*thisa*y;
                double thisRealY = thisb*y;	
		double thatX = another.getx();
                double thatY = another.gety();
		double thatRealX = a*thatX + ma*a*thatY;
		double thatRealY = b*thatY;
		Particle realParticle = new Particle(thisRealX, thisRealY);
		double result = realParticle.distanceto(new Particle(thatRealX, thatRealY));
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i != 0 || j != 0) {
					double thatImageX = thatRealX + i*a + j*ma*a;
					double thatImageY = thatRealY + j*b;
					result = Math.min(result, realParticle.distanceto(new Particle(thatImageX, thatImageY)));	
				}
			}
		}
		return result;	
	}

	/**computes if the images of the particle overlaps with a line
 * 	@param x1
 * 	@param x2
 * 	@param y1
 * 	@param y2
 * 	@param r double, radius of particle
 * 	@param sideLength
 * 	@return boolean*/
	public boolean overlapWithLine(double x1, double x2, double y1, double y2, double r, double sideLength) {
		double a = (y2 - y1)/(x2 - x1);
		double b = -1.0;
		double c = -x1*a + y1;
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				double x0 = x + i*sideLength;
				double y0 = y + j*sideLength;
				double dist = Math.abs(a*x0 + b*y0 + c)/Math.sqrt(a*a + b*b);
				double xOnTheFoot = (b*(b*x0 - a*y0) - a*c)/(a*a + b*b);
				if (dist < r && xOnTheFoot > Math.min(x1, x2) && xOnTheFoot < Math.max(x1, x2)) {
					return true;
				}
			}
		}
		return false;
	}

	public Vector3D toVector3D() {
		return new Vector3D(x, y, 0);
	}

	public static double minIn3(double a, double b, double c, double center) {
		return Math.min(Math.abs(a-center),Math.min(Math.abs(b-center),Math.abs(c-center)));
	}


	public boolean equals(Particle another) {
		return (x == another.getx() && y == another.gety());
	}

	public String toString() {
		return Double.toString(x) + " " + Double.toString(y);
	}
}
