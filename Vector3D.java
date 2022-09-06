package toShare;
import java.util.*;

/**A vector of at most 3D*/
public class Vector3D {
        //data
	private double x,y,z;//coordinates of the spheres

	//constructor
        public Vector3D(double xx, double yy, double zz) {
		x = xx;
                y = yy;
		z = zz;
        }

        //gets the coordinates
        public double getx() {
		return x;
	}
    
	public double gety() {
		return y;
	}

	public double getz() {
		return z;
	}

	//modifiers
	public void setX(double xx) {
		x = xx;
	}

	public void setY(double yy) {
		y = yy;
	}

	public void setZ(double zz) {
		z = zz;
	}

	/**decides if two vectors are the same
 *      @param another Particle, another particle */
        @Override
        public boolean equals(Object another) {
                if (!(another instanceof Vector3D)) {
                        return false;
                }
                Vector3D a = (Vector3D) another;
                return (Math.abs(x - a.getx())<1e-6 && Math.abs(y - a.gety())<1e-6 && Math.abs(z - a.getz())<1e-6);
        }

        public int hashCode() {
                return 0;
        }

	/**computes vector length
 * 	@return double, the vector length */
	public double length() {
		return Math.sqrt(x*x+y*y+z*z);
	}

	/**computes vector addition
 * 	@param another Vector3D, another 3D vector
 * 	@return Vector3D, the sum vector */
	public Vector3D add(Vector3D another) {
		double sumX = x + another.getx();
		double sumY = y + another.gety();
		double sumZ = z + another.getz();
		Vector3D result	= new Vector3D(sumX, sumY, sumZ);
		return result;
	}

	/**computes scalar multiple of the vector
 *	@param scalar double, the constant multiplying the vector
 *	@return Vector3D, the result of scalar multiplication */
	public Vector3D scalarMul(double scalar) {
		Vector3D result = new Vector3D(scalar*x,scalar*y,scalar*z);
		return result;
	}
	
	/**computes vector subtraction this - another
 * 	@param another Vector3D, another 3D vector
 * 	@return Vector3D, the difference vector */
	public Vector3D minus(Vector3D another) {
		return this.add(another.scalarMul(-1.0));
	}

	/**computes the inner product with another vector
 *	@param another Vector3D, another 3D vector
 *	@return double, the inner product */
	public double innerProd(Vector3D another) {
		double prodX = x * another.getx();
		double prodY = y * another.gety();
		double prodZ = z * another.getz();
		return prodX + prodY + prodZ;
	}

	/**decides if two vectors are the same
 * 	@param another Particle, another particle */
	public boolean equals(Vector3D another) {
		return (x == another.getx() && y == another.gety() && z == another.getz());
	}

	public String toString() {
		return Double.toString(x) + " " + Double.toString(y) + ' ' + Double.toString(z);
	}
}
