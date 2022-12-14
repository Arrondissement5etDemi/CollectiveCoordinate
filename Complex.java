package toShare;
import java.util.*;

public class Complex {
	private double re,im,norm,angle;
	//constructs a complex number from its real and imaginary parts
	public Complex(double real,double imaginary) {
		re = real;
		im = imaginary;
		setNormAngle();
	}

	//converts a real number to a complex number
	public Complex(double real) {
		re = real;
		im = 0.0;
		setNormAngle();
	}

	//constructs a complex number from its norm and angle
	public Complex(int dummy, double n, double a) {
		norm = n;
		angle = a%(2*Math.PI);
		re = n * Math.cos(angle);
		im = n * Math.sin(angle);
	}

	private void setNormAngle() {
		norm = Math.sqrt(Math.pow(re,2) + Math.pow(im,2));
		angle = Math.atan2(im,re);
	}

	//accessors
	public double getRe() {
		return re;
	}

	public double getIm() {
		return im;
	}

	public double getNorm() {
		norm = Math.sqrt(Math.pow(re,2) + Math.pow(im,2));
		return norm;
	}

	public double normSq() {
                return re*re + im*im;
        }

	public double getAngle() {
		return angle;
	}

	/**The constant i=sqrt(-1)*/
	public static final Complex i() {
		Complex result = new Complex(0.0,1.0);
		return result;
	}

	/**The constant 1*/
	public static final Complex one() {
		Complex result = new Complex(1.0,0.0);
		return result;
	}

	//arithmetic
	public Complex conjugate() {
		Complex result = new Complex(re,-im);
		return result;
	}

	public void add(Complex another) {
		re += another.getRe();
		im += another.getIm();
	}

	public Complex nonDestructiveAdd(Complex another) { 
		return new Complex(re + another.getRe(), im + another.getIm());
	}

	public Complex minus(Complex another) {
		double resultRe = re - another.getRe();
                double resultIm = im - another.getIm();
                Complex result = new Complex(resultRe, resultIm);
                return result;
	}

	public Complex times(Complex another) {
		double a = re;
		double b = im;
		double c = another.getRe();
		double d = another.getIm();
		double resultRe = a * c - b * d;
		double resultIm = a * d + b * c;
		Complex result = new Complex(resultRe, resultIm);
		return result;
	}

	public Complex inverse() {
		if (re == 0.0 && im == 0.0) {
			throw new ArithmeticException("You can't devide by 0.");
		}
		else {
			double a = re;
			double b = im;
			double resultRe = a/(a*a+b*b);
			double resultIm = -b/(a*a+b*b);
			Complex result = new Complex(resultRe, resultIm);
			return result;
		}
	}

	public Complex divide(Complex another) {
		return this.times(another.inverse());
	}

	//computes exp(this)
	public Complex expThis() {
		double n = Math.exp(re);//norm
		double a = im;//angle
		Complex result = new Complex(1,n,a);
		return result;
	}

	//computes q=exp(i*pi*this) 
	public Complex expIpithis() {
		double n = Math.exp(-1*Math.PI*im);
		double a = re*Math.PI;
		Complex result = new Complex(1,n,a);
		return result;
	}

	//computes q=exp(i*this)
        public Complex expIthis() {
                double n = Math.exp(-1*im);
                double a = re;
                Complex result = new Complex(1, n, a);
                return result;
        }

	public static Complex expI(double arg) {
		return new Complex(Math.cos(arg), Math.sin(arg));
	}

	public boolean equals(Complex another) {
		return (re==another.getRe() && im==another.getIm());
	}
	
	public String toString() {
		if (im >= 0) {
			return Double.toString(re)+" + "+Double.toString(im)+" i";
		}
		else {
			return Double.toString(re)+" "+Double.toString(im)+" i";
		}
	}
}
