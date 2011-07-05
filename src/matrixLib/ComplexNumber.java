package matrixLib;

/**
 * Represents a complex number, which has a real and imaginary part
 * Enables using matrices over the field C instead of just R
 * @author Bryan Cuccioli
 */
public class ComplexNumber implements Field, Comparable<ComplexNumber> {

	private double re;
	private double im;
	private double epsilon;

	/**
	 * Creates the complex number z= (0, 0)
	 */
	public ComplexNumber() {
		this.re = 0;
		this.im = 0;
		epsilon = 1e-15;
	}
	
	/**
	 * Creates a complex number z = (re, im)
	 * @param re The real part of z
	 * @param im The imaginary part of z
	 */
	public ComplexNumber(double re, double im) {
		
		this.re = re;
		this.im = im;
	}

	/**
	 * Returns the real part of the complex number
	 * @return The real part of the complex number
	 */
	public double Re() {
		return this.re;
	}
	
	/**
	 * Returns the imaginary part of the complex number
	 * @return The imaginary part of the complex number
	 */
	public double Im() {
		return this.im;
	}
	
	/**
	 * Sets the epsilon to use in comparing complex numbers
	 * @param eps the epsilon to use
	 */
	public void setEpsilon(double eps) {
		this.epsilon = eps;
	}
	
	/**
	 * Gets the epsilon currently in use in comparing complex numbers
	 * @return the epsilon currently in use
	 */
	public double getEpsilon() {
		return this.epsilon;
	}
	
	/**
	 * Returns the complex number created by adding this to z
	 * @param z The complex number to add to this
	 * @return The complex number created by adding this with z
	 */
	public ComplexNumber add(ComplexNumber z) {
		
		return new ComplexNumber(z.Re() + re, z.Im() + im);
	}
	
	/**
	 * Returns the complex number created by subtracting z from this
	 * @param z The complex number to subtract from this
	 * @return The complex number created by subtracting z from this
	 */
	public ComplexNumber subtract(ComplexNumber z) {
		
		return new ComplexNumber(re - z.Re(), im - z.Im());
	}
	
	/**
	 * Returns the complex number created by multiplying this with z
	 * Using the standard definition of multiplication over C
	 * @param z The complex number to multiply with this
	 * @return The complex number created by multiplying this with z
	 */
	public ComplexNumber multiply(ComplexNumber z) {
		
		return new ComplexNumber(
				re*z.Re() - im*z.Im(),
				re*z.Im() + im*z.Re());
	}
	
	/**
	 * Divides this complex number by z (through multiplying by the inverse)
	 * @param z the complex number to divide this by
	 * @return the quotient of this and the argument
	 * @throws ArithmeticException trying to divide by 0
	 */
	public ComplexNumber divide(ComplexNumber z) throws ArithmeticException {
		
		if (z.isZero()) {
			 throw new ArithmeticException();
		}
		
		// this is more numerically stable than the concise formulation
		ComplexNumber quot = this.multiply(z.conjugate());
		double recmod = (z.Re()*z.Re() + z.Im()*z.Im());
		return new ComplexNumber(quot.Re()/recmod, quot.Im()/recmod);
	}

	/**
	 * Returns the reciprocal 1/z of this complex number
	 * @return the reciprocal of this complex number
	 * @throws ArithmeticException trying to take the reciprocal of zero
	 */
	public ComplexNumber reciprocal() throws ArithmeticException {
		
		if (this.isZero()) {
			 throw new ArithmeticException();
		}
		
		// this is more numerically stable than the concise formulation
		ComplexNumber quot = this.conjugate();
		double recmod = (re*re + im*im);
		return new ComplexNumber(quot.Re()/recmod, quot.Im()/recmod);
	}
	
	/**
	 * Scales a complex number by a factor
	 * @param f the scaling factor
	 * @return the complex number scaled by f
	 */
	public ComplexNumber multiply(double f) {
		
		return new ComplexNumber(re*f, im*f);
	}
	
	/**
	 * Returns the complex conjugate a-bi of this complex number
	 * @return The complex conjugate of this complex number
	 */
	public ComplexNumber conjugate() {
		
		return new ComplexNumber(re, -im);
	}
	
	/**
	 * Returns the modulus |z| of this complex number, sqrt(Re^2+Im^2)
	 * @return The absolute value of this complex number
	 */
	public double abs() {
		
		// if the number is real or imaginary, |z| is simple:
		if (re == 0) {
			return Math.abs(im);
		}
		else if (im == 0) {
			return Math.abs(re);
		}
		// otherwise it's the more complicated modulus formula:
		else {
			return Math.sqrt(re*re + im*im);
		}
	}

	/**
	 * Returns the argument of this complex number
	 * @return the argument of this complex number
	 * @throws ArtithmeticException trying to take the argument of 0
	 */
	public double arg() throws ArithmeticException {
		
		// arg is undefined for the complex number (0,0)
		if (re == 0.0 && im == 0.0) {
			 throw new ArithmeticException();
		}
		
		if (re == 0.0) {
			if (im > 0) {
				return Math.PI/2;
			}
			else {
				return -Math.PI/2;
			}
		}
		
		else if (re > 0) {
			return Math.atan(im/re);
		}
		else {
			if (im < 0) {
				return Math.atan(im/re) - Math.PI;
			}
			else {
				return Math.atan(im/re) + Math.PI;
			}
		}
	}

	/**
	 * Returns the principle square root of this complex number
	 * @return the principle square root of this complex number
	 */
	public ComplexNumber sqrt() {
		
		double a = this.Re(), b = this.Im();
		double p = 1.0/Math.sqrt(2) * Math.sqrt(Math.sqrt(a*a+b*b)+a);
		double q = 1.0/Math.sqrt(2) * Math.sqrt(Math.sqrt(a*a+b*b)-a);
		
		if (b < 0) {
			q = -q;
		}
		
		return new ComplexNumber(p, q);
	}
	
	/**
	 * Returns the a + bi form of this complex number
	 * @return The string representation of this number
	 */
	public String toString() {

		// return only a or bi if the number is real or imaginary
		if (re == 0) {
			if (im == 0) {
				return "0";
			}
			return Double.toString(im) + "i";
		}
		else if (im == 0) {
			return Double.toString(re);
		}
		// otherwise have to return a + bi
		else {
			return re + " + " + im + "i";
		}
	}
	
	/**
	 * Tells whether two complex numbers are equal
	 * @param z The complex number to compare with this
	 * @return Whether the two complex numbers are equal
	 */
	public boolean equals(ComplexNumber z) {
		
		//return this.re == z.Re() && this.im == z.Im();
		/*if (Math.abs(re-z.Re()) > epsilon) {
			System.out.println("---------------");
			System.out.println(re);
			System.out.println(z.Re());
		}*/
		
		return Math.abs(re-z.Re()) < epsilon
			&& Math.abs(im-z.Im()) < epsilon;
		//return re == z.Re() && im == z.Im();
	}
	
	/**
	 * Tells whether this is the 0 vector in C
	 * @return whether this is the complex number equal to 0
	 */
	public boolean isZero() {
		
		return Math.abs(this.re)<epsilon && Math.abs(this.im)<epsilon;
		//return this.re == 0 && this.im == 0;
	}

	/**
	 * Compares this to z by comparing the modulus of this with |z|
	 */
	public int compareTo(ComplexNumber z) {
		
		if (this.abs() > z.abs()) {
			return -1;
		}
		else if (this.abs() < z.abs()) {
			return 1;
		}
		else {
			return 0;
		}
	}
}
