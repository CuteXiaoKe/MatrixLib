package matrixLib;

/**
 * Represents a complex number, which has a real and imaginary part
 * Enables using matrices over the field C instead of just R
 * @author Bryan Cuccioli
 */
public class ComplexNumber {

	private double re;
	private double im;
	public final double pi = 3.1415926536; // used in arg()

	/**
	 * Creates the complex number z= (0, 0)
	 */
	public ComplexNumber() {
		this.re = 0;
		this.im = 0;
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
	 * Scales a complex number by a factor
	 * @param f the scaling factor
	 * @return the complex number scaled by f
	 */
	public ComplexNumber multiply(float f) {
		
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
	 * Returns the modulus |z| of this complex number
	 * @return The absolute value of this complex number
	 */
	public double abs() {
		
		// if the number is real or imaginary, |z| is simple:
		if (re == 0.0) {
			return Math.abs(im);
		}
		else if (im == 0.0) {
			return Math.abs(re);
		}
		// otherwise it's the more complicated modulus formula:
		else {
			return Math.sqrt(re*re + im*im);
		}
	}

	/**
	 * Returns the argument of this complex number
	 * according to the definition at
	 * http://upload.wikimedia.org/math/c/e/6/ce60fee0f163f054266695db496b8c34.png
	 * @return the argument of this complex number
	 */
	public double arg() throws ArithmeticException {
		
		// arg is undefined for the complex number (0,0)
		if (re == 0.0 && im == 0.0) {
			throw new ArithmeticException();
		}
		
		if (re == 0.0) {
			if (im > 0) {
				return pi/2;
			}
			else {
				return -pi/2;
			}
		}
		
		else if (re > 0) {
			return Math.atan(im/re);
		}
		else {
			if (im < 0) {
				return Math.atan(im/re) - pi;
			}
			else {
				return Math.atan(im/re) + pi;
			}
		}
	}
	
	/**
	 * Returns the a + bi form of this complex number
	 * @return The string representation of this number
	 */
	public String toString() {

		return re + " + " + im + "i";
	}
	
	/**
	 * Tells whether two complex numbers are equal
	 * @param z The complex number to compare with this
	 * @return Whether the two complex numbers are equal
	 */
	public boolean equals(ComplexNumber z) {
		
		return this.re == z.Re() && this.im == z.Im();
	}
}
