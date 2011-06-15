package matrixLib;

/**
 * Implementation of vectors (column matrices)
 * @author Bryan Cuccioli
 */

public class Vector extends Matrix {
	
	/**
	 * Constructs an empty vector of size n
	 * @param n the size of the vector
	 */
	public Vector(int n) {
		
		super(n, 1);
	}

	/**
	 * Constructs a wrapping vector for entries
	 * @param entries the underlying data for the new vector
	 */
	public Vector(ComplexNumber[] entries) {
		
		super(parse(entries));
	}
	
	/**
	 * Constructs a wrapping vector for entries
	 * @param entries the underlying data for the new vector
	 */
	public Vector(float[] entries) {
	
		super(parse(entries));
	}
	
	// wraps a 1d array into a 2d column array
	// used as a helper to the Vector(ComplexNumber[]) constructor
	private static ComplexNumber[][] parse(ComplexNumber[] array) {
		
		ComplexNumber[][] ret = new ComplexNumber[array.length][1];
		for (int i = 0; i < array.length; i++) {
			ret[i][0] = array[i];
		}
		
		return ret;
	}
	
	// wraps a 1d array into a 2d column array
	// used as a helper to the Vector(float[]) constructor
	private static float[][] parse(float[] array) {
		
		float[][] ret = new float[array.length][1];
		for (int i = 0; i < array.length; i++) {
			ret[i][0] = array[i];
		}
		
		return ret;
	}
	
	/**
	 * Gets a value at a certain coordinate
	 * @param i the coordinate whose value we wish to retrieve
	 * @return the value at the ith coordinate in the vector
	 */
	public ComplexNumber getAt(int i) {
		
		return super.getAt(i, 0);
	}
	
	/**
	 * Sets the element at a given coordinate to something
	 * @param coord the coordinate to set the value at
	 * @param val the value to set at the given coordinate
	 */
	public void set(int coord, ComplexNumber val) {
		
		super.set(coord, 0, val);
	}
	
	/**
	 * Returns the dimension of the vector space that the vector resides in
	 * @return the number of coordinates of the vector
	 */
	public int dim() {
		
		return super.rows();
	}
	
	/**
	 * Adds two vectors
	 * @param v the vector to add to this one
	 * @return the vector sum of this vector and v
	 * @throws DimensionMismatchException
	 */
	public Vector add(Vector v) throws DimensionMismatchException {
		
		return super.add((Matrix)v).toVector();
	}

	/**
	 * Subtracts a vector from this one
	 * @param v the vector to subtract
	 * @return the vector difference of this vector and v
	 * @throws DimensionMismatchException
	 */
	public Vector subtract(Vector v) throws DimensionMismatchException {
		
		return super.add((Matrix)v.scale(-1)).toVector();
	}
	
	/**
	 * Compute the standard inner product (dot product) of two vectors
	 * @param v the vector to dot against this
	 * @return the inner product <this, v> of this vector with v 
	 */
	public ComplexNumber dot(Vector v) throws DimensionMismatchException {
		
		if (this.dim() != v.dim()) {
			throw new DimensionMismatchException();
		}
		
		//return this.conjugateTranspose().multiply(v).getAt(0, 0);
		
		ComplexNumber dotprod = new ComplexNumber(0, 0);
		for (int i = 0; i < dim(); i++) {
			dotprod = dotprod.add(this.getAt(i).multiply(v.getAt(i)));
		}
		return dotprod;
	}
	
	/**
	 * Compute the cross product between two vectors
	 * @param v The vector to cross this with
	 * @return The cross product between this and v, this x v
	 * @throws DimensionMismatchException
	 */
	public Vector cross(Vector v) throws DimensionMismatchException {
		
		// cross product is only defined for vectors in F^3, e.g. R^3
		if (this.dim() != 3 || v.dim() != 3) {
			throw new DimensionMismatchException();
		}
		
		ComplexNumber[] coords = new ComplexNumber[3];
		
		coords[0] = getAt(2).multiply(v.getAt(3)).subtract(getAt(3).multiply(v.getAt(2)));
		coords[1] = getAt(3).multiply(v.getAt(1)).subtract(getAt(1).multiply(v.getAt(3)));
		coords[2] = getAt(1).multiply(v.getAt(2)).subtract(getAt(2).multiply(v.getAt(1)));
		
		return new Vector(coords);
	}
	
	/**
	 * Return the vector that is the projection of this onto a certain vector
	 * @param onto the vector onto which we project this
	 * @return the projection proj_onto(this)
	 * @throws DimensionMismatchException
	 */
	public Vector proj(Vector onto) throws DimensionMismatchException {
		
		if (this.dim() != onto.dim()) {
			throw new DimensionMismatchException();
		}
		
		return onto.scale(onto.dot(this).divide(onto.dot(onto))).toVector();
		//return onto.scale(onto.dot(this).multiply(1/Math.pow(onto.length(),2))).toVector();

	}
	
	/**
	 * Returns the unit vector pointing in the same direction as this vector
	 * @return the unit vector pointing in the same direction as this vector
	 */
	public Vector normalize() {
		
		ComplexNumber factor = new ComplexNumber(0, 0);
		
		for (int i = 0; i < dim(); i++) {
			factor = factor.add(getAt(i).multiply(getAt(i)));
		}
		
		factor = factor.sqrt();
		System.out.println("asdf: " + factor);
		
		ComplexNumber[] entries = new ComplexNumber[dim()];
		for (int i = 0; i < dim(); i++) {
			entries[i] = getAt(i).divide(factor);
		}
		
		return new Vector(entries);
	}
	
	/**
	 * Returns a string representation of this vector
	 * @return a string representation <a, b, ...>
	 */
	public String toString() {
		
		String str = "<";
		for (int i = 0; i < this.dim(); i++) {
			str += this.getAt(i);
			if (i != this.dim() - 1) {
				str += ", ";
			}
			
		}
		str += ">";
		
		return str;
	}
}
