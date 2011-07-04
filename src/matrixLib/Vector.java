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
	public Vector(double[] entries) {
	
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
	// used as a helper to the Vector(double[]) constructor
	private static double[][] parse(double[] array) {
		
		double[][] ret = new double[array.length][1];
		for (int i = 0; i < array.length; i++) {
			ret[i][0] = array[i];
		}
		
		return ret;
	}
	
	/**
	 * Gets a value at a certain coordinate
	 * @param i the coordinate whose value we wish to retrieve
	 * @return the value at the ith coordinate in the vector
	 * @throws ArrayIndexOutOfBoundsException tried to access an invalid coordinate
	 */
	public ComplexNumber getAt(int i) throws IndexOutOfBoundsException {
		
		if (i < 0 || i >= dim()) {
			throw new IndexOutOfBoundsException();
		}
		
		return super.getAt(i, 0);
	}
	
	/**
	 * Sets the element at a given coordinate to something
	 * @param coord the coordinate to set the value at
	 * @param val the value to set at the given coordinate
	 * @throws ArrayIndexOutOfBoundsException tried to access an invalid coordinate
	 */
	public void set(int coord, ComplexNumber val) throws IndexOutOfBoundsException {
		
		if (coord < 0 || coord >= dim()) {
			throw new IndexOutOfBoundsException();
		}
		
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
	 * @throws DimensionMismatchException the vectors to add do not have matching dimension
	 */
	public Vector add(Vector v) throws DimensionMismatchException {
		
		return super.add((Matrix)v).toVector();
	}

	/**
	 * Subtracts a vector from this one
	 * @param v the vector to subtract
	 * @return the vector difference of this vector and v
	 * @throws DimensionMismatchException the vectors do not have matching dimension
	 */
	public Vector subtract(Vector v) throws DimensionMismatchException {
		
		return super.add((Matrix)v.scale(-1)).toVector();
	}
	
	/**
	 * Compute the Hermetian inner product (dot product) of two vectors
	 * @param v the vector to dot against this
	 * @return the inner product <this, v> of this vector with v
	 * @throws DimensionMismatchException the vectors do not have matching dimension 
	 */
	public ComplexNumber dot(Vector v) throws DimensionMismatchException {
		
		if (this.dim() != v.dim()) {
			 new DimensionMismatchException();
		}
		
		ComplexNumber dotprod = new ComplexNumber(0, 0);
		for (int i = 0; i < dim(); i++) {
			dotprod = dotprod.add(this.getAt(i).multiply(v.getAt(i).conjugate()));
		}
		return dotprod;
	}
	
	/**
	 * Compute the cross product between two vectors
	 * @param v The vector to cross this with
	 * @return The cross product between this and v, this x v
	 * @throws DimensionMismatchException the vectors do not have three coordinates
	 */
	public Vector cross(Vector v) throws DimensionMismatchException {
		
		// cross product is only defined for vectors in F^3, e.g. R^3
		if (this.dim() != 3 || v.dim() != 3) {
			 new DimensionMismatchException();
		}
		
		ComplexNumber[] coords = new ComplexNumber[3];
		
		coords[0] = getAt(1).multiply(v.getAt(2)).subtract(getAt(2).multiply(v.getAt(1)));
		coords[1] = getAt(2).multiply(v.getAt(0)).subtract(getAt(0).multiply(v.getAt(2)));
		coords[2] = getAt(0).multiply(v.getAt(1)).subtract(getAt(1).multiply(v.getAt(0)));
		
		return new Vector(coords);
	}
	
	/**
	 * Return the vector that is the projection of this onto a certain vector
	 * @param onto the vector onto which we project this
	 * @return the projection proj_onto(this)
	 * @throws DimensionMismatchException the vectors do not have matching dimension
	 */
	public Vector proj(Vector onto) throws DimensionMismatchException {
		
		if (this.dim() != onto.dim()) {
			 new DimensionMismatchException();
		}
		
		ComplexNumber factor = onto.dot(this).divide(onto.dot(onto));
		ComplexNumber[] entries = new ComplexNumber[dim()];
		
		for (int i = 0; i < dim(); i++) {
			entries[i] = onto.getAt(i).multiply(factor);
		}
		
		return new Vector(entries);
	}
	
	/**
	 * Returns the unit vector pointing in the same direction as this vector
	 * @return the unit vector pointing in the same direction as this vector
	 */
	public Vector normalize() {
		
		ComplexNumber norm = new ComplexNumber(0, 0);
		for (int i = 0; i < dim(); i++) {
			norm = norm.add(getAt(i).multiply(getAt(i).conjugate()));
		}
		norm = norm.sqrt();
		
		ComplexNumber[] entries = new ComplexNumber[dim()];
		for (int i = 0; i < dim(); i++) {
			entries[i] = getAt(i).divide(norm);
		}
		
		return new Vector(entries);
	}
	
	/**
	 * Tells whether the absolute value of each component of this vector
	 * is within epsilon of the corresponding component of the given vector
	 * @param v the vector to compare this one too
	 * @return whether this vector is "almost" v
	 */
	public boolean isAlmost(Vector v) {
		
		double epsilon = 1e-6;
		
		for (int i = 0; i < dim(); i++) {
			if (Math.abs(getAt(i).abs() - v.getAt(i).abs()) > epsilon) {
				return false;
			}
		}
		
		return true;
	}

	/**
	 * Generates a unitary matrix whose first column is this vector
	 * @return a unitary matrix whose first column is this vector 
	 */
	public Matrix generateUnitaryMatrix() {
		
		Vector[] overspan = new Vector[this.dim()+1];
		overspan[0] = this.normalize();
		ComplexNumber[] oscontent = new ComplexNumber[this.dim()];
		oscontent[0] = new ComplexNumber(1,0); // debug - should be (1,1)
		
		for (int i = 1; i < this.dim(); i++) {
			oscontent[i] = new ComplexNumber(0,0);
		}
		
		for (int i = 1; i <= this.dim(); i++) {
			overspan[i] = new Vector(oscontent);
			if (i != this.dim()) {
				oscontent[i] = oscontent[i-1];
			}
			oscontent[i-1] = new ComplexNumber(0,0);
		}
		
		// choose a maximal subset of independent vectors
		System.out.println(new Matrix(overspan));
		Matrix allvectors = new Matrix(overspan);
		Matrix rref = allvectors.rref();
		System.out.println(rref);
		Vector[] basis = new Vector[this.dim()];
		int vec_pos = 0;
		for (int i = 0; i < this.dim(); i++) {
			for (int j = 0; j <= this.dim(); j++) {
				if (!getAt(i,j).isZero()) {
					basis[vec_pos++] = allvectors.getVector(j);
					System.out.println("added vector " + j + ": " + allvectors.getVector(j));
					break;
				}
			}
		}
		System.out.println(new Matrix(basis));
		
		// apply Gram-Schmidt to orthogonalize the vectors
		return (new Matrix(basis)).orthonormalize();
	}
	
	/**
	 * Determine the elementary reflector associated with this vector
	 * @return the elementary reflector associated with this vector 
	 */
	public Matrix reflector() {
		
		double factor = 2.0/Math.pow(Norm.pnorm(this, 2), 2);
		ComplexNumber[][] ref = new ComplexNumber[this.dim()][this.dim()];
		
		for (int i = 0; i < this.dim(); i++) {
			for (int j = 0; j < this.dim(); j++) {
				ref[i][j] = (i==j) ? new ComplexNumber(1,0) : new ComplexNumber(0,0);
				ref[i][j] = ref[i][j].subtract(getAt(i).multiply(getAt(j)).multiply(factor));
			}
		}
		
		return new Matrix(ref);
	}
	
	/**
	 * Returns a string representation of this vector <a, b, ...>
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
