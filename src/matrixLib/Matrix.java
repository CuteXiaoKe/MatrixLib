package matrixLib;

/**
 * Represents a matrix over either R or C
 * @author Bryan Cuccioli
 *
 */
public class Matrix {

	private ComplexNumber[][] matrix;
	private int rows, cols;
	
	/**
	 * Blank constructor
	 */
	public Matrix() {
		
	}
	
	/**
	 * Constructs the nxn identity matrix
	 * @param n the number of rows and columns in this identity matrix
	 */
	public Matrix(int n) {
		
		matrix = new RealNumber[n][n];
		rows = n;
		cols = n;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				matrix[i][j] = new RealNumber((i == j) ? 1 : 0);
			}
		}
	}

	/**
	 * Constructs the matrix that wraps the array mat[][]
	 * @param mat the data to go in the matrix
	 */
	public Matrix(float[][] mat) {
		
		matrix = new RealNumber[mat.length][mat[0].length];
		
		for (int i = 0; i < mat.length; i++) {
			for (int j = 0; j < mat[0].length; j++) {
				matrix[i][j] = new RealNumber(mat[i][j]);
			}
		}
		
		rows = mat.length;
		cols = mat[0].length; // the number of columns in the matrix
	}
	
	/**
	 * Constructs the matrix that wraps the array mat[][]
	 * @param mat the data to go in the matrix
	 */
	public Matrix(ComplexNumber[][] mat) {
		
		matrix = new ComplexNumber[mat.length][mat[0].length];
		
		for (int i = 0; i < mat.length; i++) {
			for (int j = 0; j < mat[0].length; j++) {
				matrix[i][j] = mat[i][j];
			}
		}
		
		rows = mat.length;
		cols = mat[0].length; // the number of columns in the matrix
	}
	
	/**
	 * Create a new unpopulated rxc matrix
	 * @param r the number of rows in the matrix
	 * @param c the number of columns in the matrix
	 */
	public Matrix(int r, int c) {
		
		matrix = new ComplexNumber[r][c];
		this.rows = r;
		this.cols = c;
		
		// build the zero matrix
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				matrix[i][j] = new ComplexNumber(0, 0);
			}
		}
	}

	/**
	 * Constructs a matrix from a collection of vectors
	 * @param vectors the vectors to create the matrix from
	 * @throws DimensionMismatchException
	 */
	public Matrix(Vector[] vectors) throws DimensionMismatchException {
		
		this.rows = vectors[0].dim();
		this.cols = vectors.length;
		matrix = new ComplexNumber[rows][cols];
		
		for (int j = 0; j < cols; j++) {
			if (j != 0 && vectors[j].dim() != vectors[j-1].dim()) {
				// the vectors did not all have the same number of coordinates
				throw new DimensionMismatchException();
			}
			for (int i = 0; i < rows; i++) {
				matrix[i][j] = vectors[j].getAt(i);
			}
		}
	}
	
	/**
	 * Retrieves the nth vector from the matrix
	 * @param n the column/vector to get from the matrix
	 * @return the nth vector from the matrix
	 */
	public Vector getVector(int n) {
		
		ComplexNumber[] entries = new ComplexNumber[rows()];
		
		for (int i = 0; i < rows(); i++) {
			entries[i] = matrix[i][n];
		}
		
		return new Vector(entries);
	}
	
	/**
	 * Gets the element at location (r,c) in the matrix
	 * @param r the row to retrieve the element from
	 * @param c the column to retrieve the element from
	 * @return the element in the matrix at (r,c)
	 */
	public ComplexNumber getAt(int r, int c) {
		return matrix[r][c];
	}
	
	/**
	 * Sets the element at (r, c) to something
	 * @param r the row of the element to set
	 * @param c the column of the element to set
	 * @param val the value to set (r, c) to
	 */
	public void set(int r, int c, ComplexNumber val) {
		
		matrix[r][c] = val;
	}
	
	/**
	 * Returns the underlying data array
	 * @return the data array underlying the matrix
	 */
	public ComplexNumber[][] getMatrix() {
		return matrix;
	}

	/**
	 * Computes the length (Frobenius norm) of the matrix,
	 * which is the square root of the sum of each element squared
	 * @return the length (Frobenius norm) of the matrix
	 */
	public double length() {
		
		float sum = 0;
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				sum += Math.pow(matrix[i][j].Re(), 2) + Math.pow(matrix[i][j].Im(), 2); 
			}
		}
		
		return Math.sqrt(sum);
	}

	/**
	 * Returns the matrix transpose of this matrix
	 * @return the transpsoe of this matrix
	 */
	public Matrix transpose() {
		
		ComplexNumber[][] trans = new ComplexNumber[cols][rows];
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				trans[i][j] = matrix[j][i];
			}
		}
		
		return new Matrix(trans);
	}
	
	/**
	 * Returns the conjugate transpose of this matrix
	 * @return the matrix that is the conjugate transpose of this matrix
	 */
	public Matrix conjugateTranspose() {
		
		ComplexNumber[][] ct = new ComplexNumber[cols][rows];
		
		for (int i = 0; i < cols; i++) {
			for (int j = 0; j < rows; j++) {
				ct[i][j] = matrix[j][i].conjugate();
			}
		}
		
		return new Matrix(ct);
	}
	
	/**
	 * Returns the number of rows in the matrix
	 * @return the number of rows
	 */
	public int rows() {
		return rows;
	}
	
	/**
	 * Returns the number of columns in the matrix
	 * @return the number of columns
	 */
	public int cols() {
		return cols;
	}

	/**
	 * Returns the result of multiplying this matrix by m
	 * @param m the matrix to multiply this one by
	 * @return the product this*m
	 * @throws DimensionMismatchException
	 */
	public Matrix multiply(Matrix m) throws DimensionMismatchException {
		
		// can only multiply matrices of dimension nxm by mxp
		if (cols != m.rows()) {
			throw new DimensionMismatchException();
		}
		
		Matrix prod = new Matrix(rows, m.cols());
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < m.cols(); j++) {
				for (int k = 0; k < cols; k++) {
					prod.set(i, j, prod.getAt(i, j).add(matrix[i][k].multiply(m.getAt(k, j))));
				}
			}
		}
		
		//return new Matrix(prod);
		return prod;
	}
	
	/**
	 * Adds this matrix to m
	 * @param m the matrix to add to this one
	 * @return the sum this+m
	 * @throws DimensionMismatchException
	 */
	public Matrix add(Matrix m) throws DimensionMismatchException {
		
		if (m.rows() != rows || m.cols() != cols) {
			throw new DimensionMismatchException();
		}
		
		Matrix sum = new Matrix(rows, cols);
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				sum.set(i, j, matrix[i][j].add(m.getAt(i, j)));
			}
		}
		
		return sum;
	}

	/**
	 * Multiplies each element of the matrix by a complex number
	 * @param factor the scalar to multiply the matrix by
	 * @return the matrix factor*this
	 */
	public Matrix scale(ComplexNumber factor) {
		
		Matrix scaled = new Matrix(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				scaled.set(i, j, matrix[i][j].multiply(factor));
			}
		}
		
		return scaled;
	}

	/**
	 * Multiplies each element of the matrix by a float
	 * @param factor the scalar to multiply the matrix by
	 * @return the matrix factor*this
	 */
	public Matrix scale(float factor) {
		
		Matrix scaled = new Matrix(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				scaled.set(i, j, matrix[i][j].multiply(factor));
			}
		}
		
		return scaled;
	}

	/**
	 * Performs a QR decomposition on this matrix
	 * @return an ordered pair {Q, R}
	 * @throws DimensionMismatchException 
	 */
	public Matrix[] QRDecompose() throws DimensionMismatchException {
		
		Matrix[] qr = new Matrix[2];
		Vector[] u = new Vector[cols()];
		Vector[] e = new Vector[cols()];
		
		for (int i = 0; i < cols(); i++) {
			
			u[i] = getVector(i);
			System.out.println("u["+i+"]: " + u[i]);
			for (int j = 0; j < i; j++) {
				u[i] = u[i].subtract(getVector(i).proj(e[j]));
			}
			System.out.println(u[i]);
			e[i] = u[i].normalize();
			System.out.println("e["+i+"]: " + e[i]);
		}
		
		qr[0] = new Matrix(e);
		qr[1] = qr[0].transpose().multiply(this);
		
		return qr;
	}
	
	/**
	 * Performs row reduction (Gauss-Jordan elimination) on this matrix
	 * @return the reduced row echelon form of this matrix
	 */
	public Matrix rowReduce() {
		
		return this; // placeholder
	}
	
	/**
	 * Gets the underlying data array for this matrix
	 * @return the underlying data array for this matrix
	 */
	protected ComplexNumber[][] getData() {
		
		return this.matrix;
	}
	
	/**
	 * Converts this matrix to a vector through the standard isomorphism
	 * @return a vector with the elements from this matrix
	 */
	public Vector toVector() {
		
		Vector v = new Vector(rows * cols);
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				v.set(i*cols + j, getAt(i, j));
			}
		}
		
		return v;
	}
	
	public String toString() {
		
		String mstr = "[";

		for (int i = 0; i < rows; i++) {
			if (i != 0) mstr += " "; // align left margin horizontally
			mstr += "[";

			for (int j = 0; j < cols; j++) {
				mstr += getAt(i, j).toString();
				if (j != cols - 1) mstr += ", ";
				else {
					mstr += "]";
					if (i == rows - 1) mstr += "]";
					else mstr += "\n";
				}
			}
		}
		
		return mstr;
	}
	
	/**
	 * Tells if two matrices are equal (if all of their corresponding elements are equal)
	 * @param m the matrix to compare this one too
	 * @return whether the two matrices are equal
	 */
	public boolean equals(Matrix m) {
		
		// have to have matching dimension to be equal
		if (rows != m.rows() || cols != m.cols()) {
			return false;
		}
		
		// each element has to match for them to be equal
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (!m.getAt(i, j).equals(matrix[i][j])) {
					return false;
				}
			}
		}
		
		return true; // if it got this far, they are equal
	}
}
