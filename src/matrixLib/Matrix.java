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
	 * Computes the length of the matrix,
	 * which is the square root of the sum of each element squared
	 * @return the length of the matrix
	 */
	public double getLength() {
		
		float sum = 0;
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				sum += Math.pow(matrix[i][j].abs(), 2);
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
		//ComplexNumber[][] prod = new ComplexNumber[rows][m.cols()];
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				for (int k = 0; k < cols; k++) {
					prod.set(i, j, prod.getAt(i, j).add(matrix[i][k].multiply(m.getAt(k, j))));
					//prod[i][j] = prod[i][j].add(matrix[i][k].multiply(m.getAt(k, j)));
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
	 * Multiplies each element of the matrix by a field scalar
	 * @param factor the scalar to multiply the matrix by
	 * @return the matrix factor*this
	 */
	public Matrix scale(float factor) {
		
		Matrix scaled = new Matrix(rows, cols);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				scaled.set(i, j, matrix[i][j].multiply(factor));
			}
		}
		
		return scaled;
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
	 * Gets the underlying data array for this matrix
	 * @return the underlying data array for this matrix
	 */
	protected ComplexNumber[][] getData() {
		
		return this.matrix;
	}
	
	public boolean equals(Matrix m) {
		
		// have to have matching dimension to be equal
		if (rows != m.rows() || cols != m.cols()) {
			return false;
		}
		
		// each element has to match for them to be equal
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (m.getAt(i, j) != matrix[i][j]) {
					return false;
				}
			}
		}
		
		return true; // if it got this far, they are equal
	}
}
