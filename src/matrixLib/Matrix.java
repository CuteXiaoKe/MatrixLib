package matrixLib;

/**
 * Represents a matrix over either R or C
 * @author bryan
 *
 */
public class Matrix {

	private float[][] matrix;
	private int rows, cols;
	
	/**
	 * Constructs the nxn identity matrix
	 * @param n the number of rows and columns in this identity matrix
	 */
	public Matrix(int n) {
		
		matrix = new float[n][n];
		rows = n;
		cols = n;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				matrix[i][j] = (i == j) ? 1 : 0;
			}
		}
	}

	/**
	 * Constructs the matrix that wraps the array mat[][]
	 * @param mat the data to go in the matrix
	 */
	public Matrix(float[][] mat) {
		
		matrix = new float[mat.length][mat[0].length];
		
		for (int i = 0; i < mat.length; i++) {
			for (int j = 0; j < mat[0].length; j++) {
				matrix[i][j] = mat[i][j];
			}
		}
		
		rows = mat.length;
		cols = mat[0].length; // the number of columns in the matrix
	}
	
	/**
	 * Gets the element at location (r,c) in the matrix
	 * @param r the row to retrieve the element from
	 * @param c the column to retrieve the element from
	 * @return the element in the matrix at (r,c)
	 */
	public float getAt(int r, int c) {
		return matrix[r][c];
	}
	
	/**
	 * Returns the underlying data array
	 * @return the data array underlying the matrix
	 */
	public float[][] getMatrix() {
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
				sum += matrix[i][j]*matrix[i][j];
			}
		}
		
		return Math.sqrt(sum);
	}

	/**
	 * Returns the matrix transpose of this matrix
	 * @return the transpsoe of this matrix
	 */
	public Matrix getTranspose() {
		
		float[][] trans = new float[cols][rows];
		
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
	public int getRows() {
		return rows;
	}
	
	/**
	 * Returns the number of columns in the matrix
	 * @return the number of columns
	 */
	public int getCols() {
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
		if (cols != m.getRows()) {
			throw new DimensionMismatchException();
		}
		
		float[][] prod = new float[rows][m.getCols()];
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {

				float sum = 0;
				for (int k = 0; k < cols; k++) {
					sum += matrix[i][k] * m.getAt(k, j);
				}
				prod[i][j] = sum;
			}
		}
		
		return new Matrix(prod);
	}
	
	/**
	 * Adds this matrix to m
	 * @param m the matrix to add to this one
	 * @return the sum this+m
	 * @throws DimensionMismatchException
	 */
	public Matrix add(Matrix m) throws DimensionMismatchException {
		
		if (m.getRows() != rows || m.getCols() != cols) {
			throw new DimensionMismatchException();
		}
		
		float[][] sum = new float[rows][cols];
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				sum[i][j] = matrix[i][j] + m.getAt(i, j);
			}
		}
		
		return new Matrix(sum);
	}

	/**
	 * Multiplies each element of the matrix by a field scalar
	 * @param factor the scalar to multiply the matrix by
	 * @return the matrix factor*this
	 */
	public Matrix scale(float factor) {
		
		float[][] scaled = new float[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				scaled[i][j] = factor * matrix[i][j];
			}
		}
		
		return new Matrix(scaled);
	}

	public String toString() {
		
		String mstr = "[";

		for (int i = 0; i < rows; i++) {
			if (i != 0) mstr += " "; // align left margin horizontally
			mstr += "[";

			for (int j = 0; j < cols; j++) {
				mstr += Float.toString(getAt(i, j));
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
	protected float[][] getData() {
		
		return this.matrix;
	}
	
	public boolean equals(Matrix m) {
		
		// have to have matching dimension to be equal
		if (rows != m.getRows() || cols != m.getCols()) {
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
