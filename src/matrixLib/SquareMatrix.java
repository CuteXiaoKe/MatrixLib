package matrixLib;

/**
 * Represents a square matrix, on which
 * many more operations can be performed
 * @author Bryan Cuccioli
 *
 */
public class SquareMatrix extends Matrix {

	/**
	 * Construct the nxn identity matrix
	 * @param n the dimension of the identity matrix
	 */
	public SquareMatrix(int n) {
		super(n);
	}

	/**
	 * Construct the matrix with specified underlying data
	 * @param mat the underlying data of the matrix
	 * @throws NotSquareException
	 */
	public SquareMatrix(ComplexNumber[][] mat) throws NotSquareException {
		
		super(mat);
		
		if (mat.length != mat[0].length) {
			throw new NotSquareException();
		}
	}

	/**
	 * Construct the matrix with specified underlying data
	 * @param mat the underlying data of the matrix
	 * @throws NotSquareException
	 */
	public SquareMatrix(float[][] mat) throws NotSquareException {
		
		super(mat);
		
		if (mat.length != mat[0].length) {
			throw new NotSquareException();
		}
	}
	
	/**
	 * Tells whether the matrix is upper triangular
	 * @return whether the matrix is upper triangular or not
	 */
	public boolean isUpperTriangular() {
		
		for (int i = 1; i < this.rows(); i++) {
			for (int j = 0; j < i; j++) {
				if (!getAt(i,j).equals(new ComplexNumber(0, 0))) {
					return false;
				}
			}
		}
		
		return true;
	}

	/**
	 * Tells whether the matrix is lower triangular
	 * @return whether the matrix is lower triangular or not
	 */
	public boolean isLowerTriangular() {

		for (int i = 0; i < this.rows(); i++) {
			for (int j = i+1; j < this.cols(); j++) {
				if (!getAt(i,j).equals(new ComplexNumber(0, 0))) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Tells if the matrix is diagonal (upper and lower triangular) or not
	 * @return whether the matrix is diagonal
	 */
	public boolean isDiagonal() {
		
		// a matrix is diagonal if it is upper and lower triangular
		
		return isUpperTriangular() && isLowerTriangular();
	}
	
	/**
	 * Returns the inverse of this matrix
	 * @return the corresponding inverse matrix
	 * @throws SingularMatrixException
	 */
	public Matrix getInverse() throws SingularMatrixException {
				
		if (determinant().equals(new ComplexNumber(0, 0))) {
			// if the matrix has det 0 it is not invertible
			throw new SingularMatrixException();
		}
		
		return this; // placeholder
	}

	/**
	 * Returns the matrix missing row r and column c
	 * @param r the row to remove from the matrix
	 * @param c the column to remove from the matrix
	 * @return the matrix missing row r and column c
	 * @throws Exception
	 */
	protected ComplexNumber[][] getReducted(int r, int c) throws Exception {
		// returns a matrix without row r or column c

		if (this.rows() == 1 || this.cols() == 1) {
			throw new Exception("This matrix is too small to reduce any further.");
		}
		
		int targetRow = 0, targetCol = 0;
		ComplexNumber[][] reduced = new ComplexNumber[this.rows()-1][this.cols()-1];
		
		for (int sourceRow = 0; sourceRow < this.rows(); sourceRow++) {
			// if we're at the row to skip, jump back to the top
			if (sourceRow == r) {
				continue;
			}
			
			for (int sourceCol = 0; sourceCol < this.cols(); sourceCol++) {
				if (sourceCol == c) {
					continue;
				}
				reduced[targetRow][targetCol] = getAt(sourceRow, sourceCol);
				targetCol++;
			}
			// only update targetRow if we haven't jumped over the row to cut
			targetRow++;
			targetCol = 0; // already 0 if we skipped the row
		}
		
		return reduced;
	}

	/**
	 * Returns the determinant of this matrix in O(n!)
	 * @return the determinant of the matrix
	 */
	public ComplexNumber determinant() {

		// if triangular just take the product along the diagonal
		if (isUpperTriangular() || isLowerTriangular()) {
			ComplexNumber prod = new ComplexNumber(1, 0);
			for (int i = 0; i < rows(); i++) {
				prod = prod.multiply(getAt(i, i));
			}
			return prod;
		}
		else {
			// otherwise we have to go through the usual algorithm
			return determinant(this.getData());
		}
	}

	/**
	 * Tells whether the matrix is symmetric or not
	 * @return whether the matrix is symmetric
	 */
	public boolean isSymmetric() {
		
		// a symmetric matrix equals its transpose
		return (this.equals(this.transpose()));
	}
	
	/**
	 * Tells whether the matrix is anti-(skew-)symmetric or not
	 * @return whether the mattrix is anti-symmetric
	 */
	public boolean isAntiSymmetric() {
				
		// an antisymmetric matrix equals the negative of its transpose
		return (this.equals(this.transpose().scale(-1)));
	}

	/**
	 * Tells whether the matrix is the identity matrix
	 * @return whether the matrix is the identity matrix
	 */
	public boolean isIdentity() {
		
		for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols(); j++) {
				if ((i != j && !getAt(i, j).equals(new ComplexNumber(0, 0)))
					|| (i == j && !getAt(i, j).equals(new ComplexNumber(1, 0)))) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Tells whether the matrix is orthogonal (A^T A = I) or not
	 * @return whether the matrix is orthogonal
	 * @throws DimensionMismatchException 
	 */
	public boolean isOrthogonal() throws DimensionMismatchException {
		
		return ((SquareMatrix)this.multiply(this.transpose())).isIdentity();
	}
	
	/**
	 * Computes and returns the LU decomposition of the matrix
	 * @return the result of the LU decomposition {L,D,U}, or null if no LU decomposition is admitted
	 */
	public Matrix[] luDecompose() {
		
		// returns {null, null, null} if the matrix does not admit an LU-decomposition
		
		Matrix[] ldu = new Matrix[3];
		
		
		
		
		return ldu;
	}
	
	// determinant helper method
	private ComplexNumber determinant(ComplexNumber[][] mat) {
		
		// precondition: this is a square matrix
		
		// base case is 2x2 matrix:
		if (mat.length == 2) {
			return mat[0][0].multiply(mat[1][1]).subtract(mat[0][1].multiply(mat[1][0]));
		}
		
		int negator = 1;
		ComplexNumber det = new ComplexNumber(0, 0);
		
		for (int i = 0; i < cols(); i++) {
			try {
				det = det.add(getAt(0, i).multiply(negator).multiply(determinant(getReducted(0, i))));
			} catch (Exception e) {
				// the matrix ran out of rows or columns to remove - shouldn't happen
				e.printStackTrace();
			}
			negator *= -1;
		}
		
		return det;
	}
	
	/**
	 * Returns the trace of the matrix, the sum of the diagonal elements
	 * @return the trace of the matrix
	 */
	public ComplexNumber trace() {
				
		ComplexNumber tr = new ComplexNumber(0, 0);
		
		for (int i = 0; i < rows(); i++) {
			tr = tr.add(getAt(i, i));
		}
		
		return tr;
	}

}
