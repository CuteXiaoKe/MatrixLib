package matrixLib;

/**
 * A library of common operations that are restricted to square matrices
 * @author Bryan Cuccioli
 */

public class SquareMatrixOps {

	// inverts a lower triangular matrix; called by inverse()
	private static Matrix inverse_lt(Matrix tri) {
		ComplexNumber[][] tri_inv = new ComplexNumber[tri.rows()][tri.rows()];

		for (int j = 0; j < tri.rows(); j++) {
			tri_inv[0][j] = new ComplexNumber(0,0);
		}
		
		// compute the inverse of the triangular matrix, which is easier
		tri_inv[0][0] = tri.getAt(0,0).reciprocal();
		for (int i = 1; i < tri.rows(); i++) {
			// just compute the reciprocal of the diagonal elements
			if (tri.getAt(i,i).isZero()) {
				return null; // if 0 is on the diagonal, matrix is singular
			}
			tri_inv[i][i] = tri.getAt(i,i).reciprocal();
			
			for (int j = 0; j < i; j++) {
				// found by solving L L^-1 = I
				ComplexNumber sum = new ComplexNumber(0,0);
				for (int k = j; k < i; k++) {
					sum = sum.add(tri.getAt(i,k).multiply(tri_inv[k][j]));
				}
				tri_inv[i][j] = sum.divide(tri.getAt(i,i)).negative();
			}
			// inverse of lower triangular is lower triangular, so fill the rest with 0's
			for (int j = i+1; j < tri.rows(); j++) {
				tri_inv[i][j] = new ComplexNumber(0,0);
			}
		}
		
		return new Matrix(tri_inv);
	}
	
	/**
	 * Computes the inverse of the matrix via Gauss-Jordan elimination
	 * @param m the matrix of which the inverse is computed
	 * @return the corresponding inverse matrix, or null if the matrix is not invertible
	 * @throws NotSquareException the matrix is not square
	 */
	public static Matrix inverse(Matrix m) throws NotSquareException {
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}

		// test for common easy cases
		if (Pattern.isUpperTriangular(m)) {
			return inverse_lt(m.transpose());
		}
		else if (Pattern.isLowerTriangular(m)) {
			return inverse_lt(m);
		}
		else if (Pattern.isHermetian(m)) {
			Matrix tri = Factorization.choleskyDecompose(m);
			Matrix l_inv = inverse_lt(tri);
			return l_inv.conjugateTranspose().multiply(l_inv);
		}
		
		// otherwise we have to use the generic algorithm
		ComplexNumber[][] augmented = new ComplexNumber[m.rows()][m.cols()*2];
		
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.cols()*2; j++) {
				if (j < m.cols()) {
					augmented[i][j] = m.getAt(i, j);
				}
				else {
					augmented[i][j] = new ComplexNumber((i==j-m.cols())?1:0,0);
				}
			}
		}
		
		Matrix aug_rref = (new Matrix(augmented)).rref();
		
		// if the left wasn't reduced to Id, inverse doesn't exist
		for (int i = 0; i < m.rows(); i++) {
			boolean all_zero = true;
			for (int j = 0; j < m.cols(); j++) {
				if (!aug_rref.getAt(i, j).isZero()) {
					all_zero = false;
					break;
				}
			}
			if (all_zero) {
				// there was an all 0 row; not identity
				return null;
			}
		}
		
		ComplexNumber[][] inv = new ComplexNumber[m.rows()][m.cols()];
		
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.cols(); j++) {
				inv[i][j] = aug_rref.getAt(i, j+m.cols());
			}
		}
		
		return new Matrix(inv);
	}
	
	/**
	 * Returns the determinant of this matrix in O(n^3) time complexity
	 * @param m the matrix whose determinant is computed
	 * @return the determinant of the matrix
	 * @throws NotSquareException the matrix is not square
	 */
	public static ComplexNumber determinant(Matrix m) throws NotSquareException {

		// try shortcuts for common sizes first
		if (m.rows() == 1) {
			return m.getAt(0,0);
		}
		else if (m.rows() == 2) {
			return m.getAt(0,0).multiply(m.getAt(1,1)).subtract(m.getAt(1,0).multiply(m.getAt(0,1)));
		}
		else {
			return m.rref(1).getAt(0, 0);
		}
	}
	
	/**
	 * Returns the trace of the matrix, the sum of the diagonal elements
	 * @param m the matrix whose trace to compute
	 * @return the trace of the matrix
	 * @throws NotSquareException the given matrix is not square
	 */
	public static ComplexNumber trace(Matrix m) throws NotSquareException {
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		ComplexNumber tr = new ComplexNumber(0, 0);
		
		for (int i = 0; i < m.rows(); i++) {
			tr = tr.add(m.getAt(i, i));
		}
		
		return tr;
	}
	
	public static ComplexNumber[] fast_eigenvalues(Matrix m) {
		
		return null;
	}
	
	/**
	 * Returns a list of the eigenvalues of the matrix
	 * @param m the matrix whose eigenvalues we compute
	 * @return an array containing the eigenvalues of the matrix
	 * @throws NotSquareException the given matrix is not square 
	 */
	public static ComplexNumber[] eigenvalues(Matrix m) {
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}

		// max amount of eigenvalues is the dimension of the matrix
		ComplexNumber[] evals = new ComplexNumber[m.rows()];
		double[][] shift_arr = {{2.8,0,0},{0,2.8,0},{0,0,2.8}};
		Matrix shift = new Matrix (shift_arr);
		
		Matrix temp = Pattern.hessenberg(m);
		//Matrix temp = new Matrix(m.getData());
		
		ComplexNumber.setEpsilon(1e-8);
		int count = 0;

		System.out.println("loop");
		while (!Pattern.isUpperTriangular(temp)) {
			count++;
			//System.out.println("QR decomposing...");
			
			//Matrix[] qr = Factorization.QRDecompose(temp.subtract(shift));
			//temp = qr[1].multiply(qr[0].transpose()).add(shift);
			
			Matrix[] qr = Factorization.QRDecompose(temp);
			temp = qr[1].multiply(qr[0]);
			
			if(count>650000)break;
		}
		System.out.printf("Completed in %d iterations.\n", count);
		
		// temp is now upper triangular, so the evals are on the diagonal
		for (int i = 0; i < m.rows(); i++) {
			evals[i] = temp.getAt(i, i);
		}
		
		return evals;
	}

	/**
	 * Computes the normalized eigenvectors
	 * @param m the matrix whose eigenvectors are computed
	 * @return an array of normalized eigenvectors for the matrix
	 * @throws NotSquareException the given matrix is not square
	 */
	public static Vector[] eigenvectors(Matrix m) throws NotSquareException {

		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		return eigenvectors(m, SquareMatrixOps.eigenvalues(m));
	}
	
	/**
	 * Computes the normalized eigenvectors of the matrix corresponding to certain eigenvalues
	 * @param m the matrix whose eigenvectors are being computed 
	 * @param evals the list of eigenvalues to compute the associated eigenvectors of
	 * @return an array of normalized eigenvectors for the matrix
	 */
	public static Vector[] eigenvectors(Matrix m, ComplexNumber[] evals) {
		
		Vector[] evecs = new Vector[evals.length];
		int vec_pos = 0;
		
		ComplexNumber[] test_vec = new ComplexNumber[evals.length];
		for (int i = 0; i < evals.length; i++) {
			test_vec[i] = new ComplexNumber(1, 0);
		}
		
		for (ComplexNumber ev : evals) {
			
			Matrix diag = new Matrix(m.getData());
			for (int i = 0; i < evals.length; i++) {
				diag.set(i, i, diag.getAt(i, i).subtract(ev));
			}
			diag = SquareMatrixOps.inverse(diag);
			
			Vector prev = new Vector(test_vec);
			Vector next = diag.multiply(prev).normalize();
			
			do {
				prev = next;
				next = diag.multiply(prev).normalize();
			//} while (!next.isAlmost(prev));
			} while (!next.equals(prev));
			evecs[vec_pos++] = next.normalize();
		}
		
		return evecs;
	}	
	
	/**
	 * Computes the value of the matrix raised to a particular power
	 * @param m the matrix that is raised to the specified power
	 * @param power the power to which the matrix is raised
	 * @return the matrix raised to the specified power
	 * @throws NotSquareException the supplied matrix is not square
	 */
	public static Matrix pow(Matrix m, int power) throws NotSquareException {
		
		if (m.rows() != m.cols()) {
			// A^p is only really defined for square A
			throw new NotSquareException();
		}
		
		if (power == 0) {
			// a matrix to the zero power is the identity matrix
			return new Matrix(m.rows());
		}
		else if (power < 0) {
			// A^-p = (A^-1)^p
			return pow(SquareMatrixOps.inverse(m), -power);
		}
		else {
			// multiply the matrix by itself a number of times
			Matrix temp = new Matrix(m.getData());
			for (int i = 0; i < power-1; i++) {
				temp = temp.multiply(m);
			}
			return temp;
		}
	}
}
