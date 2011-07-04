package matrixLib;

/**
 * A library of common operations that are restricted to square matrices
 * @author Bryan Cuccioli
 */

public class SquareMatrixOps {

	/**
	 * Tells whether the matrix is upper triangular
	 * @param m the matrix to check
	 * @return whether the matrix is upper triangular or not
	 */
	public static boolean isUpperTriangular(Matrix m) {
		
		// the matrix must be square to be triangular
		if (m.rows() != m.cols()) {
			return false;
		}
		
		for (int i = 1; i < m.rows(); i++) {
			for (int j = 0; j < i; j++) {
				if (!m.getAt(i,j).equals(new ComplexNumber(0, 0))) {
					return false;
				}
			}
		}
		
		return true;
	}

	/**
	 * Tells whether the matrix is lower triangular
	 * @param m the matrix to check
	 * @return whether the matrix is lower triangular or not
	 */
	public static boolean isLowerTriangular(Matrix m) {
		
		// matrix must be square to be triangular
		if (m.rows() != m.cols()) {
			return false;
		}

		for (int i = 0; i < m.rows(); i++) {
			for (int j = i+1; j < m.cols(); j++) {
				if (!m.getAt(i,j).equals(new ComplexNumber(0, 0))) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Tells if the matrix is diagonal (upper and lower triangular) or not
	 * @param m the matrix to check
	 * @return whether the matrix is diagonal
	 */
	public static boolean isDiagonal(Matrix m) {
		
		// a matrix is diagonal if it is upper and lower triangular
		
		return isUpperTriangular(m) && isLowerTriangular(m);
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
		
		ComplexNumber[][] augmented = new ComplexNumber[m.rows()][m.cols()*2];
		
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.cols()*2; j++) {
				if (j < m.cols()) {
					augmented[i][j] = m.getAt(i, j);
				}
				else {
					if (i == j - m.cols()) {
						augmented[i][j] = new ComplexNumber(1,0);
					}
					else {
						augmented[i][j] = new ComplexNumber(0,0);
					}
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
	 * Computes the tridiagonalized form of the matrix via Householder reflections
	 * @param m the matrix to tridiagonalize
	 * @return the tridiagonalized form of the matrix
	 * @throws NotSquareException the supplied matrix is not square
	 */
	public static Matrix tridiagonalized(Matrix m) throws NotSquareException {
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		Matrix hess = new Matrix(m.getData());
		
		for (int k = 0; k < m.rows(); k++) {
			ComplexNumber sum = new ComplexNumber(0,0);
			for (int j = k+1; j < m.rows(); j++) {
				sum = sum.add(hess.getAt(j, k).multiply(hess.getAt(j, k)));
			}
			ComplexNumber alpha = sum.sqrt().multiply(-Math.signum(hess.getAt(k+1,k).Re()));
			ComplexNumber r = alpha.multiply(alpha).subtract(hess.getAt(k+1,k).multiply(alpha)).multiply(.5).sqrt().reciprocal();
			
			ComplexNumber[] vec = new ComplexNumber[m.rows()];
			for (int j = 0; j <= k; j++) {
				vec[j] = new ComplexNumber(0,0);
			}
			vec[k+1] = hess.getAt(k+1,k).subtract(alpha).multiply(0.5*r.Re());
			for (int j = k+2; j < m.rows(); j++) {
				vec[j] = hess.getAt(j,k).multiply(0.5*r.Re());
			}
			
			ComplexNumber[][] p = new ComplexNumber[m.rows()][m.rows()];
			for (int i = 0; i < m.rows(); i++) {
				for (int j = 0; j < m.rows(); j++) {
					p[i][j] = (i==j) ? new ComplexNumber(1,0) : new ComplexNumber(0,0);
					p[i][j] = p[i][j].subtract(vec[i].multiply(vec[j]).multiply(2.0));
				}
			}
			
			Matrix hh = new Matrix(p);
			hess = hh.multiply(hess).multiply(hh);
		}
		
		return hess;
	}
	
	/**
	 * Returns the upper Hessenberg form of this matrix,
	 * which has zero entries below the first subdiagonal
	 * @param m the matrix whose upper Hessenberg form is computed
	 * @return the upper Hessenberg form of this matrix
	 * @throws NotSquareException the given matrix is not square
	 */
	public static Matrix hessenbergForm(Matrix m) throws NotSquareException {
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		Matrix prev = new Matrix(m.getData());
		
		for (int k = 0; k < m.rows()-2; k++) {
			// x is the kth column, restricted to below the diagonal
			ComplexNumber[] vec = new ComplexNumber[m.rows()-k-1];
			for (int i = k+1; i < m.rows(); i++) {
				vec[i-k-1] = prev.getAt(i, k);
			}
			Vector x = new Vector(vec);
			x.set(0, x.getAt(0).add(new ComplexNumber(x.length(),0)));
			
			Matrix ref = x.reflector();
			ComplexNumber[][] p = new ComplexNumber[k*2][k*2];
			for (int i = 0; i < k*2; i++) {
				for (int j = 0; j < k*2; j++) {
					if (i >= k && j >= k) {
						p[i][j] = ref.getAt(i-k, j-k);
					}
					else if (i == j) {
						p[i][j] = new ComplexNumber(1,0);
					}
					else {
						p[i][j] = new ComplexNumber(0,0);
					}
				}
			}
			Matrix block = new Matrix(p);
			prev = block.transpose().multiply(prev).multiply(block);
		}
		
		return prev;
	}
	
	/**
	 * Tells whether the matrix is symmetric or not
	 * @param m the matrix to check
	 * @return whether the matrix is symmetric
	 */
	public static boolean isSymmetric(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // only square matrices can be symmetric
		}
		
		// a symmetric matrix equals its transpose
		return (m.equals(m.transpose()));
		
		/*for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.rows(); j++) {
				
				if
			}
		}*/
	}
	
	/**
	 * Tells whether the matrix is anti-(skew-)symmetric or not,
	 * that is whether reflecting it over its diagonal is its negative
	 * @param m the matrix to check
	 * @return whether the mattrix is anti-symmetric
	 */
	public static boolean isAntiSymmetric(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // the matrix must be square to be antisymmetric
		}
		
		// an antisymmetric matrix equals the negative of its transpose
		return (m.equals(m.transpose().scale(-1)));
	}
	
	/**
	 * Tells whether the matrix is Hermetian/self-adjoint (equal to its own conjugate transpose)
	 * @param m the matrix to check
	 * @return whether the matrix is Hermetian
	 */
	public static boolean isHermetian(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // must be square to be Hermetian
		}
		
		return m.conjugateTranspose().equals(m);
	}

	/**
	 * Tells whether the matrix is orthogonal (A^T A = I) or not
	 * @param m the matrix to check
	 * @return whether the matrix is orthogonal 
	 */
	public static boolean isOrthogonal(Matrix m) {

		if (m.rows() != m.cols()) {
			return false; // must be square to be orthogonal
		}
		
		return SquareMatrixOps.isIdentity(m.multiply(m.transpose()));
	}
	
	/**
	 * Tells whether the matrix is unitary (its inverse is its conjugate transpose)
	 * @param m the matrix to check
	 * @return whether the matrix is unitary or not
	 */
	public static boolean isUnitary(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // must be square to be unitary
		}
		
		return SquareMatrixOps.isIdentity(m.multiply(m.conjugateTranspose()));
	}
	
	/**
	 * Tells whether the matrix is upper Hessenberg
	 * (whether it has zeros below the first subdiagonal)
	 * @param m the matrix to check
	 * @return whether the matrix is upper Hessenberg
	 */
	public static boolean isUpperHessenberg(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // must be square to be upper Hessenberg
		}
		
		for (int i = 2; i < m.rows(); i++) {
			for (int j = 0; j < i-1; j++) {
				if (!m.getAt(i,j).isZero()) {
					return false;
				}
			}
		}
		
		return true;
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

	private static boolean isAlmostUpperTriangular(Matrix m) {
		
		double epsilon = 1E-6;
		
		for (int i = 1; i < m.rows(); i++) {
			for (int j = 0; j < i; j++) {
				if (m.getAt(i, j).abs() > epsilon) {
					return false;
				}
			}
		}
		
		return true;
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
		
		Matrix temp = m;
		while (!isAlmostUpperTriangular(temp)) {
			//System.out.println("QR decomposing...");
			Matrix[] qr = Factorization.QRDecompose(temp);
			temp = qr[1].multiply(qr[0]);
			//System.out.println(temp);
		}
		
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
			} while (!next.isAlmost(prev));
			evecs[vec_pos++] = next.normalize();
		}
		
		return evecs;
	}
	
	/**
	 * Tells whether the matrix is the identity matrix
	 * @param m the matrix to check
	 * @return whether the matrix is the identity matrix
	 */
	public static boolean isIdentity(Matrix m) {
		
		if (m.rows() != m.cols()) {
			return false; // must be square to be the identity
		}
		
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 0; j < m.cols(); j++) {
				if ((i != j && !m.getAt(i, j).equals(new ComplexNumber(0, 0)))
					|| (i == j && !m.getAt(i, j).equals(new ComplexNumber(1, 0)))) {
					return false;
				}
			}
		}
		
		return true;
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
		else if (power == 1) {
			// A^1 = A
			return m;
		}
		else {
			// multiply the matrix by itself a number of times
			Matrix temp = new Matrix(m.getData());
			for (int i = 0; i < power; i++) {
				temp = temp.multiply(m);
			}
			return temp;
		}
	}
}
