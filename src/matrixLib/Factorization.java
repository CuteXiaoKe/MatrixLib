package matrixLib;

/**
 * A library of routines for performing different factorizations on matrices
 * @author Bryan Cuccioli
 */

public class Factorization {

	/**
	 * Performs a QR decomposition on the given matrix, writing a matrix as a product of
	 * a unitary matrix Q and upper triangular matrix R
	 * @param m the matrix whose QR factorization we wish to compute
	 * @return an ordered pair {Q, R} 
	 */
	public static Matrix[] QRDecompose(Matrix m) {
		
		Matrix[] qr = new Matrix[2];
		qr[0] = m.orthonormalize();
		qr[1] = qr[0].transpose().multiply(m);
		
		return qr;
	}
	
	/**
	 * Computes the singular value decomposition (SVD) of this matrix,
	 * writing it as a product of a unitary matrix, a diagonal matrix, and another unitary matrix.
	 * @param m the matrix whose SVD we wish to compute
	 * @return the matrices composing the factorization M={unitary, diagonal, unitary*}
	 */
	public static Matrix[] singularValueDecomposition(Matrix m) {
		
		Matrix[] svd = new Matrix[3];
		ComplexNumber[] singvals = m.singularValues();
		
		ComplexNumber[][] build = new ComplexNumber[singvals.length][singvals.length];
		int sv_pos = 0;
		for (int i = 0; i < build.length; i++) {
			for (int j = 0; j < build[0].length; j++) {
				build[i][j] = (i == j) ? singvals[sv_pos++] : new ComplexNumber(0, 0);
			}
		}
		Matrix diag = new Matrix(build);
		System.out.println(diag);
		System.out.println(SquareMatrixOps.inverse(diag));
		
		return null; // fix this
	}
	
	/**
	 * Computes the Schur decomposition of the matrix
	 * @param m the matrix whose Schur decomposition we are computing
	 * @return an array containing the unitary and triangular matrices
	 * @throws NotSquareException the given matrix is not square
	 */
	public static Matrix[] schurDecompose(Matrix m) throws NotSquareException {
	
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		Matrix[] ut = new Matrix[2];
		Matrix u = null, prev = new Matrix(m.getData());
		
		for (int k = 0; k < m.rows()-1; k++) {
			ComplexNumber[][] build_a = new ComplexNumber[m.rows()-k][m.rows()-k];
			
			for (int i = k; i < m.rows(); i++) {
				for (int j = k; j < m.rows(); j++) {
					build_a[i-k][j-k] = prev.getAt(i, j);
				}
			}
			Matrix temp = new Matrix(build_a);
			ComplexNumber[] eigenval = {(SquareMatrixOps.eigenvalues(temp))[0]};
			System.out.println("eigenvalue: " + eigenval[0]);
			
			Matrix unitary = SquareMatrixOps.eigenvectors(temp, eigenval)[0].generateUnitaryMatrix();
			
			if (k == 0) {
				u = unitary;
			}
			else {
				ComplexNumber[][] build_u = new ComplexNumber[k*2][k*2];
				for (int i = 0; i < k*2; i++) {
					for (int j = 0; j < k*2; j++) {
						if (i >= k && j >= k) {
							build_u[i][j] = unitary.getAt(i-k, j-k);
						}
						else if (i == j) {
							build_u[i][j] = new ComplexNumber(1,0);
						}
						else {
							build_u[i][j] = new ComplexNumber(0,0);
						}
					}
				}
				u = new Matrix(build_u);
			}
			
			prev = (Matrix) u.conjugateTranspose().multiply(prev).multiply(u);
		}
		
		ut[0] = u;
		ut[1] = prev;
		return ut;
	}
	
	/**
	 * Computes and returns the LU decomposition of the matrix
	 * @param m the matrix whose LU decomposition we are computing
	 * @return the result of the LU decomposition {L,U}, or null if no LU decomposition is admitted 
	 */
	public static Matrix[] luDecompose(Matrix m) {
		
		// returns {null, null} if the matrix does not admit an LU-decomposition
		
		Matrix[] lu = new Matrix[2];
		
		// check all leading principal minors are nonzero
		// this encompasses checking that the matrix is invertible too
		for (int i = 1; i <= m.rows(); i++) {
			ComplexNumber[][] build = new ComplexNumber[i][i];
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < i; k++) {
					build[j][k] = m.getAt(j, k);
				}
			}
			if (SquareMatrixOps.determinant(new Matrix(build)).isZero()) {
				lu[0] = lu[1] = null;
				return lu;
			}
		}
		
		Matrix[] lowers = new Matrix[m.rows()-1];
		Matrix temp = m;
		
		for (int n = 0; n < m.rows()-1; n++) {
			ComplexNumber[] ls = new ComplexNumber[m.rows()];
			ComplexNumber[][] build = new ComplexNumber[m.rows()][m.cols()];
			
			for (int i = n+1; i < m.rows(); i++) {
				ls[i] = temp.getAt(i, n).divide(temp.getAt(n, n)).multiply(-1);
			}
			
			for (int a = 0; a < m.rows(); a++) {
				for (int b = 0; b < m.cols(); b++) {
					if (a == b) {
						build[a][b] = new ComplexNumber(1,0);
					}
					else if (b == n && a > n) {
						build[a][b] = ls[a];
					}
					else {
						build[a][b] = new ComplexNumber(0,0);
					}
				}
			}
			
			lowers[n] = new Matrix(build);
			temp = lowers[n].multiply(temp);
		}
		
		lu[1] = temp;
		
		lu[0] = SquareMatrixOps.inverse(lowers[0]);
		for (int i = 1; i < m.rows()-1; i++) {
			lu[0] = lu[0].multiply(SquareMatrixOps.inverse(lowers[i]));
		}
		
		return lu;
	}
	
	/**
	 * Performs a Cholesky decomposition, which writes a 
	 * positive definite Hermetian matrix in terms of a 
	 * lower triangular matrix and its conjugate transpose
	 * @param m the matrix whose Cholesky decomposition we wish to compute
	 * @return the lower triangular matrix L in A=LL* or null if no Cholesky decomposition is admitted
	 */
	public static Matrix choleskyDecompose(Matrix m) {
		
		if (!SquareMatrixOps.isHermetian(m)) {
			return null;
		}
		
		// need to check if m is postive definite too
		
		ComplexNumber[][] L = new ComplexNumber[m.rows()][m.rows()];
		
		// this loop is for debugging
		for (int i = 0; i < m.rows(); i++) {
			for (int j = i+1; j < m.rows(); j++) {
				L[i][j] = new ComplexNumber(0, 0);
			}
		}
		
		// initialize the first column of L
		L[0][0] = m.getAt(0,0).sqrt();
		for (int i = 1; i < m.rows(); i++) {
			L[i][0] = m.getAt(i,0).divide(L[0][0]);
		}
		
		for (int j = 1; j < m.rows(); j++) {
			// the next computation is based on the previously computed column
			// each row in the below matrix is a vector; there are rows-j vectors
			ComplexNumber[][] prevcols = new ComplexNumber[m.rows()-j][j];
			// first j elements in the ith row
			for (int i = j; i < m.rows(); i++) {
				for (int k = 0; k < j; k++) {
					prevcols[i-j][k] = L[i][k];
				}
			}
			
			Vector temp_ip = new Vector(prevcols[0]);
			L[j][j] = m.getAt(j,j).subtract(temp_ip.dot(temp_ip)).sqrt();
			
			if (j != m.rows()-1) {
				Vector temp_arg = new Vector(prevcols[0]);
				for (int i = j+1; i < m.rows(); i++) {
					temp_ip = new Vector(prevcols[i-j]);
					L[i][j] = m.getAt(i,j).subtract(temp_ip.dot(temp_arg)).divide(L[j][j]);
				}
			}
		}
		
		return new Matrix(L);
	}
}
