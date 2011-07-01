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
	 * @throws NotSquareException the array to build the matrix from is not square
	 */
	public SquareMatrix(ComplexNumber[][] mat) throws NotSquareException {
		
		super(mat);
		
		if (mat.length != mat[0].length) {
			 new NotSquareException();
		}
	}

	/**
	 * Construct the matrix with specified underlying data
	 * @param mat the underlying data of the matrix
	 * @throws NotSquareException the array to build the matrix from is not square
	 */
	public SquareMatrix(double[][] mat) throws NotSquareException {
		
		super(mat);
		
		if (mat.length != mat[0].length) {
			 new NotSquareException();
		}
	}
	
	/**
	 * Construct the square matrix made up of the specified vectors
	 * @param cols the vectors that make up the matrix
	 * @throws DimensionMismatchException the vectors do not form a square matrix
	 */
	public SquareMatrix(Vector[] cols) throws DimensionMismatchException {
		
		super(cols);
		
		if (cols.length != cols[0].dim()) {
			 new DimensionMismatchException();
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
	 * Computes the inverse of the matrix via Gauss-Jordan elimination
	 * @return the corresponding inverse matrix, or null if the matrix is not invertible
	 */
	public SquareMatrix inverse() {
		
		ComplexNumber[][] augmented = new ComplexNumber[rows()][cols()*2];
		
		for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols()*2; j++) {
				if (j < cols()) {
					augmented[i][j] = getAt(i, j);
				}
				else {
					if (i == j - cols()) {
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
		for (int i = 0; i < rows(); i++) {
			boolean all_zero = true;
			for (int j = 0; j < cols(); j++) {
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
		
		ComplexNumber[][] inv = new ComplexNumber[rows()][cols()];
		
		for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols(); j++) {
				inv[i][j] = aug_rref.getAt(i, j+cols());
			}
		}
		
		return new SquareMatrix(inv);
	}
	
	/**
	 * Returns the determinant of this matrix in O(n^3) time complexity
	 * @return the determinant of the matrix
	 */
	public ComplexNumber determinant() {

		// try shortcuts for common sizes first
		if (rows() == 1) {
			return getAt(0,0);
		}
		else if (rows() == 2) {
			return getAt(0,0).multiply(getAt(1,1)).subtract(getAt(1,0).multiply(getAt(0,1)));
		}
		else {
			return rref(1).getAt(0, 0);
		}
	}

	/**
	 * Computes the tridiagonalized form of the matrix via Householder reflections
	 * @return the tridiagonalized form of the matrix
	 * @s NotSquareException 
	 */
	public SquareMatrix tridiagonalized() s NotSquareException {
		
		SquareMatrix hess = new SquareMatrix(getData());
		
		for (int k = 0; k < rows(); k++) {
			ComplexNumber sum = new ComplexNumber(0,0);
			for (int j = k+1; j < rows(); j++) {
				sum = sum.add(hess.getAt(j, k).multiply(hess.getAt(j, k)));
			}
			ComplexNumber alpha = sum.sqrt().multiply(-Math.signum(hess.getAt(k+1,k).Re()));
			ComplexNumber r = alpha.multiply(alpha).subtract(hess.getAt(k+1,k).multiply(alpha)).multiply(.5).sqrt().reciprocal();
			
			ComplexNumber[] vec = new ComplexNumber[rows()];
			for (int j = 0; j <= k; j++) {
				vec[j] = new ComplexNumber(0,0);
			}
			vec[k+1] = hess.getAt(k+1,k).subtract(alpha).multiply(0.5*r.Re());
			for (int j = k+2; j < rows(); j++) {
				vec[j] = hess.getAt(j,k).multiply(0.5*r.Re());
			}
			
			ComplexNumber[][] p = new ComplexNumber[rows()][rows()];
			for (int i = 0; i < rows(); i++) {
				for (int j = 0; j < rows(); j++) {
					p[i][j] = (i==j) ? new ComplexNumber(1,0) : new ComplexNumber(0,0);
					p[i][j] = p[i][j].subtract(vec[i].multiply(vec[j]).multiply(2.0));
				}
			}
			
			SquareMatrix hh = new SquareMatrix(p);
			hess = hh.multiply(hess).multiply(hh);
		}
		
		return hess;
	}
	
	/**
	 * Returns the upper Hessenberg form of this matrix,
	 * which has zero entries below the first subdiagonal
	 * @return the upper Hessenberg form of this matrix
	 * @s NotSquareException 
	 */
	public SquareMatrix hessenbergForm() s NotSquareException {
		
		SquareMatrix prev = new SquareMatrix(getData());
		
		for (int k = 0; k < rows()-2; k++) {
			// x is the kth column, restricted to below the diagonal
			ComplexNumber[] vec = new ComplexNumber[rows()-k-1];
			for (int i = k+1; i < rows(); i++) {
				vec[i-k-1] = prev.getAt(i, k);
			}
			Vector x = new Vector(vec);
			x.set(0, x.getAt(0).add(new ComplexNumber(x.length(),0)));
			
			SquareMatrix ref = x.reflector();
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
			SquareMatrix block = new SquareMatrix(p);
			prev = block.transpose().multiply(prev).multiply(block);
		}
		
		return prev;
	}
	
	/**
	 * Tells whether the matrix is symmetric or not
	 * @return whether the matrix is symmetric
	 */
	public boolean isSymmetric() {
		
		// a symmetric matrix equals its transpose
		return (this.equals(this.transpose()));
		
		/*for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < rows(); j++) {
				
				if
			}
		}*/
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
	 * Tells whether the matrix is Hermetian/self-adjoint (equal to its own conjugate transpose)
	 * @return whether the matrix is Hermetian
	 */
	public boolean isHermetian() {
		
		return this.conjugateTranspose().equals(this);
	}

	/**
	 * Tells whether the matrix is orthogonal (A^T A = I) or not
	 * @return whether the matrix is orthogonal 
	 */
	public boolean isOrthogonal() {
		
		return ((SquareMatrix)this.multiply(this.transpose())).isIdentity();
	}
	
	/**
	 * Tells whether the matrix is unitary (its inverse is its conjugate transpose)
	 * @return whether the matrix is unitary or not
	 * @s DimensionMismatchException
	 */
	public boolean isUnitary() {
		
		return ((SquareMatrix)this.multiply(this.conjugateTranspose())).isIdentity();
	}
	
	/**
	 * Tells whether the matrix is upper Hessenberg
	 * (whether it has zeros below the first subdiagonal)
	 * @return whether the matrix is upper Hessenberg
	 */
	public boolean isUpperHessenberg() {
		
		for (int i = 2; i < rows(); i++) {
			for (int j = 0; j < i-1; j++) {
				if (!getAt(i,j).isZero()) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	/**
	 * Computes the Schur decomposition of the matrix
	 * @return an array containing the unitary and triangular matrices
	 */
	public SquareMatrix[] schurDecompose() {
		
		SquareMatrix[] ut = new SquareMatrix[2];
		SquareMatrix u = null, prev = new SquareMatrix(getData());
		
		for (int k = 0; k < rows()-1; k++) {
			ComplexNumber[][] build_a = new ComplexNumber[rows()-k][rows()-k];
			
			for (int i = k; i < rows(); i++) {
				for (int j = k; j < rows(); j++) {
					build_a[i-k][j-k] = prev.getAt(i, j);
				}
			}
			SquareMatrix temp = new SquareMatrix(build_a);
			ComplexNumber[] eigenval = {(temp.eigenvalues())[0]};
			System.out.println("eigenvalue: " + eigenval[0]);
			
			SquareMatrix unitary = temp.eigenvectors(eigenval)[0].generateUnitaryMatrix();
			
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
				u = new SquareMatrix(build_u);
			}
			
			prev = (SquareMatrix) u.conjugateTranspose().multiply(prev).multiply(u);
		}
		
		ut[0] = u;
		ut[1] = prev;
		return ut;
	}
	
	/**
	 * Computes and returns the LU decomposition of the matrix
	 * @return the result of the LU decomposition {L,D,U}, or null if no LU decomposition is admitted 
	 */
	public Matrix[] luDecompose() {
		
		// returns {null, null} if the matrix does not admit an LU-decomposition
		
		Matrix[] lu = new Matrix[2];
		
		// check all leading principal minors are nonzero
		// this encompasses checking that the matrix is invertible too
		for (int i = 1; i <= rows(); i++) {
			ComplexNumber[][] build = new ComplexNumber[i][i];
			for (int j = 0; j < i; j++) {
				for (int k = 0; k < i; k++) {
					build[j][k] = getAt(j, k);
				}
			}
			if ((new SquareMatrix(build)).determinant().isZero()) {
				lu[0] = lu[1] = null;
				return lu;
			}
		}
		
		SquareMatrix[] lowers = new SquareMatrix[rows()-1];
		Matrix temp = this;
		
		for (int n = 0; n < rows()-1; n++) {
			ComplexNumber[] ls = new ComplexNumber[rows()];
			ComplexNumber[][] build = new ComplexNumber[rows()][cols()];
			
			for (int i = n+1; i < rows(); i++) {
				ls[i] = temp.getAt(i, n).divide(temp.getAt(n, n)).multiply(-1);
			}
			
			for (int a = 0; a < rows(); a++) {
				for (int b = 0; b < cols(); b++) {
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
			
			lowers[n] = new SquareMatrix(build);
			temp = lowers[n].multiply(temp);
		}
		
		lu[1] = temp;
		
		lu[0] = lowers[0].inverse();
		for (int i = 1; i < rows()-1; i++) {
			lu[0] = lu[0].multiply(lowers[i].inverse());
		}
		
		return lu;
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

	private boolean isAlmostUpperTriangular(Matrix m) {
		
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
	 * @return an array containing the eigenvalues of the matrix 
	 */
	public ComplexNumber[] eigenvalues() {
		
		// max amount of eigenvalues is the dimension of the matrix
		ComplexNumber[] evals = new ComplexNumber[rows()];
		
		Matrix temp = this;
		while (!isAlmostUpperTriangular(temp)) {
			//System.out.println("QR decomposing...");
			Matrix[] qr = temp.QRDecompose();
			temp = qr[1].multiply(qr[0]);
			//System.out.println(temp);
		}
		
		// temp is now upper triangular, so the evals are on the diagonal
		for (int i = 0; i < rows(); i++) {
			evals[i] = temp.getAt(i, i);
		}
		
		return evals;
	}

	/**
	 * Computes the normalized eigenvectors
	 * @return an array of normalized eigenvectors for the matrix
	 */
	public Vector[] eigenvectors() {
		
		return eigenvectors(eigenvalues());
	}
	
	/**
	 * Computes the normalized eigenvectors of the matrix corresponding to certain eigenvalues
	 * @param evals The list of eigenvalues to compute the associated eigenvectors of
	 * @return an array of normalized eigenvectors for the matrix
	 * @s DimensionMismatchException 
	 * @s NotSquareException 
	 */
	public Vector[] eigenvectors(ComplexNumber[] evals) {
		
		Vector[] evecs = new Vector[rows()];
		int vec_pos = 0;
		
		ComplexNumber[] test_vec = new ComplexNumber[rows()];
		for (int i = 0; i < rows(); i++) {
			test_vec[i] = new ComplexNumber(1, 0);
		}
		
		for (ComplexNumber ev : evals) {
			
			SquareMatrix diag = new SquareMatrix(getData());
			for (int i = 0; i < rows(); i++) {
				diag.set(i, i, diag.getAt(i, i).subtract(ev));
			}
			diag = diag.inverse();
			
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
	 * Performs a Cholesky decomposition, which writes a 
	 * positive definite Hermetian matrix in terms of a 
	 * lower triangular matrix and its conjugate transpose
	 * @return the lower triangular matrix L in A=LL* or null if no Cholesky decomposition is admitted
	 */
	public SquareMatrix choleskyDecompose() {
		
		if (!isHermetian()) {
			return null;
		}
		
		ComplexNumber[][] L = new ComplexNumber[rows()][rows()];
		
		for (int i = 0; i < cols(); i++) {
			for (int j = 0; j < rows(); j++) {
				if (j < i) {
					L[j][i] = new ComplexNumber(0, 0);
				}
				else {
					ComplexNumber entry = getAt(j, i);
					for (int k = 0; k < i-1; k++) {
						//System.out.printf("(%d, %d)\n", j, k);
						//System.out.printf("(%d, %d)\n", i, k);
						entry = entry.subtract(L[j][k].multiply(L[i][k].conjugate()));
					}
					if (i == j) {
						entry = entry.sqrt();
					}
					else {
						//System.out.printf("(%d, %d)\n", i, i);
						entry = entry.divide(L[i][i]);
					}
					//System.out.printf("adding (%d, %d)\n", j, i);
					L[j][i] = new ComplexNumber(entry.Re(), entry.Im());
				}
			}
		}
		
		return new SquareMatrix(L);
	}
}