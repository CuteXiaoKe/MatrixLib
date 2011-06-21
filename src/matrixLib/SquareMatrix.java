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
	 * @return the corresponding inverse matrix, or null if the matrix is not invertible
	 * @throws DimensionMismatchException 
	 * @throws NotSquareException 
	 */
	public SquareMatrix inverse() throws DimensionMismatchException, NotSquareException {
		
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
		
		// if the left wasn't reduced to the identity, inverse doesn't exist
		ComplexNumber[][] id = new ComplexNumber[rows()][cols()];
		for (int i = 0; i < rows(); i++) {
			for (int j = 0; j < cols(); j++) {
				id[i][j] = aug_rref.getAt(i, j);
			}
		}
		
		if (!(new SquareMatrix(id)).isIdentity()) {
			return null;
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
	 * @throws DimensionMismatchException 
	 */
	public boolean isOrthogonal() throws DimensionMismatchException {
		
		return ((SquareMatrix)this.multiply(this.transpose())).isIdentity();
	}
	
	/**
	 * Tells whether the matrix is unitary (its inverse is its conjugate transpose)
	 * @return whether the matrix is unitary or not
	 * @throws DimensionMismatchException
	 */
	public boolean isUnitary() throws DimensionMismatchException {
		
		return ((SquareMatrix)this.multiply(this.conjugateTranspose())).isIdentity();
	}
	
	/**
	 * Computes and returns the LU decomposition of the matrix
	 * @return the result of the LU decomposition {L,D,U}, or null if no LU decomposition is admitted
	 * @throws DimensionMismatchException 
	 * @throws SingularMatrixException 
	 * @throws NotSquareException 
	 */
	public Matrix[] luDecompose() throws DimensionMismatchException, SingularMatrixException, NotSquareException {
		
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
	 * @throws DimensionMismatchException 
	 */
	public ComplexNumber[] eigenvalues() throws DimensionMismatchException {
		
		// max amount of eigenvalues is the dimension of the matrix
		ComplexNumber[] evals = new ComplexNumber[rows()];
		
		Matrix temp = this;
		while (!isAlmostUpperTriangular(temp)) {
			System.out.println("QR decomposing...");
			Matrix[] qr = temp.QRDecompose();
			temp = qr[1].multiply(qr[0]);
			System.out.println(temp);
		}
		
		// temp is now upper triangular, so the evals are on the diagonal
		for (int i = 0; i < rows(); i++) {
			evals[i] = temp.getAt(i, i);
		}
		
		return evals;
	}
	
	/**
	 * Performs a Cholesky decomposition, which writes a 
	 * positive definite Hermetian matrix in terms of a 
	 * lower triangular matrix and its conjugate transpose
	 * @return the lower triangular matrix L in A=LL* or null if no Cholesky decomposition is admitted
	 * @throws NotSquareException 
	 */
	public SquareMatrix choleskyDecompose() throws NotSquareException {
		
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