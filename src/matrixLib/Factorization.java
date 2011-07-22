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
		
		Matrix curr = new Matrix(m.getData()); // iterate on this matrix
		Matrix tri = new Matrix(m.getData()), uni = new Matrix(m.rows());
		boolean real = m.isReal(); // whether the matrix is real or complex
		
		for (int k = 0; k < m.cols()-1; k++) {
			
			Vector column = curr.getVector(0), original;
			int first = 0;
			while (column.isZero()) {
				column = curr.getVector(++first);
			}
			original = curr.getVector(first);
			
			ComplexNumber factor = new ComplexNumber(1,0);
			if (real) {
				// find the pivot element
				int pivot = column.dim()-1;
				while (column.getAt(pivot).isZero()) {
					pivot--;
					if (pivot == 0) {
						break; // 0 vector
					}
				}
				ComplexNumber piv_el = column.getAt(pivot);
				factor = new ComplexNumber(piv_el.Re() < 0 ? -1 : 1, 0);
			}
			else {
				// have to compute exp(i*arg(x)) instead of just +/- 1
				// uses the fact that exp(ix)=cos(x)+isin(x) and the cos(arctan(x)) formula
				ComplexNumber coord = original.getAt(0);
				if (coord.Re() == 0) {
					factor = new ComplexNumber(coord.Im() > 0 ? -1 : 1, 0);
				}
				else {
					double ratio = coord.Im() / coord.Re();
					factor = new ComplexNumber(1.0/Math.sqrt(ratio*ratio+1),ratio/Math.sqrt(ratio*ratio+1)).negative();
					if (coord.Re() < 0) {
						// arg adds/subs pi to the result of atan, negating cos/sin
						factor = factor.negative();
					}
				}
			}
			
			column.set(0, column.getAt(0).add((new ComplexNumber(Norm.pnorm(column,2),0)).multiply(factor)));
			
			// compute the 2-norm squared (avoiding sqrt)
			double normsqr = 0;
			for (int j = 0; j < column.dim(); j++) {
				ComplexNumber coord = column.getAt(j);
				normsqr += coord.Re()*coord.Re()+coord.Im()*coord.Im();
			}
			
			// build the householder matrix, Q = I - 2uu^T/||u||^2
			ComplexNumber[][] hh_arr = new ComplexNumber[curr.rows()][curr.rows()]; // suspect
			// cfactor is used in computing the householder matrix when using complex numbers
			// remember complex inner product does not commute
			ComplexNumber cfactor;
			if (real) {
				cfactor = new ComplexNumber(2, 0);
			}
			else {
				cfactor = column.dot(original).divide(original.dot(column)).add(new ComplexNumber(1,0));
			}
			
			// for debugging:
			ComplexNumber[][] vv = new ComplexNumber[column.dim()][column.dim()];
			for (int i = 0; i < column.dim(); i++) {
				for (int j = 0; j < column.dim(); j++) {
					vv[i][j] = column.getAt(i).multiply(column.getAt(j).conjugate());
				}
			}
			
			for (int i = 0; i < curr.rows(); i++) {
				for (int j = 0; j < curr.rows(); j++) {
					hh_arr[i][j] = column.getAt(i).multiply(column.getAt(j).conjugate()).multiply(cfactor).multiply(1.0/normsqr).negative();
					// account for subtracting the above from the identity matrix
					if (i == j) {
						hh_arr[i][j] = hh_arr[i][j].add(new ComplexNumber(1,0));
					}
				}
			}
			Matrix householder = new Matrix(hh_arr);
			
			if (k > 0) {
				// append 1's in upper left diagonal so we can multiply them all together
				householder = Pattern.blockDiagonal(householder, k);
			}
			
			uni = uni.multiply(householder.conjugateTranspose()); // compute running Q
			tri = householder.multiply(tri); // compute running R
			
			// now operate on the (1,1)-minor of the current matrix
			ComplexNumber[][] minor = new ComplexNumber[curr.rows()-1][curr.cols()-1];
			for (int i = 0; i < minor.length; i++) {
				for (int j = 0; j < minor[0].length; j++) {
					minor[i][j] = tri.getAt(i+1,j+1);
				}
			}
			curr = new Matrix(minor);
		}
		
		Matrix[] qr = {uni, tri};
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
		
		Matrix diag = Pattern.diag(singvals);
		
		System.out.println(diag);
		System.out.println(SquareMatrixOps.inverse(diag));
		
		return null; // fix this
	}
	
	/**
	 * Computes the Schur decomposition of the matrix
	 * @param m the matrix whose Schur decomposition we are computing
	 * @return an array containing the unitary matrix that can transform it to triangular form and the triangular form
	 * @throws NotSquareException the given matrix is not square
	 */
	public static Matrix[] schurDecompose(Matrix m) throws NotSquareException {
	
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		Matrix curr = new Matrix(m.getData()), minor = new Matrix(m.getData()); // start with the given matrix
		Matrix unitary_prod = new Matrix(m.rows());
		
		for (int k = 0; k < m.rows()-1; k++) {
			
			ComplexNumber[][] newmat = new ComplexNumber[m.rows()-k][m.rows()-k];
			for (int i = 0; i < newmat.length; i++) {
				for (int j = 0; j < newmat[0].length; j++) {
					newmat[i][j] = curr.getAt(i+k,j+k);
				}
			}
			minor = new Matrix(newmat);
			
			ComplexNumber[] eval = {SquareMatrixOps.eigenvalues(minor)[0]};
			Vector evec = SquareMatrixOps.eigenvectors(minor, eval)[0];
			
			Matrix unitary = evec.generateUnitaryMatrix();
			Matrix uni_ref = (k == 0) ? unitary : Pattern.blockDiagonal(unitary, k);
			
			unitary_prod = unitary_prod.multiply(uni_ref); // running unitary matrix
			curr = uni_ref.conjugateTranspose().multiply(curr).multiply(uni_ref);
		}
		
		Matrix[] ut = {unitary_prod, curr};
		return ut;
	}
	
	/**
	 * Computes and returns the LU decomposition of the matrix
	 * @param m the matrix whose LU decomposition we are computing
	 * @return the result of the LU decomposition {L,U}, or null if no LU decomposition is admitted
	 * @throws NotSquareException the supplied matrix is not square 
	 */
	public static Matrix[] luDecompose(Matrix m) throws NotSquareException {
		
		// returns {null, null} if the matrix does not admit an LU-decomposition
		
		if (m.rows() != m.cols()) {
			throw new NotSquareException();
		}
		
		ComplexNumber[][] L = new ComplexNumber[m.rows()][m.rows()];
		ComplexNumber[][] U = new ComplexNumber[m.rows()][m.rows()];
		
		// initialize the L and U matrices
		if (m.getAt(0,0).isZero()) {
			return null; // no factorization if (0,0) is 0
		}
		
		// initialize the first row of U
		for (int j = 0; j < m.cols(); j++) {
			U[0][j] = m.getAt(0,j).divide(m.getAt(0,0));
		}
		// initialize the first column of L
		for (int j = 0; j < m.rows(); j++) {
			L[j][0] = m.getAt(j,0);
		}
		// initialize the remainder of L and U to zero
		for (int i = 0; i < m.rows(); i++) {
			for (int j = 1; j < m.cols(); j++) {
				L[i][j] = new ComplexNumber(0,0);
				U[j][i] = new ComplexNumber(0,0);
			}
		}
		
		for (int n = 1; n < m.rows(); n++) {
			
			// next computations are based on the previously determined rows
			// there are rows-n vectors of n coordinates to consider
			ComplexNumber[][] l_prev = new ComplexNumber[m.rows()-n][n];
			ComplexNumber[][] u_prev = new ComplexNumber[m.rows()-n][n];
			
			for (int i = n; i < m.rows(); i++) {
				for (int j = 0; j < n; j++) {
					l_prev[i-n][j] = L[i][j];
					u_prev[i-n][j] = U[j][i];
				}
			}
			
			// compute the nth column of L
			for (int i = n; i < m.rows(); i++) {
				ComplexNumber sum = new ComplexNumber(0,0);
				for (int k = 0; k < n; k++) {
					sum = sum.add(l_prev[i-n][k].multiply(u_prev[0][k]));
				}
				L[i][n] = m.getAt(i,n).subtract(sum);
			}
			if (L[n][n].isZero() && n != m.rows()-1) {
				// the factorization is not possible
				return null;
			}
			
			U[n][n] = new ComplexNumber(1,0);
			if (n != m.rows() - 1) {
				// compute the nth row of U, right of the diagonal
				for (int j = n+1; j < m.cols(); j++) {
					ComplexNumber sum = new ComplexNumber(0,0);
					for (int k = 0; k < n; k++) {
						sum = sum.add(l_prev[0][k].multiply(u_prev[j-n][k]));
					}
					U[n][j] = m.getAt(n,j).subtract(sum).divide(L[n][n]);
				}
			}
		}
		
		// return the matrices {L, U}
		Matrix[] lu = new Matrix[2];
		lu[0] = new Matrix(L);
		lu[1] = new Matrix(U);
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
		
		if (!Pattern.isHermetian(m)) {
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
