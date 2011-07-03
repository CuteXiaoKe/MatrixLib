package matrixLib;

public class Demo {

	public static void test_eigen() {
		double[][] f = {{1,2},{2,1}};
		SquareMatrix m = new SquareMatrix(f);
		ComplexNumber[] evals = m.eigenvalues();
		for (ComplexNumber val : evals) {
			System.out.println(val);
		}
		Vector[] eigenvecs = m.eigenvectors(evals);
		for (Vector ev: eigenvecs) {
			System.out.println(ev);
		}
	}
	
	public static void test_proj() {
		ComplexNumber[] a = {new ComplexNumber(2,1), new ComplexNumber(3, -1)};
		ComplexNumber[] b = {new ComplexNumber(1,1), new ComplexNumber(4, -1)};
		Vector v1 = new Vector(a);
		Vector v2 = new Vector(b);

		System.out.println("dot: " + v1.dot(v2));
		System.out.println("dot: " + v2.dot(v2));
		
		// projection of v1 onto v2
		System.out.println(v1.proj(v2));
	}
	
	public static void test_dot() {
		
		ComplexNumber[] a = {new ComplexNumber(1,1),new ComplexNumber(2,1)};
		ComplexNumber[] b = {new ComplexNumber(3,-1), new ComplexNumber(4,1)};
		
		System.out.println((new Vector(a)).dot(new Vector(b)));
	}
	
	public static void test_qr() {
		double[][] c = {{1,2},{3,4}};
		SquareMatrix m = new SquareMatrix(c);
		
		//Matrix[] qr1 = m.QRDecompose();
		//System.out.println(qr1[0]);
		//System.out.println(qr1[1]);
		
		//ComplexNumber[][] z = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		ComplexNumber[][] z = {{new ComplexNumber(0,1), new ComplexNumber(2, 0)},{new ComplexNumber(1,0),new ComplexNumber(1,1)}};
		SquareMatrix zm = new SquareMatrix(z);
		
		Matrix[] qr = zm.QRDecompose();
		System.out.println(qr[0]);
		System.out.println(qr[1]);
	}
	
	public static void test_cholesky() {
		
		double[][] f = {{2,-1,0},{-1,2,-1},{0,-1,2}};
		SquareMatrix m = new SquareMatrix(f);
		
		SquareMatrix L = m.choleskyDecompose();
		System.out.println(L);
		System.out.println(L.multiply(L.conjugateTranspose()));
	}
	
	public static void test_det() {
		double[][] f = {{1,2,3},{4,5,6},{7,8,7}};
		double[][] g = {{1,2},{3,4}};
		double[][] h = {{3,0,6,-3},{0,2,3,0},{-4,-7,2,0},{2,0,1,10}};
		SquareMatrix m = new SquareMatrix(f);
		SquareMatrix t = new SquareMatrix(g);
		SquareMatrix u = new SquareMatrix(h);
		
		System.out.println("det(f) = " + m.determinant());
		System.out.println("det(g) = " + t.determinant());
		System.out.println("det(I_9) = " + (new SquareMatrix(9)).determinant());
		System.out.println("det(h) = " + u.determinant());
	}
	
	public static void test_lu() {
		double[][] a = {{1,2},{3,4}};
		SquareMatrix m = new SquareMatrix(a);
		
		Matrix[] lu = m.luDecompose();
		
		System.out.println(lu[0]);
		System.out.println(lu[1]);
	}
	
	public static void test_norm() {
		ComplexNumber[] a = {new ComplexNumber(1,0),new ComplexNumber(0,1)};
		Vector v = new Vector(a);
		System.out.println(v.normalize());
	}
	
	public static void test_rref() {
		double[][] f = {{1,2,3},{4,5,6},{7,8,9}};
		//double[][] f = {{1,2},{2,4}};
		Matrix m = new Matrix(f);
		System.out.println(m.rref());
		//System.out.println("rank: " + m.rank());
		//System.out.println("nullity: " + m.nullity());
	}
	
	public static void test_inverse() {
		double[][] f = {{5,19},{1,4}};
		System.out.println((new SquareMatrix(f)).inverse());
		//double[][] f = {{2,(double) 2.09999},{(double) 2.09999,2}};
		//double[][] data = {{1,2},{2,1}};
		//SquareMatrix m = new SquareMatrix(f);
		//System.out.println(m.inverse());
		
		/*ComplexNumber ev = new ComplexNumber(-1,0);
		ev = ev.multiply(1.1);
		
		SquareMatrix diag = new SquareMatrix(data);
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				if (i == j) {
					diag.set(i, j, diag.getAt(i, j).subtract(ev));
				}
			}
		}
		System.out.println(diag);
		diag = diag.inverse();
		System.out.println(diag);*/
	}
	
	public static void test_multiply() {
		//ComplexNumber[][] z = {{new ComplexNumber(1,1),new ComplexNumber(2,0)},{new ComplexNumber(1,0),new ComplexNumber(0,1)}};
		//Matrix m = new Matrix(z);
		//double[][] f = {{1,2,3},{4,5,6}};
		//double[] g = {7, 8, 9};
		double[][] f= {{1,2},{3,4}};
		SquareMatrix m1 = new SquareMatrix(f);
		SquareMatrix m2 = new SquareMatrix(f);
		SquareMatrix res = m1.multiply(m2);
		System.out.println("blah");
		System.out.println(res);
	}
	
	public static void test_bases() {
		double[][] f = {{0,0,0},{1,2,3},{2,4,7}};
		Matrix m = new Matrix(f);
		System.out.println("Basis for image: " + m.imageBasis());
	}
	
	public static void test_svd() {
		double[][] f = {{4,0},{3,-5}};
		Matrix m = new Matrix(f);
		System.out.println(m.singularValueDecomposition());
	}
	
	public static void test_schur() {
		double[][] f = {{4,0,1},{1,3,-1},{-1,0,2}};
		SquareMatrix m = new SquareMatrix(f);
		SquareMatrix[] schur = m.schurDecompose();
		System.out.println(schur[0]);
		System.out.println(schur[1]);
	}
	
	public static void test_gs() {
		double[][] f = {{.7000, .70711},{.70001, .70711}};
		SquareMatrix m = new SquareMatrix(f);
		System.out.println(m.orthonormalize());
	}
	
	public static void test_reflector() {
		double[] f = {1,2};
		Vector v = new Vector(f);
		System.out.println(v.reflector());
	}
	
	public static void main(String[] args) {
		
		//test_reflector();
		
		//test_gs();
		//test_schur();
		
		//test_svd();
		
		//test_bases();
		
		test_multiply();
		
		//test_det();
		//test_rref();
		//test_lu();
		//test_eigen();
		//test_inverse();
		//test_cholesky();
		//test_qr();
		//test_proj();
		//test_dot();
		
	}
}
