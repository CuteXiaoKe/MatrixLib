package matrixLib;

/**
 * Library of tests/demos for the other libraries; uses JUnit 4
 * @author Bryan Cuccioli
 */

import static org.junit.Assert.*;
import org.junit.Test;

public class Demo {
	
	@Test public void test_det() {
		// set up matrices for a testing a variety of determinants
		double[][] f = {{1,2,3},{4,5,6},{7,8,7}};
		double[][] g = {{1,2},{3,4}};
		double[][] h = {{3,0,6,-3},{0,2,3,0},{-4,-7,2,0},{2,0,1,10}};
		ComplexNumber[][] i = {{new ComplexNumber(2,1),new ComplexNumber(3,-1),new ComplexNumber(4,-3)},
				{new ComplexNumber(4,0),new ComplexNumber(6,-1),new ComplexNumber(2,5)},
				{new ComplexNumber(0,3),new ComplexNumber(2,-1),new ComplexNumber(1,3)}};
		Matrix m = new Matrix(f);
		Matrix t = new Matrix(g);
		Matrix u = new Matrix(h);
		Matrix z = new Matrix(i);
		
		System.out.println(SquareMatrixOps.determinant(m));
		
		assertTrue(SquareMatrixOps.determinant(m).equals(new ComplexNumber(6, 0)));
		assertTrue(SquareMatrixOps.determinant(t).equals(new ComplexNumber(-2, 0)));
		assertTrue(SquareMatrixOps.determinant(new Matrix(9)).equals(new ComplexNumber(1, 0)));
		assertTrue(SquareMatrixOps.determinant(u).equals(new ComplexNumber(1404, 0)));
		assertTrue(SquareMatrixOps.determinant(z).equals(new ComplexNumber(-118,-84)));
	}
	
	public static void test_eigen() {
		double[][] f = {{1,2},{2,1}};
		Matrix m = new Matrix(f);
		ComplexNumber[] evals = SquareMatrixOps.eigenvalues(m);
		for (ComplexNumber val : evals) {
			System.out.println(val);
		}
		Vector[] eigenvecs = SquareMatrixOps.eigenvectors(m, evals);
		for (Vector ev: eigenvecs) {
			System.out.println(ev);
		}
	}
	
	@Test public void test_proj() {
		ComplexNumber[] a = {new ComplexNumber(3,1), new ComplexNumber(3, -1)};
		ComplexNumber[] b = {new ComplexNumber(2,1), new ComplexNumber(2, -1)};
		Vector v1 = new Vector(a);
		Vector v2 = new Vector(b);

		assertTrue(v1.dot(v2) == new ComplexNumber(14,0));
		//assertTrue(v2.dot(v2))
		
		System.out.println("dot: " + v1.dot(v2));
		System.out.println("dot: " + v2.dot(v2));
		
		// projection of v1 onto v2
		System.out.println(v1.proj(v2));
		ComplexNumber[] expected={new ComplexNumber(2.8,1.4),new ComplexNumber(2.8,-1.4)};

		assertTrue(v1.proj(v2).equals(new Vector(expected)));
	}
	
	@Test public void test_dot() {
		// dot product over complex numbers multiplies
		ComplexNumber[] a = {new ComplexNumber(1,1),new ComplexNumber(2,1)};
		ComplexNumber[] b = {new ComplexNumber(3,-1), new ComplexNumber(4,1)};
		
		assertTrue((new Vector(a)).dot(new Vector(b)).equals(new ComplexNumber(14, -6)));
		
		System.out.println((new Vector(a)).dot(new Vector(b)));
	}
	
	public static void test_qr() {
		double[][] c = {{1,2},{3,4}};
		Matrix m = new Matrix(c);
		
		//Matrix[] qr1 = m.QRDecompose();
		//System.out.println(qr1[0]);
		//System.out.println(qr1[1]);
		
		//ComplexNumber[][] z = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		ComplexNumber[][] z = {{new ComplexNumber(0,1), new ComplexNumber(2, 0)},{new ComplexNumber(1,0),new ComplexNumber(1,1)}};
		Matrix zm = new Matrix(z);
		
		Matrix[] qr = Factorization.QRDecompose(zm);
		System.out.println(qr[0]);
		System.out.println(qr[1]);
	}
	
	public static void test_cholesky() {
		
		double[][] f = {{2,-1,0},{-1,2,-1},{0,-1,2}};
		Matrix m = new Matrix(f);
		
		Matrix L = Factorization.choleskyDecompose(m);
		System.out.println(L);
		System.out.println(L.multiply(L.conjugateTranspose()));
	}
	
	public static void test_lu() {
		double[][] a = {{1,2},{3,4}};
		Matrix m = new Matrix(a);
		
		Matrix[] lu = Factorization.luDecompose(m);
		
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
		System.out.println(SquareMatrixOps.inverse(new Matrix(f)));
		//double[][] f = {{2,(double) 2.09999},{(double) 2.09999,2}};
		//double[][] data = {{1,2},{2,1}};
		//Matrix m = new Matrix(f);
		//System.out.println(m.inverse());
		
		/*ComplexNumber ev = new ComplexNumber(-1,0);
		ev = ev.multiply(1.1);
		
		Matrix diag = new Matrix(data);
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
		Matrix m1 = new Matrix(f);
		Matrix m2 = new Matrix(f);
		Matrix res = m1.multiply(m2);
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
		System.out.println(Factorization.singularValueDecomposition(m));
	}
	
	public static void test_schur() {
		double[][] f = {{4,0,1},{1,3,-1},{-1,0,2}};
		Matrix m = new Matrix(f);
		Matrix[] schur = Factorization.schurDecompose(m);
		System.out.println(schur[0]);
		System.out.println(schur[1]);
	}
	
	public static void test_gs() {
		double[][] f = {{.7000, .70711},{.70001, .70711}};
		Matrix m = new Matrix(f);
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
		test_schur();
		
		//test_svd();
		
		//test_bases();
		
		//test_multiply();
		
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
