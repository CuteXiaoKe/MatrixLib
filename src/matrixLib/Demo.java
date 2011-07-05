package matrixLib;

/**
 * Library of tests/demos for the other libraries; uses JUnit 4
 * @author Bryan Cuccioli
 */

import static org.junit.Assert.*;
import org.junit.Test;

public class Demo {
	
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
