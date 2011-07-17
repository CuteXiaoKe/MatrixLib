package matrixLib.test;

import matrixLib.ComplexNumber;
import matrixLib.Factorization;
import matrixLib.Matrix;
import matrixLib.SquareMatrixOps;
import matrixLib.Pattern;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Library of routines for testing the Factorization library
 * @author Bryan Cuccioli
 */

public class FactorizationTest {

	@Test public void QRDecompose() {
		double[][] m1 = {{12,-51,4},{6,167,-68},{-4,24,-41}};
		Factorization.QRDecompose(new Matrix(m1));
		
		/*double[][] m1 = {{1.92,2.32},{.56,2.76}};
		double[][] q1 = {{.96, -.28},{.28, .96}};
		double[][] r1 = {{2,3},{0,2}};
		
		Matrix[] qr1 = Factorization.QRDecompose(new Matrix(m1));
		assertTrue(qr1[0].equals(new Matrix(q1)));
		assertTrue(qr1[1].equals(new Matrix(r1)));*/
		
		//double[][] q2 = {{1,2},{3,4}};
		//System.out.println(Factorization.QRDecompose(new Matrix(q2))[1]);
		
		//Matrix[] qr1 = m.QRDecompose();
		//System.out.println(qr1[0]);
		//System.out.println(qr1[1]);
		
		//ComplexNumber[][] z = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		/*ComplexNumber[][] z = {{new ComplexNumber(0,1), new ComplexNumber(2, 0)},{new ComplexNumber(1,0),new ComplexNumber(1,1)}};
		Matrix zm = new Matrix(z);
		
		Matrix[] qr = Factorization.QRDecompose(zm);
		System.out.println(qr[0]);
		System.out.println(qr[1]);*/
	}
	
	@Test public void luDecompose() {
		double[][] mat1 = {{1,2},{3,4}};
		double[][] mat2 = {{1,2,-2,3},{-1,1,0,2},{3,-3,4,1},{2,1,1,-2}};
		ComplexNumber[][] mat3 = {{new ComplexNumber(0,1),new ComplexNumber(2,1)},{new ComplexNumber(3,-1),new ComplexNumber(5,3)}};
		double[][] mat4 = {{1,2,-2,3},{-1,-2,0,2},{3,-3,0,1},{2,1,1,-2}};
		
		Matrix[] lu = Factorization.luDecompose(new Matrix(mat1));
		assertTrue(Pattern.isLowerTriangular(lu[0]) && Pattern.isUpperTriangular(lu[1])
				&& lu[0].multiply(lu[1]).equals(new Matrix(mat1)));
		
		lu = Factorization.luDecompose(new Matrix(mat2));
		assertTrue(Pattern.isLowerTriangular(lu[0]) && Pattern.isUpperTriangular(lu[1])
				&& lu[0].multiply(lu[1]).equals(new Matrix(mat2)));
		
		lu = Factorization.luDecompose(new Matrix(mat3));
		assertTrue(Pattern.isLowerTriangular(lu[0]) && Pattern.isUpperTriangular(lu[1])
				&& lu[0].multiply(lu[1]).equals(new Matrix(mat3)));
		
		lu = Factorization.luDecompose(new Matrix(mat4));
		assertTrue(lu == null); // this one does not admit an LU factorization
	}
	
	@Test public void choleskyDecompose() {

		double[][] diag = {{1,0,0,0},{0,25,0,0},{0,0,64,0},{0,0,0,100}};
		double[][] d_exp = {{1,0,0,0},{0,5,0,0},{0,0,8,0},{0,0,0,10}};;
		double[][] m = {{2,-1,0},{-1,2,-1},{0,-1,2}};
		double[][] m_exp = {{Math.sqrt(2),0,0},{-.5*Math.sqrt(2),.5*Math.sqrt(6),0},{0,-Math.sqrt(6)/3.0,2.0/Math.sqrt(3)}};
		//double[][] big = {{16,-3,5,-8},{-3,16,-5,-8},{5,-5,24,0},{-8,-8,0,21}};
		ComplexNumber[][] c = {{new ComplexNumber(4,0),new ComplexNumber(0,2),new ComplexNumber(0,-1)},{new ComplexNumber(0,-2),new ComplexNumber(10,0),new ComplexNumber(1,0)},{new ComplexNumber(0,1),new ComplexNumber(1,0),new ComplexNumber(9,0)}};
		ComplexNumber[][] c_exp = {{new ComplexNumber(2,0),new ComplexNumber(0,0),new ComplexNumber(0,0)},{new ComplexNumber(0,-1),new ComplexNumber(3,0),new ComplexNumber(0,0)},{new ComplexNumber(0,.5),new ComplexNumber(.5,0),new ComplexNumber(Math.sqrt(8.5),0)}};
		
		assertTrue(Factorization.choleskyDecompose(new Matrix(diag)).equals(new Matrix(d_exp)));
		assertTrue(Factorization.choleskyDecompose(new Matrix(m)).equals(new Matrix(m_exp)));
		assertTrue(Factorization.choleskyDecompose(new Matrix(c)).equals(new Matrix(c_exp)));
	}
	
	@Test public void test_svd() {
		/*double[][] f = {{4,0},{3,-5}};
		Matrix m = new Matrix(f);
		System.out.println(Factorization.singularValueDecomposition(m));*/
	}
	
	@Test public void test_schur() {
		double[][] f = {{4,0,1},{1,3,-1},{-1,0,2}};
		//Matrix schur = Factorization.schurDecompose(new Matrix(f));
		
	}
}
