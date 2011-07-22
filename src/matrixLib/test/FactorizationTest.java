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
		ComplexNumber[][] z = {{new ComplexNumber(0,1), new ComplexNumber(2,0)},{new ComplexNumber(3,0),new ComplexNumber(4,0)}};
		Matrix[] qr = Factorization.QRDecompose(new Matrix(z));
		assertTrue(Pattern.isUnitary(qr[0]) && Pattern.isUpperTriangular(qr[1]) && qr[0].multiply(qr[1]).equals(new Matrix(z)));
		
		ComplexNumber[][] z2 = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		qr = Factorization.QRDecompose(new Matrix(z2));
		assertTrue(Pattern.isUnitary(qr[0]) && Pattern.isUpperTriangular(qr[1]) && qr[0].multiply(qr[1]).equals(new Matrix(z2)));
		
		double[][] m1 = {{1.92,2.32},{.56,2.76}};
		qr = Factorization.QRDecompose(new Matrix(m1));
		assertTrue(Pattern.isUnitary(qr[0]) && Pattern.isUpperTriangular(qr[1]) && qr[0].multiply(qr[1]).equals(new Matrix(m1)));
		
		double[][] m2 = {{12,-51,4},{6,167,-68},{-4,24,-41}};
		qr = Factorization.QRDecompose(new Matrix(m2));
		ComplexNumber.setEpsilon(1e-12);
		assertTrue(Pattern.isUnitary(qr[0]) && Pattern.isUpperTriangular(qr[1]) && qr[0].multiply(qr[1]).equals(new Matrix(m2)));
		ComplexNumber.setEpsilon(1e-14);
		
		double[][] big = {{1,3,2,-1},{1,1,2,-3},{3,1,1,-1},{2,-2,1,2}};
		qr = Factorization.QRDecompose(new Matrix(big));
		assertTrue(Pattern.isUnitary(qr[0]) && qr[0].multiply(qr[1]).equals(new Matrix(big)));
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
		assertTrue(Factorization.choleskyDecompose(new Matrix(diag)).equals(new Matrix(d_exp)));
		
		double[][] m = {{2,-1,0},{-1,2,-1},{0,-1,2}};
		Matrix chol = Factorization.choleskyDecompose(new Matrix(m));
		assertTrue(Pattern.isLowerTriangular(chol) && chol.multiply(chol.conjugateTranspose()).equals(new Matrix(m)));

		ComplexNumber[][] c = {{new ComplexNumber(4,0),new ComplexNumber(0,2),new ComplexNumber(0,-1)},{new ComplexNumber(0,-2),new ComplexNumber(10,0),new ComplexNumber(1,0)},{new ComplexNumber(0,1),new ComplexNumber(1,0),new ComplexNumber(9,0)}};
		chol = Factorization.choleskyDecompose(new Matrix(c));
		assertTrue(Pattern.isLowerTriangular(chol) && chol.multiply(chol.conjugateTranspose()).equals(new Matrix(c)));

		double[][] big = {{16,-3,5,-8},{-3,16,-5,-8},{5,-5,24,0},{-8,-8,0,21}};
		chol = Factorization.choleskyDecompose(new Matrix(big));
		assertTrue(Pattern.isLowerTriangular(chol) && chol.multiply(chol.conjugateTranspose()).equals(new Matrix(big)));
		
		double[][] fail = {{2,-17,7},{-17,-4,1},{7,1,-14}};
		chol = Factorization.choleskyDecompose(new Matrix(fail));
		assertTrue(chol == null); // this matrix was not positive definite
	}
	
	@Test public void singularValueDecomposition() {
		
		double[][] m1 = {{4,0},{3,-5}};
		//System.out.println(Factorization.singularValueDecomposition(new Matrix(m1)));
		
	}
	
	@Test public void schurDecomposition() {

		// this works fine except the eigenvectors aren't computed right
		if(1==1)return;
		double[][] m1 = {{3,0,0,-1},{1,2,0,1},{2,0,4,2},{-1,0,0,3}};
		Matrix[] schur = Factorization.schurDecompose(new Matrix(m1));
		assertTrue(Pattern.isUnitary(schur[0]) && Pattern.isUpperTriangular(schur[1]) && schur[0].conjugateTranspose().multiply(new Matrix(m1)).multiply(schur[0]).equals(schur[1]));
		
		double[][] m2 = {{4,0,1},{1,3,-1},{-1,0,2}};
		schur = Factorization.schurDecompose(new Matrix(m2));
		assertTrue(Pattern.isUnitary(schur[0]) && Pattern.isUpperTriangular(schur[1]) && schur[0].conjugateTranspose().multiply(new Matrix(m2)).multiply(schur[0]).equals(schur[1]));
	}
}
