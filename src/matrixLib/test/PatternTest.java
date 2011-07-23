package matrixLib.test;

import matrixLib.ComplexNumber;
import matrixLib.Pattern;
import matrixLib.Matrix;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Library of routines for testing the Pattern class
 * @author Bryan Cuccioli
 */

public class PatternTest {

	@Test public void hessenberg() {
		
		double[][] hess = {{1,2,3,4},{5,6,7,8},{0,9,10,0},{0,0,11,12}};
		double[][] mat = {{2,1,-2},{-3,1,0},{4,3,1}};
	
		assertTrue(Pattern.isUpperHessenberg(new Matrix(hess)));
		assertTrue(Pattern.isUpperHessenberg(Pattern.hessenberg(new Matrix(mat))));
	}

	@Test public void triangular() {
		
		double[][] ut = {{1,2,3,4},{0,5,6,7},{0,0,8,9},{0,0,0,10}};
		double[][] lt = {{1,0,0,0},{2,3,0,0},{4,5,6,0},{7,8,9,10}};
		double[][] diag = {{1,0,0,0},{0,2,0,0},{0,0,3,0},{0,0,0,4}};
		
		assertTrue(Pattern.isUpperTriangular(new Matrix(ut)));
		assertTrue(Pattern.isLowerTriangular(new Matrix(lt)));
		assertTrue(Pattern.isDiagonal(new Matrix(diag)));
	}
	
	@Test public void diagonal() {
		
		double[][] mat1 = {{1,2},{3,4}};
		double[][] exp1 = {{1,0,0,0},{0,1,0,0},{0,0,1,2},{0,0,3,4}};
		double[][] diag1 = {{3,0,0},{0,3,0},{0,0,3}};
		double[][] diag2 = {{3,0,0},{0,4,0},{0,0,7}};
		ComplexNumber[] diag_els = {new ComplexNumber(3,0),new ComplexNumber(4,0),new ComplexNumber(7,0)};
		
		assertTrue(Pattern.blockDiagonal(new Matrix(mat1),2).equals(new Matrix(exp1)));
		assertTrue(Pattern.isIdentity(Pattern.diag(new ComplexNumber(2,1),3).multiply((new ComplexNumber(2,1)).reciprocal())));
		assertTrue(Pattern.diag(new ComplexNumber(3,0), 3).equals(new Matrix(diag1)));
		assertTrue(Pattern.diag(diag_els).equals(new Matrix(diag2)));
	}
	
	@Test public void symmetry() {
		
		double[][] sym = {{1,7,3},{7,4,-5},{3,-5,6}};
		double[][] asym = {{0,2,-1},{-2,0,-4},{1,4,0}};
		ComplexNumber[][] herm = {{new ComplexNumber(2,0),new ComplexNumber(1,3),new ComplexNumber(4,2)},
				{new ComplexNumber(1,-3),new ComplexNumber(6,0),new ComplexNumber(5,3)},
				{new ComplexNumber(4,-2),new ComplexNumber(5,-3),new ComplexNumber(7,0)}};
		
		assertTrue(Pattern.isSymmetric(new Matrix(sym)));
		assertTrue(Pattern.isAntiSymmetric(new Matrix(asym)));
		assertTrue(Pattern.isHermetian(new Matrix(herm)));
	}
	
	@Test public void other() {
		
		double[][] ortho = {{0,-.8,-.6},{.8,-.36,.48},{.6,.48,-.64}};
		ComplexNumber[][] uni = {{new ComplexNumber(0,1),new ComplexNumber(0,1)},
				{new ComplexNumber(0,1),new ComplexNumber(0,-1)}};
		double[][] id = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
		
		assertTrue(Pattern.isOrthogonal(new Matrix(ortho)));
		assertTrue(Pattern.isUnitary((new Matrix(uni)).multiply(1.0/Math.sqrt(2))));
		assertTrue(Pattern.isIdentity(new Matrix(id)));
	}
}
