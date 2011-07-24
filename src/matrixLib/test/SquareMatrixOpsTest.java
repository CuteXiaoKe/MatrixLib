package matrixLib.test;

import matrixLib.ComplexNumber;
import matrixLib.Matrix;
import matrixLib.SquareMatrixOps;
import matrixLib.Pattern;
import matrixLib.Vector;
import matrixLib.exception.SingularMatrixException;

import org.junit.Test;
import static org.junit.Assert.assertTrue;

/**
 * Library of test routines for the SquareMatrixOps class
 * @author Bryan Cuccioli
 */

public class SquareMatrixOpsTest {
	
	@Test public void det() {
		// set up matrices for a testing a variety of determinants
		double[][] f = {{1,2,3},{4,5,6},{7,8,7}};
		double[][] g = {{1,2},{3,4}};
		double[][] h = {{3,0,6,-3},{0,2,3,0},{-4,-7,2,0},{2,0,1,10}};
		double[][] singular = {{0,0,0},{1,0,1},{2,-1,-2}};
		ComplexNumber[][] i = {{new ComplexNumber(2,1),new ComplexNumber(3,-1),new ComplexNumber(4,-3)},
				{new ComplexNumber(4,0),new ComplexNumber(6,-1),new ComplexNumber(2,5)},
				{new ComplexNumber(0,3),new ComplexNumber(2,-1),new ComplexNumber(1,3)}};
		
		assertTrue(SquareMatrixOps.determinant(new Matrix(singular)).isZero());
		assertTrue(SquareMatrixOps.determinant(new Matrix(f)).equals(new ComplexNumber(6, 0)));
		assertTrue(SquareMatrixOps.determinant(new Matrix(g)).equals(new ComplexNumber(-2, 0)));
		assertTrue(SquareMatrixOps.determinant(new Matrix(9)).equals(new ComplexNumber(1, 0)));
		assertTrue(SquareMatrixOps.determinant(new Matrix(h)).equals(new ComplexNumber(1404, 0)));
		assertTrue(SquareMatrixOps.determinant(new Matrix(i)).equals(new ComplexNumber(-118,-84)));
	}
	
	@Test public void inverse() {
		double[][] mat = {{5,19},{1,4}};
		double[][] sing = {{1,2,3},{4,5,6},{7,8,9}}; // not invertible
		double[][] inv = {{4,-19},{-1,5}};
		
		double[][] tri = {{1,0,0},{3,2,0},{4,6,5}}, tri_exp = {{1,0,0},{-1.5,.5,0},{1,-.6,.2}};
		double[][] herm = {{6,2,-2},{2,6,-2},{-2,-2,10}}, herm_exp = {{7.0/36,-1.0/18,1.0/36},{-1.0/18,7.0/36,1.0/36},{1.0/36,1.0/36,1.0/9}};
		
		ComplexNumber[][] c = {{new ComplexNumber(4,0),new ComplexNumber(0,2),new ComplexNumber(0,-1)},{new ComplexNumber(0,-2),new ComplexNumber(10,0),new ComplexNumber(1,0)},{new ComplexNumber(0,1),new ComplexNumber(1,0),new ComplexNumber(9,0)}};
		ComplexNumber[][] c_exp = {{new ComplexNumber(4,0),new ComplexNumber(0,2),new ComplexNumber(0,-1)},{new ComplexNumber(0,-2),new ComplexNumber(10,0),new ComplexNumber(1,0)},{new ComplexNumber(0,1),new ComplexNumber(1,0),new ComplexNumber(9,0)}};
		
		assertTrue(SquareMatrixOps.inverse(new Matrix(tri)).equals(new Matrix(tri_exp)));
		assertTrue(SquareMatrixOps.inverse(new Matrix(herm)).equals(new Matrix(herm_exp)));
		assertTrue(SquareMatrixOps.inverse(new Matrix(mat)).equals(new Matrix(inv)));
		
		boolean failed = false;
		try {
			SquareMatrixOps.inverse(new Matrix(sing));
		}
		catch (SingularMatrixException e) {
			failed = true;
		}
		assertTrue(failed);
	}
	
	@Test public void eigenvalues() {
		double[][] mat1 = {{3,0,0},{1,3,1},{2,-1,1}};
		Matrix m = new Matrix(mat1);
		ComplexNumber[] evs = SquareMatrixOps.eigenvalues(m);
		for (ComplexNumber ev : evs) {
			assertTrue(SquareMatrixOps.determinant(m.subtract(Pattern.diag(ev, m.rows()))).isZero());
		}
		Vector[] vecs = SquareMatrixOps.eigenvectors(m, evs);
		for (Vector v : vecs) {
			assertTrue(m.multiply(v).normalize().equals(v.normalize()));
		}
		
		double[][] mat2 = {{0,-1},{1,0}};
		m = new Matrix(mat2);
		evs = SquareMatrixOps.eigenvalues(m);
		for (ComplexNumber ev : evs) {
			assertTrue(SquareMatrixOps.determinant(m.subtract(Pattern.diag(ev, m.rows()))).isZero());
		}
		
		double[][] mat3 = {{1,3,2,-1},{1,1,2,-3},{3,1,1,-1},{2,-2,1,2}};
		m = new Matrix(mat3);
		evs = SquareMatrixOps.eigenvalues(m);
		for (ComplexNumber ev : evs) {
			assertTrue(SquareMatrixOps.determinant(m.subtract(Pattern.diag(ev, m.rows()))).isZero());
		}
	}
	
	@Test public void pow() {
		double[][] mat_arr = {{1,2},{3,4}}, matcube = {{37,54},{81,118}};
		Matrix mat = new Matrix(mat_arr);
		
		assertTrue(SquareMatrixOps.pow(mat,1).equals(mat));
		assertTrue(SquareMatrixOps.pow(mat,3).equals(new Matrix(matcube)));
		assertTrue(Pattern.isIdentity(SquareMatrixOps.pow(mat,0)));
	}
}
