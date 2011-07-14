package matrixLib.test;

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
		
		double[][] mat1 = {{2,1,-2},{-3,1,0},{4,3,1}};
		double[][] exp1 = {{2,11.0/5,2.0/5},{-5,-11.0/25,48.0/25},{0,-27.0/25,61.0/25}};
	
		//System.out.println(SquareMatrixOps.hessenberg(new Matrix(mat1)));
		
		//assertTrue(SquareMatrixOps.hessenbergForm(new Matrix(mat1)).equals(new Matrix(exp1)));
	}

	@Test public void triangular() {
		double[][] mat = {{1,2,3,4},{0,5,6,7},{0,0,8,9},{0,0,0,10}};
		Matrix m = new Matrix(mat);
		assertTrue(Pattern.isUpperTriangular(m));
		assertTrue(Pattern.isLowerTriangular(m.transpose()));
	}
}
