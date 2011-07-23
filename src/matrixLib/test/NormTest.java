package matrixLib.test;

import matrixLib.Matrix;
import matrixLib.Norm;
import matrixLib.Vector;
import matrixLib.ComplexNumber;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Collection of routines for testing the Norm library
 * @author Bryan Cuccioli
 */

public class NormTest {

	@Test public void vecnorm() {
		double[] vec = {1,2,3};
		
		assertTrue(Norm.frobeniusNorm(new Vector(vec)) == Math.sqrt(14));
		assertTrue(Norm.pnorm(new Vector(vec), 1) == 6);
		assertTrue(Norm.infinityNorm(new Vector(vec)) == 3);
		assertTrue(Norm.pnorm(new Vector(vec), 5) == Math.pow(276, 0.2));
	}
	
	@Test public void matnorm() {
		ComplexNumber[][] mat = {{new ComplexNumber(4,3),new ComplexNumber(-7,0)},{new ComplexNumber(3,0),new ComplexNumber(0,4)}};
		Matrix m = new Matrix(mat);
		
		assertTrue(Norm.frobeniusNorm(m) == Math.sqrt(99));
		assertTrue(Norm.pnorm(m, 1) == 11);
		assertTrue(Norm.infinityNorm(m) == 12);
		
		double[][] spectest = {{1,2,3},{4,5,6},{7,8,7}};
		System.out.println(Norm.spectralNorm(new Matrix(spectest)));
		ComplexNumber.setEpsilon(1e-3);
		assertTrue(Norm.spectralNorm(new Matrix(spectest)).equals(new ComplexNumber(15.8101,0)));
	}
}
