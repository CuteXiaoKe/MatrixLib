package matrixLib.test;

/**
 * Library of JUnit 4 tests for the Matrix class
 * @author Bryan Cuccioli
 */

import static org.junit.Assert.*;
import org.junit.Test;

import matrixLib.ComplexNumber;
import matrixLib.Matrix;
import matrixLib.Vector;

public class MatrixTest {

	@Test public void constructors() {
		// test all of the constructors
		double[][] zero = {{0,0,0},{0,0,0},{0,0,0}};
		double[][] id = {{1,0,0},{0,1,0},{0,0,1}};
		ComplexNumber[][] vm = {{new ComplexNumber(1,2),new ComplexNumber(5,6)},{new ComplexNumber(3,4),new ComplexNumber(7,8)}};
		ComplexNumber[][] vecs = {{new ComplexNumber(1,2),new ComplexNumber(3,4)},{new ComplexNumber(5,6),new ComplexNumber(7,8)}};
		Vector[] vectors = {new Vector(vecs[0]), new Vector(vecs[1])};
		
		assertTrue((new Matrix(3,3)).equals(new Matrix(zero)));
		assertTrue((new Matrix(3)).equals(new Matrix(id)));
		assertTrue((new Matrix(vectors)).equals(new Matrix(vm)));
	}
	
	@Test public void getVector() {
		double[][] mat = {{1,2,3},{4,5,6},{7,8,9}};
		double[] vec = {2,5,8};
		assertTrue((new Matrix(mat)).getVector(1).equals(new Vector(vec)));
	}
	
	@Test public void transposes() {
		ComplexNumber[][] mat = {{new ComplexNumber(2,1),new ComplexNumber(3,-1),new ComplexNumber(4,-3)},
								{new ComplexNumber(4,0),new ComplexNumber(6,-1),new ComplexNumber(2,5)},
								{new ComplexNumber(0,3),new ComplexNumber(2,-1),new ComplexNumber(1,3)}};
		ComplexNumber[][] trans = {{new ComplexNumber(2,1),new ComplexNumber(4,0),new ComplexNumber(0,3)},
				{new ComplexNumber(3,-1),new ComplexNumber(6,-1),new ComplexNumber(2,-1)},
				{new ComplexNumber(4,-3),new ComplexNumber(2,5),new ComplexNumber(1,3)}};
		ComplexNumber[][] herm = {{new ComplexNumber(2,-1),new ComplexNumber(4,0),new ComplexNumber(0,-3)},
				{new ComplexNumber(3,1),new ComplexNumber(6,1),new ComplexNumber(2,1)},
				{new ComplexNumber(4,3),new ComplexNumber(2,-5),new ComplexNumber(1,-3)}};
		
		// test transpose and conjugate transpose operations
		assertTrue((new Matrix(mat)).transpose().equals(new Matrix(trans)));
		assertTrue((new Matrix(mat)).conjugateTranspose().equals(new Matrix(herm)));
	}
	
	@Test public void multiplications() {
		
		ComplexNumber[][] z1 = {{new ComplexNumber(1,1),new ComplexNumber(2,3)},{new ComplexNumber(1,2),new ComplexNumber(3,1)}};
		ComplexNumber[][] z2 = {{new ComplexNumber(2,1),new ComplexNumber(2,5)},{new ComplexNumber(1,1),new ComplexNumber(-5,1)}};
		ComplexNumber[][] res = {{new ComplexNumber(0,8),new ComplexNumber(-16,-6)},{new ComplexNumber(2,9),new ComplexNumber(-24,7)}};
		ComplexNumber[] vec = {new ComplexNumber(4,2),new ComplexNumber(5,3)};
		ComplexNumber[] vres = {new ComplexNumber(1,39),new ComplexNumber(-26,-4)};
		
		assertTrue((new Matrix(z1)).multiply(new Matrix(z2)).equals(new Matrix(res)));
		assertTrue((new Matrix(z2)).multiply(new Vector(vec)).equals(new Vector(vres)));
	}
	
	@Test public void add() {
		double[][] f = {{1,2},{3,4}};
		double[][] g = {{5,6},{7,8}};
		double[][] h = {{6,8},{10,12}};
		assertTrue((new Matrix(f)).add(new Matrix(g)).equals(new Matrix(h)));
	}
	
	@Test public void orthonormalize() {
		double[][] original = {{1,-1,4},{1,4,-2},{1,4,2},{1,-1,0}};
		double[][] basis = {{.5,-.5,.5},{.5,.5,-.5},{.5,.5,.5},{.5,-.5,-.5}};
		assertTrue((new Matrix(original)).orthonormalize().equals(new Matrix(basis)));
	}
}
