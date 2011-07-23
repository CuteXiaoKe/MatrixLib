package matrixLib.test;

import static org.junit.Assert.*;
import org.junit.Test;

import matrixLib.ComplexNumber;
import matrixLib.Matrix;
import matrixLib.Vector;
import java.util.LinkedList;

/**
 * Library of JUnit 4 tests for the Matrix class
 * @author Bryan Cuccioli
 */

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
	
	@Test public void add_sub_scale() {
		double[][] f = {{1,2},{3,4}};
		double[][] g = {{5,6},{7,8}};
		double[][] sum = {{6,8},{10,12}};
		double[][] doubled = {{2,4},{6,8}};
		
		assertTrue((new Matrix(f)).add(new Matrix(g)).equals(new Matrix(sum)));
		assertTrue((new Matrix(sum)).subtract(new Matrix(f)).equals(new Matrix(g)));
		assertTrue((new Matrix(f)).multiply(2).equals(new Matrix(doubled)));
	}
	
	@Test public void orthonormalize() {
		double[][] original = {{1,-1,4},{1,4,-2},{1,4,2},{1,-1,0}};
		double[][] basis = {{.5,-.5,.5},{.5,.5,-.5},{.5,.5,.5},{.5,-.5,-.5}};
		assertTrue((new Matrix(original)).orthonormalize().equals(new Matrix(basis)));
	}
	
	@Test public void rref() {
		double[][] mat1 = {{1,2},{3,4}};
		double[][] ref1 = {{1,0},{0,1}};
		double[][] mat2 = {{1,2,3},{4,5,6},{7,8,9}};
		double[][] ref2 = {{1,0,-1},{0,1,2},{0,0,0}};
		double[][] mat3 = {{1,2,3},{4,5,6}};
		double[][] ref3 = {{1,0,-1},{0,1,2}};
		
		// test row reductions
		assertTrue((new Matrix(mat1)).rref().equals(new Matrix(ref1)));
		assertTrue((new Matrix(mat2)).rref().equals(new Matrix(ref2)));
		assertTrue((new Matrix(mat3)).rref().equals(new Matrix(ref3)));
		
		// test rank and nullity computations
		assertTrue((new Matrix(mat1)).rank() == 2);
		assertTrue((new Matrix(mat1)).nullity() == 0);
		assertTrue((new Matrix(mat2)).rank() == 2);
		assertTrue((new Matrix(mat2)).nullity() == 1);
	}
	
	@Test public void imageBasis() {
		double[][] mat = {{1,0,1,3,0},{0,1,1,2,0},{1,1,2,5,1},{0,0,0,0,0}};
		LinkedList<Vector> basis = (new Matrix(mat)).imageBasis();
		
		// build the list of expected vectors
		LinkedList<Vector> expected = new LinkedList<Vector>();
		double[] vec1 = {1,0,1,0}; double[] vec2 = {0,1,1,0}; double[] vec3 = {0,0,1,0};
		Vector v1 = new Vector(vec1), v2 = new Vector(vec2), v3 = new Vector(vec3);
		expected.add(v1); expected.add(v2); expected.add(v3);
		
		for (Vector v : basis) {
			// make sure each computed vector is in the expecteds list
			boolean contains = false;
			for (Vector exp : expected) {
				if (v.equals(exp)) {
					contains = true;
				}
			}
			assertTrue(contains);
		}
	}
}
