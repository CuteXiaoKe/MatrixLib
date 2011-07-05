package matrixLib.test;

/**
 * Library of JUnit 4 tests for the Vector class
 * @author Bryan Cuccioli
 */

import static org.junit.Assert.*;

import org.junit.Test;

import matrixLib.ComplexNumber;
import matrixLib.Matrix;
import matrixLib.Vector;

public class VectorTest {
	
	@Test public void add() {
		ComplexNumber[] a = {new ComplexNumber(1,3),new ComplexNumber(4,2),new ComplexNumber(-5.5,-6)};
		ComplexNumber[] b = {new ComplexNumber(2,-3),new ComplexNumber(-2,-2),new ComplexNumber(6,-3)};
		ComplexNumber[] sum = {new ComplexNumber(3,0),new ComplexNumber(2,0),new ComplexNumber(.5,-9)};
		
		assertTrue((new Vector(a)).add(new Vector(b)).equals(new Vector(sum)));
	}
	
	@Test public void subtract() {
		ComplexNumber[] a = {new ComplexNumber(1,3),new ComplexNumber(4,2),new ComplexNumber(-5.5,-6)};
		ComplexNumber[] b = {new ComplexNumber(2,-3),new ComplexNumber(-2,-2),new ComplexNumber(6,-3)};
		ComplexNumber[] sum = {new ComplexNumber(-1,6),new ComplexNumber(6,4),new ComplexNumber(-11.5,-3)};
		
		assertTrue((new Vector(a)).subtract(new Vector(b)).equals(new Vector(sum)));
	}
	
	@Test public void dot() {
		ComplexNumber[] a = {new ComplexNumber(1,2),new ComplexNumber(3,4)};
		ComplexNumber[] b = {new ComplexNumber(3,1),new ComplexNumber(4,5)};
		ComplexNumber c = new ComplexNumber(37,6);
		
		ComplexNumber[] x = {new ComplexNumber(2,1),new ComplexNumber(0,0),new ComplexNumber(4,-5)};
		ComplexNumber[] y = {new ComplexNumber(1,1),new ComplexNumber(2,1),new ComplexNumber(0,0)};
		ComplexNumber z = new ComplexNumber(3,-1);
		
		assertTrue((new Vector(a)).dot(new Vector(b)).equals(c));
		assertTrue((new Vector(x)).dot(new Vector(y)).equals(z));
	}
	
	@Test public void cross() {
		double[] a = {3,-3,1};
		double[] b = {4,9,2};
		double[] cp = {-15,-2,39};
		
		assertTrue((new Vector(a)).cross(new Vector(b)).equals(new Vector(cp)));
	}
	
	@Test public void proj() {
		ComplexNumber[] a = {new ComplexNumber(3,1), new ComplexNumber(3, -1)};
		ComplexNumber[] b = {new ComplexNumber(2,1), new ComplexNumber(2, -1)};
		Vector v1 = new Vector(a);
		Vector v2 = new Vector(b);

		assertTrue(v1.dot(v2).equals(new ComplexNumber(14,0)));
		assertTrue(v2.dot(v2).equals(new ComplexNumber(10,0)));
		
		// projection of v1 onto v2
		ComplexNumber[] expected={new ComplexNumber(2.8,1.4),new ComplexNumber(2.8,-1.4)};

		assertTrue(v1.proj(v2).equals(new Vector(expected)));
	}
	
	@Test public void normalize() {
		ComplexNumber[] a = {new ComplexNumber(3,1), new ComplexNumber(3, -1)};
		Vector normal = (new Vector(a)).normalize();
		
		// make sure the normal vector has magnitude 1
		double magnitude = matrixLib.Norm.pnorm(normal, 2);
		assertTrue(Math.abs(magnitude-1) < 1e-15);
	}
	
	@Test public void reflector() {
		double[] vec1 = {1,2};
		double[] vec2 = {9,3,-6};
		double[][] ref1 = {{.6,-.8},{-.8,-.6}};
		double[][] ref2 = {{-2.0/7,-3.0/7,6.0/7},{-3.0/7,6.0/7,2.0/7},{6.0/7,2.0/7,3.0/7}};
		
		assertTrue((new Vector(vec1)).reflector().equals(new Matrix(ref1)));
		assertTrue((new Vector(vec2)).reflector().equals(new Matrix(ref2)));
	}
	
	@Test public void generateUnitaryMatrix() {
		double[] firstcol = {0,1,0};
		double[][] mat = {{0,1,0},{1,0,0},{0,0,1}};
		assertTrue((new Vector(firstcol)).generateUnitaryMatrix().equals(new Matrix(mat)));
	}
}
