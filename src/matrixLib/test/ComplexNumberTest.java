package matrixLib.test;

/**
 * Library of JUnit 4 tests for the ComplexNumber class
 * @author Bryan Cuccioli
 */

import static org.junit.Assert.*;
import org.junit.Test;
import matrixLib.ComplexNumber;

public class ComplexNumberTest {

	@Test public void equals() {
		ComplexNumber a = new ComplexNumber(3, 4);
		ComplexNumber b = new ComplexNumber(3, 4);
		assertTrue(a.equals(b));
	}
	
	@Test public void add() {
		ComplexNumber a = new ComplexNumber(4,3);
		ComplexNumber b = new ComplexNumber(2,-6.5);
		assertTrue(a.add(b).equals(new ComplexNumber(6, -3.5)));
	}
	
	@Test public void subtract() {
		ComplexNumber a = new ComplexNumber(4,3);
		ComplexNumber b = new ComplexNumber(2,-6.5);
		assertTrue(a.subtract(b).equals(new ComplexNumber(2, 9.5)));
	}
	
	@Test public void multiply() {
		ComplexNumber a = new ComplexNumber(4,3);
		ComplexNumber b = new ComplexNumber(2,-6.5);
		assertTrue(a.multiply(b).equals(new ComplexNumber(27.5, -20)));
	}
	
	@Test public void divide() {
		ComplexNumber a = new ComplexNumber(1,2);
		ComplexNumber b = new ComplexNumber(2,-4);
		assertTrue(a.divide(b).equals(new ComplexNumber(-.3, .4)));
	}
	
	@Test public void reciprocal() {
		ComplexNumber a = new ComplexNumber(4,2);
		assertTrue(a.reciprocal().equals(new ComplexNumber(.2,-.1)));
	}
	
	@Test public void abs() {
		ComplexNumber a = new ComplexNumber(5,-3);
		assertTrue(a.abs()*a.abs() == 34);
	}
	
	@Test public void arg() {
		ComplexNumber quot = new ComplexNumber(-1,-1).divide(new ComplexNumber(0,1));
		assertTrue(quot.arg() == 3*Math.PI/4);
	}
	
	@Test public void sqrt() {
		ComplexNumber c = new ComplexNumber(12, 16);
		assertTrue(c.sqrt().equals(new ComplexNumber(4,2)));
	}
	
	@Test public void epsilon() {
		ComplexNumber.setEpsilon(1337);
		assertTrue(ComplexNumber.getEpsilon() == 1337);
		ComplexNumber.setEpsilon(1e-15);
	}
}
