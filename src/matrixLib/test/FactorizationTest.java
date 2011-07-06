package matrixLib.test;

import matrixLib.ComplexNumber;
import matrixLib.Factorization;
import matrixLib.Matrix;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Library of routines for testing the Factorization library
 * @author Bryan Cuccioli
 */

public class FactorizationTest {

	@Test public void QRDecompose() {
		double[][] m1 = {{1.92,2.32},{.56,2.76}};
		double[][] q1 = {{.96, -.28},{.28, .96}};
		double[][] r1 = {{2,3},{0,2}};
		
		Matrix[] qr1 = Factorization.QRDecompose(new Matrix(m1));
		assertTrue(qr1[0].equals(new Matrix(q1)));
		assertTrue(qr1[1].equals(new Matrix(r1)));
		
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
	
	@Test public void choleskyDecompose() {

		double[][] diag = {{1,0,0,0},{0,25,0,0},{0,0,64,0},{0,0,0,100}};
		double[][] d_exp = {{1,0,0,0},{0,5,0,0},{0,0,8,0},{0,0,0,10}};;
		double[][] m = {{2,-1,0},{-1,2,-1},{0,-1,2}};
		double[][] m_exp = {{Math.sqrt(2),0,0},{-.5*Math.sqrt(2),.5*Math.sqrt(6),0},{0,-Math.sqrt(6)/3.0,2.0/Math.sqrt(3)}};
		//double[][] big = {{16,-3,5,-8},{-3,16,-5,-8},{5,-5,24,0},{-8,-8,0,21}};
		
		assertTrue(Factorization.choleskyDecompose(new Matrix(diag)).equals(new Matrix(d_exp)));
		assertTrue(Factorization.choleskyDecompose(new Matrix(m)).equals(new Matrix(m_exp)));
	}
}
