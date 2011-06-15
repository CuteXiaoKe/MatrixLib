package matrixLib;

public class Demo {

	public static void test_eigen() throws NotSquareException, DimensionMismatchException {
		ComplexNumber[][] z = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		SquareMatrix zm = new SquareMatrix(z);
		ComplexNumber[] evals = zm.eigenvalues();
		System.out.println(evals[0] + " " + evals[1]);
	}
	
	public static void test_proj() throws DimensionMismatchException {
		ComplexNumber[] a = {new ComplexNumber(2,1), new ComplexNumber(3, -1)};
		ComplexNumber[] b = {new ComplexNumber(1,1), new ComplexNumber(4, -1)};
		Vector v1 = new Vector(a);
		Vector v2 = new Vector(b);
		
		// projection of v1 onto v2
		System.out.println(v1.proj(v2));
	}
	
	public static void test_qr() throws NotSquareException, DimensionMismatchException {
		ComplexNumber[][] z = {{new ComplexNumber(1,0), new ComplexNumber(1,1)},{new ComplexNumber(2,-1),new ComplexNumber(3,0)}};
		SquareMatrix zm = new SquareMatrix(z);
		
		Matrix[] qr = zm.QRDecompose();
		System.out.println(qr[0]);
		System.out.println(qr[1]);
	}
	
	public static void test_det() throws NotSquareException {
		float[][] f = {{1,2,3},{4,5,6},{7,8,7}};
		float[][] g = {{1,2},{3,4}};
		SquareMatrix m = new SquareMatrix(f);
		SquareMatrix t = new SquareMatrix(g);	
		
		System.out.println("det(f) = " + m.determinant());
		System.out.println("det(g) = " + t.determinant());
		System.out.println("det(I_5) = " + (new SquareMatrix(4)).determinant());
	}
	
	public static void main(String[] args) throws Exception {
		
		test_proj();
	}

}
