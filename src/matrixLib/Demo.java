package matrixLib;

public class Demo {

	public static void main(String[] args) throws Exception {
		
		float[][] f = {{1,2,3},{4,5,6},{7,8,7}};
		float[][] g = {{1,2},{3,4}};
		SquareMatrix m = new SquareMatrix(f);
		SquareMatrix t = new SquareMatrix(g);
		
		System.out.println("det(f) = " + m.determinant());
		System.out.println("det(g) = " + t.determinant());
		System.out.println("det(I_5) = " + (new SquareMatrix(4)).determinant());
		
		System.out.println(m.toString());
		
		float[][] h = {{1,0,0},{2,3,0},{4,5,6}};
		SquareMatrix blah = new SquareMatrix(h);
		System.out.println("h is lower triangular: " + blah.isLowerTriangular());
		System.out.println(blah.toString());
	}

}
