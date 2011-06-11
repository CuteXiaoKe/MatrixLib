package matrixLib;

/**
 * Implementation of vectors (column matrices)
 * @author Bryan Cuccioli
 */

public class Vector extends Matrix {
	
	public Vector(int n) {
		
		
	}
	
	public Vector(float[][] entries) throws DimensionMismatchException {
		
		//float[][] data = new float[entries.length][1];
		
		super(entries);
		
		if (entries[0].length > 1) {
			throw new DimensionMismatchException();
		}
	}

}
