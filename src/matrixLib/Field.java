package matrixLib;

/**
 * Represents an abstract field;
 * either ComplexNumber or RealNumber will instantiate this
 * @author Bryan Cuccioli
 */
public abstract class Field {
	
	public abstract String toString();
	public abstract boolean equals(Field f);
}
