package matrixLib;

/**
 * Specifies common and necessary operations on elements of a field
 * @author Bryan Cuccioli
 */

public interface Field {

	/**
	 * Adds two elements of a field
	 * @param f the field element to add to this
	 * @return the sum of this and f
	 */
	public Field add(Field f);
	
	/**
	 * Subtracts one element of a field from another
	 * @param f the field element to subtract from this
	 * @return the difference between this and f
	 */
	public Field subtract(Field f);
	
	/**
	 * Multiplies two elements of a field
	 * @param f the field element to multiply with this
	 * @return the product of this and f
	 */
	public Field multiply(Field f);
	
	/**
	 * Divides this field element by f
	 * @param f the field element by which to divide this
	 * @return the quotient of this and f
	 * @throws ArithmeticException tried to divide by zero
	 */
	public Field divide(Field f) throws ArithmeticException;
	
	/**
	 * Computes the absolute value/magnitude/modulus of this field element
	 * @return the magnitude of this field element
	 */
	public double abs();
	
	/**
	 * Determines whether two field elements are equal
	 * @param f the field element to which this is compared
	 * @return whether the two field elements are equal
	 */
	public boolean equals(Field f);
	
	/**
	 * Computes the string representation of this field element
	 * @return the string representation of this field element
	 */
	public String toString();
}
