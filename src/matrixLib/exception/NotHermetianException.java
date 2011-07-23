package matrixLib.exception;

/**
 * Represents the exception that the given matrix is not Hermetian,
 * which is required for the current operation
 * @author Bryan Cuccioli
 */

public class NotHermetianException extends RuntimeException {
	
	private static final long serialVersionUID = 2960821471544223636L;

	/**
	 * Constructs a NotHermetianException with the default message
	 */
	public NotHermetianException() {
		super("The matrix is not Hermetian.");
	}

	/**
	 * Constructs a NotHermetianException with a custom message
	 * @param message the custom message
	 */
	public NotHermetianException(String message) {
		super(message);
	}

	/**
	 * Constructs a NotHermetianException with a custom cause
	 * @param cause the custom cause
	 */
	public NotHermetianException(Throwable cause) {
		super(cause);
	}

	/**
	 * Constructs a NotHermetianException with a custom message and cause
	 * @param message the custom message
	 * @param cause the custom cause
	 */
	public NotHermetianException(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * Computes a string representation of this exception
	 */
	public String toString() {
		return "NotHermetianException - " + this.getMessage();
	}
}
