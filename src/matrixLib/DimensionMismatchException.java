package matrixLib;

/**
 * Represents the exception that the dimensions of the matrices
 * do not correspond to the current operation
 * @author Bryan Cuccioli
 *
 */

public class DimensionMismatchException extends Exception {

	private static final long serialVersionUID = -7526472151302870197L;

	public DimensionMismatchException() {
		super("The matrices do not have the correct dimensions.");
	}

	public DimensionMismatchException(String message) {
		super(message);
	}

	public DimensionMismatchException(Throwable cause) {
		super(cause);
	}

	public DimensionMismatchException(String message, Throwable cause) {
		super(message, cause);
	}

	public String toString() {
		return "DimensionMismatchException - the dimensions are not right for the current operation.";
	}
}
