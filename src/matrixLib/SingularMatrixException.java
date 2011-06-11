package matrixLib;

public class SingularMatrixException extends Exception {

	private static final long serialVersionUID = 5248977736388665535L;

	public SingularMatrixException() {
		super("The matrix is singular.");
	}

	public SingularMatrixException(String message) {
		super(message);
	}

	public SingularMatrixException(Throwable cause) {
		super(cause);
	}

	public SingularMatrixException(String message, Throwable cause) {
		super(message, cause);
	}

	public String toString() {
		return "SingularMatrixException - matrix is singular.";
	}
}
