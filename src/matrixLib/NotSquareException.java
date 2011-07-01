package matrixLib;

public class NotSquareException extends RuntimeException {

	private static final long serialVersionUID = -1623777304308609186L;

	public NotSquareException() {
		super("The matrix is not square.");
	}

	public NotSquareException(String arg0) {
		super(arg0);
	}

	public NotSquareException(Throwable arg0) {
		super(arg0);
	}

	public NotSquareException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

	public String toString() {
		return "NotSquareExpception - matrix was not square.";
	}
}
