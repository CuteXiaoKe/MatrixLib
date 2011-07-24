package matrixLib.test;

/**
 * Suite of JUnit tests to test the entire project at once
 * @author Bryan Cuccioli
 */

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

// specify a runner class: Suite.class
@RunWith(Suite.class)

// specify an array of test classes
@Suite.SuiteClasses({
	ComplexNumberTest.class,
	FactorizationTest.class,
	MatrixTest.class,
	NormTest.class,
	PatternTest.class,
	SquareMatrixOpsTest.class,
	VectorTest.class}
)

// empty class just used to make the compiler happy
public class GlobalTest {}