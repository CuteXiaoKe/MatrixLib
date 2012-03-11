<?require("head.php");?>

<p>This is the SourceForge page for MatrixLib, a matrix library written in Java. MatrixLib is an efficient, numerically precise library of Java classes and routines that can be incorporated into larger projects. It is adept at handling very large matrices, and the algorithms it employs are designed to be as efficient as possible. For example, the determinant algorithm works via a modification of Gauss-Jordan reduction, achieving an asymptotic time complexity of O(n<sup>3</sup>) instead of O(n!), as per the standard formula.</p>

<p>You can download the .jar file that can be imported into your projects <a href="https://sourceforge.net/projects/jmatrixlib/files/">here</a>, and you can view the source tree on GitHub <a href="https://github.com/smango/MatrixLib">here</a>. Additionally, you can view the development blog, wherein issues related to design and implementation are discussed, as a subsection of my blog <a href="http://mebrah.net/?cat=3">here</a>.</p>

<h3>Features</h3>
<ul>
<li>Support for complex numbers and common operations over the complex numbers, including arithmetic and square root</li><br />
<li>Support for matrices over complex numbers</li><br />
<li>Matrix class that handles common matrix operations, including pivoted Gauss-Jordan reduction</li><br />
<li>Vector class that supports common operations, including cross product, Hermitian inner product (a generalization of the dot product) and projection</li><br />
<li>Computation of determinant and matrix inverse in O(n<sup>3</sup>) asymptotic time complexity</li><br />
<li>Routines to compute numerically stable pivoted Cholesky factorization, LU decomposition, QR decomposition, and Schur decomposition for real and complex arbitrarily large matrices.</li><br />
<li>Computation of eigenvalues and singular values via explicit shifted QR algorithm and eigenvectors via inverse iteration</li><br />
<li>Computation of standard norms, including the general p-norm, Frobenius norm, infinity norm and spectral norm</li><br />
<li>Orthonormalization of matrices via modified Gram-Schmidt procedure and conversion of matrices to Hessenberg form via Householder transformations</li><br />
</ul>

<h3>Known issues</h3>
<ul>
<li>Computation of eigenvectors produces incorrect results for some matrices. This causes their Schur decomposition to be computed incorrectly.</li>
</ul>

<?require("foot.php");?>
