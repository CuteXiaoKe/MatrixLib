# MatrixLib

MatrixLib is an efficient, numerically precise library of Java classes and routines that can be incorporated into larger projects. It is adept at handling very large matrices, and the algorithms it employs are designed to be as efficient as possible. For example, the determinant algorithm works via a modification of Gauss-Jordan reduction, achieving an asymptotic time complexity of O(n^3) instead of O(n!), as per the standard formula.

## Important Links

* Homepage: [https://bcuccioli.github.com/MatrixLib](http://bcuccioli.github.com/MatrixLib)
* Source tree: [https://github.com/bcuccioli/MatrixLib](https://github.com/bcuccioli/MatrixLib)
* Dev. blog: [http://bcuccioli.com/blog](http://bcuccioli.com/blog)

## Features

* Support for complex numbers and common operations over the complex numbers, including arithmetic and square root

* Support for matrices over complex numbers

* Matrix class that handles common matrix operations, including pivoted Gauss-Jordan reduction

* Vector class that supports common operations, including cross product, Hermitian inner product (a generalization of the dot product) and projection

* Computation of determinant and matrix inverse in O(n^3) asymptotic time complexity

* Routines to compute numerically stable pivoted Cholesky factorization, LU decomposition, QR decomposition, and Schur decomposition for real and complex arbitrarily large matrices.

* Computation of eigenvalues and singular values via explicit shifted QR algorithm and eigenvectors via inverse iteration

* Computation of standard norms, including the general p-norm, Frobenius norm, infinity norm and spectral norm

* Orthonormalization of matrices via modified Gram-Schmidt procedure and conversion of matrices to Hessenberg form via Householder transformations

## Known Issues

* Computation of eigenvectors produces incorrect results for some matrices. This causes their Schur decomposition to be computed incorrectly.

