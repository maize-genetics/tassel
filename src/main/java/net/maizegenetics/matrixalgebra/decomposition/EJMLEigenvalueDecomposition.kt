package net.maizegenetics.matrixalgebra.decomposition

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.EJMLDoubleMatrix
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.row.CommonOps_DDRM
import org.ejml.dense.row.EigenOps_DDRM
import org.ejml.dense.row.NormOps_DDRM
import org.ejml.dense.row.factory.DecompositionFactory_DDRM
import org.ejml.interfaces.decomposition.EigenDecomposition_F64

class EJMLEigenvalueDecomposition(matrix: DMatrixRMaj, val isSymmetric: Boolean = true): EigenvalueDecomposition {
    private val myDecomposition: EigenDecomposition_F64<DMatrixRMaj>
    private val myShift: Double

    init {
        var shift = 0.0
        var decomposition = decompose(matrix)
        if (decomposition == null) {
            // EJML's symmetric QR algorithm can fail to converge on a matrix with a large null space, such as
            // the compressed kinship matrix ZKZ' in MLM. Shifting the diagonal by at least the spectral radius
            // moves those eigenvalues away from zero and leaves the eigenvectors unchanged.
            shift = NormOps_DDRM.normF(matrix)
            val shifted = matrix.copy()
            for (i in 0 until shifted.numRows) shifted.add(i, i, shift)
            decomposition = decompose(shifted)
        }
        myDecomposition = checkNotNull(decomposition) {"EJML Eigenvalue decomposition failed"}
        myShift = shift
    }

    private fun decompose(matrix: DMatrixRMaj): EigenDecomposition_F64<DMatrixRMaj>? {
        val decomposition = DecompositionFactory_DDRM.eig(matrix.numRows, true,true)
        val inputMatrix = if (decomposition.inputModified()) matrix.copy() else matrix
        return if (decomposition.decompose(inputMatrix)) decomposition else null
    }

    override fun getEigenvalues(): DoubleArray {
        val numberOfValues = myDecomposition.numberOfEigenvalues
        return DoubleArray(numberOfValues) { ndx -> getEigenvalue(ndx) }
    }

    override fun getEigenvalue(i: Int): Double {
        return myDecomposition.getEigenvalue(i).real - myShift
    }

    override fun getEigenvectors(): DoubleMatrix {
        //create a matrix whose rows are eigenvectors
        return EJMLDoubleMatrix(EigenOps_DDRM.createMatrixV(myDecomposition))
    }

    override fun getEigenvalueMatrix(): DoubleMatrix {
        return EJMLDoubleMatrix(CommonOps_DDRM.diag(*eigenvalues))
    }

}
