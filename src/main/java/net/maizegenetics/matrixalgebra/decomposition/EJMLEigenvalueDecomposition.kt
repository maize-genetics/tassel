package net.maizegenetics.matrixalgebra.decomposition

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.EJMLDoubleMatrix
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.row.CommonOps_DDRM
import org.ejml.dense.row.EigenOps_DDRM
import org.ejml.dense.row.factory.DecompositionFactory_DDRM
import org.ejml.interfaces.decomposition.EigenDecomposition_F64

class EJMLEigenvalueDecomposition(matrix: DMatrixRMaj, val isSymmetric: Boolean = true): EigenvalueDecomposition {
    val myDecomposition: EigenDecomposition_F64<DMatrixRMaj>

    init {
        myDecomposition = DecompositionFactory_DDRM.eig(matrix.numRows, true,true)
        val inputMatrix = if (myDecomposition.inputModified()) matrix.copy() else matrix
        val success = myDecomposition.decompose(inputMatrix)
        check(success) {"EJML Eigenvalue decomposition failed"}
    }

    override fun getEigenvalues(): DoubleArray {
        val numberOfValues = myDecomposition.numberOfEigenvalues
        return DoubleArray(numberOfValues) { ndx -> myDecomposition.getEigenvalue(ndx).real}
    }

    override fun getEigenvalue(i: Int): Double {
        return myDecomposition.getEigenvalue(i).real
    }

    override fun getEigenvectors(): DoubleMatrix {
        //create a matrix whose rows are eigenvectors
        return EJMLDoubleMatrix(EigenOps_DDRM.createMatrixV(myDecomposition))
    }

    override fun getEigenvalueMatrix(): DoubleMatrix {
        return EJMLDoubleMatrix(CommonOps_DDRM.diag(*eigenvalues))
    }

}