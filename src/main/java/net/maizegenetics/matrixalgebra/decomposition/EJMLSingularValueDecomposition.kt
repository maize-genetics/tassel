package net.maizegenetics.matrixalgebra.decomposition

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.EJMLDoubleMatrix
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.row.factory.DecompositionFactory_DDRM

class EJMLSingularValueDecomposition(val matrix: DMatrixRMaj): SingularValueDecomposition {
    val mySvd = DecompositionFactory_DDRM.svd(matrix.numRows, matrix.numCols, true, true, true)
    val myMatrix : DMatrixRMaj
    val tol = 1e-10

    init {
        myMatrix = if (mySvd.inputModified()) matrix.copy() else matrix
        check(mySvd.decompose(matrix)) {"unable to create singular value decomposition"}
    }

    override fun getU(transpose: Boolean): DoubleMatrix {
        return EJMLDoubleMatrix(mySvd.getU(null, transpose))
    }

    override fun getV(transpose: Boolean): DoubleMatrix {
        return EJMLDoubleMatrix(mySvd.getV(null, transpose))
    }

    override fun getS(): DoubleMatrix {
        return EJMLDoubleMatrix(mySvd.getW(null))
    }

    override fun getSingularValues(): DoubleArray {
        return mySvd.singularValues
    }

    override fun getRank(): Int {
        return singularValues.filter { it >= tol }.count()
    }

}