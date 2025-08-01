package net.maizegenetics.matrixalgebra.Matrix

import net.maizegenetics.matrixalgebra.decomposition.*
import net.maizegenetics.taxa.distance.DistanceMatrix
import org.ejml.data.DMatrixRMaj
import org.ejml.dense.row.CommonOps_DDRM
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM

class EJMLDoubleMatrix: DoubleMatrix {
    val tol = 1e-10
    var myMatrix: DMatrixRMaj

    constructor(matrix: DMatrixRMaj) {
        myMatrix = matrix
    }

    constructor(row: Int, col: Int) {
        myMatrix = DMatrixRMaj(row, col)
    }

    constructor(row: Int, col: Int, values: DoubleArray) {
        myMatrix = DMatrixRMaj(row, col, true, *values)
    }

    constructor(row: Int, col: Int, value: Double) {
        myMatrix = DMatrixRMaj(row, col)
        myMatrix.fill(value)
    }

    constructor(values: Array<DoubleArray>) {
        myMatrix = DMatrixRMaj(values)
    }

    constructor(dist: DistanceMatrix) {
        myMatrix = DMatrixRMaj(dist.distances)
    }

    constructor(size: Int) {
        myMatrix = CommonOps_DDRM.identity(size)
    }

    constructor(diagonal: DoubleArray) {
        myMatrix = CommonOps_DDRM.diag(*diagonal)
    }

    override fun get(row: Int, col: Int): Double {
        return myMatrix.unsafe_get(row, col)
    }

    override fun getChecked(row: Int, col: Int): Double {
        return myMatrix[row,col]
    }

    override fun set(row: Int, col: Int, value: Double) {
        myMatrix.unsafe_set(row, col, value)
    }

    override fun setChecked(row: Int, col: Int, value: Double) {
        myMatrix[row, col] = value
    }

    override fun transpose(): DoubleMatrix {
        return EJMLDoubleMatrix(CommonOps_DDRM.transpose(myMatrix, DMatrixRMaj(myMatrix.numRows, myMatrix.numCols)))
    }

    override fun mult(dm: DoubleMatrix, transpose: Boolean, transposedm: Boolean): DoubleMatrix {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        val resultMatrix = DMatrixRMaj(1,1)
        when {
            transpose && transposedm -> CommonOps_DDRM.multTransAB(myMatrix, otherMatrix, resultMatrix)
            transpose && !transposedm -> CommonOps_DDRM.multTransA(myMatrix, otherMatrix, resultMatrix)
            !transpose && transposedm -> CommonOps_DDRM.multTransB(myMatrix, otherMatrix, resultMatrix)
            else -> CommonOps_DDRM.mult(myMatrix, otherMatrix, resultMatrix)
        }
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun mult(dm: DoubleMatrix): DoubleMatrix {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.mult(myMatrix, otherMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun multadd(
        A: DoubleMatrix,
        B: DoubleMatrix?,
        alpha: Double,
        beta: Double,
        transpose: Boolean,
        transposeA: Boolean,
    ): DoubleMatrix {
        //result = alpha * myMatrix * A + beta * B
        val Amatrix = (A as EJMLDoubleMatrix).myMatrix
        val Bmatrix = if (B == null) {
            when {
                transpose && transposeA -> DMatrixRMaj(myMatrix.numCols, Amatrix.numRows)
                transpose && !transposeA -> DMatrixRMaj(myMatrix.numCols, Amatrix.numCols)
                !transpose && transposeA -> DMatrixRMaj(myMatrix.numRows, Amatrix.numRows)
                else -> DMatrixRMaj(myMatrix.numRows, Amatrix.numCols)
            }
        }
        else if (beta == 1.0) (B as EJMLDoubleMatrix).myMatrix
        else {
            val scaledB = DMatrixRMaj(1,1)
            CommonOps_DDRM.scale(beta, (B as EJMLDoubleMatrix).myMatrix, scaledB)
            scaledB
        }

        when {
            transpose && transposeA -> CommonOps_DDRM.multAddTransAB(alpha, myMatrix, Amatrix, Bmatrix)
            transpose && !transposeA -> CommonOps_DDRM.multAddTransA(alpha, myMatrix, Amatrix, Bmatrix)
            !transpose && transposeA -> CommonOps_DDRM.multAddTransB(alpha, myMatrix, Amatrix, Bmatrix)
            else -> CommonOps_DDRM.multAdd(alpha, myMatrix, Amatrix, Bmatrix)
        }
        return EJMLDoubleMatrix(Bmatrix)
    }

    override fun crossproduct(): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.multTransA(myMatrix, myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun crossproduct(dm: DoubleMatrix): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.multTransA(myMatrix, (dm as EJMLDoubleMatrix).myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun tcrossproduct(): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.multTransB(myMatrix, myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun tcrossproduct(dm: DoubleMatrix): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.multTransB(myMatrix, (dm as EJMLDoubleMatrix).myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun concatenate(dm: DoubleMatrix, rows: Boolean): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        if (rows) CommonOps_DDRM.concatRows(myMatrix, (dm as EJMLDoubleMatrix).myMatrix, resultMatrix)
        else CommonOps_DDRM.concatColumns(myMatrix, (dm as EJMLDoubleMatrix).myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun inverse(): DoubleMatrix? {
        val resultMatrix = DMatrixRMaj(1,1)
        return if (CommonOps_DDRM.invert(myMatrix, resultMatrix)) EJMLDoubleMatrix(resultMatrix) else null
    }

    override fun invert(): Boolean {
        return CommonOps_DDRM.invert(myMatrix)
    }

    override fun generalizedInverse(): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(1,1)
        CommonOps_DDRM.pinv(myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun generalizedInverseWithRank(rank: IntArray): DoubleMatrix {
        val svd = EJMLSingularValueDecomposition(myMatrix)
        rank[0] = svd.rank

        //get S (W) and invert it
        val invS = svd.s
        for (ndx in 0 until invS.numberOfRows()) {
            invS[ndx,ndx] = if (invS[ndx,ndx] < tol) 0.0 else 1.0 / invS[ndx,ndx]
        }

        val V = svd.getV(false)
        val UT = svd.getU(true)
        return V.mult(invS).mult(UT)
    }

    override fun solve(Y: DoubleMatrix?): DoubleMatrix {
        //solve AX = B for X
        val solver = LinearSolverFactory_DDRM.leastSquares(myMatrix.numRows, myMatrix.numCols)
        val Amatrix = if (solver.modifiesA()) myMatrix.copy() else myMatrix
        check(solver.setA(Amatrix)) {"solver unable to set A"}
        val Bmatrix = if (solver.modifiesB()) (Y as EJMLDoubleMatrix).myMatrix.copy() else (Y as EJMLDoubleMatrix).myMatrix
        val Xmatrix = DMatrixRMaj(Amatrix.numCols, Bmatrix.numCols)
        solver.solve(Bmatrix, Xmatrix)
        return EJMLDoubleMatrix(Xmatrix)
    }

    override fun numberOfRows(): Int {
        return myMatrix.numRows
    }

    override fun numberOfColumns(): Int {
        return myMatrix.numCols
    }

    override fun row(i: Int): DoubleMatrix {
        val out = CommonOps_DDRM.extractRow(myMatrix, i, null)
        out.reshape(myMatrix.numCols, 1)
        return EJMLDoubleMatrix(out)
    }

    override fun column(j: Int): DoubleMatrix {
        return EJMLDoubleMatrix(CommonOps_DDRM.extractColumn(myMatrix, j, null))
    }

    override fun getXtXGM(): Array<DoubleMatrix> {
        val result0 = DMatrixRMaj(myMatrix.numCols, myMatrix.numCols) //X'X
        CommonOps_DDRM.multTransA(myMatrix, myMatrix, result0)

        //try inverting using faster code first before calculating pseudo inverse
        var result1 = result0.copy() //invX'X=G
        if (!CommonOps_DDRM.invert(result1)) {
            CommonOps_DDRM.pinv(result0, result1)
        }

        val tmp1 = DMatrixRMaj(myMatrix.numRows, result1.numCols)
        CommonOps_DDRM.mult(myMatrix, result1, tmp1) //tmp1 = XG
        val tmp2 = DMatrixRMaj(tmp1.numRows, myMatrix.numRows)
        CommonOps_DDRM.multTransB(tmp1,myMatrix,tmp2) //tmp2 = XGX'
        val result2 = CommonOps_DDRM.identity(tmp2.numRows)  //I - XGX'
        CommonOps_DDRM.subtractEquals(result2, tmp2)

        return arrayOf(EJMLDoubleMatrix(result0), EJMLDoubleMatrix(result1), EJMLDoubleMatrix(result2))
    }

    override fun copy(): DoubleMatrix {
        return EJMLDoubleMatrix(myMatrix.copy())
    }

    override fun getEigenvalueDecomposition(): EigenvalueDecomposition {
        return EJMLEigenvalueDecomposition(myMatrix)
    }

    override fun getSingularValueDecomposition(): SingularValueDecomposition {
        return EJMLSingularValueDecomposition(myMatrix)
    }

    override fun getQRDecomposition(): QRDecomposition {
        TODO("Not yet implemented")
    }

    override fun minus(dm: DoubleMatrix): DoubleMatrix {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        check(myMatrix.numRows == otherMatrix.numRows && myMatrix.numCols == otherMatrix.numCols) {"Attempted to subtract unequal size matrices"}
        val resultMatrix = DMatrixRMaj(myMatrix.numRows, myMatrix.numCols)
        CommonOps_DDRM.subtract(myMatrix, otherMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun minusEquals(dm: DoubleMatrix) {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        check(myMatrix.numRows == otherMatrix.numRows && myMatrix.numCols == otherMatrix.numCols) {"Attempted to subtract unequal size matrices"}
        CommonOps_DDRM.subtractEquals(myMatrix, otherMatrix)
    }

    override fun plus(dm: DoubleMatrix): DoubleMatrix {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        check(myMatrix.numRows == otherMatrix.numRows && myMatrix.numCols == otherMatrix.numCols) {"Attempted to add unequal size matrices"}
        val resultMatrix = DMatrixRMaj(myMatrix.numRows, myMatrix.numCols)
        CommonOps_DDRM.add(myMatrix, otherMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun plusEquals(dm: DoubleMatrix) {
        val otherMatrix = (dm as EJMLDoubleMatrix).myMatrix
        check(myMatrix.numRows == otherMatrix.numRows && myMatrix.numCols == otherMatrix.numCols) {"Attempted to add unequal size matrices"}
        CommonOps_DDRM.addEquals(myMatrix, otherMatrix)
    }

    override fun scalarAdd(s: Double): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(myMatrix.numRows, myMatrix.numCols)
        CommonOps_DDRM.add(myMatrix, s, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun scalarAddEquals(s: Double) {
        CommonOps_DDRM.add(myMatrix, s)
    }

    override fun scalarMult(s: Double): DoubleMatrix {
        val resultMatrix = DMatrixRMaj(myMatrix.numRows, myMatrix.numCols)
        CommonOps_DDRM.scale(s, myMatrix, resultMatrix)
        return EJMLDoubleMatrix(resultMatrix)
    }

    override fun scalarMultEquals(s: Double) {
        CommonOps_DDRM.scale(s, myMatrix)
    }

    override fun getSelection(rows: IntArray?, columns: IntArray?): DoubleMatrix {
        val selectedRows = if (rows == null) IntArray(myMatrix.numRows) { it } else rows
        val selectedColumns = if (columns == null) IntArray(myMatrix.numCols) { it } else columns
        return EJMLDoubleMatrix(CommonOps_DDRM.extract(myMatrix, selectedRows, selectedRows.size, selectedColumns, selectedColumns.size, null))
    }

    override fun rowSum(row: Int): Double {
        return (0 until myMatrix.numCols).map { myMatrix[row, it] }.sum()
    }

    override fun columnSum(column: Int): Double {
        return (0 until myMatrix.numRows).map { myMatrix[it, column] }.sum()
    }

    override fun columnRank(): Int {
        return EJMLSingularValueDecomposition(myMatrix).rank
    }

    override fun to1DArray(): DoubleArray {
        return myMatrix.data
    }

    override fun toArray(): Array<DoubleArray> {
        return Array(myMatrix.numRows) { rowIndex ->
            DoubleArray(myMatrix.numCols) { colIndex -> myMatrix[rowIndex, colIndex]}
        }
    }


}