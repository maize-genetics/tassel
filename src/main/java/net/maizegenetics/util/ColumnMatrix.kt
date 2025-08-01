package net.maizegenetics.util

/**
 * @author Terry Casstevens
 * Created December 06, 2018
 *
 * This class optimizes read-only column access of a byte matrix.
 */

class ColumnMatrix private constructor(private val matrix: SuperByteMatrix) {

    val numColumns = matrix.numRows
    val numRows = matrix.numColumns

    fun column(index: Int) = matrix.getAllColumns(index)

    class Builder(numRows: Int, numColumns: Int) {

        private val matrix = SuperByteMatrixBuilder.getInstance(numColumns, numRows)

        fun set(row: Int, column: Int, value: Byte) = matrix.set(column, row, value)

        fun build() = ColumnMatrix(matrix)
    }

}
