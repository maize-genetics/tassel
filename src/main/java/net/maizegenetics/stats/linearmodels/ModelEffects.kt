package net.maizegenetics.stats.linearmodels

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory

class NestedFactorModelEffect(val nestedEffect : FactorModelEffect, val outerEffect : FactorModelEffect, val id : Any) : ModelEffect {
    //so that the designMatrix is only calculated once
    var designMatrix : DoubleMatrix = DoubleMatrixFactory.DEFAULT.make(1,1,0.0);

    init {
        if (nestedEffect.size != outerEffect.size) throw IllegalArgumentException("nested and outer effects must have the same size (number of observations)")
    }

    override fun getID(): Any {
        return id
    }

    override fun setID(id: Any?) {
        throw NotImplementedError("setId is not implemented and must be set in the constructor.")
    }

    override fun getSize(): Int {
        return nestedEffect.getSize()
    }

    override fun getX(): DoubleMatrix {

        if (designMatrix.numberOfRows() == 1) {
            //do element-wise multiplication of columns of nestedEffect with columns of outerEffect
            val nestedDM = nestedEffect.x
            val outerDM = outerEffect.x
            val nrows = nestedDM.numberOfRows()
            val ncols = nestedDM.numberOfColumns() * outerDM.numberOfColumns();

            designMatrix = DoubleMatrixFactory.DEFAULT.make(nrows, ncols)

            var designCol = 0
            for (nestedCol in 0 until nestedDM.numberOfColumns()) {
                for (outerCol in 0 until outerDM.numberOfColumns()) {
                    for (row in 0 until nrows) {
                        designMatrix.set(row, designCol, nestedDM.get(row, nestedCol) * outerDM.get(row, outerCol))
                    }
                    designCol++
                }
            }

        }
        return designMatrix

    }

    override fun getXtX(): DoubleMatrix {
        return getX().crossproduct()
    }

    override fun getXty(y: DoubleArray?): DoubleMatrix {
        if ( y!= null) return getX().mult(DoubleMatrixFactory.DEFAULT.make(y.size, 1, y))
        else throw java.lang.IllegalArgumentException("null agrument to getXty not allowed")
    }

    override fun getyhat(beta: DoubleMatrix?): DoubleMatrix {
        throw NotImplementedError("getyhat not implemented for NestFactorModelEffects")
    }

    override fun getyhat(beta: DoubleArray?): DoubleMatrix {
        throw NotImplementedError("getyhat not implemented for NestFactorModelEffects")
    }

    override fun getLevelCounts(): IntArray {
        throw NotImplementedError("getyhat not implemented for NestFactorModelEffects")
    }

    override fun getNumberOfLevels(): Int {
        throw NotImplementedError("getLevelCounts not implemented for NestFactorModelEffects")
    }

    override fun getEffectSize(): Int {
        throw NotImplementedError("getEffectSize not implemented for NestFactorModelEffects")
    }

    override fun getCopy(): ModelEffect {
        throw NotImplementedError("getCopy not implemented for NestFactorModelEffects")
    }

    override fun getSubSample(sample: IntArray?): ModelEffect {
        throw NotImplementedError("getSubSample not implemented for NestFactorModelEffects")
    }
}