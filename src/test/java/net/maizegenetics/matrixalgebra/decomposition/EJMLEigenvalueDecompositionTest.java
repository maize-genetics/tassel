package net.maizegenetics.matrixalgebra.decomposition;

import static org.junit.Assert.assertEquals;

import java.util.Random;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.EJMLDoubleMatrix;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.junit.Test;

public class EJMLEigenvalueDecompositionTest {

    /**
     * Builds ZKZ' for a random 3 x 3 kinship matrix K and an incidence matrix Z that assigns 100 taxa to 3 groups.
     * This is the form of the compressed kinship matrix that MLM decomposes. It has rank 3, and EJML's symmetric
     * QR algorithm does not converge on it unless the diagonal is shifted first (#39).
     */
    private static DMatrixRMaj compressedKinship() {
        int ntaxa = 100;
        int ngroups = 3;
        Random random = new Random(2);
        DMatrixRMaj B = new DMatrixRMaj(ngroups, ngroups);
        for (int i = 0; i < B.data.length; i++) B.data[i] = random.nextGaussian();
        DMatrixRMaj K = new DMatrixRMaj(ngroups, ngroups);
        CommonOps_DDRM.multTransB(B, B, K);

        DMatrixRMaj Z = new DMatrixRMaj(ntaxa, ngroups);
        for (int t = 0; t < ntaxa; t++) Z.set(t, t < ngroups ? t : random.nextInt(ngroups), 1.0);

        DMatrixRMaj ZK = new DMatrixRMaj(ntaxa, ngroups);
        CommonOps_DDRM.mult(Z, K, ZK);
        DMatrixRMaj ZKZ = new DMatrixRMaj(ntaxa, ntaxa);
        CommonOps_DDRM.multTransB(ZK, Z, ZKZ);
        return ZKZ;
    }

    @Test
    public void testRankDeficientMatrix() {
        DMatrixRMaj matrix = compressedKinship();
        int n = matrix.numRows;
        DoubleMatrix A = new EJMLDoubleMatrix(matrix);
        EigenvalueDecomposition eig = A.getEigenvalueDecomposition();
        double[] eigenvalues = eig.getEigenvalues();
        DoubleMatrix V = eig.getEigenvectors();
        assertEquals(n, eigenvalues.length);

        DoubleMatrix VtV = V.crossproduct();
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                assertEquals("V'V at " + r + ", " + c, r == c ? 1.0 : 0.0, VtV.get(r, c), 1e-10);
            }
        }

        double tol = 1e-10 * NormOps_DDRM.normF(matrix);
        DoubleMatrix AV = A.mult(V);
        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                assertEquals("AV at " + r + ", " + c, V.get(r, c) * eigenvalues[c], AV.get(r, c), tol);
            }
        }
    }
}
