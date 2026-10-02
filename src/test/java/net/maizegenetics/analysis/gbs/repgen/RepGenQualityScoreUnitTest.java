/**
 *
 */
package net.maizegenetics.analysis.gbs.repgen;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * Pure-unit coverage for the averaged quality-score algorithm used by RepGen loading, extracted from
 * {@code RepGenLoadSeqToDBPluginTest} (which is excluded because its other tests require hardcoded
 * dev-machine paths and external data). This test is fully self-contained and runs in CI.
 *
 * @author lcj34 (original), extracted for CI unit coverage
 */
public class RepGenQualityScoreUnitTest {

    @Test
    public void testCreateQualityScore() {
        String qualSFFShort1 = "ABD;->DFABF-../9";
        String qualSFFShort2 = "BBD;->DEABF-../9";
        String qualSFFShort3 = "GB-;->DFABF-../G";
        int[] qualShortNum = new int[qualSFFShort1.length()];

        int qualityScoreBase = 33;

        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = (qualSFFShort1.charAt(idx) - qualityScoreBase);
        }
        StringBuilder sb = new StringBuilder();
        for (int idx = 0; idx < qualShortNum.length; idx++) {
            char ch = (char) (qualShortNum[idx] + qualityScoreBase);
            sb.append(ch);
        }
        // Round-tripping the numeric values back to characters must reproduce the original string.
        assertEquals(qualSFFShort1, sb.toString());

        // Triple each score, then divide by three: the recreated string is unchanged.
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] += (qualSFFShort1.charAt(idx) - qualityScoreBase);
            qualShortNum[idx] += (qualSFFShort1.charAt(idx) - qualityScoreBase);
        }
        sb.setLength(0);
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = qualShortNum[idx] / 3;
            char ch = (char) (qualShortNum[idx] + qualityScoreBase);
            sb.append(ch);
        }
        assertEquals("Averaging three identical quality strings must reproduce the original",
                qualSFFShort1, sb.toString());

        // Averaging three different strings must not throw and must yield a same-length result.
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] += (qualSFFShort2.charAt(idx) - qualityScoreBase);
            qualShortNum[idx] += (qualSFFShort3.charAt(idx) - qualityScoreBase);
        }
        sb.setLength(0);
        for (int idx = 0; idx < qualSFFShort1.length(); idx++) {
            qualShortNum[idx] = qualShortNum[idx] / 3;
            char ch = (char) (qualShortNum[idx] + qualityScoreBase);
            sb.append(ch);
        }
        assertEquals(qualSFFShort1.length(), sb.length());
    }
}
