package net.maizegenetics.analysis.popgen;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

/**
 * Container class for reporting LD results.
 * This class could still be expanded to include the two sites with various sort approaches.  The reason for this
 * functionality is so that classes can calculate LD for billions of site pairs, but only retain the most significant.
 *
 * @author Ed Buckler
 */
public final class LDResult {

    private final float r2;
    private final float dprime;
    private final float p;
    private final int n;
    private final int site1;
    private final int site2;

    public static Ordering<LDResult> byR2Ordering = new Ordering<LDResult>() {
        public int compare(LDResult left, LDResult right) {
            return ComparisonChain.start()
                    .compareTrueFirst(Float.isNaN(left.r2()),Float.isNaN(right.r2()))
                    .compare(left.r2(),right.r2())
                    .compare(left.n(),right.n())
                    .result();
        }
    };

    public static Ordering<LDResult> byPValueOrdering = new Ordering<LDResult>() {
        public int compare(LDResult left, LDResult right) {
            return ComparisonChain.start()
                    .compareTrueFirst(Float.isNaN(left.p()),Float.isNaN(right.p()))
                    .compare(right.p(),left.p())
                    .compare(left.r2(), right.r2())
                    .compare(left.n(),right.n())
                    .result();
        }
    };

    public LDResult(int site1, int site2, float r2, float dprime, float p, int n) {
        this.site1=site1;
        this.site2=site2;
        this.r2=r2;
        this.dprime=dprime;
        this.p=p;
        this.n=n;
    }

    public float r2() {
        return r2;
    }

    public float dPrime() {
        return dprime;
    }

    public float p() {
        return p;
    }

    public int n() {
        return n;
    }

    public int site1() {
        return site1;
    }

    public int site2() {
        return site2;
    }

    public static class Builder {
        private float r2=Float.NaN;
        private float dprime=Float.NaN;
        private float p=Float.NaN;
        private int n=0;
        private final int site1;
        private final int site2;

        public Builder(int site1, int site2) {
            this.site1=site1;
            this.site2=site2;
        }

        public LDResult build() {
            return new LDResult(site1,site2,r2, dprime, p, n);
        }

        public Builder r2(float value) {
            r2=value;
            return this;
        }

        public Builder dprime(float value) {
            dprime=value;
            return this;
        }

        public Builder p(float value) {
            p=value;
            return this;
        }

        public Builder n(int value) {
            n=value;
            return this;
        }

    }

    @Override
    public String toString() {
        return "LDResult{"+
                "site1="+site1+
                ", site2="+site2+
                ", r2="+r2+
                ", dprime="+dprime+
                ", p="+p+
                ", n="+n+
                '}';
    }
}
