/**
 * Data structures for holding raw sequence data (only Fastq currently), including Illumina, PacBio, Nanapore, etc.
 * This package provides methods for quality control and reformatting of sequence.
 *<p>
 * Definitions:
 * <li>Read:  a single sequence from a sequencer
 *<p>
 * All sequences are stored in string. Large chunk of reads do not go into memory
 */
package net.maizegenetics.dna.read;
