/*
* By downloading the PROGRAM you agree to the following terms of use:
* 
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
* 
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
* 
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
* 
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
* 
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
* 
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
* 
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
* 
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
* 
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
* 
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
* 
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.gatk.utils.IndexRange;
import org.broadinstitute.gatk.utils.MathUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Collection of allele counts for a genotype. It encompasses what alleles are present in the genotype and in what number.</p>
 *
 * <p>Alleles are represented herein by their indices running from <b>0</b> to <b>N-1</b> where <i>N</i> is the number of alleles.</p>
 *
 * <p>Each allele present in a genotype (count != 0) has a <i>rank</i>, that is the 0-based ordinal of
 * that allele amongst the ones present in the genotype as sorted by their index.</p>
 *
 * <p>For example:</p>
 *
 * <p><b>0/0/2/2</b> has two alleles with indices <b>0</b> and <b>2</b>, both with count 2.
 * The rank of <b>0</b> is <i>0</i> whereas the rank of <b>2</b> is <i>1</i>.</p>
 *
 * <p><b>2/4/4/7</b> has three alleles with indices <b>2</b>, <b>4</b> and <b>7</b>. <b>2</b> and <b>7</b> have count 1 whereas <b>4</b> has count 2.
 * The rank of <b>2</b> is <i>0</i>, the rank of <b>4</b> is <i>1</i>. and the rank of <b>7</b> is <i>2</i>.</p>
 *
 * <p>In contrast, in both examples above both <b>3</b> and <b>10</b> (and many others) are absent thus they have no rank (represented by <i>-1</i> whenever applies).</p>
 *
 * <p>{@link GenotypeAlleleCounts} instances have themselves their own index (returned by {@link #index() index()}, that indicate their 0-based ordinal within the possible genotype combinations with the same ploidy.</p>
 *
 * <p>For example, for ploidy 3:</p>
 *
 * <table>
 *     <th>Index</th><th>Genotype</th>
 *     <tr><td>0</td><td><b>0/0/0</b></td></tr>
 *     <tr><td>1</td><td><b>0/0/1</b></td></tr>
 *     <tr><td>2</td><td><b>0/1/1</b></td></tr>
 *     <tr><td>3</td><td><b>1/1/1</b></td></tr>
 *     <tr><td>4</td><td><b>0/0/2</b></td></tr>
 *     <tr><td>6</td><td><b>0/1/2</b></td></tr>
 *     <tr><td>7</td><td><b>1/1/2</b></td></tr>
 *     <tr><td>8</td><td><b>0/2/2</b></td></tr>
 *     <tr><td>9</td><td><b>1/2/2</b></td></tr>
 *     <tr><td>10</td><td><b>2/2/2</b></td></tr>
 *     <tr><td>11</td><td><b>0/0/3</b></td></tr>
 *     <tr><td>12</td><td><b>0/1/3</b></td></tr>
 *     <tr><td>13</td><td><b>1/1/3</b></td></tr>
 *     <tr><td>14</td><td><b>0/2/3</b></td></tr>
 *     <tr><td>15</td><td><b>1/2/3</b></td></tr>
 *     <tr><td>16</td><td><b>2/2/3</b></td></tr>
 *     <tr><td>17</td><td><b>0/3/3</b></td></tr>
 *     <tr><td>...</td><td>...</td></tr>
 * </table>
 *
 * The total number of possible genotypes is only bounded by the maximum allele index.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class GenotypeAlleleCounts implements Comparable<GenotypeAlleleCounts>, Cloneable {

    private static final double UNCOMPUTED_LOG_10_COMBINATION_COUNT = -1;

    /**
     * The log10 number of phased genotypes corresponding to this unphased genotype.  For example,
     * [0, 1, 1, 1] = AB:  log10(2)
     * [0, 2] = AA:  log10(1)
     * [0, 1, 1, 1, 2, 1] = ABC: log10(6)
     * [0, 2, 1, 2] = AABB: log10(4!/(2!2!))
     * This is evaluated lazily i.e. it is initialized to {@link GenotypeAlleleCounts::UNCOMPUTED_LOG_10_COMBINATION_COUNT}
     * and only calculated if its getter is invoked.
     */
    private double log10CombinationCount = UNCOMPUTED_LOG_10_COMBINATION_COUNT;

    /**
     * The ploidy of the genotype.
     */
    private final int ploidy;

    /**
     * Sorted array of integer pairs as described in {@link #GenotypeAlleleCounts(int, int, int...)}.
     */
    private int[] sortedAlleleCounts;

    /**
     * Number of different alleles in the genotype.
     */
    private int distinctAlleleCount;

    /**
     * Index of this genotype within genotypes of the same ploidy.
     */
    private int index;

    /**
     * Creates a new unphased genotype.
     *
     * <p>This method assumes that the invoker is passing a well formatted and sorted allele frequency array.
     * Not checks are done for the sake of performance.</p>
     *
     * <p>
     *     The input argument {@code sortedAlleleCounts} list the index of alleles included in the unphased genotype
     *     and their frequency in the genotype in a single array using consecutive pairs:<br/>
     *
     *     <pre> [allele_1,freq_1,allele_2,freq_2, ... , allele_i, freq_i, ... , allele_n, freq_n]</pre>
     *
     *     <br/>
     *     No entry can have frequency == 0 (these must be omitted) and entries are sorted by allele index without
     *     any repetitions so that if <i>i < j</i> then <i>allele_i < allele_j</i>.
     *
     * </p>
     *
     * <p>
     *     The {@code ploidy} provided must be equal to the sum of all frequencies in {@code sortedAlleleCounts}
     * </p>
     * @param ploidy the genotype ploidy.
     * @param sortedAlleleCounts sorted allele counts following the restrictions above.
     * @param index the genotype index.
     */
    private GenotypeAlleleCounts(final int ploidy, final int index, final int... sortedAlleleCounts) {
        this(ploidy, index, sortedAlleleCounts, sortedAlleleCounts.length >> 1);
    }

    private GenotypeAlleleCounts(final int ploidy, final int index, final int[] sortedAlleleCounts, final int distinctAlleleCount){
        this.ploidy = ploidy;
        this.index = index;
        this.sortedAlleleCounts = sortedAlleleCounts;
        this.distinctAlleleCount = distinctAlleleCount;
    }

    /**
     * Gets the log10 combination count, computing it if uninitialized.  Note that the invoked MathUtils method uses fast cached
     * log10 values of integers for any reasonable ploidy.
     *
     * This method should be invoked on instances of {@link GenotypeAlleleCounts} cached in {@link GenotypeLikelihoodCalculators::genotypeTableByPloidy}.
     * Such usage allows the result of this computation to be cached once for an entire run of HaplotypeCaller.
     * @return
     */
    public double log10CombinationCount() {
        if (log10CombinationCount == UNCOMPUTED_LOG_10_COMBINATION_COUNT) {
            log10CombinationCount = MathUtils.log10Factorial(ploidy)
                    - new IndexRange(0, distinctAlleleCount).sum(n -> MathUtils.log10Factorial(sortedAlleleCounts[2*n+1]));
        }
        return log10CombinationCount;
    }

    /**
     * Returns the genotype's ploidy.
     * @return 0 or greater.
     */
    public int ploidy() {
        return ploidy;
    }

    /**
     * Increases the allele counts a number of times.
     *
     * <p>
     *     This method must not be invoked on cached genotype-allele-counts that are meant to remain constant,
     *     such as the ones contained in {@link GenotypeLikelihoodCalculators#genotypeTableByPloidy}.
     * </p>
     *
     * @param times the number of times to increase.
     *
     * @throws IllegalArgumentException if {@code times} is negative.
     */
    protected void increase(final int times) {
        for (int i = 0; i < times; i++)
            increase();
    }

    /**
     * Updates the genotype counts to match the next genotype.
     *
     * <p>
     *     This method must not be invoked on cached genotype-allele-counts that are meant to remain constant,
     *     such as the ones contained in {@link GenotypeLikelihoodCalculators#genotypeTableByPloidy}
     * </p>
     */
    protected void increase() {
        // if the ploidy is zero there is only one possible genotype.
        if (distinctAlleleCount == 0)
            return;

        // Worth make this case faster.
        if (distinctAlleleCount == 1) {
            if (ploidy == 1) {
                sortedAlleleCounts[0]++;
            } else {
                if (sortedAlleleCounts.length < 4)
                    sortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts,4);
                sortedAlleleCounts[2] = sortedAlleleCounts[0] + 1;
                sortedAlleleCounts[3] = 1;
                sortedAlleleCounts[0] = 0;
                sortedAlleleCounts[1] = ploidy - 1;
                distinctAlleleCount = 2;
            }
        } else {
            // Now, all the following ifs are just the way to avoid working with dynamically sizing List<int[]>
            // as the final size of the resulting new sorted-allele-counts array varies depending on the situation.
            // this is considerably faster and the logic complexity would not be that different actually so it is worth
            // the if indentations.
            //
            // Notice that at this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
            //
            // We only need to look at the two lowest allele indices to decide what to do.

            final int allele0 = sortedAlleleCounts[0];
            final int freq0 = sortedAlleleCounts[1];
            final int allele1 = sortedAlleleCounts[2];
            final int allele0Plus1 = allele0 + 1;
            final boolean allele0And1AreConsecutive = allele0Plus1 == allele1;
            final int[] newSortedAlleleCounts;
            // The rest of the sorted allele counts array contains junk
            final int sortedAlleleCountsLength = distinctAlleleCount << 1;


            if (freq0 == 1) {   // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
                if (allele0And1AreConsecutive) {  // need just to remove the first allele and add 1 to the frequency of the second (freq1 += 1).
                    System.arraycopy(sortedAlleleCounts, 2, sortedAlleleCounts, 0, sortedAlleleCountsLength - 2); // shift left the first component away.
                    sortedAlleleCounts[1]++; // freq1 has become freq0.
                    distinctAlleleCount--;
                } else  // just need to mutate allele0 to allele0 + 1.
                    sortedAlleleCounts[0] = allele0Plus1;
            } else { // && freq0 > 1 as per sortedAlleleCounts format restrictions. In this case allele0 will mutated to '0' with frequency decreased by 1.
                if (allele0And1AreConsecutive) { // we don't need to add a component for allele0 + 1 since it already exists.
                    sortedAlleleCounts[0] = 0;
                    sortedAlleleCounts[1] = freq0 - 1;
                    sortedAlleleCounts[3]++;
                } else { // we need to insert allele0 + 1 in the sorted-allele-counts array and give it frequency 1.
                    if (sortedAlleleCounts.length < sortedAlleleCountsLength + 2) // make room for the new component.
                        sortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts,sortedAlleleCountsLength + 2);
                    System.arraycopy(sortedAlleleCounts, 2, sortedAlleleCounts, 4, sortedAlleleCountsLength - 2);
                    sortedAlleleCounts[0] = 0;
                    sortedAlleleCounts[1] = freq0 - 1;
                    sortedAlleleCounts[2] = allele0Plus1;
                    sortedAlleleCounts[3] = 1;
                    distinctAlleleCount++;
                }
            }
        }
        index++;
        log10CombinationCount = -1;
    }

    /**
     * Calculates the next genotype in likelihood indexing order.
     * @return never null.
     */
    protected GenotypeAlleleCounts next() {
        // if the ploidy is zero there is only one possible genotype.
        if (distinctAlleleCount == 0)
            return this;

        // Worth make this case faster.
        if (distinctAlleleCount == 1) {
            if (ploidy == 1)  // A -> B , D -> E etc...
                return new GenotypeAlleleCounts(1, index + 1, sortedAlleleCounts[0] + 1, 1);
            else  // AAAAA -> AAAAB, DDD -> AAE etc...
                return new GenotypeAlleleCounts(ploidy, index + 1, 0, ploidy - 1, sortedAlleleCounts[0] + 1, 1);
        }
        // Now, all the following ifs are just the way to avoid working with dynamically sizing List<int[]>
        // as the final size of the resulting new sorted-allele-counts array varies depending on the situation.
        // this is considerably faster and the logic complexity would not be that different actually so it is worth
        // the if indentations.
        //
        // Notice that at this point distinctAlleleCount >= 2 thus sortedAlleleCounts.length >= 4.
        //
        // We only need to look at the two lowest allele indices to decide what to do.

        final int allele0 = sortedAlleleCounts[0];
        final int freq0 = sortedAlleleCounts[1];
        final int allele1 = sortedAlleleCounts[2];
        final int allele0Plus1 = allele0 + 1;
        final boolean allele0And1AreConsecutive = allele0Plus1 == allele1;
        final int[] newSortedAlleleCounts;
        // The rest of the sorted allele counts array contains junk
        final int sortedAlleleCountsLength = distinctAlleleCount << 1;

        if (freq0 == 1) {   // in this case allele0 wont be present in the result and all is frequency should go to allele0 + 1.
            if (allele0And1AreConsecutive) {  // need just to remove the first allele and 1 to the frequency of the second (freq1 += 1).
                newSortedAlleleCounts = Arrays.copyOfRange(sortedAlleleCounts,2,sortedAlleleCountsLength);
                newSortedAlleleCounts[1]++;
            } else {  // just need to mutate allele0 to allele0 + 1.
                newSortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts,sortedAlleleCountsLength);
                newSortedAlleleCounts[0] = allele0Plus1;
                // newSortedAlleleCounts[1] = 1; // :) no need to do it because it is already the case (freq0 == 1).
            }
        } else { // && freq0 > 1 as per sortedAlleleCounts format restrictions. In this case allele0 will muttated to '0' with frequency decreased by 1.
            if (allele0And1AreConsecutive) { // we don't need to add a component for allele0 + 1 since it already exists.
                newSortedAlleleCounts = sortedAlleleCounts.clone();
                newSortedAlleleCounts[0] = 0;
                newSortedAlleleCounts[1] = freq0 - 1;
                newSortedAlleleCounts[3]++;
            } else { // we need to insert allele0 + 1 in the sorted-allele-counts array.
                newSortedAlleleCounts = new int[sortedAlleleCountsLength + 2];
                newSortedAlleleCounts[0] = 0;
                newSortedAlleleCounts[1] = freq0 - 1;
                newSortedAlleleCounts[2] = allele0Plus1;
                newSortedAlleleCounts[3]++; // = 1 as the array was freshly created with 0s.
                System.arraycopy(sortedAlleleCounts,2,newSortedAlleleCounts,4,sortedAlleleCountsLength - 2);
            }
        }
        return new GenotypeAlleleCounts(ploidy, index + 1, newSortedAlleleCounts);
    }

    /**
     * Returns the number of different alleles that participate in the genotype.
     *
     * @return 0 or greater.
     */
    public int distinctAlleleCount() {
        return distinctAlleleCount;
    }

    /**
     * Returns the index of the allele from its rank in the genotype.
     *
     * @param rank the query rank.
     *
     * @throws IllegalArgumentException if the {@code rank} provided is outside the valid range [0,{@link #distinctAlleleCount()}).
     *
     * @return 0 or greater.
     */
    public int alleleIndexAt(final int rank) {
        if (rank < 0 || rank >= distinctAlleleCount)
            throw new IllegalArgumentException("the requested rank " + rank + " is out of range [0," + distinctAlleleCount + ")");
        return sortedAlleleCounts[rank << 1];
    }

    /**
     * Returns the rank of an allele in the genotype by its index.
     *
     * @param index the target index.
     *
     * @throws IllegalArgumentException if {@code index} is less that 0. Indices can be arbitrarily large.
     *
     * @return -1 or less if the allele index is not present in the genotype, 0 to {@link #distinctAlleleCount()} - 1 otherwise.
     *   If negative, the absolute value can be used to determine where would be that index inserted within {@code [0,{@link #distinctAlleleCount()}]} as
     *   {@code - result - 1}.
     *
     */
    public int alleleRankFor(final int index) {
        if (index < 0)
            throw new IllegalArgumentException("the index must be 0 or greater");
        return alleleIndexToRank(index, 0, distinctAlleleCount);
    }

    /**
     * Generates a string that would represent the unphased genotype with this allele counts.
     *
     * <p>
     *     In this string allele calls appear in alleleIndex order with as many repeats as copies of each allele. So
     *     for example:<br/>
     *     <pre>
     *         0         # haploid reference.
     *         0/0       # typical diploid calls
     *         0/1
     *         1/1
     *         0/0/1/3/3 # pentaploid with to ref, one first alt. and 2 third alt. allele
     *     </pre>
     *
     * </p>
     *
     * @return never {@code null}.
     */
    public String toUnphasedGenotypeString() {
        if (ploidy == 0) return "";
        final StringBuilder sb = new StringBuilder(distinctAlleleCount * 3);
        for (int i = 0; i < distinctAlleleCount; i += 2) {
            final int alleleIndex = sortedAlleleCounts[i];
            final int alleleCount = sortedAlleleCounts[i + 1];
            for (int j = 0; j < alleleCount; j++)
                sb.append(alleleIndex).append('/');

        }
        sb.setLength(sb.length() - 1);
        return sb.toString();
    }

    @Override
    public String toString() {
        // Perhaps we should change in the future, but the unphased genotype representation seems to be
        // a good one.
        return toUnphasedGenotypeString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(final Object o) {
        if (o instanceof GenotypeAlleleCounts)
            return equals((GenotypeAlleleCounts)o);
        else
            return false;
    }

    /**
     * Compares with another genotype.
     * @param o the other genotype.
     * @return never {@code null}.
     */
    public boolean equals(final GenotypeAlleleCounts o) {
        if (o == this)
            return true;
        if (o == null)
            return false;
        if (ploidy != o.ploidy)
            return false;
        return Arrays.equals(sortedAlleleCounts, o.sortedAlleleCounts);
    }

    /**
     * Returns the index of this genotype allele count within all possible genotypes with the same ploidy.
     *
     * @return 0 or greater.
     */
    public int index() {
        return index;
    }

    /**
     * Compares to genotypes.
     *
     * <p>A genotype with larger ploidy is considered greater than one with a lower ploidy. If both genotypes have
     * the same ploidy, then the genotype with the largest allele index or largest count if these are the same</p>.
     *
     * @param other genotype to compare to.
     *
     * @throws IllegalArgumentException if {@code other} is {@code null}.
     *
     * @return 0 if both genotypes are equivalent, < 0 if this genotype is less than {@code other} and > 0
     * if this genotype is greater than {@code other}.
     */
    @Override
    public int compareTo(final GenotypeAlleleCounts other) {
        if (other == this)
            return 0;
        if (other == null)
            throw new IllegalArgumentException("input genotype cannot be null");
        if (other.ploidy == ploidy)
            return index - other.index;
        else
            return ploidy - other.ploidy;
    }

    @Override
    public int hashCode() {
        return ((31 + ploidy) * 31 ) + index;
    }

    /**
     * Implements binary search across allele indexes.
     * @param index the target index.
     * @param from first inclusive possible rank.
     * @param to last exclusive possible rank.
     * @return -1 or less if the allele index is not in the genotype false otherwise. You can obtain
     *  the potential insertion point (within the interval [from,to]) as {@code -result - 1}
     */
    private int alleleIndexToRank(final int index,final int from, final int to) {
        if (to <= from)
            return - from - 1;
        if (from == to - 1) {
            final int onlyIndex = sortedAlleleCounts[from << 1];
            return onlyIndex == index ? from : (onlyIndex > index) ? -from - 1 : -to - 1;
        }

        final int mid = (to + from) >> 1;
        final int midIndex = sortedAlleleCounts[mid << 1];
        if (midIndex == index)
            return mid;
        else if (midIndex < index)
            return alleleIndexToRank(index,mid + 1,to);
        else
            return alleleIndexToRank(index,0,mid);
    }

    /**
     * Returns the count of an allele in the genotype given is rank in the genotype (not the allele index itself).
     *
     * @param rank of the requested allele within the genotype.
     *
     * @throws IllegalArgumentException if {@code rank} is out the the valid range [0,{@link #distinctAlleleCount})
     *
     * @return 1 or greater.
     */
    public int alleleCountAt(final int rank) {
        if (rank < 0 || rank >= distinctAlleleCount)
            throw new IllegalArgumentException("the rank is out of range");
        return sortedAlleleCounts[(rank << 1) + 1];
    }

    /**
     * Checks whether this genotype contain at least one call on a particular allele index.
     *
     * @param index the target allele.
     *
     * @throws IllegalArgumentException if {@code index} is negative.
     *
     * @return {@code true} iff the genotype contains that allele index.
     */
    public boolean containsAllele(final int index) {
        return alleleRankFor(index) >= 0;
    }

    /**
     * Returns the count of an allele in the genotype given it index.
     *
     * @return 0 if the allele is not present in the genotype, 1 or more otherwise.
     */
    public int alleleCountFor(final int index) {
        final int rank = alleleRankFor(index);
        return rank < 0 ? 0 : alleleCountAt(rank);
    }

    /**
     * Returns the allele counts for each allele index to maximum.
     * @param maximumAlleleIndex the maximum allele index required.
     * @throws IllegalArgumentException if {@code maximumAlleleIndex} is less than 0.
     * @return never {@code null}, an array of exactly {@code maximumAlleleIndex + 1} positions with the counts
     * of each allele where the position in the array is equal to its index.
     */
    public int[] alleleCountsByIndex(final int maximumAlleleIndex) {
        if (maximumAlleleIndex < 0)
            throw new IllegalArgumentException("the requested allele count cannot be less than 0");
        final int[] result = new int[maximumAlleleIndex + 1];
        copyAlleleCountsByIndex(result, 0, 0, maximumAlleleIndex);
        return result;
    }


    private void copyAlleleCountsByIndex(final int[] dest, final int offset, final int minimumAlleleIndex, final int maximumAlleleIndex) {

        // First we determine what section of the sortedAlleleCounts array contains the counts of interest,
        // By the present allele rank range of interest.
        final int minimumAlleleRank = alleleRankFor(minimumAlleleIndex);
        final int maximumAlleleRank = alleleRankFor(maximumAlleleIndex);

        // If the min or max allele index are absent (returned rank < 0) we note where the would be inserted; that
        // way we avoid going through the rest of positions in the sortedAlleleCounts array.
        // The range of interest is then [startRank,endRank].
        final int startRank = minimumAlleleRank < 0 ? - minimumAlleleRank - 1 : minimumAlleleRank;
        final int endRank = maximumAlleleRank < 0 ? - maximumAlleleRank - 2 : maximumAlleleRank;

        // Iteration variables:
        int nextIndex = minimumAlleleIndex; // next index that we want to output the count for.
        int nextRank = startRank; // next rank to query in sortedAlleleCounts.
        int nextSortedAlleleCountsOffset = nextRank << 1; // offset in sortedAlleleCounts where the info is present for the next rank.
        int nextDestOffset = offset; // next offset in destination array where to set the count for the nextIndex.

        while (nextRank++ <= endRank) {
            final int alleleIndex = sortedAlleleCounts[nextSortedAlleleCountsOffset++];
            // fill non-present allele counts with 0s.
            while (alleleIndex > nextIndex) {
                dest[nextDestOffset++] = 0;
                nextIndex++;
            }
            // It is guaranteed that at this point alleleIndex == nextIndex
            // thanks to the condition of the enclosing while: there must be at least one index of interest that
            // is present in the remaining (nextRank,endRank] interval as otherwise endRank would be less than nextRank.
            dest[nextDestOffset++] = sortedAlleleCounts[nextSortedAlleleCountsOffset++];
            nextIndex++;
        }
        // Finally we take care of trailing requested allele indices.
        while (nextIndex++ <= maximumAlleleIndex)
            dest[nextDestOffset++] = 0;
    }

    /**
     * Copies the sorted allele counts into an array.
     *
     * <p>
     *     Sorted allele counts are disposed as an even-sized array where even positions indicate the allele index and
     *     the following odd positions the number of copies of that allele in this genotype allele count:
     * </p>
     * <p><pre>
     *     [ allele_0, freq_0, allele_1, freq_1 ... ]
     * </pre></p>
     *
     * <p>
     *     With {@code offset} you can indicate an alternative first position in the destination array.
     * </p>
     *
     * @param dest where to copy the counts.
     * @param offset starting position.
     *
     * @throws IllegalArgumentException if {@code dest} is {@code null}, {@code offset} is less than 0
     *   or {@code dest} is not large enough considering the number of alleles present in this genotype
     *   allele counts and the {@code offset} provided. A total of
     *   <code>{@link #distinctAlleleCount()} * 2 positions</code>
     *   are required for the job.
     */
    public void copyAlleleCounts(final int[] dest, final int offset) {
        if (dest == null)
            throw new IllegalArgumentException("the destination cannot be null");
        if (offset < 0)
            throw new IllegalArgumentException("the offset cannot be negative");
        final int sortedAlleleCountsLength = distinctAlleleCount << 1;
        if (offset + sortedAlleleCountsLength > dest.length)
            throw new IllegalArgumentException("the input array does not have enough capacity");
        System.arraycopy(sortedAlleleCounts,0,dest,offset,sortedAlleleCountsLength);
    }

    /**
     * Instantiates the first genotype possible provided a total ploidy.
     * @param ploidy the ploidy of the genotype.
     *
     * @throws java.lang.IllegalArgumentException if ploidy is less than 0.
     *
     * @return never {@code null}.
     */
    protected static GenotypeAlleleCounts first(final int ploidy) {
        if (ploidy < 0)
            throw new IllegalArgumentException("the ploidy must be 0 or greater");
        else if (ploidy == 0)
            return new GenotypeAlleleCounts(0,0);
        else
            return new GenotypeAlleleCounts(ploidy, 0, 0, ploidy);
    }

    /**
     * Makes the next genotype in likelihood indexing order.
     *
     * @param g the original genotype.
     *
     * @throws IllegalArgumentException if {@code g} is {@code null}.
     *
     * @return never {@code null}.
     */
    public static GenotypeAlleleCounts makeNextGenotype(final GenotypeAlleleCounts g) {
        if (g == null)
            throw new IllegalArgumentException("the next genotype");
        return g.next();
    }

    /**
     * Returns the largest allele index present in the genotype.
     *
     * @return -1 if there is no alleles (ploidy == 0), 0 or greater otherwise.
     */
    public int maximumAlleleIndex() {
        if (distinctAlleleCount == 0)
            return -1;
        else
            return sortedAlleleCounts[(distinctAlleleCount - 1) << 1];
    }

    /**
     * Returns the smallest allele index present in the genotype.
     *
     * @return -1 if there is no allele (ploidy == 0), 0 or greater otherwise.
     */
    public int minimumAlleleIndex() {
        if (distinctAlleleCount == 0)
            return -1;
        else
            return sortedAlleleCounts[0];
    }

    /**
     * Creates an independent copy of this genotype.
     * @return never {@code null}.
     */
    @Override
    protected GenotypeAlleleCounts clone() {
        final GenotypeAlleleCounts result;
        try {
            result = (GenotypeAlleleCounts) super.clone();
        } catch (final CloneNotSupportedException e) {
            throw new IllegalStateException(e);
        }
        result.sortedAlleleCounts = Arrays.copyOf(sortedAlleleCounts,distinctAlleleCount << 1);
        return result;
    }

    /**
     * Composes a list with the alleles.
     *
     * @param allelesToUse alleles to use.
     *
     * @throws IllegalArgumentException if {@code allelesToUse} is {@code null},
     *          or does not contain enough elements to accommodate the maximum allele index in this allele-counts.
     *
     * @return never null, but it might be restricted (unmodifiable or non-expandable).
     */
    public <T extends Allele> List<T> asAlleleList(final List<T> allelesToUse) {
        if (allelesToUse == null)
            throw new IllegalArgumentException("the input allele list cannot be null");
        if (allelesToUse.size() < maximumAlleleIndex())
            throw new IllegalArgumentException("the provided alleles to use does not contain an element for the maximum allele ");
        if (distinctAlleleCount == 1 ) {
            if (ploidy == 1)
                return Collections.singletonList(allelesToUse.get(sortedAlleleCounts[0]));
            else
                return Collections.nCopies(ploidy,allelesToUse.get(sortedAlleleCounts[0]));
        } else {
            final Allele[] myAlleles = new Allele[ploidy];
            int nextIndex = 0;
            for (int i = 0, ii = 0; i < distinctAlleleCount; i++) {
                final Allele allele = allelesToUse.get(sortedAlleleCounts[ii++]);
                final int repeats = sortedAlleleCounts[ii++];
                for (int j = 0; j < repeats; j++)
                    myAlleles[nextIndex++] = allele;
            }
            return (List<T>) Arrays.asList(myAlleles);
        }
    }

    /**
     * Returns an array with the allele indices repeated based on the number of occurrences in the genotype.
     *
     * <p>
     *     indices are sorted from the smallest to the greatest.
     * </p>
     *
     * <p>
     *     If a sufficiently large array is provided as {@code dest}, this is used as the destination. Unnecessary
     *     positions at the back of the array are left untouched.
     * </p>
     *
     * <p>
     *     However if {@code dest} is {@code null} or it does not have enough space, a new array of with length equal to
     *     the ploidy will be used and returned instead.
     * </p>
     *
     * @param dest destination array. Can be {@code null} or not have sufficient positions (ploidy); in that case a new
     *             one is created.
     * @return never {@code null}, {@code dest} if sufficiently large otherwise an array of ploidy positions.
     */
    public int[] toAlleleIndicesArray(final int[] dest) {
        final int[] result = dest == null || dest.length < ploidy ? new int[ploidy] : dest;
        int k = 0;
        for (int i = 0,ii = 0; i < distinctAlleleCount; i++) {
            final int index = sortedAlleleCounts[ii++];
            final int repeats = sortedAlleleCounts[ii++];
            for (int j = 0; j < repeats; j++)
                result[k++] = index;
        }
        return result;
    }

    @FunctionalInterface
    public interface IntBiConsumer {
        void accept(final int alleleIndex, final int alleleCount);
    }

    @FunctionalInterface
    public interface IntToDoubleBiFunction {
        double apply(final int alleleIndex, final int alleleCount);
    }

    public void forEachAlleleIndexAndCount(final IntBiConsumer action) {
        new IndexRange(0, distinctAlleleCount).forEach(n -> action.accept(sortedAlleleCounts[2*n], sortedAlleleCounts[2*n+1]));
    }

    public double sumOverAlleleIndicesAndCounts(final IntToDoubleBiFunction func) {
        return new IndexRange(0, distinctAlleleCount).sum(n -> func.apply(sortedAlleleCounts[2*n], sortedAlleleCounts[2*n+1]));
    }
}
