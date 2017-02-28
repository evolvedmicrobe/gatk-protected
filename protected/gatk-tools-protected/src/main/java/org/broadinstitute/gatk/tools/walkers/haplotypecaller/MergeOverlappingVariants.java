package org.broadinstitute.gatk.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.gatk.utils.haplotype.Haplotype;

import java.util.*;

/**
 * Created by Nigel Delaney on 2/27/17.
 */
public class MergeOverlappingVariants extends MergeVariantsAcrossHaplotypes {

    // An accessory class for storing and sorting variants
    class SortableVariant implements Comparable<SortableVariant> {
        long start;
        long end;
        Haplotype haplotype;
        VariantContext variant;

        SortableVariant(Haplotype h, VariantContext vc) {
            haplotype = h;
            start = vc.getStart();
            end = vc.getEnd();
            // Note: Could consider not caching start/end in this class
            // and instead always get it from the variant context
            variant = vc;
        }

        boolean overlaps(long start, long end) {
           return this.start <= end && this.end >= start;
        }
        /* Compare first by start location, than longest location given a start. */
        @Override
        public int compareTo(SortableVariant sv) {
            if (start == sv.start) {
                if(this.end > sv.end) {
                    return -1;
                } else if (this.end == sv.end) {
                    return 0;
                } else {
                    return 1;
                }
            } else {
                return (int)(this.start - sv.start);
            }
        }
    }
    private final static Logger logger = Logger.getLogger(MergeOverlappingVariants.class);
    public static String MERGED_SOURCE_NAME = "merged";
    public MergeOverlappingVariants() {
        super();
    }

    /**
     * Merge as events among Haplotypes that overlap with each others start/end positions
     *
     * @param haplotypes a list of haplotypes whose events we want to merge
     * @param readLikelihoods map from sample name -> read likelihoods for each haplotype
     * @param startPosKeySet a set of starting positions of all events among the haplotypes
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     */
    @Override
    public boolean merge( final List<Haplotype> haplotypes,
                          final ReadLikelihoods<Haplotype> readLikelihoods,
                          final TreeSet<Integer> startPosKeySet,
                          final byte[] ref,
                          final GenomeLoc refLoc ) {
        if ( haplotypes == null ) throw new IllegalArgumentException("haplotypes cannot be null");
        if ( readLikelihoods == null ) throw new IllegalArgumentException("readLikelihoods cannot be null");
        if ( startPosKeySet == null ) throw new IllegalArgumentException("startPosKeySet cannot be null");
        if ( ref == null ) throw new IllegalArgumentException("ref cannot be null");
        if ( refLoc == null ) throw new IllegalArgumentException("refLoc cannot be null");
        if ( refLoc.size() != ref.length ) throw new IllegalArgumentException("refLoc size " + refLoc.size() + " != ref.length " + ref.length + " at " + refLoc);
        if( startPosKeySet.size() <= 1 ) { return false; }

        // To find variants in need of merging, we first load up all the variants and sort them
        List<SortableVariant> variants = new ArrayList<SortableVariant>(startPosKeySet.size() * 2);
        for(Haplotype h : haplotypes) {
            for (Map.Entry<Integer, VariantContext> entry : h.getEventMap().entrySet()) {
                VariantContext vc = entry.getValue();
                variants.add(new SortableVariant(h, entry.getValue()));
                if (vc.getStart() != entry.getKey()) { throw new IllegalArgumentException("Positions in mapping were not consistent."); }
            }
        }
        Collections.sort(variants);

        // Now pass through the sorted list and find any overlapping variants to merge
        long curStart = variants.get(0).start;
        long curEnd = variants.get(0).end;
        int intervalStart = 0;
        boolean mergeEventFound = false;
        boolean somethingWasMerged = false;
        for(int i = 1; i < variants.size(); i++) {
            SortableVariant sv = variants.get(i);
            if (sv.start == curStart) {continue;}
            else {
                if (sv.overlaps(curStart, curEnd)) {
                    mergeEventFound = true;
                    curEnd = Math.max(curEnd, sv.end);
                } else {
                    // No overlap, so we either combine the set of variants
                    // found or continue on if we didn't find any.
                    if(mergeEventFound) {
                        somethingWasMerged = true;
                        mergeRange(variants, intervalStart, i, curStart, curEnd, ref, refLoc);
                    }
                    // Now reset the whole process
                    curStart = sv.start;
                    curEnd = sv.end;
                    intervalStart = i;
                    mergeEventFound = false;
                }
            }
        }
        // Now handle a merge than went through end of list
        if (mergeEventFound) {
            somethingWasMerged = true;
            mergeRange(variants, intervalStart, variants.size(), curStart, curEnd, ref, refLoc);
        }
        return somethingWasMerged;
    }

    /**
     * Given a range of overlapping variants, go through and create a new all-encompassing variant by them.
     *
     * To accomplish this, we change all their start positions to the same start location, and pad the reference
     * and alt alleles to completely cover the range given, adding bases as needed.
     *
     * @param variants a list of sorted variants with a range to be merged
     * @param intervalStart start of range in list to merge (inclusive)
     * @param intervalStart end of range in list to merge (exclusive)
     * @param curStart start position of merged variant
     * @param curEnd end position of merged variant
     * @param ref the reference bases
     * @param refLoc the span of the reference bases
     */
    private void mergeRange(final List<SortableVariant> variants,
                            final int intervalStart,
                            final int intervalEnd,
                            final long curStart,
                            final long curEnd,
                            final byte[] ref,
                            final GenomeLoc refLoc ) {
        List<Allele> alleleSet = new ArrayList<>(2);
        // Note that for variant contexts, the end coordinate is inclusive!
        // Create the reference
        byte[] refSection = Arrays.copyOfRange(ref, (int) (curStart - refLoc.getStart()), (int) (curEnd - refLoc.getStart() + 1));
        alleleSet.add(Allele.create(refSection, true));

        // Now go through and update all the variants with incorrect contexts
        // TODO: Could cache this so that am not recreating for identical variants.
        for (int j = intervalStart; j < intervalEnd; j++) {
            // Let's pad the start/end to make  a new variant
            VariantContext vc = variants.get(j).variant;
            long leadNeeded = vc.getStart() - curStart;
            long tailNeeded = curEnd - vc.getEnd();
            // Should always be true
            if (!vc.isBiallelic()) {
                throw new IllegalArgumentException("Allele from haplotype eventmap had more than 2 types");
            }
            Allele altAllele = vc.getAlternateAllele(0);
            byte[] altbases = altAllele.getBases();
            byte[] altseq = new byte[(int) (leadNeeded + tailNeeded + altbases.length)];
            System.arraycopy(refSection, 0, altseq, 0, (int) leadNeeded);
            System.arraycopy(altbases, 0, altseq, (int) leadNeeded, altbases.length);
            System.arraycopy(refSection, (int) (refSection.length - tailNeeded), altseq, (int) (altbases.length + leadNeeded), (int) tailNeeded);
            alleleSet.add(Allele.create(altseq, false));

            // Now to replace variants in map
            VariantContext newVariant = (new VariantContextBuilder(MERGED_SOURCE_NAME, vc.getContig(), curStart, curEnd, alleleSet)).make();
            Haplotype h = variants.get(j).haplotype;
            h.getEventMap().remove(vc.getStart());
            h.getEventMap().addVC(newVariant, false);
            // Clear the allele we added
            alleleSet.remove(1);
        }
    }
}
