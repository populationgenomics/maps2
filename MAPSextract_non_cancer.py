# LT 21/10 version 1
## including all the needed functions in same file
## runs very slow.

# version 2 
# extract controls only, as a proxy to subsampling

# v3
# extract subsets in a generic way using dictionary

import hail as hl;
from typing import Optional, Union;

## The functions below are from Konrad's original code
#
# handles filtering for callrate > 80%
#

AN_ADJ_FILTER = 0.8
data_type_sex_counts = {
    'exomes': {'female': 57787, 'male': 67961},
    'genomes': {'female': 6967, 'male': 8741}
}


def get_an_adj_criteria(mt, samples_by_sex: Optional[dict] = None, an_adj_proportion_cutoff: float = AN_ADJ_FILTER):
    if samples_by_sex is None:
        samples_by_sex = mt.aggregate_cols(hl.agg.counter(mt.meta.sex))
    return (hl.case()
            .when(mt.locus.in_autosome_or_par(), mt.freq[0].AN >= an_adj_proportion_cutoff * 2 * sum(samples_by_sex.values()))
            .when(mt.locus.in_x_nonpar(),
                  mt.freq[0].AN >= an_adj_proportion_cutoff * (samples_by_sex['male'] + samples_by_sex['female'] * 2))
            .when(mt.locus.in_y_nonpar(), mt.freq[0].AN >= an_adj_proportion_cutoff * samples_by_sex['male'])
            .or_missing())

## The functions below are from Konrad's original code
#
# TODO use the Broad Institute gnomAD utils (vep.py)
# https://github.com/broadinstitute/gnomad_methods
#

def filter_vep_to_canonical_transcripts(
    mt: Union[hl.MatrixTable, hl.Table], vep_root: str = "vep"
) -> Union[hl.MatrixTable, hl.Table]:
    """Filter VEP transcript consequences to those in the canonical transcript."""
    canonical = mt[vep_root].transcript_consequences.filter(
        lambda csq: csq.canonical == 1
    )
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return (
        mt.annotate_rows(**{vep_root: vep_data})
        if isinstance(mt, hl.MatrixTable)
        else mt.annotate(**{vep_root: vep_data})
    )

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant",
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
]

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

CSQ_ORDER = (
    CSQ_CODING_HIGH_IMPACT
    + CSQ_CODING_MEDIUM_IMPACT
    + CSQ_CODING_LOW_IMPACT
    + CSQ_NON_CODING
)

def get_worst_consequence_with_non_coding(ht):
    def get_worst_csq(csq_list: hl.expr.ArrayExpression, protein_coding: bool) -> hl.struct:
        lof = hl.null(hl.tstr)
        no_lof_flags = hl.null(hl.tbool)
        # lof_filters = hl.null(hl.tstr)
        # lof_flags = hl.null(hl.tstr)
        if protein_coding:
            all_lofs = csq_list.map(lambda x: x.lof)
            lof = hl.literal(['HC', 'OS', 'LC']).find(lambda x: all_lofs.contains(x))
            csq_list = hl.cond(hl.is_defined(lof), csq_list.filter(lambda x: x.lof == lof), csq_list)
            no_lof_flags = hl.or_missing(hl.is_defined(lof),
                                         csq_list.any(lambda x: (x.lof == lof) & hl.is_missing(x.lof_flags)))
            # lof_filters = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_filter).filter(lambda x: hl.is_defined(x))), '|')
            # lof_flags = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_flags).filter(lambda x: hl.is_defined(x))), '|')
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        worst_csq = hl.literal(CSQ_ORDER).find(lambda x: all_csq_terms.contains(x))
        return hl.struct(worst_csq=worst_csq, protein_coding=protein_coding, lof=lof, no_lof_flags=no_lof_flags,
                         # lof_filters=lof_filters, lof_flags=lof_flags
                         )

    protein_coding = ht.vep.transcript_consequences.filter(lambda x: x.biotype == 'protein_coding')
    return ht.annotate(**hl.case(missing_false=True)
                       .when(hl.len(protein_coding) > 0, get_worst_csq(protein_coding, True))
                       .when(hl.len(ht.vep.transcript_consequences) > 0, get_worst_csq(ht.vep.transcript_consequences, False))
                       .when(hl.len(ht.vep.regulatory_feature_consequences) > 0, get_worst_csq(ht.vep.regulatory_feature_consequences, False))
                       .when(hl.len(ht.vep.motif_feature_consequences) > 0, get_worst_csq(ht.vep.motif_feature_consequences, False))
                       .default(get_worst_csq(ht.vep.intergenic_consequences, False)))

## The functions below are from Konrad's original code
#
# handles the symetry of the strands, recoding the context so the focal nucleotide is either A or C
#

def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))

def collapse_strand(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    collapse_expr = {
        'ref': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.ref), ht.ref),
        'alt': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.alt), ht.alt),
        'context': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                           reverse_complement_bases(ht.context), ht.context),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T')
    }
    return ht.annotate(**collapse_expr) if isinstance(ht, hl.Table) else ht.annotate_rows(**collapse_expr)


#
# https://github.com/macarthur-lab/gnomad_lof/blob/49ebacbcb42d0bd14d8601761ce532433682e33f/constraint_utils/generic.py

def trimer_from_heptamer(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
    return t.annotate_rows(context=trimer_expr) if isinstance(t, hl.MatrixTable) else t.annotate(context=trimer_expr)

def annotate_variant_types(t: Union[hl.MatrixTable, hl.Table],
                           heptamers: bool = False) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (((t.ref == "A") & (t.alt == "G")) | ((t.ref == "G") & (t.alt == "A")) |
                       ((t.ref == "T") & (t.alt == "C")) | ((t.ref == "C") & (t.alt == "T")))
    cpg_expr = (((t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1:mid_index] == 'C')) |
                ((t.ref == "C") & (t.alt == "T") & (t.context[mid_index + 1:mid_index + 2] == 'G')))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (hl.case()
                         .when(t.cpg, 'CpG')
                         .when(t.transition, 'non-CpG transition')
                         .default('transversion'))
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)
    else:
        return t.annotate(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)

#
# https://github.com/macarthur-lab/gnomad_lof/blob/49ebacbcb42d0bd14d8601761ce532433682e33f/constraint_utils/constraint_basics.py
#

def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        'methylation_level': hl.case().when(
            ht.cpg & (ht.methylation.MEAN > 0.6), 2
        ).when(
            ht.cpg & (ht.methylation.MEAN > 0.2), 1
        ).default(0)
    }
    if annotate_coverage:
        annotation['exome_coverage'] = ht.coverage.exomes.median
    return ht.annotate(**annotation) if isinstance(ht, hl.Table) else ht.annotate_rows(**annotation)


# Locations of ressources

gnomad_ht_path = 'gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht/'

context_ht_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'
# this other location seems to contain the same file
# gs://gcp-public-data--gnomad/resources/context/grch37_context_vep_annotated.ht/

mutation_rate_ht_path = 'gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.ht'

context_ht = hl.read_table(context_ht_path)
mutation_ht = hl.read_table(mutation_rate_ht_path)

ht = hl.read_table(gnomad_ht_path)

ht = ht.filter((hl.len(ht.filters) == 0) & get_an_adj_criteria(ht, {'female': 57787, 'male':67961}))
ht = filter_vep_to_canonical_transcripts(ht)
ht = get_worst_consequence_with_non_coding(ht)

# for MAPS data extractionm filter on synonymous as early as possible for performance
ht = ht.filter(ht.worst_csq == 'synonymous_variant')

context = context_ht[ht.key]
ht = prepare_ht(ht.annotate(context=context.context, methylation=context.methylation), True, False)

print('Counting singletons')

# count variants
# 214 is index for gnomAD controls
# redo this in a more generic way
subset = 'non_cancer'
d = hl.eval(ht.globals.freq_index_dict)
isubset = d[subset]

print(f'in subset {subset}')

ht = ht.group_by('context', 'ref', 'alt', 'methylation_level', 'protein_coding', 'worst_csq').aggregate(variant_count=hl.agg.count_where(ht.freq[isubset].AC > 0), singleton_count=hl.agg.count_where(ht.freq[isubset].AC == 1))

# annotate with mutability and proportion of singleton 
ht = ht.annotate(mu=mutation_ht[
  hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt, methylation_level=ht.methylation_level)].mu_snp)


## calibrate maps on synonymous variants
syn_ps_ht = ht

# group by mutability
syn_ps_ht = syn_ps_ht.group_by(syn_ps_ht.mu).aggregate(singleton_count=hl.agg.sum(syn_ps_ht.singleton_count),
variant_count=hl.agg.sum(syn_ps_ht.variant_count))

# proportion of singleton per mutability
syn_ps_ht = syn_ps_ht.annotate(ps=syn_ps_ht.singleton_count / syn_ps_ht.variant_count)

syn_ps_ht.export(f'gs://20211021-maps/syn_ps_{subset}.tsv')