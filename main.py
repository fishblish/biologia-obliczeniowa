import argparse
import os
import firstmodule as fm

base_path = '/home/julia/Desktop/uni/projekt_warianty/all/brain_reg_var-main/'
programs_path = '/home/julia/Desktop/uni/projekt_warianty/programs'
data_path = base_path + '/data'

#połączyć dane hic z plikiem gtf aby uzyskać dane jak teraz o kontaktach
parser = argparse.ArgumentParser()
parser.add_argument("--gatk_path", type=str, nargs='?', default= programs_path + \
                    "/gatk/gatk-4.4.0.0/gatk", help="Path to directory containing GATK program.")
parser.add_argument("--annovar_path", type=str, nargs='?', default= programs_path + \
                    "/annovar/", help="Path to directory containing ANNOVAR program.")
parser.add_argument("--input_vcf", type=str, nargs='?', default= data_path + "/test_variants_chr16.vcf.gz", \
                    help="Path to directory containing input VCF file. We provide default example file, but it sholud be replaced by your file")
parser.add_argument("--promoter_regions", type=str, nargs='?', default= data_path + "/brain_promoters_active.bed", \
                    help="Path to directory containing bed file with promoter regions. Last column should contain gene names, comma separated if promoters of several genes overlap.")
parser.add_argument("--enhancer_regions", type=str, nargs='?', default= data_path + "/brain_enhancers_active.bed", \
                    help="Path to directory containing bed file with enhancer regions. It should have 4 columns: chr, start, end, gene.\
                        Last column should contain gene names if enhancer is located inside gene, \
                        comma separated, '.' for intergenic enhancers = no gene overlaps.")
parser.add_argument("--output", type=str, nargs='?', default= "output", help="Name of directory to save outputs in.")
parser.add_argument("--gtex", type=str, nargs='?', default= data_path + "/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", \
                    help="Path to directory containing GTEX file (.gct extension).")
parser.add_argument("--enhancer_activity", type=str, nargs='?', default= data_path + "/h3k27ac_coverage_quantile_normalized.csv", \
                    help="Path to directory containing csv file with enhancer activity.")
parser.add_argument("--gene_expression", type=str, nargs='?', default= data_path + "/transcrtpts_rnaseq_quantile_normalized.csv", \
                    help="Path to directory containing csv file with gene expression.")
parser.add_argument("--chromatin_contacts", type=str, nargs='?', default= data_path + "/predicted_contacts.bed", \
                    help="Path to directory containing bed file with chromatin contacts.")
parser.add_argument("--genes_info", type=str, nargs='?', default= data_path + '/hg38_full.genes.gtf', \
                    help="Path to directory containing file with information about genes.")

parser.add_argument("--freq_filter_target", type=str, choices=['r', 'c'], default='r', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to select rare or common variants considering frequency in population")
parser.add_argument("--freq_filter_cutoff", type=float, default=0.01, choices=range(0,1), nargs='?', \
                    help="Set cut-off point for frequency filter. It expresses how rare or how common variant you what to select (rare/common depends on freq_filter_target argument)")
parser.add_argument("--population", type=str, nargs='*', choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH'], \
                    help="Choose population which input data will be copared with. You can pass multiple arguments. If you won't pass any argument, all populations will be used.")
parser.add_argument("--freq_filter_missing", type=str, choices=['r', 'c'], default='c', nargs='?', \
                    help="Choose 'r' (rare) or 'c' (common). It expresses if you want to treat variant with missing frequency data as rare or common variant.")
parser.add_argument("--reference_population", type=str, choices=['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH'], default='ALL', nargs='?', \
                    help="Choose reference population for binomial test.")
parser.add_argument("--bh_alpha", type=float, default=0.01, choices=range(0,1), nargs='?', \
                    help="Set cut-off point for Benjamini-Hochberg correction.")


args = parser.parse_args()

ANNOVAR = args.annovar_path
GATK = args.gatk_path
INPUT_VCF = args.input_vcf
PROMOTER_REGIONS = args.promoter_regions
ENHANCER_REGIONS = args.enhancer_regions
OUTPUT = args.output
ANNOTATED_PROMOTER_SNPs = OUTPUT + "/annotated_promoter_snps.csv"
ANNOTATED_ENHANCER_SNPs = OUTPUT + "/annotated_enhancer_snps.csv"
GTEX = args.gtex
ENHANCER_ACTIVITY = args.enhancer_activity
GENE_EXPRESSION = args.gene_expression
CHROMATIN_CONTACTS = args.chromatin_contacts
GENES_INFO = args.genes_info
freq_filter_target = args.freq_filter_target
freq_filter_cutoff = args.freq_filter_cutoff
population = args.population
freq_filter_missing = args.freq_filter_missing
reference_population = args.reference_population
bh_alpha = args.bh_alpha

if population and reference_population!='ALL':
    assert reference_population in population, "Reference population must be in population list"
if not population:
        population = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
population = population + ['ALL']

try:
    os.mkdir(OUTPUT)
except:

    print('Directory ', OUTPUT, ' already exist. Output files will be saved in this localization.')

#check if ANNOVAR and GATK are installed in provided localization
flags = fm.check_programs(GATK, ANNOVAR) 
if flags!= [0,0]:
    from itertools import compress
    raise TypeError('You should provide correct path to program(s)', list(compress(['GATK', 'ANNOVAR'], flags)))


# print how many variants are in input file
fm.count_variants(GATK, INPUT_VCF) 

biallelic_variants_in_regulatory_regions_path = fm.select_biallelic_inside_regulatory_regions(GATK, INPUT_VCF, PROMOTER_REGIONS, ENHANCER_REGIONS, OUTPUT)
print(biallelic_variants_in_regulatory_regions_path)
if fm.prepare_annovar_database(ANNOVAR)==0:
    raise TypeError(f'ANNOVAR database cound not be installed. You can install it manually from annovar website. It should be placed in {ANNOVAR}/humandb directory.')

annotated_variants_path = fm.annotate_freq(biallelic_variants_in_regulatory_regions_path,OUTPUT, ANNOVAR)
annotated_variants_path = {'promoter': 'output/promoter_SNPs.hg38_multianno.vcf', 'enhancer': 'output/enhancer_SNPs.hg38_multianno.vcf'}
filtered_by_freq_variants_file = fm.select_snps_by_freq(annotated_variants_path, OUTPUT,GATK, target=freq_filter_target, cutoff=freq_filter_cutoff, population=population, treating_gaps=freq_filter_missing)

filtered_by_freq_variants_file = {'promoter': 'output/promoter_SNPs.hg38_multianno.nomissing.gnomad_below_0.01.vcf', 'enhancer': 'output/enhancer_SNPs.hg38_multianno.nomissing.gnomad_below_0.01.vcf'}
rare_enriched_promoter_snps_df, rare_enriched_enhancer_snps_df = fm.select_enriched_snps(filtered_by_freq_variants_file, GATK, OUTPUT, bh_alpha=bh_alpha)

rare_enriched_promoter_snps_gene = fm.assign_genes_to_promoter_snps(rare_enriched_promoter_snps_df, PROMOTER_REGIONS)

rare_enriched_enhancer_snps_gene = fm.assign_genes_intronic_enhancer_snps(rare_enriched_enhancer_snps_df, ENHANCER_REGIONS)

genes_info_prepared = fm.prepare_genes_info(GENES_INFO)
rare_enriched_enhancer_snps_gene_closest = \
    fm.assign_closest_gene_to_enhancers(rare_enriched_enhancer_snps_gene, genes_info_prepared)
rare_enriched_enhancer_snps_gene_closest_contacts = \
    fm.assign_chromatin_contacting_gene_to_enhancer(rare_enriched_enhancer_snps_gene_closest, genes_info_prepared, CHROMATIN_CONTACTS)
rare_enriched_enhancer_snps_genes_collected = \
    fm.reformat_target_genes_enh(rare_enriched_enhancer_snps_gene_closest_contacts)

rare_enriched_enhancer_snps_genes_collected_corelations = \
    fm.check_signal_gene_expression_correlation_enhancer(rare_enriched_enhancer_snps_genes_collected, ENHANCER_ACTIVITY, GENE_EXPRESSION)
print(rare_enriched_enhancer_snps_genes_collected_corelations.iloc[:, -4:])


rare_enriched_promoter_snps_gene, rare_enriched_enhancer_snps_genes_collected_corelations = \
    fm.assign_median_tpm_gtex(rare_enriched_promoter_snps_gene, rare_enriched_enhancer_snps_genes_collected_corelations, GTEX)

promoter_snps_path = OUTPUT+'/final_rare_enriched_promoter_snps_gene.csv'
enhancer_snps_path = OUTPUT+'/final_rare_enriched_enhancer_snps_genes_collected_corelations.csv'
rare_enriched_promoter_snps_gene.to_csv(promoter_snps_path, sep = "\t", index = False)
rare_enriched_enhancer_snps_genes_collected_corelations.to_csv(enhancer_snps_path, sep = "\t", index = False)

promoter_snps_sample_lvl, enhancer_snps_sample_lvl = fm.import_vcf_sample_level(INPUT_VCF, OUTPUT, GATK, promoter_snps_path, enhancer_snps_path)


#funkcje związane z motywami

# snps_bed_files = fm.snps_to_bed_file(rare_enriched_promoter_snps, rare_enriched_enhancer_snps, OUTPUT)
# fm.search_motifs()
# fm.score_motifs(snps_bed_files)
# rare_enriched_promoter_snps_motif, rare_enriched_enhancer_snps_motif = \
#     fm.select_results(OUTPUT, rare_enriched_promoter_snps, rare_enriched_enhancer_snps)

#saving final results
#rare_enriched_promoter_snps_motif_gene.to_csv(ANNOTATED_PROMOTER_SNPs, sep = "\t", index = False)
#rare_enriched_enhancer_snps_motif_genes_collected_corelations.to_csv(ANNOTATED_ENHANCER_SNPs, sep = "\t", index = False)
