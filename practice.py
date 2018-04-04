# import the algorithm
from fastlmm.association import single_snp

# set up data
##############################
bed_fn = "/home/conrad/Data/primary_project_3/gwas/H19_indels_plink/indels_binary"
pheno_fn = "/home/conrad/Data/primary_project_3/gwas/H19_indels_plink/indels.phe"

# run gwas
###################################################################
results_df = single_snp(bed_fn,  pheno_fn, leave_out_one_chrom=False, count_A1=True)

# manhattan plot
import fastlmm.util.util as flutil
flutil.manhattan_plot(results_df.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-5,xaxis_unit_bp=False)

# # qq plot
from fastlmm.util.stats import plotp
plotp.qqplot(results_df["PValue"].values, xlim=[0,5], ylim=[0,5], fileout="test_qq")

# from pysnptools.snpreader import Ped
# from pysnptools.snpreader import Pheno
# from pysnptools.snpreader import wrap_plink_parser

# # Load snp data:
# print "Loading variant data..."
# ped_file = Ped(bed_fn)
# print "Loading phenotype data..."
# pheno_fn = Pheno(pheno_fn)

# # Run basic association test:
# print "Running FaST-LMM single_snp test..."
# results_df = single_snp(test_snps=ped_file, pheno=pheno_fn, leave_out_one_chrom=True)

# print results_df