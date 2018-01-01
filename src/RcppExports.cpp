// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// create_block_cpp
Rcpp::List create_block_cpp(const Rcpp::DataFrame ann);
RcppExport SEXP _SGNB_create_block_cpp(SEXP annSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame >::type ann(annSEXP);
    rcpp_result_gen = Rcpp::wrap(create_block_cpp(ann));
    return rcpp_result_gen;
END_RCPP
}
// create_read_type_cpp
Rcpp::DataFrame create_read_type_cpp(const std::string& input_sam_path, const Rcpp::List block_ann, const Rcpp::List gene_range, int min_overlap);
RcppExport SEXP _SGNB_create_read_type_cpp(SEXP input_sam_pathSEXP, SEXP block_annSEXP, SEXP gene_rangeSEXP, SEXP min_overlapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type input_sam_path(input_sam_pathSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type block_ann(block_annSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type gene_range(gene_rangeSEXP);
    Rcpp::traits::input_parameter< int >::type min_overlap(min_overlapSEXP);
    rcpp_result_gen = Rcpp::wrap(create_read_type_cpp(input_sam_path, block_ann, gene_range, min_overlap));
    return rcpp_result_gen;
END_RCPP
}
// create_read_type_group_cpp
Rcpp::DataFrame create_read_type_group_cpp(const std::vector<std::string>& read_gene_unique_vec, const std::vector<std::string>& read_gene_vec, const std::vector<std::string>& read_type_vec, double min_reduce);
RcppExport SEXP _SGNB_create_read_type_group_cpp(SEXP read_gene_unique_vecSEXP, SEXP read_gene_vecSEXP, SEXP read_type_vecSEXP, SEXP min_reduceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_gene_unique_vec(read_gene_unique_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_gene_vec(read_gene_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_type_vec(read_type_vecSEXP);
    Rcpp::traits::input_parameter< double >::type min_reduce(min_reduceSEXP);
    rcpp_result_gen = Rcpp::wrap(create_read_type_group_cpp(read_gene_unique_vec, read_gene_vec, read_type_vec, min_reduce));
    return rcpp_result_gen;
END_RCPP
}
// fit_SGNB_cpp
Rcpp::DataFrame fit_SGNB_cpp(const std::vector<std::string>& read_gene_unique_vec, const std::vector<std::string>& read_gene_vec, Rcpp::NumericMatrix data_matrix, Rcpp::NumericVector lib_size_norm, Rcpp::IntegerVector group_sample_num, Rcpp::List gene_size_ls, double tol, int times);
RcppExport SEXP _SGNB_fit_SGNB_cpp(SEXP read_gene_unique_vecSEXP, SEXP read_gene_vecSEXP, SEXP data_matrixSEXP, SEXP lib_size_normSEXP, SEXP group_sample_numSEXP, SEXP gene_size_lsSEXP, SEXP tolSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_gene_unique_vec(read_gene_unique_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_gene_vec(read_gene_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data_matrix(data_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lib_size_norm(lib_size_normSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type group_sample_num(group_sample_numSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type gene_size_ls(gene_size_lsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_SGNB_cpp(read_gene_unique_vec, read_gene_vec, data_matrix, lib_size_norm, group_sample_num, gene_size_ls, tol, times));
    return rcpp_result_gen;
END_RCPP
}
// fit_SGNB_exact_cpp
Rcpp::List fit_SGNB_exact_cpp(const std::vector<std::string>& read_gene_vec, const std::vector<int>& read_type_vec, Rcpp::NumericMatrix data_matrix, Rcpp::NumericVector lib_size_norm, Rcpp::IntegerVector group_sample_num, double tol, int times);
RcppExport SEXP _SGNB_fit_SGNB_exact_cpp(SEXP read_gene_vecSEXP, SEXP read_type_vecSEXP, SEXP data_matrixSEXP, SEXP lib_size_normSEXP, SEXP group_sample_numSEXP, SEXP tolSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type read_gene_vec(read_gene_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type read_type_vec(read_type_vecSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data_matrix(data_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lib_size_norm(lib_size_normSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type group_sample_num(group_sample_numSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(fit_SGNB_exact_cpp(read_gene_vec, read_type_vec, data_matrix, lib_size_norm, group_sample_num, tol, times));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SGNB_create_block_cpp", (DL_FUNC) &_SGNB_create_block_cpp, 1},
    {"_SGNB_create_read_type_cpp", (DL_FUNC) &_SGNB_create_read_type_cpp, 4},
    {"_SGNB_create_read_type_group_cpp", (DL_FUNC) &_SGNB_create_read_type_group_cpp, 4},
    {"_SGNB_fit_SGNB_cpp", (DL_FUNC) &_SGNB_fit_SGNB_cpp, 8},
    {"_SGNB_fit_SGNB_exact_cpp", (DL_FUNC) &_SGNB_fit_SGNB_exact_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SGNB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
