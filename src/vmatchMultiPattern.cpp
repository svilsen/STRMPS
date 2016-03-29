#include <Rcpp.h>

#include <IRanges_interface.h>
#include <Biostrings_interface.h>
#include <XVector_interface.h>

extern "C" {
#include "Biostrings.h"
}

//' @title
//' vmatchMultiPattern
//' @description
//' Wrapper function for \code{XStringSet_vmatch_pattern} allowing for multiple patterns.
//' Alternate version of vmatchMultiPattern.
//'
//' @param pattern an XStringSet-object.
//'
//' @param subject an XStringSet-object.
//'
//' @param max_mismatch the maximum allowed number of mismatching bases between the pattern and the subject.
//'
//' @param min_mismatch the minimum allowed number of mismatching bases between the pattern and the subject.
//'
//' @details
//' \code{vmatchMultiPattern}...
//'
//' @export
// [[Rcpp::export]]
SEXP vmatchMultiPattern(SEXP pattern, SEXP subject, SEXP max_mismatch, SEXP min_mismatch,
                        SEXP with_indels, SEXP fixed, SEXP algorithm, SEXP matches_as, SEXP envir) {
    XStringSet_holder Pset;
    XStringSet_holder Sset;

    Pset = _hold_XStringSet(pattern);
    Sset = _hold_XStringSet(subject);

    int Pset_length, i;
    int Sset_length, j;

    Chars_holder P, P_elt;
    Chars_holder S, S_elt;

    const char *algo, *ms_mode;

    algo = CHAR(STRING_ELT(algorithm, 0));
    ms_mode = CHAR(STRING_ELT(matches_as, 0));

    Pset_length = _get_length_from_XStringSet_holder(&Pset);
    Sset_length = _get_length_from_XStringSet_holder(&Sset);

    _init_match_reporting(ms_mode, Pset_length * Sset_length);
    for (i = 0; i < Pset_length; i++) {
        P_elt = _get_elt_from_XStringSet_holder(&Pset, i);

        for (j = 0; j < Sset_length; j++) {
            S_elt = _get_elt_from_XStringSet_holder(&Sset, j);
            _set_active_PSpair(i*Sset_length + j);

            _match_pattern_XString(&P_elt, &S_elt,
                                   max_mismatch, min_mismatch, with_indels, fixed,
                                   algo);
        }
    }

    return _MatchBuf_as_SEXP(_get_internal_match_buf(), envir);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

SEXP XStringSet_vmatch_pattern(Chars_holder &P, XStringSet_holder &S, int &Sset_length,
                               SEXP max_mismatch, SEXP min_mismatch, SEXP with_indels, SEXP fixed,
                               SEXP algorithm, SEXP ms_mode) {
    Chars_holder S_elt;
    int j;
    const char *algo;

    algo = CHAR(STRING_ELT(algorithm, 0));
    _init_match_reporting(CHAR(STRING_ELT(ms_mode, 0)), Sset_length);

    for (j = 0; j < Sset_length; j++) {
        S_elt = _get_elt_from_XStringSet_holder(&S, j);
        _set_active_PSpair(j);
        _match_pattern_XString(&P, &S_elt, max_mismatch, min_mismatch, with_indels, fixed,
                               algo);
    }
    return _MatchBuf_as_SEXP(_get_internal_match_buf(), R_NilValue);
}

//' @title
//' vmatchMultiPatternAlternate
//' @description
//' Wrapper function for \code{XStringSet_vmatch_pattern} allowing for multiple patterns.
//'
//' @param pattern an XStringSet-object.
//'
//' @param subject an XStringSet-object.
//'
//' @param max_mismatch the maximum allowed number of mismatching bases between the pattern and the subject.
//'
//' @param min_mismatch the minimum allowed number of mismatching bases between the pattern and the subject.
//'
//' @details
//' \code{vmatchMultiPatternAlternate}...
//'
//' @export
// [[Rcpp::export]]
Rcpp::List vmatchMultiPatternAlternate(SEXP pattern, SEXP subject, SEXP max_mismatch, SEXP min_mismatch,
                                       SEXP with_indels, SEXP fixed, SEXP algorithm, SEXP matches_as, SEXP envir) {
    XStringSet_holder Pset, Sset;
    Pset = _hold_XStringSet(pattern);
    Sset = _hold_XStringSet(subject);
    int Pset_length, Sset_length, i;
    Chars_holder P_elt;

    Pset_length = _get_length_from_XStringSet_holder(&Pset);
    Sset_length = _get_XStringSet_length(subject);

    Rcpp::List res(Pset_length);
    for (i = 0; i < Pset_length; i++) {
        P_elt = _get_elt_from_XStringSet_holder(&Pset, i);
        res[i] = XStringSet_vmatch_pattern(P_elt, Sset, Sset_length, max_mismatch, min_mismatch,
                                           with_indels, fixed, algorithm, matches_as);
    }
    return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////
#include <seqan/find.h>
using namespace seqan;

//' @title
//' vmatchMultiPatternSeqAn
//' @description
//' Function achieving the same as the wrapper functions for \code{XStringSet_vmatch_pattern} allowing for multiple patterns.
//' This alternate method is based on the external SeqAn C++ library.
//'
//' @param pattern a vector of strings.
//'
//' @param subject a vector of strings.
//'
//' @param max_mismatch the maximum allowed number of mismatching bases between the pattern and the subject.
//'
//'
//' @details
//' \code{vmatchMultiPatternSeqAn}.
//'
//' @export
//[[Rcpp::export]]
Rcpp::List vmatchMultiPatternSeqAn(SEXP pattern, SEXP subject, int max_mismatch) {
    std::vector<std::string> P = Rcpp::as<std::vector<std::string> >(pattern);
    std::vector<std::string> S = Rcpp::as<std::vector<std::string> >(subject);
    size_t P_length = P.size();
    size_t S_length = S.size();

    size_t i, j;
    Rcpp::List res(P_length);
    for (i = 0; i < P_length; i++) {
        String<char> P_elt = P[i].c_str();
        int P_length = P[i].length();
        Pattern<String<char>, Myers<FindInfix> > p_pattern(P_elt);

        Rcpp::List res_p(S_length);
        for (j = 0; j < S_length; j++) {
            String<char> S_elt = S[j].c_str();
            Finder<String<char> > s_finder(S_elt);

            std::vector<int> ends;
            while (find(s_finder, p_pattern, -max_mismatch))
            {
                while(findBegin(s_finder, p_pattern, -max_mismatch)) {
                    if((endPosition(s_finder) - beginPosition(s_finder)) == (P_length - 1)) {
                        ends.push_back(endPosition(s_finder) + 1);
                    }
                }
            }

            res_p[j] = ends;
        }

        res[i] = res_p;
    }

    return res;
}
