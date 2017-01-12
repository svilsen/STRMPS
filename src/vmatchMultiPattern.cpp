#include <Rcpp.h>
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
