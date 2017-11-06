#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <seqan/find.h>
using namespace seqan;

//' @title
//' vmatchMultiPatternSeqAn
//' @description
//' Function achieving the same as the wrapper functions for \code{XStringSet_vmatch_pattern} allowing for multiple patterns.
//' This alternate method is based on the external SeqAn C++ library.
//'
//' @param forwardPattern a vector of strings.
//' @param reversePattern a vector of strings.
//' @param subject a vector of strings.
//' @param idVector a vector of numbers used for sequences identification.
//' @param max_mismatch the maximum allowed number of mismatching bases between the pattern and the subject.
//'
//' @details
//' \code{vmatchMultiPatternSeqAn}.
//[[Rcpp::export()]]
Eigen::MatrixXd vmatchMultiPatternSeqAn(SEXP forwardPattern, SEXP reversePattern, SEXP subject, SEXP idVector, int max_mismatch, int indels = 0) {
    std::vector<std::string> FP = Rcpp::as<std::vector<std::string> >(forwardPattern);
    std::vector<std::string> RP = Rcpp::as<std::vector<std::string> >(reversePattern);
    std::vector<std::string> S = Rcpp::as<std::vector<std::string> >(subject);
    std::vector<double> id = Rcpp::as<std::vector<double> >(idVector);

    std::size_t P_length = FP.size();
    std::size_t S_length = S.size();

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(S_length, 4);
    for (std::size_t i = 0; i < S_length; i++)
    {
        CharString S_elt = S[i].c_str();

        std::vector<int> forwardEnds_s;
        std::vector<int> reverseEnds_s;

        int markerIdentified = 0;
        int j = 0, h = 0;
        bool condition = ((j < P_length) & (markerIdentified <= 1));
        while (condition)
        {
            CharString FP_j = FP[j].c_str();
            Finder<CharString> s_finder_forward(S_elt);
            // Pattern<String<char>, Myers<FindInfix> > fp_pattern(FP[j].c_str());
            Pattern<CharString, DPSearch<SimpleScore> > fp_pattern(FP_j, SimpleScore(0, -1, -(2 - indels)));

            std::vector<int> fp_ends;
            while (find(s_finder_forward, fp_pattern, -max_mismatch))
            {
                while(findBegin(s_finder_forward, fp_pattern, -max_mismatch)) {
                    if((endPosition(s_finder_forward) - beginPosition(s_finder_forward)) == (FP[j].length() - 1)) {
                        fp_ends.push_back(endPosition(s_finder_forward) + 1);
                    }
                }
            }


            if (fp_ends.size() > 0)
            {
                CharString RP_j = RP[j].c_str();
                Finder<CharString> s_finder_reverse(S_elt);
                // Pattern<String<char>, Myers<FindInfix> > rp_pattern(RP[j].c_str());
                Pattern<CharString, DPSearch<SimpleScore> > rp_pattern(RP_j, SimpleScore(0, -1, -(2 - indels)));

                std::vector<int> rp_ends;
                while (find(s_finder_reverse, rp_pattern, -max_mismatch))
                {
                    while(findBegin(s_finder_reverse, rp_pattern, -max_mismatch)) {
                        if((endPosition(s_finder_reverse) - beginPosition(s_finder_reverse)) == (RP[j].length() - 1)) {
                            rp_ends.push_back(endPosition(s_finder_reverse) + 1);
                        }
                    }
                }

                if (rp_ends.size() > 0)
                {
                    markerIdentified++;

                    h = j + 1;
                    forwardEnds_s = fp_ends;
                    reverseEnds_s = rp_ends;
                }
            }

            j++;

            condition = ((j < P_length) & (markerIdentified <= 1));
        }

        res(i, 0) = id[i];
        if (markerIdentified == 1)
        {
            res(i, 1) = h;
            res(i, 2) = forwardEnds_s[0];
            res(i, 3) = reverseEnds_s[reverseEnds_s.size() - 1];
        }
    }

    return res;
}

