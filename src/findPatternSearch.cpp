#include <Rcpp.h>
#include <seqan/align.h>
using namespace seqan;

//' @title
//' getAlignmentNeighbourScore
//'
//' @description
//' Gets the alignment score of specified neighbours.
//'
//' @param parent the parent string with which the neighbours are matched.
//'
//' @param neighbours potential variations of the parent string.
//'
//' @param gapOpening the penalty for opening a new gap (positive value).
//'
//' @param gapExtension the penalty for extenting a gap (positive value).
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector getAlignmentNeighbourScore(SEXP parent, Rcpp::StringVector neighbours, int mismatch, int gapOpening, int gapExtension) {
    String<char> pattern = Rcpp::as<std::string>(parent);
    std::vector<std::string> subject = Rcpp::as<std::vector<std::string> >(neighbours);

    size_t pattern_length = length(pattern);
    size_t subject_length = subject.size();

    Rcpp::NumericVector res(subject_length);
    for (size_t i; i < subject_length; i++) {
        String<char> subject_i = subject[i];
        Align<String<char> > align;
        resize(rows(align), 2);
        assignSource(row(align, 0), subject_i);
        assignSource(row(align, 1), pattern);

        int alignmentScore = globalAlignment(align, Score<int, Simple>(1, -mismatch, -gapExtension, -gapOpening));
        res[i] = alignmentScore - gapExtension;
    }

    return res;
}
