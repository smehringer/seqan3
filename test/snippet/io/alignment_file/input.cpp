#include <iostream>
#include <fstream>

#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/single_pass_input.hpp>

//! [my_traits]
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

using namespace seqan3;

struct my_traits : alignment_file_input_default_traits_sam
{
    using sequence_alphabet = dna4;                        // instead of dna5

    template <typename alph>
    using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
};

// ...

// alignment_file_input<my_traits> fin{"/tmp/my.sam"};

// ...
//! [my_traits]

auto sam_file_raw = R"//![sam_file](
@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
r002    0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA    *
r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
r004    0 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC       *
r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9,+,5S6M,30,1;
r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1
)//![sam_file]";

auto write_sam_file = [] (auto & sam_file_string)
    {
        std::ofstream file_stream{"/tmp/my.sam"};
        std::ranges::ostreambuf_iterator stream_it{file_stream};

        std::string sam_file{sam_file_string};
        auto sam_view = sam_file | view::single_pass_input;

        ++begin(sam_view); // skip first new line from raw string
        ranges::copy(sam_view | view::take_until(is_char<' '>), stream_it);

        while(begin(sam_view) != end(sam_view))
        {
            stream_it = '\t';
            detail::consume(sam_view | view::take_until(!is_char<' '>));
            ranges::copy(sam_view | view::take_until(is_char<' '>), stream_it);
        }
    };

int main()
{

write_sam_file(sam_file_raw);
std::vector<std::string> ref_ids{{"ref"}};
std::vector<dna5_vector> refs{{"ATCGAGAGCTAGAGAGCGAGAGCGAGCAGAGCGACGAGCGAGCAG"_dna5}};

{
//![get_header]
alignment_file_input fin{"/tmp/my.sam", ref_ids, refs}; // TODO change this to not take the ref info once the gap_decorator is in

// access the header information
debug_stream << fin.header().format_version << std::endl; // 1.6
debug_stream << fin.header().ref_dict << std::endl;       // [(ref,(45,))] (this only works with our debug_stream!)
//![get_header]
}

}
