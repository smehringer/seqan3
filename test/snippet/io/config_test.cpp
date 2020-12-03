#include <utility>
#include <fstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::configuration config = seqan3::io_cfg::select_fields<seqan3::field::seq> |
                                   seqan3::io_cfg::select_formats<seqan3::format_fasta>;

    seqan3::sequence_file_input fin{"/tmp/test.fa", config};
    seqan3::debug_stream << *fin.begin() << '\n';

    std::ifstream stream{"/tmp/test.fa"};
    seqan3::sequence_file_input fin_stream{stream, config};
    seqan3::debug_stream << *fin_stream.begin() << '\n';
}
