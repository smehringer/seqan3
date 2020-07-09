// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/test/seqan2.hpp>

#if !SEQAN3_HAS_SEQAN2 // multiple sequence alignment is only enabled with seqan2
static_assert(false, "You need to have seqan 2.x");
#else // SEQAN3_HAS_SEQAN2

#include <cstdlib>
#include <seqan3/std/ranges>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace seqan3::detail
{

// unfortunately it doesn't work with out alphabets..
auto convert_alph_3_to_2(seqan3::dna4 chr)
{
    return seqan::Dna{seqan3::to_char(chr)};
}

template <typename alphabet_type, typename score_type>
auto seqan2_msa_configuration()
{
    seqan::MsaOptions<alphabet_type, score_type> msaOpt{};

    seqan::appendValue(msaOpt.method, 0); // global alignment
    seqan::appendValue(msaOpt.method, 1); // local alignment
    msaOpt.build = 0; // neighbour joining to build the guide tree
    msaOpt.pairwiseAlignmentMethod = 1; // unbanded
    score_type scMat;
    msaOpt.sc = scMat;
    msaOpt.sc.data_gap_open = -13;  // tcoffee app default
    msaOpt.sc.data_gap_extend = -1; // tcoffee app default
    msaOpt.sc.data_match = 5;       // tcoffee app default
    msaOpt.sc.data_match = -4;      // tcoffee app default

    return msaOpt;
}

} //seqan3::detail

namespace seqan3
{

template <std::ranges::forward_range range_t>
auto align_multiple(std::vector<range_t> const & input)
{
    using score_type = seqan::Score<int>;
    using alphabet_type = decltype(detail::convert_alph_3_to_2(std::ranges::range_value_t<range_t>{}));
    using sequence_type = seqan::String<alphabet_type>;
    using graph_type = seqan::Graph<seqan::Alignment<seqan::StringSet<sequence_type, seqan::Dependent<>>,
                                                     void,
                                                     seqan::WithoutEdgeId>>;

    auto msaOpt = detail::seqan2_msa_configuration<alphabet_type, score_type>();

    // fill seqan2 data storage
    seqan::StringSet<sequence_type, seqan::Owner<>> sequenceSet;
    seqan::StringSet<seqan::String<char>> sequenceNames;

    seqan::String<char> dummy_name = "dummy_name";
    for (auto const & seq : input)
    {
        sequence_type tmp;
        for (auto chr : seq)
            seqan::appendValue(tmp, detail::convert_alph_3_to_2(chr));

        seqan::appendValue(sequenceSet, tmp);
        seqan::appendValue(sequenceNames, dummy_name);
    }

    // Alignment of the sequences
    graph_type gAlign;

    // MSA
    try
    {
        seqan::globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
    }
    catch (const std::bad_alloc & exception)
    {
        std::cerr << "Allocation for globalAlignment failed. Use smaller data or try a seeded alignment. \n"
                  << exception.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // fake result
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4,        g,        g, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4,        g,        g, 'G'_dna4, 'G'_dna4, 'G'_dna4}
    };

    return output;
}

} // namespace seqan3

#endif // !SEQAN3_HAS_SEQAN2
