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

#include <seqan3/std/algorithm>
#include <cstdlib>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/alignment/multiple/detail.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/chunk.hpp>

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

    std::string mat;
    seqan::convertAlignment(gAlign, mat);

    using gapped_alphabet_type = seqan3::gapped<seqan3::dna4>;
    std::vector<std::vector<gapped_alphabet_type>> output{seqan::length(seqan::stringSet(gAlign))};
    size_t ali_len = mat.size() / output.size();

    auto iter = output.begin();
    for (auto && gapped_seq : mat | seqan3::views::char_to<gapped_alphabet_type> | seqan3::views::chunk(ali_len))
    {
        iter->reserve(ali_len);
        std::ranges::copy(gapped_seq, std::cpp20::back_inserter(*iter++));
    }

    return output;
}

} // namespace seqan3

#endif // !SEQAN3_HAS_SEQAN2
