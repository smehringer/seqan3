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

namespace seqan3
{

template <std::ranges::forward_range range_t>
auto align_multiple(std::vector<range_t> const & /*input*/)
{
    // do something

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
