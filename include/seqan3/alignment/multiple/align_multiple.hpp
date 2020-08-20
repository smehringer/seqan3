// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \author Simon Sasse <simon.sasse AT fu-berlin.de>
 */

#pragma once

#include <seqan3/test/seqan2.hpp>

#if !SEQAN3_HAS_SEQAN2 // multiple sequence alignment is only enabled with seqan2
static_assert(false, "You need to have seqan 2.x");
#else // SEQAN3_HAS_SEQAN2

#include <seqan3/std/ranges>
#include <vector>

#include <seqan/graph_msa.h>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/multiple/detail/align_multiple_seqan2_adaptation.hpp>

namespace seqan3::align_cfg
{

constexpr configuration msa_default_configuration = gap{gap_scheme{gap_score{-1}, gap_open_score{-13}}};

}

namespace seqan3
{

template <std::ranges::forward_range range_t, typename config_t = decltype(align_cfg::msa_default_configuration)>
auto align_multiple(std::vector<range_t> const & input, config_t config = align_cfg::msa_default_configuration)
{
    using seqan3_alphabet_type = std::ranges::range_value_t<range_t>;
    using seqan2_adaptation_type = detail::align_multiple_seqan2_adaptation<seqan3_alphabet_type>;
    using graph_type = typename seqan2_adaptation_type::graph_type;

    seqan2_adaptation_type seqan2_adaptation{};

    auto msaOpt = seqan2_adaptation.create_msa_configuration(config);
    auto && [sequenceSet, sequenceNames] = seqan2_adaptation.convert_sequences(input);

    graph_type gAlign;

    seqan::globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);

    auto && output = seqan2_adaptation.create_output(gAlign);

    return output;
}

} // namespace seqan3

#endif // !SEQAN3_HAS_SEQAN2
