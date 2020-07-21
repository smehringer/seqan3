// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides some detail functionality to the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace seqan3::detail
{

// unfortunately it doesn't work with out alphabets..
auto convert_alph_3_to_2(seqan3::dna4 chr)
{
    return seqan::Dna{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::aa27 chr)
{
    return seqan::AminoAcid{seqan3::to_char(chr)};
}

template <typename alphabet_type, typename score_type, typename seqan3_configuration_t>
auto seqan2_msa_configuration(seqan3_configuration_t const & config)
{
    seqan::MsaOptions<alphabet_type, score_type> msaOpt{};

    seqan::appendValue(msaOpt.method, 0); // global alignment
    seqan::appendValue(msaOpt.method, 1); // local alignment
    msaOpt.build = 0; // neighbour joining to build the guide tree

    if constexpr (config.template exists<seqan3::align_cfg::band_fixed_size>())
    {
        msaOpt.pairwiseAlignmentMethod = 2; // banded
        auto const & band = get<seqan3::align_cfg::band_fixed_size>(config);
        msaOpt.bandWidth = band.upper_diagonal.get() - band.lower_diagonal.get();
    }
    else
    {
        msaOpt.pairwiseAlignmentMethod = 1; // unbanded
    }

    score_type scMat;
    msaOpt.sc = scMat;

    // seqan2 tcoffee app default: gap -1, gap open -13
    auto const & gaps = config.get_or(align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-13}}}).value;
    msaOpt.sc.data_gap_open = gaps.get_gap_open_score();
    msaOpt.sc.data_gap_extend = gaps.get_gap_score();
    if constexpr (std::is_same_v<score_type, seqan::Score<int>>)
    {
        msaOpt.sc.data_match = 5;       // tcoffee app default
        msaOpt.sc.data_mismatch = -4;   // tcoffee app default
    }

    return msaOpt;
}

} //seqan3::detail
