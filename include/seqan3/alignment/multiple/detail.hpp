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
    msaOpt.sc.data_mismatch = -4;   // tcoffee app default

    return msaOpt;
}

} //seqan3::detail
