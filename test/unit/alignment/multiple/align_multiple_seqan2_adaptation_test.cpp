// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>
#include <vector>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/multiple/detail/align_multiple_seqan2_adaptation.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>

TEST(detail_initialise_scoring_scheme, no_scoring_configuration)
{
    // no scoring information
    seqan3::configuration config = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-2},
                                                                             seqan3::gap_open_score{-8}}};
    auto msaOpt = seqan3::detail::initialise_scoring_scheme<seqan::Dna>(config);

    EXPECT_TRUE((std::same_as<decltype(msaOpt.sc), seqan::Score<int>>));
}

TEST(detail_initialise_scoring_scheme, blosum62)
{
    // no scoring information
    seqan3::aminoacid_scoring_scheme scheme{seqan3::aminoacid_similarity_matrix::BLOSUM62};
    seqan3::configuration config = seqan3::align_cfg::scoring{scheme};
    auto msaOpt = seqan3::detail::initialise_scoring_scheme<seqan::AminoAcid>(config);

    // compare matrix size and abort test if not equal so the value comparison does not segfault
    ASSERT_EQ(static_cast<size_t>(decltype(msaOpt.sc)::TAB_SIZE), static_cast<size_t>(seqan::Blosum62::TAB_SIZE));

    seqan::Blosum62 expected_matrix{};
    for (size_t i = 0; i < static_cast<size_t>(seqan::Blosum62::TAB_SIZE); ++i)
        EXPECT_EQ(msaOpt.sc.data_tab[i], expected_matrix.data_tab[i]);
}

TEST(align_multiple_test, gap_score_conversion_test)
{
    // with go = -1, g = -1
    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                                    seqan3::gap_open_score{-1}}};
    // seqan2 does not add a gap score for the first gap character but just the gap open score
    // seqan3 does add the gap extension score additionally to the open score for the first gap character

    detail::align_multiple_seqan2_adaptation seqan2_adaptation{};

    auto msaOpt = seqan2_adaptation.msa_configuration<seqan3::dna4>(cfg);

    EXPECT_EQ(msaOpt.sc.data_gap_extend, -1);
    EXPECT_EQ(msaOpt.sc.data_gap_open, -2);
}
