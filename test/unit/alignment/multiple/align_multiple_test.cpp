// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/multiple/align_multiple.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna4;

TEST(align_multiple_test, the_first_test)
{
    std::vector<seqan3::dna4_vector> input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4,        g,        g, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4,        g,        g, 'G'_dna4, 'G'_dna4, 'G'_dna4}
    };

    auto result = seqan3::align_multiple(input);

    EXPECT_RANGE_EQ(output, result);
}
