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
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/multiple/align_multiple.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna15;
using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_rna5;
using seqan3::operator""_rna4;

TEST(align_multiple_test, the_first_dna4_test)
{
    std::vector<seqan3::dna4_vector> input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {       g,        g, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {       g,        g, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4}
    };

    auto result = seqan3::align_multiple(input);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_first_banded_test)
{
    std::vector<seqan3::dna4_vector> input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {       g,        g, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {       g,        g, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4}
    };

    auto cfg = seqan3::align_cfg::msa_default_configuration |
               seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                                  seqan3::align_cfg::upper_diagonal{4}};

    auto result = seqan3::align_multiple(input, cfg);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_first_aminoacid_test)
{
    // sequences taken from seqan/apps/seqan_tcoffee/tests/1aab.fa
    using seqan3::operator""_aa27;
    std::vector<seqan3::aa27_vector> input
    {
        "KKDSNAPKRAMTSFMFFSSDFRSKHSDLSIVEMSKAAGAAWKELGPEERKVYEEMAEKDKERYKREM"_aa27,
        "KPKRPRSAYNIYVSESFQEAKDDSAQGKLKLVNEAWKNLSPEEKQAYIQLAKDDRIRYDNEMKSWEEQMAE"_aa27,
        "ADKPKRPLSAYMLWLNSARESIKRENPDFKVTEVAKKGGELWRGLKDKSEWEAKAATAKQNYIRALQEYERNGG"_aa27,
        "DPNKPKRAPSAFFVFMGEFREEFKQKNPKNKSVAAVGKAAGERWKSLSESEKAPYVAKANKLKGEYNKAIAAYNKGESA"_aa27
    };

    // alignment taken from seqan/apps/seqan_tcoffee/tests/1aab.fasta
    std::vector<std::vector<seqan3::gapped<seqan3::aa27>>> output{4};
    using std::ranges::copy;
    copy(std::string_view{"KKDSNAPKRAMTSFMFFSSDFRSKHSDLS-----IVEMSKAAGAAWKELGPEERKVYEEMAEKDKERYKREM---------"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::aa27>>, std::cpp20::back_inserter(output[0]));
    copy(std::string_view{"-----KPKRPRSAYNIYVSESFQEAKDDS-----AQGKLKLVNEAWKNLSPEEKQAYIQLAKDDRIRYDNEMKSWEEQMAE"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::aa27>>, std::cpp20::back_inserter(output[1]));
    copy(std::string_view{"---ADKPKRPLSAYMLWLNSARESIKRENPDFK-VTEVAKKGGELWRGL--KDKSEWEAKAATAKQNYIRALQEYER-NGG"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::aa27>>, std::cpp20::back_inserter(output[2]));
    copy(std::string_view{"--DPNKPKRAPSAFFVFMGEFREEFKQKNPKNKSVAAVGKAAGERWKSLSESEKAPYVAKANKLKGEYNKAIAAYNKGESA"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::aa27>>, std::cpp20::back_inserter(output[3]));

    seqan3::aminoacid_scoring_scheme scheme{seqan3::aminoacid_similarity_matrix::BLOSUM62};
    seqan3::configuration config = seqan3::align_cfg::scoring{scheme};

    auto result = seqan3::align_multiple(input, config);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_first_rna5_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/rna5.fa -a rna -o data/out.fa
    std::vector<seqan3::rna5_vector> input{"UUUNCCCGGG"_rna5, "UUCCCGGG"_rna5, "UUUNCGGG"_rna5};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::rna5>>> output
    {
        {'U'_rna5, 'U'_rna5, 'U'_rna5, 'N'_rna5, 'C'_rna5, 'C'_rna5, 'C'_rna5, 'G'_rna5, 'G'_rna5, 'G'_rna5},
        {'U'_rna5, 'U'_rna5,        g,        g, 'C'_rna5, 'C'_rna5, 'C'_rna5, 'G'_rna5, 'G'_rna5, 'G'_rna5},
        {'U'_rna5, 'U'_rna5,        g,        g, 'U'_rna5, 'N'_rna5, 'C'_rna5, 'G'_rna5, 'G'_rna5, 'G'_rna5}
    };

    auto result = seqan3::align_multiple(input);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_first_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna1.fa -a dna -g -10 -e -2 -o data/out_g1.fa
    std::vector<seqan3::dna4_vector> input{"ACGGTGG"_dna4, "ACCGTGCC"_dna4, "GCCGGTGCC"_dna4};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4,        g, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'T'_dna4, 'G'_dna4, 'G'_dna4,        g},
        {'A'_dna4,        g, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'G'_dna4, 'C'_dna4, 'C'_dna4},
        {'G'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'T'_dna4, 'G'_dna4, 'C'_dna4, 'C'_dna4}
    };

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-2},
                                                                                    seqan3::gap_open_score{-8}}};

    auto result = seqan3::align_multiple(input, cfg);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_second_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna2.fa -a dna -g -2 -e -2 -o data/out_g2.fa
    std::vector<seqan3::dna5_vector> input{"NNTGTNN"_dna5, "GGTNTNNGT"_dna5, "NGTNTGGG"_dna5};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna5>>> output
    {
        {'N'_dna5, 'N'_dna5, 'T'_dna5, 'G'_dna5, 'T'_dna5, 'N'_dna5,        g, 'N'_dna5,        g,        g,        g,        g,        g},
        {'G'_dna5,        g,        g, 'G'_dna5, 'T'_dna5, 'N'_dna5, 'T'_dna5, 'N'_dna5, 'N'_dna5, 'G'_dna5, 'T'_dna5,        g,        g},
        {       g, 'N'_dna5,        g, 'G'_dna5, 'T'_dna5, 'N'_dna5, 'T'_dna5,        g,        g, 'G'_dna5,        g, 'G'_dna5, 'G'_dna5}
    };

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-2},
                                                                                    seqan3::gap_open_score{0}}};

    auto result = seqan3::align_multiple(input, cfg);

    EXPECT_RANGE_EQ(output, result);
}

TEST(align_multiple_test, the_third_gap_score_test)
{
    // alignment generated with seqan_tcoffee app: ./seqan_tcoffee -s data/dna3.fa -a dna -g -16 -e -4 -o data/out_g3.fa
    std::vector<seqan3::dna15_vector> input{"GGGTGGYTG"_dna15, "KTGTGGYTYTG"_dna15, "KTGTYYYTG"_dna15};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna15>>> output
    {
        {'G'_dna15, 'G'_dna15, 'G'_dna15, 'T'_dna15, 'G'_dna15, 'G'_dna15,         g,         g, 'Y'_dna15, 'T'_dna15, 'G'_dna15},
        {'K'_dna15, 'T'_dna15, 'G'_dna15, 'T'_dna15, 'G'_dna15, 'G'_dna15, 'Y'_dna15, 'T'_dna15, 'Y'_dna15, 'T'_dna15, 'G'_dna15},
        {'K'_dna15, 'T'_dna15, 'G'_dna15, 'T'_dna15, 'Y'_dna15, 'Y'_dna15,         g,         g, 'Y'_dna15, 'T'_dna15, 'G'_dna15}
    };

    constexpr seqan3::configuration cfg = seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-4},
                                                                                    seqan3::gap_open_score{-12}}};

    auto result = seqan3::align_multiple(input, cfg);

    EXPECT_RANGE_EQ(output, result);
}
