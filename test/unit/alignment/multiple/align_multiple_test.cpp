// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/iterator>
#include <string_view>
#include <vector>

#include <seqan3/alignment/multiple/align_multiple.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna4;

TEST(align_multiple_test, the_first_test)
{
    std::vector<seqan3::dna4_vector> input{"AAAACCCGGG"_dna4, "AACCCGGG"_dna4, "AAAACGGG"_dna4};
    seqan3::gap g{};
    std::vector<std::vector<seqan3::gapped<seqan3::dna4>>> output
    {
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {       g, 'A'_dna4, 'A'_dna4, 'C'_dna4,        g, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4,        g,        g, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4}
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
        {       g, 'A'_dna4, 'A'_dna4, 'C'_dna4,        g, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4},
        {'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4,        g,        g, 'C'_dna4, 'G'_dna4, 'G'_dna4, 'G'_dna4}
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

    auto result = seqan3::align_multiple(input);

    EXPECT_RANGE_EQ(output, result);
}
