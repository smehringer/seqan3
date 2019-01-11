// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/view/filter.hpp>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(alignment_file_input_iterator, concepts)
{
    using it_t = typename alignment_file_input<>::iterator;
    using sen_t = typename alignment_file_input<>::sentinel;

    EXPECT_TRUE((std::InputIterator<it_t>));
    EXPECT_TRUE((std::Sentinel<sen_t, it_t>));
}

struct alignment_file_input_f : public ::testing::Test
{
    std::string input
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    std::vector<dna5_vector> seq_comp
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> id_comp
    {
        "read1",
        "read2",
        "read3"
    };

    std::vector<std::vector<phred42>> qual_comp
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };
};

TEST_F(alignment_file_input_f, concepts)
{
    using t = alignment_file_input<>;
    EXPECT_TRUE((std::ranges::InputRange<t>));

    using ct = alignment_file_input<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::InputRange<ct>));
}

TEST_F(alignment_file_input_f, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            // must contain at least one record
            filecreator << "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n";
        }

        EXPECT_NO_THROW( alignment_file_input<>{filename.get_path()} );
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        test::tmp_filename filename{"alignment_file_input_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( alignment_file_input<>{filename.get_path()} ,
                      unhandled_extension_error );
    }

    /* non-existent file*/
    {
        EXPECT_THROW( alignment_file_input<>{"/dev/nonexistant/foobarOOO"},
                      file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            // must contain at least one record
            filecreator << "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n";
        }

        EXPECT_NO_THROW(( alignment_file_input<alignment_file_input_default_traits_sam,
                                             fields<field::SEQ>,
                                             type_list<alignment_file_format_sam>,
                                             char>{filename.get_path(), fields<field::SEQ>{}} ));
    }
}

TEST_F(alignment_file_input_f, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( alignment_file_input<alignment_file_input_default_traits_sam,
                                         fields<field::SEQ, field::ID, field::QUAL>,
                                         type_list<alignment_file_format_sam>,
                                         char>{std::istringstream{input},
                                               alignment_file_format_sam{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( alignment_file_input<alignment_file_input_default_traits_sam,
                                          fields<field::SEQ, field::ID, field::QUAL>,
                                          type_list<alignment_file_format_sam>,
                                          char>{std::istringstream{input},
                                                alignment_file_format_sam{},
                                                fields<field::SEQ, field::ID, field::QUAL>{}} ));
}

TEST_F(alignment_file_input_f, default_template_args_and_deduction_guides)
{
    using comp0 = alignment_file_input_default_traits_sam;
    using comp1 = fields<field::SEQ, field::ID, field::OFFSET, field::REF_SEQ,
                         field::REF_ID, field::REF_OFFSET, field::ALIGNMENT,
                         field::MAPQ, field::QUAL, field::FLAG, field::MATE,
                         field::TAGS, field::EVALUE, field::BIT_SCORE, field::HEADER_PTR>;
    using comp2 = type_list<alignment_file_format_sam>;
    using comp3 = char;

    /* default template args */
    {
        using t = alignment_file_input<>;
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"; // must contain at least one record
        }

        alignment_file_input fin{filename.get_path()};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"; // must contain at least one record
        }

        alignment_file_input fin{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::istringstream ext{input};
        alignment_file_input fin{ext, alignment_file_format_sam{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        alignment_file_input fin{std::istringstream{input}, alignment_file_format_sam{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor + custom fields + different stream char type */
    {
        auto winput = input | view::convert<wchar_t>;
        std::wistringstream ext{winput};
        alignment_file_input fin{ext, alignment_file_format_sam{}, fields<field::SEQ>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));
    }

    /* guided stream temporary constructor + custom fields + different stream char type */
    {
        auto winput = input | view::convert<wchar_t>;
        alignment_file_input fin{std::wistringstream{winput}, alignment_file_format_sam{}, fields<field::SEQ>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));
    }
}

TEST_F(alignment_file_input_f, record_reading)
{
    /* record based reading */
    alignment_file_input fin{std::istringstream{input}, alignment_file_format_sam{}};

    size_t counter = 0;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::ID>(rec), id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::QUAL>(rec), qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_f, record_reading_custom_fields)
{
    /* record based reading */
    alignment_file_input fin{std::istringstream{input},
                         alignment_file_format_sam{},
                         fields<field::ID, field::SEQ>{}};

    size_t counter = 0;
    for (auto & [ id, seq ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(qual, qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_f, file_view)
{
    alignment_file_input fin{std::istringstream{input}, alignment_file_format_sam{}};

    auto minimum_length_filter = view::filter([] (auto const & rec)
    {
        return size(get<field::SEQ>(rec)) >= 5;
    });

    size_t counter = 1; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::ID>(rec), id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::QUAL>(rec), qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

// ----------------------------------------------------------------------------
// decompression
// ----------------------------------------------------------------------------

template <typename fixture_t, typename input_file_t>
void decompression_impl(fixture_t & fix, input_file_t & fin)
{
    size_t counter = 0;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::SEQ>(rec), fix.seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::ID>(rec),  fix.id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::QUAL>(rec), fix.qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

#ifdef SEQAN3_HAS_ZLIB

std::string input_gz
{
    '\x1F','\x8B','\x08','\x08','\x9D','\x5B','\x38','\x5C','\x00','\x03','\x74','\x65','\x73','\x74','\x2E','\x73',
    '\x61','\x6D','\x00','\x6D','\xCE','\xBF','\x0A','\xC2','\x30','\x10','\xC7','\xF1','\xF9','\xD7','\xB7','\x08',
    '\x15','\xFF','\xD4','\xA8','\xB9','\x24','\xB6','\x90','\x2D','\x56','\xB8','\x29','\x5D','\x92','\x17','\x28',
    '\x58','\xC1','\x35','\x93','\xBE','\xBD','\x21','\xBA','\x08','\x4E','\x07','\xC7','\x7D','\xBE','\x5C','\x5E',
    '\xE6','\x1B','\xC1','\x12','\xF2','\x72','\x07','\xA1','\x27','\x50','\xA4','\x40','\x57','\x1D','\x3E','\x1B',
    '\x05','\xA3','\x14','\xFC','\xC8','\x09','\xA2','\x6D','\x57','\xF0','\xD1','\x3D','\x9C','\xC6','\x14','\xCA',
    '\x18','\x9A','\x5C','\xB4','\x86','\xD5','\xF5','\x56','\xA3','\xD7','\x18','\x8A','\x2D','\x3E','\xFE','\x68',
    '\xE6','\x31','\xF1','\xE4','\xB9','\x26','\xD6','\x9B','\xED','\xAE','\xC3','\xF3','\xE5','\x2E','\x2E','\x4A',
    '\x23','\xAD','\x3C','\xD7','\x8C','\x81','\x35','\x15','\x19','\xF4','\xE6','\xFB','\x84','\xFD','\x13','\x63',
    '\xF6','\x9C','\x7C','\xF2','\x10','\xA2','\xDB','\xCB','\xC3','\xF1','\xD4','\xBC','\x01','\xDB','\x85','\xA3',
    '\xD3','\xC3','\x00','\x00','\x00'
};

TEST_F(alignment_file_input_f, decompression_by_filename_gz)
{
    test::tmp_filename filename{"alignment_file_output_test.sam.gz"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(begin(input_gz), end(input_gz), std::ostreambuf_iterator<char>{of});
    }

    alignment_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, decompression_by_stream_gz)
{
    alignment_file_input fin{std::istringstream{input_gz}, alignment_file_format_sam{}};

    decompression_impl(*this, fin);
}
#endif

#ifdef SEQAN3_HAS_BZIP2
std::string input_bz2
{
    '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\x7B','\xE2','\xE1','\x92','\x00','\x00',
    '\x5C','\x5F','\x80','\x00','\x30','\x2D','\xFF','\xFF','\x90','\x3C','\x83','\x0C','\x00','\x27','\x20','\x10',
    '\x60','\x20','\x00','\x8A','\x86','\x82','\x4D','\x4D','\xA6','\xA6','\x9A','\x60','\xD4','\xC8','\xC9','\x99',
    '\x35','\x34','\x06','\x44','\x9B','\x51','\xA0','\x83','\x4C','\x06','\x88','\x31','\x19','\x32','\xDF','\x59',
    '\x81','\x84','\x10','\x62','\x4B','\x06','\x22','\x21','\xA8','\xEA','\x68','\xCD','\xA2','\x15','\xB7','\xE5',
    '\xA7','\xAB','\x0A','\xD2','\xB8','\x0A','\xEF','\xC3','\x18','\x35','\xFE','\x2C','\xE9','\x1C','\x72','\x8D',
    '\xA6','\xE2','\xC7','\x3D','\xBC','\x41','\x0E','\x00','\x50','\x3E','\x05','\x0E','\x0F','\x46','\xF5','\x2B',
    '\x39','\xEF','\x92','\xA5','\x28','\x85','\xEA','\xA5','\x93','\xE0','\xFD','\x27','\xBF','\x76','\xCC','\xE2',
    '\x6A','\xE9','\x32','\xE0','\x11','\x05','\x09','\x44','\xD2','\x51','\xB6','\x90','\x2A','\x73','\x94','\x54',
    '\x62','\x96','\x19','\xBB','\x92','\xE4','\xB8','\x20','\x28','\x32','\x8E','\x0C','\x09','\xE7','\xF8','\xBB',
    '\x92','\x29','\xC2','\x84','\x83','\xDF','\x17','\x0C','\x90'
};

TEST_F(alignment_file_input_f, decompression_by_filename_bz2)
{
    test::tmp_filename filename{"alignment_file_output_test.sam.bz2"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(begin(input_bz2), end(input_bz2), std::ostreambuf_iterator<char>{of});
    }

    alignment_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, decompression_by_stream_bz2)
{
    alignment_file_input fin{std::istringstream{input_bz2}, alignment_file_format_sam{}};

    decompression_impl(*this, fin);
}
#endif

// ----------------------------------------------------------------------------
// SAM format specificities
// ----------------------------------------------------------------------------

struct alignment_file_input_sam_format_f : public alignment_file_input_f
{

    std::vector<unsigned> offsets
    {
        1,
        0,
        1
    };

    std::vector<dna5_vector> ref_seqs = {"ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5};

    std::vector<gapped<dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped2 = {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5,
                                                 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5,
                                                 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped3 = {'T'_dna5, 'G'_dna5, 'A'_dna5, 'T'_dna5,
                                                 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5,};

    std::vector<std::string> ref_ids = {"ref"};

    std::vector<unsigned> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<gapped<dna5>>{'C'_dna5, gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<gapped<dna5>>{'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                    'G'_dna5, 'N'_dna5, gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<gapped<dna5>>{'G'_dna5, gap{}, 'A'_dna5, 'G'_dna5,
                                                    'T'_dna5, 'A'_dna5, gap{}, 'T'_dna5}}
    };
};

TEST_F(alignment_file_input_sam_format_f, construct_by_filename_and_read_alignments)
{
    test::tmp_filename filename{"alignment_file_input_constructor.sam"};
    {
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        // must contain at least one record
        filecreator << input;
    }

    alignment_file_input fin{filename.get_path(), ref_ids, ref_seqs, fields<field::ALIGNMENT>{}};

    size_t counter = 0;
    for (auto & [ alignment ] : fin)
    {
        EXPECT_EQ(get<0>(alignment), get<0>(alignments[counter]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[counter]));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_sam_format_f, construct_from_stream_and_read_alignments)
{
    // stream + format_tag + fields
    alignment_file_input<alignment_file_input_default_traits_sam,
                                          fields<field::ALIGNMENT>,
                                          type_list<alignment_file_format_sam>,
                                          char> fin{std::istringstream{input},
                                                ref_ids,
                                                ref_seqs,
                                                alignment_file_format_sam{},
                                                fields<field::ALIGNMENT>{}};

    size_t counter = 0;
    for (auto & [ alignment ] : fin)
    {
        EXPECT_EQ(get<0>(alignment), get<0>(alignments[counter]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[counter]));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}
