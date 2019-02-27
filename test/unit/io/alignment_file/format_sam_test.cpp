// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((AlignmentFileOutputFormat<alignment_file_format_sam>));
    EXPECT_TRUE((AlignmentFileInputFormat<alignment_file_format_sam>));
}

struct alignment_data : public ::testing::Test
{
    std::vector<dna5_vector> seqs
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> ids
    {
        "read1",
        "read2",
        "read3"
    };

    std::vector<std::vector<phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };

    std::vector<unsigned> offsets
    {
        1,
        0,
        1
    };

    dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<gapped<dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped2 = {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5,
                                                 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5,
                                                 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped3 = {'T'_dna5, 'G'_dna5, 'A'_dna5, 'T'_dna5,
                                                 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5,};

    std::string ref_id = "ref";

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

    std::vector<unsigned> flags
    {
        41,
        42,
        43
    };

    std::vector<unsigned> mapqs
    {
        61,
        62,
        63
    };

    std::vector<std::tuple<std::string, unsigned, unsigned>> mates
    {
        {"ref", 10, 300},
        {"ref", 10, 300},
        {"ref", 10, 300}
    };

    std::vector<sam_tag_dictionary> tag_dicts
    {
        sam_tag_dictionary{},
        sam_tag_dictionary{},
        sam_tag_dictionary{}
    };
};

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

TEST(read_header, sucess)
{
    std::string header_str
    {
        "@HD\tVN:1.0\tSO:coordinate\tSS:coordinate:queryname\tGO:none\n"
        "@PG\tID:qc\tPN:quality_control\tCL:qc -f file1\tDS:trim reads with low qual\tVN:1.0.0\n"
        "@PG\tID:novoalign\tPN:novoalign\tVN:V3.02.07\tCL:novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz -F STDFQ"
        " --Q2Off -t 400 -o SAM -c 10\tPP:qc\n"
        "@SQ\tSN:1\tLN:249250621\n"
        "@SQ\tSN:2\tLN:243199373\tAS:hs37d5\n"
        "@RG\tID:U0a_A2_L1\tPL:illumina\tPU:1\tLB:1\tSM:NA12878\n"
        "@RG\tID:U0a_A2_L2\tPL:illumina\tSM:NA12878\tPU:1\tLB:1\n"
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
    };

    std::istringstream istream(header_str);

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header{});

    alignment_file_format_sam format;
    alignment_file_input_options<dna5> options;

    ASSERT_NO_THROW((format.read(istream, options, std::ignore, std::ignore,
                                 header_ptr,  std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore)));

    EXPECT_EQ(header_ptr->format_version, "1.0");
    EXPECT_EQ(header_ptr->sorting, "coordinate");
    EXPECT_EQ(header_ptr->subsorting, "coordinate:queryname");
    EXPECT_EQ(header_ptr->grouping, "none");

    EXPECT_EQ(header_ptr->program_infos[0].id, "qc");
    EXPECT_EQ(header_ptr->program_infos[0].name, "quality_control");
    EXPECT_EQ(header_ptr->program_infos[0].version, "1.0.0");
    EXPECT_EQ(header_ptr->program_infos[0].description, "trim reads with low qual");
    EXPECT_EQ(header_ptr->program_infos[0].previous, "");
    EXPECT_EQ(header_ptr->program_infos[0].command_line_call, "qc -f file1");
    EXPECT_EQ(header_ptr->program_infos[1].id, "novoalign");
    EXPECT_EQ(header_ptr->program_infos[1].name, "novoalign");
    EXPECT_EQ(header_ptr->program_infos[1].version, "V3.02.07");
    EXPECT_EQ(header_ptr->program_infos[1].description, "");
    EXPECT_EQ(header_ptr->program_infos[1].previous, "qc");
    EXPECT_EQ(header_ptr->program_infos[1].command_line_call, "novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz -F"
                                                              " STDFQ --Q2Off -t 400 -o SAM -c 10");

    EXPECT_EQ((header_ptr->ref_dict["1"]), (std::tuple<uint32_t, std::string>{249250621u, ""}));
    EXPECT_EQ((header_ptr->ref_dict["2"]), (std::tuple<uint32_t, std::string>{243199373u, "AS:hs37d5"}));

    EXPECT_EQ(header_ptr->read_groups[0],
              (std::pair<std::string, std::string>{"U0a_A2_L1", "PL:illumina\tPU:1\tLB:1\tSM:NA12878"}));
    EXPECT_EQ(header_ptr->read_groups[1],
              (std::pair<std::string, std::string>{"U0a_A2_L2", "PL:illumina\tSM:NA12878\tPU:1\tLB:1"}));
}

TEST(read_header, errors)
{
    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header{});

    alignment_file_format_sam format;
    alignment_file_input_options<dna5> options;

    {
        std::string header_str
        {
            "@HD\tVN:1.0\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW((format.read(istream, options, std::ignore, std::ignore,
                                  header_ptr,  std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore)), parse_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\tSI:this is not a valid tag starting with S\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW((format.read(istream, options, std::ignore, std::ignore,
                                  header_ptr,  std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore)), parse_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@TT\tthis is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW((format.read(istream, options, std::ignore, std::ignore,
                                  header_ptr,  std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore)), parse_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@PG\tID:prog\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW((format.read(istream, options, std::ignore, std::ignore,
                                  header_ptr,  std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore,
                                  std::ignore, std::ignore, std::ignore)), parse_error);
    }
}

TEST_F(alignment_data, read_in_all_data)
{
    std::string file_in_str
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\taa:A:c\tAS:i:2\tff:f:3.1\tzz:Z:str\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\tbc:B:c,-3\tbC:B:C,3,200\tbs:B:s,-3,200,-300"
                                                                           "\tbS:B:S,300,40,500\tbi:B:i,-3,200,-66000"
                                                                           "\tbI:B:I,294967296\tbf:B:f,3.5,0.1,43.8\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[0]["aa"_tag] = 'c';
    tag_dicts[0]["ff"_tag] = 3.1f;
    tag_dicts[0]["zz"_tag] = "str";
    tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
    tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
    tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
    tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
    tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
    tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
    tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};

    std::unordered_map<std::string, size_t> ref_id_to_position{{ref_id, 0}};
    std::vector<dna5_vector> ref_sequences{ref_seq};
    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header{});
    alignment_file_header header;
    dna5_vector seq;
    std::string id;
    std::vector<phred42> qual;
    unsigned offset;
    std::string ref_id_in;
    unsigned ref_offset;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    unsigned flag;
    unsigned mapq;
    std::tuple<std::string, unsigned, unsigned> mate;
    sam_tag_dictionary tag_dict;

    std::istringstream istream(file_in_str);

    alignment_file_format_sam format;
    alignment_file_input_options<dna5> options;

    for (unsigned i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW((format.read(istream, options, ref_id_to_position,
                                     ref_sequences, header_ptr,
                                     seq, qual, id, offset, std::ignore,
                                     ref_id_in, ref_offset, alignment,
                                     flag, mapq, mate,
                                     tag_dict, std::ignore, std::ignore)));

        EXPECT_EQ(seq, seqs[i]);
        EXPECT_EQ(id, ids[i]);
        EXPECT_EQ(qual, quals[i]);
        EXPECT_EQ(offset, offsets[i]);
        EXPECT_EQ(ref_id_in, ref_id);
        EXPECT_EQ(ref_offset, ref_offsets[i]);
        EXPECT_EQ(get<0>(alignment), get<0>(alignments[i]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[i]));
        EXPECT_EQ(flag, flags[i]);
        EXPECT_EQ(mapq, mapqs[i]);
        EXPECT_EQ(mate, mates[i]);
        EXPECT_EQ(tag_dict, tag_dicts[i]);

        seq.clear();
        id.clear();
        qual.clear();
        offset = 0;
        ref_id_in.clear();
        ref_offset = 0;
        alignment = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{};
        flag = 0;
        mapq = 0;
        mate = std::tuple<std::string, unsigned, unsigned>{};
        tag_dict.clear();
    }
}

TEST_F(alignment_data, read_in_nothing)
{
    std::string file_in_str
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    std::istringstream istream(file_in_str);

    alignment_file_format_sam format;
    alignment_file_input_options<dna5> options;
    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header{});

    for (unsigned i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW((format.read(istream, options, std::ignore,
                                     std::ignore, header_ptr, std::ignore,
                                     std::ignore, std::ignore, std::ignore, std::ignore,
                                     std::ignore, std::ignore, std::ignore,
                                     std::ignore, std::ignore, std::ignore,
                                     std::ignore, std::ignore, std::ignore)));

    }
}

TEST_F(alignment_data, read_in_alignment_only)
{
    std::string file_in_str
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header{});
    std::unordered_map<std::string, size_t> ref_id_to_position{{ref_id, 0}};
    std::vector<dna5_vector> ref_sequences{ref_seq};
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    std::string ref_id_in;

    std::istringstream istream(file_in_str);

    alignment_file_format_sam format;
    alignment_file_input_options<dna5> options;

    /*with reference information*/
    for (unsigned i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW((format.read(istream, options, ref_id_to_position, ref_sequences, header_ptr,
                                     std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                     ref_id_in,
                                     std::ignore,
                                     alignment,
                                     std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore)));

        EXPECT_EQ(get<0>(alignment), get<0>(alignments[i]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[i]));

        alignment = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{};
        ref_id_in.clear();

    }

    /*no reference information*/

    gap_decorator_anchor_set dec{ranges::view::repeat_n(dna5{}, size_t{3}) | view::transform(detail::restrictor)};
    dec.insert_gap(dec.begin(), 2);

    using dummy_type = decltype(dec);
        // gap_decorator_anchor_set<decltype(ranges::view::repeat_n(dna5{}, size_t{}) | detail::restrict_access_view)>;
    std::pair<dummy_type, std::vector<gapped<dna5>>> alignment2;

    // auto dummy = ranges::view::repeat_n(dna5{}, size_t{}) | detail::restrict_access_view;

    debug_stream << alignment2 << std::endl;
// debug_stream << dummy.size() << std::endl;

    // istream = std::istringstream(file_in_str);
    // for (unsigned i = 0; i < 3; ++i)
    // {
    //     ASSERT_NO_THROW((format.read(istream, options,
    //                                  std::ignore, std::ignore,
    //                                  header_ptr,
    //                                  std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
    //                                  ref_id_in,
    //                                  std::ignore,
    //                                  alignment2,
    //                                  std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore)));

    //     EXPECT_EQ(get<1>(alignment2), get<1>(alignments[i]));

    //     alignment2 = std::pair<dummy_type, std::vector<gapped<dna5>>>{};
    //     ref_id_in.clear();

    // }
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

TEST_F(alignment_data, default_options_all_members_specified)
{
    alignment_file_format_sam format;

    alignment_file_output_options options;

    std::ostringstream ostream;

    std::unique_ptr<alignment_file_header> header_ptr(nullptr);

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    for (unsigned i = 0; i < 3; ++i)
        ASSERT_NO_THROW(( format.write(ostream, options, header_ptr,
                                       seqs[i], quals[i], ids[i], offsets[i], std::string{},
                                       ref_id, ref_offsets[i], alignments[i],
                                       flags[i], mapqs[i], mates[i],
                                       tag_dicts[i], 0u, 0u)));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(alignment_data, with_header)
{
    alignment_file_format_sam format;

    alignment_file_output_options options;

    std::ostringstream ostream;

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header);
    header_ptr->sorting = "unknown";
    header_ptr->grouping = "none";
    header_ptr->ref_dict["ref"] = {ref_seq.size(), ""};
    header_ptr->program_infos.push_back({"prog1", "cool_program", "", "", "", ""});
    header_ptr->comments.push_back("This is a comment.");

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp
    {
        "@HD\tVN:1.6\tSO:unknown\tGO:none\n"
        "@SQ\tSN:ref\tLN:34\n"
        "@PG\tID:prog1\tPN:cool_program\n"
        "@CO\tThis is a comment.\n"
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    for (unsigned i = 0; i < 3; ++i)
        ASSERT_NO_THROW(( format.write(ostream, options, header_ptr,
                                       seqs[i], quals[i], ids[i], offsets[i], std::string{},
                                       ref_id, ref_offsets[i], alignments[i],
                                       flags[i], mapqs[i], mates[i],
                                       tag_dicts[i], 0u, 0u)));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(alignment_data, format_errors)
{
    alignment_file_format_sam format;
    alignment_file_output_options options;

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header);
    header_ptr->ref_dict["ref"] = {ref_seq.size(), ""};

    std::ostringstream ostream;

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(( format.write(ostream, options, header_ptr,
                                   seqs[0], quals[0], ids[0], offsets[0], std::string{},
                                   "ref_id_that_does_not_exist", ref_offsets[0], alignments[0],
                                   flags[0], mapqs[0], mates[0],
                                   tag_dicts[0], 0u, 0u)),
                 format_error);
}
