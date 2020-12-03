// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::io_cfg::default_configuration configuration object.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/io/configuration/all.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>

namespace seqan3
{

/*!\brief The default traits for seqan3::sequence_file_input
 * \implements sequence_file_input_traits
 * \ingroup sequence
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_trait_overwrite.cpp
 */
struct sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::dna5.
    using sequence_alphabet                 = dna5;

    //!\brief The legal sequence alphabet for parsing is seqan3::dna15.
    using sequence_legal_alphabet           = dna15;

    //!\brief The type of a DNA sequence is std::vector.
    template <typename _sequence_alphabet>
    using sequence_container                = std::vector<_sequence_alphabet>;

    //!\brief The alphabet for an identifier string is char.
    using id_alphabet                       = char;

    //!\brief The string type for an identifier is std::basic_string.
    template <typename _id_alphabet>
    using id_container                      = std::basic_string<_id_alphabet>;

    //!\brief The alphabet for a quality annotation is seqan3::phred42.
    using quality_alphabet                  = phred42;

    //!\brief The string type for a quality annotation is std::vector.
    template <typename _quality_alphabet>
    using quality_container                 = std::vector<_quality_alphabet>;

    //!\}
};

//!\brief A traits type that specifies input as amino acids.
//!\ingroup sequence
struct sequence_file_input_default_traits_aa : sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    //!\brief The sequence alphabet is seqan3::aa27.
    using sequence_alphabet = aa27;

    //!\brief The legal sequence alphabet for parsing is seqan3::aa27.
    using sequence_legal_alphabet = aa27;
    //!\}
};

} // namespace seqan3

namespace seqan3::io_cfg
{

configuration sequence_file_default_configuration = select_fields<field::seq, field::id, field::qual>
                                                  | select_traits<sequence_file_input_default_traits_dna>
                                                  | select_formats<format_embl,
                                                                   format_fasta,
                                                                   format_fastq,
                                                                   format_genbank,
                                                                   format_sam>;

} // namespace seqan3::io_cfg
