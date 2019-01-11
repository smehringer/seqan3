// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

/*!\file
 * \brief Provides seqan3::alignment_file_input and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

// ----------------------------------------------------------------------------
// alignment_file_input_traits_concept
// ----------------------------------------------------------------------------

/*!\interface seqan3::alignment_file_input_traits_concept <>
 * \brief The requirements a traits_type for seqan3::alignment_file_input must meet.
 * \ingroup alignment
 */
/*!\name Requirements for seqan3::alignment_file_input_traits_concept
 * \brief You can expect these **member types** of all types that satisfy seqan3::alignment_file_input_traits_concept.
 * \memberof seqan3::alignment_file_input_traits_concept
 *
 * \details
 *
 * Note that the alphabet type of the seqan3::field::SEQ_QUAL cannot be specified directly, it is always
 * seqan3::qualified<sequence_alphabet, quality_alphabet> and the container type templates for
 * the field are those of seqan3::field::SEQ.
 *
 * \{
 */
/*!\typedef using sequence_alphabet
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::SEQ; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using sequence_legal_alphabet
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Intermediate alphabet for seqan3::field::SEQ; must satisfy seqan3::alphabet_concept and be convertible to
 * `sequence_alphabet`.
 *
 * \details
 *
 * This alphabet can be a superset of `sequence_alphabet` to allow conversion of some characters
 * without producing an error, e.g. if this is set to seqan3::dna15 and `sequence_alphabet` is set to seqan3::dna5,
 * 'M' will be an accepted character and automatically converted to 'N', while 'Z' will still be an illegal
 * character and produce an error.
 */
/*!\typedef using sequence_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of the seqan3::field::SEQ, a container template over `sequence_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using sequence_container_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::SEQ, a container template that can hold multiple
 * `sequence_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_alphabet
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::ID; must satisfy seqan3::alphabet_concept.
 */
/*!\typedef using id_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of the seqan3::field::ID, a container template over `id_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using id_container_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::ID, a container template that can hold multiple
 * `id_container`; must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using quality_alphabet
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Alphabet of the characters for the seqan3::field::QUAL; must satisfy seqan3::quality_concept.
 */
/*!\typedef using quality_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of the seqan3::field::QUAL, a container template over `quality_alphabet`;
 * must satisfy seqan3::sequence_container_concept.
 */
/*!\typedef using quality_container_container
 * \memberof seqan3::alignment_file_input_traits_concept
 * \brief Type template of a column of seqan3::field::QUAL, a container template that can hold multiple
 * `quality_container`; must satisfy seqan3::sequence_container_concept.
 */
//!\}
//!\cond
template <typename t>
concept alignment_file_input_traits_concept = requires (t v)
{
    // field::SEQ
    requires alphabet_concept<typename t::sequence_alphabet>;
    requires alphabet_concept<typename t::sequence_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
    requires sequence_container_concept<typename t::template sequence_container<typename t::sequence_alphabet>>;

    // field::ID
    requires alphabet_concept<typename t::id_alphabet>;
    requires sequence_container_concept<typename t::template id_container<typename t::id_alphabet>>;

    // field::QUAL
    requires quality_concept<typename t::quality_alphabet>;
    requires sequence_container_concept<typename t::template quality_container<typename t::quality_alphabet>>;

    // field::REF_SEQ
    requires alphabet_concept<typename t::ref_sequence_alphabet>;
    requires alphabet_concept<typename t::ref_sequence_legal_alphabet>;
    requires explicitly_convertible_to_concept<typename t::ref_sequence_legal_alphabet, typename t::ref_sequence_alphabet>;
    requires sequence_container_concept<typename t::template ref_sequence_container<typename t::ref_sequence_alphabet>>;

    // field::REF_ID
    requires alphabet_concept<typename t::ref_id_alphabet>;
    requires sequence_container_concept<typename t::template ref_id_container<typename t::ref_id_alphabet>>;

    // field::OFFSET
    requires std::UnsignedIntegral<typename t::offset_type>;

    // field::REF_OFFSET
    requires std::UnsignedIntegral<typename t::ref_offset_type>;

    // field::FLAG
    requires std::UnsignedIntegral<typename t::flag_type>;

    // field::MAPQ
    requires std::UnsignedIntegral<typename t::mapq_type>;

    // field::EVALUE
    requires arithmetic_concept<typename t::e_value_type>;

    // field::BITSCORE
    requires arithmetic_concept<typename t::bitscore_type>;

    // field::ALIGNMENT
    requires tuple_like_concept<typename t::template alignment_type<typename t::alignment_ref_sequence,
                                                                    typename t::alignment_query_sequence>>;
    requires std::tuple_size_v<typename t::template alignment_type<typename t::alignment_ref_sequence,
                                                                   typename t::alignment_query_sequence>> == 2;

    // field::MATE
    // the mate type cannot be configured as it will be set to
    // std::tuple<ref_id_container<ref_id_alphabet>, ref_offset_type, uint32_t>
};
//!\endcond

// ----------------------------------------------------------------------------
// alignment_file_input_default_traits
// ----------------------------------------------------------------------------

/*!\brief The default traits for seqan3::alignment_file_input
 * \implements alignment_file_input_traits_concept
 * \ingroup sequence
 *
 * \details
 *
 * If you wish to change a single or a few types from the default, just inherit from this class and
 * "overwrite" the respective type definitions.
 *
 * This example will make the file read into a smaller alphabet and a compressed container:
 *
 * ```cpp
 * struct my_traits : alignment_file_input_default_traits_dna
 * {
 *     using sequence_alphabet = dna4;                        // instead of dna5
 *
 *     template <typename alph>
 *     using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
 * };
 *
 * alignment_file_input<my_traits> fin{"/tmp/my.fasta"};
 *
 * //...
 * ```
 */
struct alignment_file_input_default_traits_sam
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::alignment_file_input_traits_concept.
     * \{
     */
    using sequence_alphabet                     = dna5;
    using sequence_legal_alphabet               = dna15;
    template <typename _sequence_alphabet>
    using sequence_container                    = std::vector<_sequence_alphabet>;

    using id_alphabet                           = char;
    template <typename _id_alphabet>
    using id_container                          = std::basic_string<_id_alphabet>;

    using quality_alphabet                      = phred42;
    template <typename _quality_alphabet>
    using quality_container                     = std::vector<_quality_alphabet>;

    using ref_sequence_alphabet                 = dna5;
    using ref_sequence_legal_alphabet           = dna15;
    template <typename _ref_sequence_alphabet>
    using ref_sequence_container                = std::vector<_ref_sequence_alphabet>;

    using ref_id_alphabet                       = char;
    template <typename _ref_id_alphabet>
    using ref_id_container                      = std::basic_string<_ref_id_alphabet>;

    using offset_type                           = uint32_t;
    using ref_offset_type                       = uint32_t;
    using flag_type                             = uint16_t;
    using mapq_type                             = uint32_t;
    using e_value_type                          = double;
    using bitscore_type                         = double;

    using alignment_ref_sequence                = std::vector<gapped<dna5>>;
    using alignment_query_sequence              = std::vector<gapped<dna5>>;
    template <typename _alignment_ref_sequence, typename _alignment_query_sequence>
    using alignment_type                        = std::pair<_alignment_ref_sequence, _alignment_query_sequence>;

    //!\}
};

// ----------------------------------------------------------------------------
// alignment_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup sequence
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 *                              seqan3::alignment_file_input_traits_concept.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::alignment_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::alignment_file_input_format_concept).
 * \tparam stream_char_type     The type of the underlying stream device(s); must model seqan3::char_concept.
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction supports reading four different fields:
 *
 *   1. seqan3::field::SEQ
 *   2. seqan3::field::ID
 *   3. seqan3::field::QUAL
 *   4. seqan3::field::SEQ_QUAL (sequence and qualities in one range)
 *
 * The first three fields are retrieved by default (and in that order). The last field may be selected to have
 * sequence and qualities directly stored in a more memory-efficient combined container. If you select the last
 * field you may not select seqan3::field::SEQ or seqan3::field::QUAL.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name and one for construction from
 * an existing stream and a known format. The first one automatically picks the format based on the extension
 * of the file name. The second can be used if you have a non-file stream, like std::cin or std::istringstream,
 * that you want to read from and/or if you cannot use file-extension based detection, but know that your input
 * file has a certain format.
 *
 * In most cases the template parameters are deduced completely automatically:
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta"}; // FastA with DNA sequences assumed, regular std::ifstream taken as stream
 * ```
 * Reading from an std::istringstream:
 * ```cpp
 * std::string input
 * {
 *     "> TEST1\n"
 *     "ACGT\n"
 *     "> Test2\n"
 *     "AGGCTGN\n"
 *     "> Test3\n"
 *     "GGAGTATAATATATATATATATAT\n"
 * };
 *
 * std::istringstream iss(input);
 *
 * alignment_file_input fin{std::move(iss), alignment_file_format_fasta{}};
 * //              ^ no need to specify the template arguments
 * ```
 *
 * Note that this is not the same as writing `alignment_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](http://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `alignment_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * ```cpp
 * alignment_file_input<alignment_file_default_traits_aa> fin{"/tmp/my.fasta"};
 * ```
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::alignment_file_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * ```cpp
 *
 *  // ... input had amino acid sequences
 * std::istringstream iss(input);
 *
 * alignment_file_input<alignment_file_default_traits_aa,
 *                  fields<field::SEQ, field::ID, field::QUAL>,
 *                  type_list<alignment_file_format_fasta>,
 *                  std::istringstream> fin{std::move(iss), alignment_file_format_fasta{}};
 * ```
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta"};
 *
 * for (auto & rec : fin)
 * {
 *     std::cout << "ID:  " << get<field::ID>(rec) << '\n';
 *     std::cout << "SEQ: " << (get<field::SEQ>(rec) | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     // a quality field also exists, but is not printed, because we know it's empty for FastA files.
 * }
 * ```
 *
 * In the above example, rec has the type \ref record_type which is a specialisation of seqan3::record and behaves
 * like an std::tuple (that's why we can access it via get). Instead of using the seqan3::field based interface on
 * the record, you could also use `std::get<0>` or even `std::get<dna4_vector>` to retrieve the sequence, but it is
 * not recommended, because it is more error-prone.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta"};
 *
 * using record_type = typename decltype(fin)::record_type;
 * std::vector<record_type> records;
 *
 * for (auto & rec : fin)
 *     records.push_back(std::move(rec));
 * ```
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta"};
 *
 * for (auto & [ seq, id, qual ] : fin)
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     std::cout << "SEQ: " << (seq | view::to_char) << '\n'; // sequence is converted to char on-the-fly
 *     // qual is empty for FastA files
 * }
 * ```
 *
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * alignment_file_input constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and QUAL (see above). Or to never actually read the QUAL, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta", fields<field::ID, field::SEQ_QUAL>{}};
 *
 * for (auto & [ id, seq_qual ] : fin) // note that the order is now different, "id" comes first, because it was specified first
 * {
 *     std::cout << "ID:  " << id << '\n';
 *     // sequence and qualities are part of the same vector, of type std::vector<dna5q>
 *     std::cout << "SEQ: "  << (seq | view::get<0> | view::to_char) << '\n'; // sequence string is extracted and converted to char  on-the-fly
 *     std::cout << "QUAL: " << (seq | view::get<1> | view::to_char) << '\n'; // quality string is extracted and converted to char  on-the-fly
 * }
 * ```
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * ```cpp
 * alignment_file_input fin{"/tmp/my.fasta"};
 *
 * auto minimum_length5_filter = view::filter([] (auto const & rec)
 * {
 *     return size(get<field::SEQ>(rec)) >= 5;
 * });
 *
 * for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
 * {
 *     // ...
 * }
 * ```
 *
 * ### End of file
 *
 * You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).
 *
 * ### Column-based reading
 *
 * The record-based interface treats the file as a range of tuples (the records), but in certain situations it
 * is desirable to read the file by field, i.e. column wise (tuple-of-ranges, instead of range-of-tuples).
 *
 * This interface is less flexible, but can save you copy operations in certain scenarios, given that
 * you have sufficient memory to load the entire file at once:
 *
 * ```cpp
 *
 * struct data_storage_t
 * {
 *     concatenated_sequences<dna5_vector>  sequences;
 *     concatenated_sequences<std::string>  ids;
 * };
 *
 * data_storage_t data_storage; // a global or globally used variable in your program
 *
 * // ... in your file reading function:
 *
 * alignment_file_input fin{"/tmp/my.fasta"};
 *
 * data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
 * data_storage.ids = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage
 * ```
 *
 * Note that for this to make sense, your storage data types need to be identical to the corresponding column types
 * of the file. If you require different column types you can specify you own traits, see
 * seqan3::alignment_file_input_traits_concept.
 *
 * ### Formats
 *
 * TODO give overview of formats, once they are all implemented
 */

template <
    alignment_file_input_traits_concept                      traits_type_        = alignment_file_input_default_traits_sam,
    detail::fields_concept                                   selected_field_ids_ = fields<field::SEQ,
                                                                                          field::ID,
                                                                                          field::OFFSET,
                                                                                          field::REF_SEQ,
                                                                                          field::REF_ID,
                                                                                          field::REF_OFFSET,
                                                                                          field::ALIGNMENT,
                                                                                          field::MAPQ,
                                                                                          field::QUAL,
                                                                                          field::FLAG,
                                                                                          field::MATE,
                                                                                          field::TAGS,
                                                                                          field::EVALUE,
                                                                                          field::BIT_SCORE,
                                                                                          field::HEADER_PTR>,
    detail::type_list_of_alignment_file_input_formats_concept valid_formats_     = type_list<alignment_file_format_sam
                                                                                         /*, ...*/>,
    char_concept                                             stream_char_type_   = char>
class alignment_file_input
{
public:
    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A traits type that defines aliases and template for storage of the fields.
    using traits_type           = traits_type_;
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids    = selected_field_ids_;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats         = valid_formats_;
    //!\brief Character type of the stream(s), usually `char`.
    using stream_char_type      = stream_char_type_;
    //!\}

    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of field::SEQ (default std::vector<seqan3::dna5>).
    using sequence_type            = typename traits_type::template sequence_container<
                                       typename traits_type::sequence_alphabet>;
    //!\brief The type of field::ID (default std::string by defaul).
    using id_type                  = typename traits_type::template id_container<
                                       typename traits_type::id_alphabet>;
    //!\brief The type of field::OFFSET (default: uint32_t).
    using offset_type              = typename traits_type::offset_type;
    //!\brief The type of field::REF_SEQ (default: std::vector<seqan3::dna5>).
    using ref_sequence_type        = typename traits_type::template ref_sequence_container<
                                       typename traits_type::ref_sequence_alphabet>;
    //!\brief The type of field::REF_ID (default: std::string).
    using ref_id_type              = typename traits_type::template ref_id_container<
                                       typename traits_type::ref_id_alphabet>;
    //!\brief The type of field::REF_OFFSET (default: uint32_t).
    using ref_offset_type          = typename traits_type::ref_offset_type;
    //!\brief The type of field::ALIGNMENT (default: std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>).
    using alignment_type           = typename traits_type::template alignment_type<
                                       typename traits_type::alignment_ref_sequence,
                                       typename traits_type::alignment_query_sequence>;
    //!\brief The type of field::MAPQ (default: uint32_t).
    using mapq_type                = typename traits_type::mapq_type;
    //!\brief The type of field::QUAL (default std::vector<seqan3::phred42>).
    using quality_type             = typename traits_type::template quality_container<
                                       typename traits_type::quality_alphabet>;
    //!\brief The type of field::FLAG (default: uint16_t).
    using flag_type                = typename traits_type::flag_type;
    //!\brief The type of field::EVALUE (default: double).
    using e_value_type             = typename traits_type::e_value_type;
    //!\brief The type of field::BITSCORE (default: double).
    using bitscore_type            = typename traits_type::bitscore_type;

    //!\brief The previously defined types aggregated in a seqan3::type_list.
    using field_types = type_list<sequence_type,
                                  id_type,
                                  offset_type,
                                  ref_sequence_type,
                                  ref_id_type,
                                  ref_offset_type,
                                  alignment_type,
                                  mapq_type,
                                  quality_type,
                                  flag_type,
                                  std::tuple<ref_id_type, ref_offset_type, uint32_t>, // mate
                                  sam_tag_dictionary, // tag type
                                  e_value_type,
                                  bitscore_type,
                                  std::unique_ptr<alignment_file_header>>;

    /*!\brief The subset of seqan3::field IDs that are valid for this file; order corresponds to the types in
     * \ref field_types.
     */
    using field_ids = fields<field::SEQ,
                             field::ID,
                             field::OFFSET,
                             field::REF_SEQ,
                             field::REF_ID,
                             field::REF_OFFSET,
                             field::ALIGNMENT,
                             field::MAPQ,
                             field::QUAL,
                             field::FLAG,
                             field::MATE,
                             field::TAGS,
                             field::EVALUE,
                             field::BIT_SCORE,
                             field::HEADER_PTR>;

    static_assert([] () constexpr
                  {
                      for (field f : selected_field_ids::as_array)
                          if (!field_ids::contains(f))
                              return false;
                      return true;
                  }(),
                  "You selected a field that is not valid for aligment files, please refer to the documentation "
                  "of alignment_file_input::field_ids for the accepted values.");

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type        = record<detail::select_types_with_ids_t<field_types, field_ids, selected_field_ids>,
                                      selected_field_ids>;
    //!\}

    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The value_type is the \ref record_type.
    using value_type        = record_type;
    //!\brief The reference type.
    using reference         = record_type &;
    //!\brief The const_reference type is void, because files are not const-iterable.
    using const_reference   = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<alignment_file_input>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = ranges::default_sentinel;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    alignment_file_input() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_input(alignment_file_input const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_input & operator=(alignment_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    alignment_file_input(alignment_file_input &&) = default;
    //!\brief Move assignment is defaulted.
    alignment_file_input & operator=(alignment_file_input &&) = default;
    //!\brief Destructor is defaulted.
    ~alignment_file_input() = default;

    /*!\brief Construct from filename.
     * \param[in] filename    Path to the file you wish to open.
     * \param[in] fields_tag  A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    alignment_file_input(filesystem::path filename,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        // open stream
        if (!primary_stream->good())
            throw file_open_error{"Could not open file for reading."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);

        // buffer first record
        read_next_record();
    }
    /* NOTE(h-2): Curiously we do not need a user-defined deduction guide for the above constructor.
     * A combination of default template parameters and auto-deduction guides works as expected,
     * independent of whether the second/optional parameter is specified or not, i.e. it is possible
     * to auto-deduct and overwrite a single template parameter out of the four if the optional parameter
     * is specified and use the default otherwise.
     */

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::alignment_file_input_format_concept.
     * \param[in] stream     The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <istream_concept2 stream_t, alignment_file_input_format_concept file_format>
    alignment_file_input(stream_t                 & stream,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}, format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    //!\overload
    template <istream_concept2 stream_t, alignment_file_input_format_concept file_format>
    alignment_file_input(stream_t                && stream,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // buffer first record
        read_next_record();
    }

    /*!\brief Construct from filename and give additional reference information.
     * \toaram    ref_ids_t       The range type that stores the reference ids.
     * \toaram    ref_sequences_t The range type that stores the reference information.
     * \param[in] filename        Path to the file you wish to open.
     * \param[in] ref_ids         A range containing the reference ids that correspond to the SAM/BAM file.
     * \param[in] ref_sequences   A range containing the reference sequences that correspond to the SAM/BAM file.
     * \param[in] fields_tag      A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name and the reference information, you may specify a custom seqan3::fields type which
     * may be easier than defining all the template parameters.
     *
     * The reference information given by the ids (names) and sequences will be used to construct a proper alignment
     * when reading in SAM or BAM files. If you are not interested in the fulll alignment you do not need to specify
     * those information.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <std::ranges::ForwardRange ref_ids_t, std::ranges::ForwardRange ref_sequences_t>
        requires std::Same<value_type_t<ref_sequences_t>, ref_sequence_type> &&
                 std::Same<value_type_t<ref_ids_t>, ref_id_type>
    alignment_file_input(filesystem::path filename,
                         ref_ids_t const & ref_ids,
                         ref_sequences_t const & ref_sequences,
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new std::ifstream{filename, std::ios_base::in | std::ios::binary}, stream_deleter_default}
    {
        // open stream
        if (!primary_stream->good())
            throw file_open_error{"Could not open file for reading."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);

        // initialize reference information
        set_references(ref_ids, ref_sequences);

        // buffer first record
        read_next_record();
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam    file_format     The format of the file in the stream,
     *                            must satisfy seqan3::alignment_file_input_format_concept.
     * \toaram    ref_ids_t       The range type that stores the reference ids.
     * \toaram    ref_sequences_t The range type that stores the reference information.
     * \param[in] stream          The stream to operate on; must be derived of std::basic_istream.
     * \param[in] ref_ids         A range containing the reference ids that correspond to the SAM/BAM file.
     * \param[in] ref_sequences   A range containing the reference sequences that correspond to the SAM/BAM file.
     * \param[in] format_tag      The file format tag.
     * \param[in] fields_tag      A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    template <istream_concept2 stream_t,
              std::ranges::ForwardRange ref_ids_t,
              std::ranges::ForwardRange ref_sequences_t,
              alignment_file_input_format_concept file_format>
        requires std::Same<value_type_t<ref_sequences_t>, ref_sequence_type> &&
                 std::Same<value_type_t<ref_ids_t>, ref_id_type>
    alignment_file_input(stream_t                 & stream,
                         ref_ids_t const & ref_ids,
                         ref_sequences_t const & ref_sequences,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{&stream, stream_deleter_noop}, format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // initialize reference information
        set_references(ref_ids, ref_sequences);

        // buffer first record
        read_next_record();
    }

    //!\overload
    template <istream_concept2 stream_t,
              std::ranges::ForwardRange ref_ids_t,
              std::ranges::ForwardRange ref_sequences_t,
              alignment_file_input_format_concept file_format>
        requires std::Same<value_type_t<ref_sequences_t>, ref_sequence_type> &&
                 std::Same<value_type_t<ref_ids_t>, ref_id_type>
    alignment_file_input(stream_t                && stream,
                         ref_ids_t const & ref_ids,
                         ref_sequences_t const & ref_sequences,
                         file_format        const & SEQAN3_DOXYGEN_ONLY(format_tag),
                         selected_field_ids const & SEQAN3_DOXYGEN_ONLY(fields_tag) = selected_field_ids{}) :
        primary_stream{new stream_t{std::move(stream)}, stream_deleter_default},
        format{file_format{}}
    {
        static_assert(meta::in<valid_formats, file_format>::value,
                      "You selected a format that is not in the valid_formats of this file.");

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // initialize reference information
        set_references(ref_ids, ref_sequences);

        // buffer first record
        read_next_record();
    }
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     *
     * Equals end() if the file is at end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {*this};
    }

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() noexcept
    {
        return {};
    }

    /*!\brief Return the record we are currently at in the file.
     * \returns A reference to the currently buffered record.
     *
     * This function returns a reference to the currently buffered record, it is identical to dereferencing begin(),
     * but begin also always points to the current record on single pass input ranges:
     *
     * ```cpp
     * alignment_file_input fin{"/tmp/my.fasta"};
     * auto it = begin(fin);
     *
     * // the following are equivalent:
     * auto & rec0 = *it;
     * auto & rec1 = fin.front();
     *
     * // both become invalid after incrementing "it"!
     * ```
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * ```cpp
     * alignment_file_input fin{"/tmp/my.fasta"};
     *
     * auto rec0 = std::move(fin.front());
     * ```
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference front() noexcept
    {
        return record_buffer;
    }
    //!\}

    //!\brief The options are public and its members can be set directly.
    alignment_file_input_options<typename traits_type::sequence_legal_alphabet> options;

    /*!\brief Access the file's header.
     *
     * \details
     *
     * ### Example
     *
     * \snippet test/snippet/io/alignment_file/output.cpp set_header
     *
     * \sa seqan3::alignment_file_header
     */
    alignment_file_header & header()
    {
        return *header_ptr;
    }

protected:
    //!\privatesection

    //!\brief The file header object.
    std::unique_ptr<alignment_file_header> header_ptr{new alignment_file_header()};

    /*!\name Data buffers
     * \{
     */
    //!\brief Buffer for a single record.
    record_type record_buffer;
    //!\}

    /*!\name Stream / file access
     * \{
     */
    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_istream<stream_char_type>,
                                         std::function<void(std::basic_istream<stream_char_type>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<stream_char_type> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<stream_char_type> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief File is at position 1 behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type = detail::transfer_template_args_onto_t<valid_formats, std::variant>;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    //!\}

    //!\brief The reference ids corresponding to the alignment records.
    std::unordered_map<ref_id_type, size_t> ref_ids_to_position;

    //!\brief The reference sequences corresponding to the alignment records.
    std::vector<ref_sequence_type> const * reference_sequences_ptr{nullptr};

    /*!\brief Set reference information (only important for SAM/BAM reading).
     *
     * \details
     *
     * The SAM format only provides semi-alignments because the reference sequence
     * is not stored explicitly. In order to be able to read in full alignment,
     * additional reference information can be given to the alignment file.
     * Note that the reference ids (names) must correspond to the exact spelling
     * in the SAM/BAM file otherwise an exception will be thrown when reading.
     *
     * ### Example
     *
     * \todo example
     */
    template <std::ranges::ForwardRange ref_ids_t, std::ranges::ForwardRange ref_sequences_t>
    void set_references(ref_ids_t const & ref_ids, ref_sequences_t const & ref_sequences)
    {
        reference_sequences_ptr = &ref_sequences;

        assert(ref_ids.size() == ref_sequences.size());

        // initialise reference map if non-empty
        for (size_t idx = 0; idx < ref_ids.size(); ++idx)
            ref_ids_to_position[ref_ids[idx]] = idx;
    }

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        if (at_end)
            return;

        // clear the record
        record_buffer.clear();

        // at end if we could not read further
        if (secondary_stream->eof())
        {
            at_end = true;
            return;
        }

        auto call_read_func = [this] (auto & map_info, auto & ref_seq_info)
            {
                std::visit([&] (alignment_file_input_format_concept & f)
                {
                    f.read(*secondary_stream,
                           options,
                           map_info,
                           ref_seq_info,
                           detail::get_or_ignore<field::HEADER_PTR>(record_buffer),
                           detail::get_or_ignore<field::SEQ>(record_buffer),
                           detail::get_or_ignore<field::QUAL>(record_buffer),
                           detail::get_or_ignore<field::ID>(record_buffer),
                           detail::get_or_ignore<field::OFFSET>(record_buffer),
                           detail::get_or_ignore<field::REF_SEQ>(record_buffer),
                           detail::get_or_ignore<field::REF_ID>(record_buffer),
                           detail::get_or_ignore<field::REF_OFFSET>(record_buffer),
                           detail::get_or_ignore<field::ALIGNMENT>(record_buffer),
                           detail::get_or_ignore<field::MAPQ>(record_buffer),
                           detail::get_or_ignore<field::FLAG>(record_buffer),
                           detail::get_or_ignore<field::MATE>(record_buffer),
                           detail::get_or_ignore<field::TAGS>(record_buffer),
                           detail::get_or_ignore<field::EVALUE>(record_buffer),
                           detail::get_or_ignore<field::BIT_SCORE>(record_buffer));

                }, format);
            }; // saves some code duplication

        assert(!format.valueless_by_exception());

        if (reference_sequences_ptr)
            call_read_func(ref_ids_to_position, *reference_sequences_ptr);
        else
            call_read_func(std::ignore, std::ignore);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_file_input
 * \{
 */
template <istream_concept2                    stream_type,
          alignment_file_input_format_concept file_format,
          detail::fields_concept              selected_field_ids>
alignment_file_input(stream_type && stream,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

template <istream_concept2                    stream_type,
          alignment_file_input_format_concept file_format,
          detail::fields_concept              selected_field_ids>
alignment_file_input(stream_type & stream,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

// template <std::ranges::ForwardRange           ref_ids_t,
//           std::ranges::ForwardRange           ref_sequences_t,
//           detail::fields_concept              selected_field_ids>
// alignment_file_input(filesystem::path,
//                      ref_ids_t const &,
//                      ref_sequences_t const &,
//                      selected_field_ids const &)
//     -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
//                             selected_field_ids,
//                             type_list<file_format>,
//                             char>;

template <istream_concept2                    stream_type,
          std::ranges::ForwardRange           ref_ids_t,
          std::ranges::ForwardRange           ref_sequences_t,
          alignment_file_input_format_concept file_format,
          detail::fields_concept              selected_field_ids>
alignment_file_input(stream_type && stream,
                     ref_ids_t const &,
                     ref_sequences_t const &,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;

template <istream_concept2                    stream_type,
          std::ranges::ForwardRange           ref_ids_t,
          std::ranges::ForwardRange           ref_sequences_t,
          alignment_file_input_format_concept file_format,
          detail::fields_concept              selected_field_ids>
alignment_file_input(stream_type & stream,
                     ref_ids_t const &,
                     ref_sequences_t const &,
                     file_format const &,
                     selected_field_ids const &)
    -> alignment_file_input<typename alignment_file_input<>::traits_type,       // actually use the default
                            selected_field_ids,
                            type_list<file_format>,
                            typename std::remove_reference_t<stream_type>::char_type>;
//!\}

} // namespace seqan3

// ------------------------------------------------------------------
// std-overloads for the tuple-like interface
// ------------------------------------------------------------------

namespace std
{
//!\brief std::tuple_size overload for column-like access. [metafunction specialisation for seqan3::alignment_file_input]
template <seqan3::alignment_file_input_traits_concept                       traits_type,
          seqan3::detail::fields_concept                                    selected_field_ids,
          seqan3::detail::type_list_of_alignment_file_input_formats_concept valid_formats,
          seqan3::char_concept                                              stream_char_t>
struct tuple_size<seqan3::alignment_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
{
    //!\brief The value equals the number of selected fields in the file.
    static constexpr size_t value = selected_field_ids::as_array.size();
};

//!\brief std::tuple_element overload for column-like access. [metafunction specialisation for seqan3::alignment_file_input]
template <size_t                                                            elem_no,
          seqan3::alignment_file_input_traits_concept                       traits_type,
          seqan3::detail::fields_concept                                    selected_field_ids,
          seqan3::detail::type_list_of_alignment_file_input_formats_concept valid_formats,
          seqan3::char_concept                                              stream_char_t>
struct tuple_element<elem_no, seqan3::alignment_file_input<traits_type, selected_field_ids, valid_formats, stream_char_t>>
    : tuple_element<elem_no, typename seqan3::alignment_file_input<traits_type,
                                                               selected_field_ids,
                                                               valid_formats,
                                                               stream_char_t>::file_as_tuple_type>
{};

} // namespace std
