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
 * \brief Provides seqan3::sequence_file_input_format_concept and auxiliary classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

namespace seqan3
{

/*!\interface seqan3::sequence_file_input_format_concept <>
 * \brief The generic concept for sequence file in formats.
 * \ingroup sequence
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept alignment_file_input_format_concept =
    requires (t                                                               & v,
              std::ifstream                                                   & stream,
              alignment_file_input_options<dna5>                              & options,
              std::unordered_map<std::string, size_t>                         & ids_map,
              std::vector<dna5_vector>                                        & ref_sequences,
              std::unique_ptr<alignment_file_header>                          & header_ptr,
              dna5_vector                                                     & seq,
              std::vector<phred42>                                            & qual,
              std::string                                                     & id,
              uint8_t                                                         & offset,
              dna5_vector                                                     & ref_seq,
              std::string                                                     & ref_id,
              uint32_t                                                        & ref_offset,
              std::pair<std::vector<gapped<dna4>>, std::vector<gapped<dna4>>> & align,
              uint16_t                                                        & flag,
              uint16_t                                                        & mapq,
              std::tuple<std::string, uint32_t, uint32_t>                     & mate, // mate_ref_id, mate_ref_offset, mate_tlen
              sam_tag_dictionary                                              & tag_dict,
              double                                                          & e_value,
              double                                                          & bit_score)
{
    t::file_extensions;
    // std::Same<decltype(t::file_extensions), std::vector<std::string>>;

    { v.read(stream,
             options,
             ids_map,
             ref_sequences,
             header_ptr,
             seq,
             qual,
             id,
             offset,
             ref_seq,
             ref_id,
             ref_offset,
             align,
             flag,
             mapq,
             mate,
             tag_dict,
             e_value,
             bit_score
             ) } -> void;
};
//!\endcond

/*!\name Requirements for seqan3::alignment_file_input_format_concept
 * \brief You can expect these **members** on all types that implement seqan3::alignment_file_input_format_concept.
 * \memberof seqan3::alignment_file_input_format_concept
 * \{
 */

/*!\fn void read(stream_type & stream, seqan3::alignment_file_input_options const & options, seq_type & alignment,
 *               id_type & id, qual_type & qualities)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \memberof seqan3::alignment_file_input_format_concept
 * \tparam stream_type      Input stream, must satisfy seqan3::istream_concept with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ input; must satisfy std::ranges::OutputRange
 * over a seqan3::alphabet_concept.
 * \tparam id_type          Type of the seqan3::field::ID input; must satisfy std::ranges::OutputRange
 * over a seqan3::alphabet_concept.
 * \tparam qual_type        Type of the seqan3::field::QUAL input; must satisfy std::ranges::OutputRange
 * over a seqan3::quality_concept.
 * \param[in,out] stream    The input stream to read from.
 * \param[in]     options   File specific options passed to the format.
 * \param[out]    alignment  The buffer for seqan3::field::SEQ input, i.e. the "alignment".
 * \param[out]    id        The buffer for seqan3::field::ID input, e.g. the header line in FastA.
 * \param[out]    qualities The buffer for seqan3::field::QUAL input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields.
 *     [this is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 */
 /*!\var static inline std::vector<std::string> seqan3::alignment_file_input_format_concept::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::alignment_file_input_format_concept [default is false].
 * \ingroup core
 * \see seqan3::type_list_of_alignment_file_input_formats_concept
 */
template <typename t>
constexpr bool is_type_list_of_alignment_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::alignment_file_input_format_concept [overload].
 * \ingroup core
  * \see seqan3::type_list_of_alignment_file_input_formats_concept
 */
template <typename ... ts>
constexpr bool is_type_list_of_alignment_file_input_formats_v<type_list<ts...>> =
    (alignment_file_input_format_concept<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::alignment_file_input_format_concept.
 * \ingroup core
 * \see seqan3::is_type_list_of_alignment_file_formats_v
 */
template <typename t>
concept type_list_of_alignment_file_input_formats_concept = is_type_list_of_alignment_file_input_formats_v<t>;

} // namespace seqan3::detail
