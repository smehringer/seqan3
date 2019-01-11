// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::alignment_file_format_sam class.
 * \author Svenja Mehringer <avenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <vector>

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/alignment_file/detail.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/charconv>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/std/view/subrange.hpp>

namespace seqan3
{

/*!\brief       The SAM format.
 * \implements  alignment_file_format_concept
 * \ingroup     alignment_file
 *
 * \details
 *
 * ### Introduction
 *
 * SAM is often used for storing alignments of several read sequences against one
 * or more reference sequences. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)) for an
 * introduction of the format or look into the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 * **SeqAn implements version 1.6 of the SAM specification**.
 *
 * ### Fields
 *
 * The SAM format provides the following fields:
 * seqan3::field::ALIGNMENT, seqan3::field::SEQ, seqan3::field::QUAL,
 * seqan3::field::ID, seqan3::field::REF_SEQ, seqan3::field::REF_ID
 * seqan3::field::REF_OSSFET, seqan3::field::OFFSET, seqan3::field::FLAG,
 * seqan3::field::MAPQ and seqan3::field::MATE.
 * In addition there is the seqan3::field::HEADER_PTR, which is usually not set but
 * needed to provide the range-based functionality of the file.
 *
 * **None of the fields are required** when writing but will be defaulted
 * to '0' for numeric fields and '*' for other fields.
 *
 * ### SAM format columns -> fields
 *
 * As many users will be accustomed to the columns of the SAM format, here is a
 * mapping of the common SAM format columns to the SeqAn3 record fields:
 *
 * | #  | SAM Column ID |  FIELD name                                       |
 * |:--:|:--------------|:--------------------------------------------------|
 * | 1  | QNAME         | seqan3::field::ID                                 |
 * | 2  | FLAG          | seqan3::field::FLAG                               |
 * | 3  | RNAME         | seqan3::field::REF_ID                             |
 * | 4  | POS           | seqan3::field::REF_OFFSET                         |
 * | 5  | MAPQ          | seqan3::field::MAPQ                               |
 * | 6  | CIGAR         | implicilty stored in seqan3::field::ALIGNMENT     |
 * | 7  | RNEXT         | seqan3::field::MATE (tuple pos 0)                 |
 * | 8  | PNEXT         | seqan3::field::MATE (tuple pos 1)                 |
 * | 9  | TLEN          | seqan3::field::MATE (tuple pos 2)                 |
 * | 10 | SEQ           | seqan3::field::SEQ                                |
 * | 11 | QUAL          | seqan3::field::QUAL                               |
 *
 * The (read sequence/query) **OFFSET** will be required to store the soft
 * clipping information at the read start (end clipping will be automatically
 * deduced by how much the read sequence length + offset is larger than the
 * alignment length).
 *
 * Note: SeqAn currently does not support hard clipping. When reading SAM,
 * hard-clipping is discarded; but the resulting alignment/sequence combination
 * is still valid.
 *
 * ### Format Check
 *
 * The format checks are implemented according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)
 * in order to ensure correct SAM file output.
 *
 * If a non-recoverable format violation is encountered on reading, or you specify
 * invalid values/combinations when writing, seqan3::format_error is thrown.
 *
 * Note: All sequence like fields in SAM (e.g. field::SEQ) are truncated at the
 *       the first white space character (see seqan3::is_space) to ensure a
 *       correct format.
 *
 * ### Header implementation
 *
 * The SAM header (if present) is read/written once in the beginning, before the
 * first record is read/written.
 */
class alignment_file_format_sam
{
public:
    /*!\name Constructors, destructor and assignment
     * \brief Rule of five explicitly defaulted.
     * \{
     */
    alignment_file_format_sam() = default;
    alignment_file_format_sam(alignment_file_format_sam const &) = delete;
    alignment_file_format_sam & operator=(alignment_file_format_sam const &) = delete;
    alignment_file_format_sam(alignment_file_format_sam &&) = default;
    alignment_file_format_sam & operator=(alignment_file_format_sam &&) = default;
    ~alignment_file_format_sam() = default;
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "sam" },
    };

    //!\copydoc AlignmentFileInputFormat::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              typename ref_id_to_pos_map_type,
              typename ref_seqs_type,
              typename header_type,
              typename seq_type,
              typename id_type,
              typename offset_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename align_type,
              typename flag_type,
              typename mapq_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void read(stream_type                                             & stream,
              alignment_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
              ref_id_to_pos_map_type                                  & ref_id_to_pos_map,
              ref_seqs_type                                           & ref_seqs,
              header_type                                             & header_ptr,
              seq_type                                                & seq,
              qual_type                                               & qual,
              id_type                                                 & id,
              offset_type                                             & offset,
              ref_seq_type                                            & SEQAN3_DOXYGEN_ONLY(ref_seq),
              ref_id_type                                             & ref_id,
              ref_offset_type                                         & ref_offset,
              align_type                                              & align,
              flag_type                                               & flag,
              mapq_type                                               & mapq,
              mate_type                                               & mate,
              tag_dict_type                                           & tag_dict,
              e_value_type                                            & SEQAN3_DOXYGEN_ONLY(e_value),
              bit_score_type                                          & SEQAN3_DOXYGEN_ONLY(bit_score))
    {
        using stream_buf_t = std::istreambuf_iterator<typename stream_type::char_type>;
        auto stream_view = view::subrange<decltype(stream_buf_t{stream}), decltype(stream_buf_t{})>
                              {stream_buf_t{stream}, stream_buf_t{}};

        auto next = view::take_until_or_throw(is_char<'\t'>);

        // these variables need to be stored to compute the ALIGNMENT
        std::conditional_t<!detail::decays_to_ignore_v<ref_offset_type>, ref_offset_type, size_t> ref_offset_tmp{};
        [[maybe_unused]] std::string ref_id_tmp{}; // in case reference information is given
        [[maybe_unused]] std::conditional_t<!detail::decays_to_ignore_v<offset_type>, offset_type, size_t> offset_tmp{};
        [[maybe_unused]] size_t soft_clipping_end{};
        [[maybe_unused]] std::vector<std::pair<char, size_t>> cigar{};
        [[maybe_unused]] size_t ref_length{0}, seq_length{0}; // length of aligned part for ref and query

        // Header
        // -------------------------------------------------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<header_type>)
        {
            if (!header_was_read && (header_ptr != nullptr) && is_char<'@'>(*begin(stream_view)))
            {
                read_header(stream_view, header_ptr);
                header_was_read = true;
            }
            else
            {
                while (is_char<'@'>(*begin(stream_view)))
                {
                    detail::consume(stream_view | view::take_line);
                    ++begin(stream_view); // skip newline
                }
            }
        }
        else
        {
            while (is_char<'@'>(*begin(stream_view)))
            {
                detail::consume(stream_view | view::take_line);
                ++begin(stream_view); // skip newline
            }
        }

        // Fields 1-5: ID FLAG REF_ID REF_OFFSET MAPQ
        // -------------------------------------------------------------------------------------------------------------
        read_field(stream_view | next, id);
        ++begin(stream_view); // skip tab

        read_field(stream_view | next, flag);
        ++begin(stream_view);

        if constexpr (detail::decays_to_ignore_v<ref_id_type> &&
                      !detail::decays_to_ignore_v<align_type> && !detail::decays_to_ignore_v<ref_id_to_pos_map_type>)
            read_field(stream_view | next, ref_id_tmp);
        else
            read_field(stream_view | next, ref_id);
        ++begin(stream_view);

        read_field(stream_view | next, ref_offset_tmp);
        ++begin(stream_view);
        assert(ref_offset_tmp > 0);
        --ref_offset_tmp;            // SAM format is 1-based but SeqAn operates 0-based.
        ref_offset = ref_offset_tmp;

        read_field(stream_view | next, mapq);
        ++begin(stream_view);

        // Field 6: CIGAR
        // -------------------------------------------------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<align_type>)
        {
            cigar = parse_cigar(offset_tmp, soft_clipping_end, stream_view | next);

            // Compute length of aligned reference and query sequence in one go.
            // The length of query sequence is needed when parsing the sequence (Field 10)
            for (auto [cigar_op, cigar_count] : cigar)
            {
                if (is_char<'M'>(cigar_op))
                {
                    ref_length += cigar_count;
                    seq_length += cigar_count;
                }
                else if (is_char<'D'>(cigar_op))
                {
                    ref_length += cigar_count;
                }
                else if (is_char<'I'>(cigar_op))
                {
                    seq_length += cigar_count;
                }
            }

            offset = offset_tmp;
        }
        else
        {
            detail::consume(stream_view | next);
        }
        ++begin(stream_view);

        // Field 7-9: (RNEXT PNEXT TLEN)=MATE
        // -------------------------------------------------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<mate_type>)
        {
            read_field(stream_view | next, get<0>(mate)); // RNEXT
            ++begin(stream_view);

            read_field(stream_view | next, get<1>(mate)); // PNEXT
            ++begin(stream_view);

            read_field(stream_view | next, get<2>(mate)); // TLEN
            ++begin(stream_view);
        }
        else
        {
            for (size_t i = 0; i < 3; ++i)
            {
                detail::consume(stream_view | next);
                ++begin(stream_view);
            }
        }

        // Field 10: Sequence
        // -------------------------------------------------------------------------------------------------------------
        auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
        auto seq_stream = stream_view | next
                                      | view::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                        {
                                            if (!is_legal_alph(c))
                                                throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                  is_legal_alph.msg.string() +
                                                                  " evaluated to false on " +
                                                                  detail::make_printable(c)};
                                            return c;
                                        });

        if constexpr (!detail::decays_to_ignore_v<align_type> && detail::decays_to_ignore_v<seq_type>)
        {
            static_assert(sequence_container_concept<std::remove_reference_t<decltype(get<1>(align))>>,
                          "If you want to read ALIGNMENT but not SEQ, the alignment"
                          " object must store a sequence container at the second position.");

            detail::consume(seq_stream | view::take_exactly_or_throw(offset_tmp));        // skip soft clipping at begin

            for (auto chr : seq_stream | view::take_exactly_or_throw(seq_length))
                get<1>(align).push_back(value_type_t<decltype(get<1>(align))>{}.assign_char(chr));
            // read_field(seq_stream | view::take_exactly_or_throw(seq_length), get<1>(align)); why does this not work??

            detail::consume(seq_stream | view::take_exactly_or_throw(soft_clipping_end)); // skip soft clipping at end
            assert(seq_stream.begin() == seq_stream.end());
        }
        else
        {
            read_field(seq_stream, seq);

            if constexpr (!detail::decays_to_ignore_v<align_type>)
            {
                // TODO instead of copying here, you can maybe gap decorate a subrange view
                for (auto it = seq.begin() + offset_tmp; it != seq.end() - soft_clipping_end; ++it)
                    get<1>(align).push_back(value_type_t<decltype(get<1>(align))>{*it});
            }
        }
        ++begin(stream_view);

        // Field 11:  Quality
        // -------------------------------------------------------------------------------------------------------------
        read_field(stream_view | view::take_until_or_throw(is_char<'\t'> || is_char<'\n'>), qual);

        // All remaining optional fields if any: SAM tags dictionary
        // -------------------------------------------------------------------------------------------------------------
        while (is_char<'\t'>(*begin(stream_view))) // read all tags if present
        {
            ++begin(stream_view);
            read_field(stream_view | view::take_until_or_throw(is_char<'\t'> || is_char<'\n'>), tag_dict);
        }
        ++begin(stream_view);                      // skip newline

        // DONE READING - wrap up
        // -------------------------------------------------------------------------------------------------------------
        // make sure "buffer at end" implies "stream at end"
        if ((stream_buf_t{stream} == stream_buf_t{}) && (!stream.eof()))
            stream.get(); // triggers error in stream and sets eof

        // Alignment object construction
        // Note that the query sequence in get<1>(align) has already been filled while reading Field 10.
        if constexpr (!detail::decays_to_ignore_v<align_type>)
        {
            if constexpr (!detail::decays_to_ignore_v<ref_id_to_pos_map_type>)
            {
                static_assert(!detail::decays_to_ignore_v<ref_seqs_type>,
                              "If you pass the id_to_ref_map you must also pass ref_seqs.");

                size_t pos; // get index for ref_seqs
                try
                {
                    if constexpr (detail::decays_to_ignore_v<ref_id_type>)
                        pos = ref_id_to_pos_map.at(ref_id_tmp);
                    else
                        pos = ref_id_to_pos_map.at(ref_id);
                }
                catch (std::out_of_range & e)
                {
                    throw parse_error(std::string("[SAM PARSE ERROR] An id in the SAM record cannot be found in the "
                                                  "given list of reference ids."));
                }

                if (ref_offset_tmp + ref_length > ref_seqs[pos].size())
                    throw parse_error("[SAM PARSE ERROR] The alignment length exceeds the reference length.");

                // copy over aligned reference sequence part
                for (auto it = ref_seqs[pos].begin() + ref_offset_tmp;
                     it != ref_seqs[pos].begin() + ref_offset_tmp + ref_length; ++it)
                    get<0>(align).push_back(value_type_t<decltype(get<0>(align))>{*it});
            }
            else
            {
                // std::cerr << "A dummy should be created now" << std::endl;
                // create a dummy sequence of length ref_length and decorate it with a gap_decorator
            }

            // insert gaps according to the cigar information
            detail::alignment_from_cigar(align, cigar);
        }
    }

    //!\copydoc AlignmentFileOutputFormat::write
    template <typename stream_type,
              typename seq_type,
              typename id_type,
              typename offset_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename align_type,
              typename flag_type,
              typename mapq_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void write(stream_type                            &  stream,
               alignment_file_output_options const    &  options,
               std::unique_ptr<alignment_file_header> &  header_ptr,
               seq_type                               && seq,
               qual_type                              && qual,
               id_type                                && id,
               offset_type                            && offset,
               ref_seq_type                           && SEQAN3_DOXYGEN_ONLY(ref_seq),
               ref_id_type                            && ref_id,
               ref_offset_type                        && ref_offset,
               align_type                             && align,
               flag_type                              && flag,
               mapq_type                              && mapq,
               mate_type                              && mate,
               tag_dict_type                          && tag_dict,
               e_value_type                           && SEQAN3_DOXYGEN_ONLY(e_value),
               bit_score_type                         && SEQAN3_DOXYGEN_ONLY(bit_score))
    {
        /* Note the following general things:
         *
         * - Given the SAM specifications, all fields may be empty
         *
         * - Arithmetic values default to 0 while all others default to '*'
         *
         * - Because of the former, arithmetic values can be directly streamed
         *   into 'stream' as operator<< is defined for all arithmetic types
         *   and the default value (0) is also the SAM default.
         *
         * - All other non-arithmetic values need to be checked for emptiness
         */

        // ---------------------------------------------------------------------
        // Type Requirements (as static asserts for user friendliness)
        // ---------------------------------------------------------------------
        static_assert((std::ranges::ForwardRange<seq_type>        &&
                      Alphabet<value_type_t<remove_cvref_t<seq_type>>>),
                      "The seq object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        static_assert((std::ranges::ForwardRange<id_type>         &&
                      Alphabet<value_type_t<remove_cvref_t<id_type>>>),
                      "The id object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        static_assert(std::UnsignedIntegral<remove_cvref_t<offset_type>>,
                      "The offset object must be a std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<ref_seq_type>    &&
                      Alphabet<value_type_t<remove_cvref_t<ref_seq_type>>>),
                      "The ref_seq object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert((std::ranges::ForwardRange<ref_id_type>     &&
                      Alphabet<value_type_t<remove_cvref_t<ref_id_type>>>),
                      "The ref_id object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert(std::Integral<remove_cvref_t<ref_offset_type>>, // -1 is given default to evaluate to 0
                      "The ref_offset object must be an std::Integral >= 0.");

        if (((ref_offset + 1) < 0))
            throw format_error("The ref_offset object must be an std::Integral >= 0.");

        static_assert(tuple_like_concept<remove_cvref_t<align_type>>,
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert((std::tuple_size_v<remove_cvref_t<align_type>> == 2 &&
                       std::EqualityComparableWith<gap, value_type_t<remove_cvref_t<decltype(std::get<0>(align))>>> &&
                       std::EqualityComparableWith<gap, value_type_t<remove_cvref_t<decltype(std::get<1>(align))>>>),
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert(std::UnsignedIntegral<remove_cvref_t<flag_type>>,
                      "The flag object must be a std::UnsignedIntegral.");

        static_assert(std::UnsignedIntegral<remove_cvref_t<mapq_type>>,
                      "The mapq object must be a std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<qual_type>       &&
                       Alphabet<value_type_t<remove_cvref_t<qual_type>>>),
                      "The qual object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert(tuple_like_concept<remove_cvref_t<mate_type>>,
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) an std::UnsignedIntegral, and"
                      "3) an std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<decltype(std::get<0>(mate))>                     &&
                      Alphabet<value_type_t<remove_cvref_t<decltype(std::get<0>(mate))>>> &&
                      std::UnsignedIntegral<remove_cvref_t<decltype(std::get<1>(mate))>>          &&
                      std::UnsignedIntegral<remove_cvref_t<decltype(std::get<2>(mate))>>),
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) an std::UnsignedIntegral, and"
                      "3) an std::UnsignedIntegral.");

        static_assert(std::Same<remove_cvref_t<tag_dict_type>, sam_tag_dictionary>,
                      "The tag_dict object must be of type seqan3::sam_tag_dictionary.");

        // ---------------------------------------------------------------------
        // logical Requirements
        // ---------------------------------------------------------------------
        if (!empty(get<1>(align)) && empty(seq))
            throw format_error("If you specify an align object you must also specify the seq object. "
                               "Hint: Check if offset needs to be set to if soft-clipping is present.");

        if (options.sam_require_header && (header_ptr != nullptr) && !empty(ref_id))
        {
            if ((header_ptr->ref_dict).count(std::string(ref_id)) == 0) // no reference id matched
                throw format_error(std::string("The ref_id '") + std::string(ref_id) +
                                   "' was not in the list of references");
        }

        // ---------------------------------------------------------------------
        // Writing the Header on first call
        // ---------------------------------------------------------------------
        if (options.sam_require_header && !written_header && (header_ptr != nullptr))
        {
            write_header(stream, options, header_ptr);
            written_header = true;
        }

        // ---------------------------------------------------------------------
        // Writing the Record
        // ---------------------------------------------------------------------
        std::ranges::ostreambuf_iterator stream_it{stream};
        char const separator{'\t'};

        write_range(stream_it, std::forward<id_type>(id));

        stream << separator;

        stream << std::forward<flag_type>(flag) << separator;

        write_range(stream_it, std::forward<ref_id_type>(ref_id));

        stream << separator;

        stream << (ref_offset + 1) << separator; // SAM is 1 based

        stream << std::forward<mapq_type>(mapq) << separator;

        if (!empty(get<1>(align)))
        {
            // compute possible distance from alignment end to sequence end
            // which indicates soft clipping at the end.
            // This should be replace by a free count_gaps function for
            // aligned sequences which is more efficient if possible.
            size_t off_end{seq.size() - offset};
            for (auto chr : get<1>(align))
                if (chr == gap{})
                    ++off_end;
            off_end -= (get<1>(align)).size();

            write_range(stream_it,
                        detail::get_cigar_string(std::forward<align_type>(align),
                                                 std::forward<offset_type>(offset),
                                                 off_end));
        }
        else
        {
            stream << '*';
        }

        stream << separator;

        write_range(stream_it, get<0>(std::forward<mate_type>(mate)));

        stream << separator;

        stream << get<1>(std::forward<mate_type>(mate)) << separator;

        stream << get<2>(std::forward<mate_type>(mate)) << separator;

        write_range(stream_it, std::forward<seq_type>(seq));

        stream << separator;

        write_range(stream_it, std::forward<qual_type>(qual));

        write_tag_fields(stream, std::forward<tag_dict_type>(tag_dict), separator);

        detail::write_eol(stream_it, options.add_carriage_return);
    }

protected:
    //!\privatesection
    //!\brief The format version string.
    static constexpr char format_version[4] = "1.6";

    char buffer[100] = "";

    //!\brief A variable that tracks whether the header_ptr as been written or not.
    bool written_header{false};

    //!\brief A variable that tracks whether the header_ptr as been read or not.
    bool header_was_read{false};

    /*!\brief Decays to detail::consume for std::ignore.
     * \tparam stream_view_type  The type of the stream as a view.
     *
     * \param[in, out] stream_view  The stream view to consume.
     * \param[in, out] target       A std::ignore placeholder.
     */
    template <typename stream_view_type, typename target_type>
        requires detail::decays_to_ignore_v<target_type>
    void read_field(stream_view_type && stream_view, target_type /*target*/)
    {
        detail::consume(stream_view);
    }

    /*!\brief Reads a range by copying from stream_view to target, converting values with seqan3::view::char_to.
     * \tparam stream_view_type  The type of the stream as a view.
     * \tparam target_range_type The type of range to parse from input; Must model std::ranges::ForwardRange.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The range to store the parsed sequence.
     */
    template <typename stream_view_type, std::ranges::ForwardRange target_range_type>
        requires requires (std::string s) { { view::char_to<value_type_t<target_range_type>>(s) }; }
    void read_field(stream_view_type && stream_view, target_range_type & target)
    {
        ranges::copy(stream_view | view::char_to<value_type_t<target_range_type>>,
                     ranges::back_inserter(target));
    }

    /*!\brief Reads arithmetic fields using std::from_chars.
     * \tparam stream_view_type The type of the stream as a view.
     * \tparam target_type      The type of value to parse from input; Must model seqan3::arithmetic_concept.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The arithmetic value object to store the parsed value.
     *
     * \throws seqan3::parse_error if the character sequence in stream_view cannot be successfully converted to a value
     *         of type target_type.
     */
    template <typename stream_view_type, seqan3::arithmetic_concept target_type>
    void read_field(stream_view_type && stream_view, target_type & target)
    {
        // unfortunately std::from_chars only accepts char const * so we need a buffer.
        size_t idx{0};
        for (auto chr : stream_view)
        {
            buffer[idx] = chr;
            idx++;
        }

        std::from_chars_result res = std::from_chars(&buffer[0], &buffer[0] + idx, target);

        if (res.ec == std::errc::invalid_argument)
            throw parse_error(std::string("[SAM file in] The string '") + std::string(buffer) +
                                          "' could not be casted into an arithmetic value.");

        if (res.ec == std::errc::result_out_of_range)
            throw parse_error(std::string("[SAM file in] Casting '") + std::string(buffer) +
                                          "' into an arithmetic value would cause an overflow.");
    }

    /*!\brief Reads a list of values separated by comma as is the case for SAM tag arrays.
     * \tparam variant_type     A specialization of std::variant that is stored in the seqan3::sam_tag_dictionary.
     * \tparam stream_view_type The type of the stream as a view.
     * \tparam value_type       The type of values to be stored in the tag array.
     *
     * \param[in, out] variant      A std::variant object to store the tag arrays in.
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in]      value        A temporary value that determines the underlying type of the tag array.
     *
     * \details
     *
     * Reading the tags is done according to the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     *
     * The function throws a seqan3::parse_error if any unknown tag type was encountered. It will also fail if the
     * format is not in a correct state (e.g. required fields are not giving), but throwing might occur downstream of
     * the actual error.
     */
    template <typename variant_type, typename stream_view_type, typename value_type>
    void read_sam_dict_vector(variant_type & variant, stream_view_type && stream_view, value_type value)
    {
        std::vector<value_type> tmp_vector;
        while (begin(stream_view) != ranges::end(stream_view)) // not fully consumed yet
        {
            read_field(stream_view | view::take_until(is_char<','>), value);
            tmp_vector.push_back(value);

            if (is_char<','>(*begin(stream_view)))
                ++begin(stream_view); // skip ','
        }
        variant = tmp_vector;
    }

    /*!\brief Reads the optional tag fields into the seqan3::sam_tag_dictionry.
     * \tparam stream_view_type   The type of the stream as a view.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The seqan3::sam_tag_dictionary to store the tag information.
     *
     * \throws seqan3::parse_error if any unexpected character or format is encountered.
     *
     * \details
     *
     * Reading the tags is done according to the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     *
     * The function throws a seqan3::parse_error if any unknown tag type was encountered. It will also fail if the
     * format is not in a correct state (e.g. required fields are not giving), but throwing might occur downstream of
     * the actual error.
     */
    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, sam_tag_dictionary & target)
    {
        uint16_t tag = static_cast<uint16_t>(*begin(stream_view)) * 256;
        ++begin(stream_view); // skip char read before
        tag += static_cast<uint16_t>(*begin(stream_view));
        ++begin(stream_view); // skip char read before
        ++begin(stream_view); // skip ':'
        char type_char = *begin(stream_view);
        ++begin(stream_view); // skip char read before
        ++begin(stream_view); // skip ':'

        switch (type_char)
        {
            case 'A' : // char
            {
                target[tag] = static_cast<char>(*begin(stream_view));
                ++begin(stream_view); // skip char that has been read
                break;
            }
            case 'i' : // int32_t
            {
                int32_t tmp;
                read_field(stream_view, tmp);
                target[tag] = tmp;
                break;
            }
            case 'f' : // float
            {
                float tmp;
                read_field(stream_view, tmp);
                target[tag] = tmp;
                break;
            }
            case 'Z' : // string
            {
                target[tag] = std::string(stream_view);
                break;
            }
            case 'H' :
            {
                // TODO
                break;
            }
            case 'B' : // Array. Value type depends on second char [cCsSiIf]
            {
                char value_type_char = *begin(stream_view);
                ++begin(stream_view); // skip char read before
                ++begin(stream_view); // skip first ','

                switch (value_type_char)
                {
                    case 'c' : // int8_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int8_t{});
                        break;
                    }
                    case 'C' : // uint8_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint8_t{});
                        break;
                    }
                    case 's' : // int16_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int16_t{});
                        break;
                    }
                    case 'S' : // uint16_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint16_t{});
                        break;
                    }
                    case 'i' : // int32_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int32_t{});
                        break;
                    }
                    case 'I' : // uint32_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint32_t{});
                        break;
                    }
                    case 'f' : // float
                    {
                        read_sam_dict_vector(target[tag], stream_view, float{});
                        break;
                    }
                    default:
                        throw parse_error(std::string("[SAM Parse Error] The first character in the numerical ") +
                                          "id of a SAM tag must be one of [cCsSiIf] but " + value_type_char +
                                          " was given.");
                }
                break;
            }
            default:
                throw parse_error(std::string("[SAM Parse Error] The second character in the numerical id of a "
                                  "SAM tag must be one of [A,i,Z,H,B,f] but ") + type_char + "was given.");
        }
    }

    /*!\brief Reads the SAM header_ptr.
     * \tparam stream_view_type     The type of the stream as a view.
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] header_ptr   The header_ptr (as a pointer) to store the parsed values.
     *
     * \throws seqan3::parse_error if any unexpected character or format is encountered.
     *
     * \details
     *
     * Reading the header format is done according to the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     *
     * The function throws a seqan3::parse_error if any unknown tag was encountered. It will also fail if the format is
     * not in a correct state (e.g. required fields are not giving), but throwing might occur downstream of the actual
     * error.
     */
    template <typename stream_view_type>
    void read_header(stream_view_type && stream_view, std::unique_ptr<alignment_file_header> & header_ptr)
    {
        alignment_file_header & hdr{*header_ptr};

        auto parse_tag_value = [&stream_view, this] (auto & value) // helper function to parse the next tag value
            {
                detail::consume(stream_view | view::take_until_or_throw(is_char<':'>)); // skip tag name
                ++begin(stream_view);                                                   // skip ':'
                read_field(stream_view | view::take_until_or_throw(is_char<'\t'> || is_char<'\n'>), value);
                ++begin(stream_view);                                                   // skip tab or newline
            };

        // @HQ line
        // -------------------------------------------------------------------------------------------------------------
        parse_tag_value(hdr.format_version); // parse required VN (version) tag

        // The SO, SS and GO tag are optional and can appear in any order
        while (!is_char<'@'>(*begin(stream_view)))
        {
            std::string * who = &hdr.grouping;

            if (is_char<'S'>(*begin(stream_view)))
            {
                ++begin(stream_view); // skip S

                if (is_char<'O'>(*begin(stream_view)))          // SO (sorting) tag
                    who = &hdr.sorting;
                else if (is_char<'S'>(*begin(stream_view)))     // SS (sub-order) tag
                    who = &hdr.subsorting;
                else
                    throw parse_error(std::string("Illegal SAM header tag: S" + (*begin(stream_view))));
            }
            else if (!is_char<'G'>(*begin(stream_view)))        // GO (grouping) tag
            {
                throw parse_error(std::string("Illegal SAM header tag in @HG starting with:" + (*begin(stream_view))));
            }

            parse_tag_value(*who);
        }

        // The rest of the header lines
        // -------------------------------------------------------------------------------------------------------------
        while (is_char<'@'>(*begin(stream_view)))
        {
            ++begin(stream_view); // skip @

            if (is_char<'S'>(*begin(stream_view)))              // SQ (sequence dictionary) tag
            {
                std::string id;
                std::tuple<uint32_t, std::string> tmp{};

                parse_tag_value(id);                            // parse required SN (sequence name) tag
                parse_tag_value(get<0>(tmp));                   // parse required LN (length) tag

                if (!is_char<'@'>(*begin(stream_view)))         // read rest of the tags
                {
                    read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), get<1>(tmp));
                    ++begin(stream_view);
                }

                hdr.ref_dict[id] = tmp;
            }
            else if (is_char<'R'>(*begin(stream_view)))         // RG (read group) tag
            {
                std::pair<std::string, std::string> tmp{};

                parse_tag_value(get<0>(tmp));                   // read required ID tag

                if (!is_char<'@'>(*begin(stream_view)))         // read rest of the tags
                {
                    read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), get<1>(tmp));
                    ++begin(stream_view);
                }

                hdr.read_groups.emplace_back(std::move(tmp));
            }
            else if (is_char<'P'>(*begin(stream_view)))         // PG (program) tag
            {
                alignment_file_header::program_info_t tmp{};

                parse_tag_value(tmp.id);                        // read required ID tag

                // The PN, CL, PP, DS, VN are optional tags and can be given in any order.
                while (!is_char<'@'>(*begin(stream_view)))
                {
                    std::string * who = &tmp.version;

                    if (is_char<'P'>(*begin(stream_view)))
                    {
                        ++begin(stream_view); // skip P

                        if (is_char<'N'>(*begin(stream_view)))  // PN (program name) tag
                            who = &tmp.name;
                        else                                    // PP (previous program) tag
                            who = &tmp.previous;
                    }
                    else if (is_char<'C'>(*begin(stream_view))) // CL (command line) tag
                    {
                        who = &tmp.command_line_call;
                    }
                    else if (is_char<'D'>(*begin(stream_view))) // DS (description) tag
                    {
                        who = &tmp.description;
                    }
                    else if (!is_char<'V'>(*begin(stream_view))) // VN (version) tag
                    {
                        throw parse_error(std::string("Illegal SAM header tag starting with:" + (*begin(stream_view))));
                    }

                    parse_tag_value(*who);
                }

                hdr.program_infos.emplace_back(std::move(tmp));
            }
            else if (is_char<'C'>(*begin(stream_view)))         // CO (comment) tag
            {
                std::string tmp;
                ++begin(stream_view); // skip C
                ++begin(stream_view); // skip O
                ++begin(stream_view); // skip :
                read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), tmp);
                ++begin(stream_view); // skip newline

                hdr.comments.emplace_back(std::move(tmp));
            }
            else
            {
                throw parse_error(std::string("Illegal SAM header tag starting with:" + (*begin(stream_view))));
            }
        }
    }

    /*!\brief Writes a field value to the stream.
     * \tparam stream_it_t The stream iterator type.
     * \tparam field_type  The type of the field value. Must model std::ranges::ForwardRange.
     *
     * \param[in,out] stream_it   The stream iterator to print to.
     * \param[in]     field_value The value to print.
     */
    template <typename stream_it_t, typename field_type>
    //!\cond
        requires std::ranges::ForwardRange<field_type>
    //!\endcond
    void write_range(stream_it_t & stream_it, field_type && field_value)
    {
        if (empty(field_value))
            stream_it = '*';
        else
            std::ranges::copy(field_value | view::to_char | view::take_until(is_space), stream_it);
    }

    /*!\brief Writes the optional fields of the seqan3::sam_tag_dictionary.
     * \tparam stream_t   The stream type.
     *
     * \param[in,out] stream    The stream to print to.
     * \param[in]     tag_dict  The tag dictionary to print.
     * \param[in]     separator The field separator to append.
     */
    template <typename stream_t>
    void write_tag_fields(stream_t & stream, sam_tag_dictionary const & tag_dict, char const separator)
    {
        auto stream_variant_fn = [&stream] (auto && arg) // helper to print an std::variant
        {
            using T = remove_cvref_t<decltype(arg)>;

            if constexpr (!container_concept<T>)
            {
                stream << arg;
            }
            else
            {
                if (arg.begin() != arg.end())
                {
                    for (auto it = arg.begin(); it != (arg.end() - 1); ++it)
                        stream << *it << ",";

                    stream << *(arg.end() - 1); // write last value without trailing ','
                }
            }
        };

        for (auto & [tag, variant] : tag_dict)
        {
            stream << separator;

            char char0 = tag / 256;
            char char1 = tag % 256;

            stream << char0 << char1 << ':' << detail::sam_tag_type_char[variant.index()] << ':';

            if (detail::sam_tag_type_char_extra[variant.index()] != '\0')
                stream << detail::sam_tag_type_char_extra[variant.index()] << ',';

            std::visit(stream_variant_fn, variant);
        }
    }

    /*!\brief Writes the SAM header_ptr.
     * \tparam stream_t   The stream type.
     *
     * \param[in,out] stream  The stream to print to.
     * \param[in]     options The options to alter printing.
     * \param[in]     header_ptr  The header_ptr (as a pointer) to print.
     *
     * \throws seqan3::format_error if the header_ptr object contains the wrong
     *         information or the contents are ill-formed.
     *
     * \details
     *
     * Before writing the header_ptr, the contents are checked for correctness
     * according to the rules of the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     */
    template <typename stream_t>
    void write_header(stream_t & stream,
                      alignment_file_output_options const & options,
                      std::unique_ptr<alignment_file_header> & header_ptr)
    {
        if (header_ptr != nullptr)
        {
            // -----------------------------------------------------------------
            // Check Header
            // -----------------------------------------------------------------

            // (@HD) Check header_ptr line
            // The format version string will be taken from the local member variable
            if (!header_ptr->sorting.empty() &&
                !(header_ptr->sorting == "unknown"   ||
                  header_ptr->sorting == "unsorted"  ||
                  header_ptr->sorting == "queryname" ||
                  header_ptr->sorting == "coordinate" ))
                throw format_error{"SAM format error: The header_ptr->sorting member must be "
                                   "one of [unknown, unsorted, queryname, coordinate]."};

            if (!header_ptr->grouping.empty() &&
                !(header_ptr->grouping == "none"   ||
                  header_ptr->grouping == "query"  ||
                  header_ptr->grouping == "reference"))
                throw format_error{"SAM format error: The header_ptr->grouping member must be "
                                   "one of [none, query, reference]."};

            // (@SQ) Check Reference Sequence Dictionary lines

            // TODO

            // - sorting order be one of ...
            // - grouping can be one of ...
            // - reference names must be unique
            // - ids of read groups must be unique
            // - program ids need to be unique
            // many more small semantic things, like fits REGEX

            // -----------------------------------------------------------------
            // Write Header
            // -----------------------------------------------------------------
            std::ranges::ostreambuf_iterator stream_it{stream};

            // (@HD) Write header_ptr line [required].
            stream << "@HD\tVN:";
            stream << format_version;

            if (!header_ptr->sorting.empty())
                stream << "\tSO:" << header_ptr->sorting;

            if (!header_ptr->subsorting.empty())
                stream << "\tSS:" << header_ptr->subsorting;

            if (!header_ptr->grouping.empty())
                stream << "\tGO:" << header_ptr->grouping;

            detail::write_eol(stream_it, options.add_carriage_return);

            // (@SQ) Write Reference Sequence Dictionary lines [required].
            for (auto const & [ref_name, ref_info] : header_ptr->ref_dict)
            {
                stream << "@SQ"
                       << "\tSN:" << ref_name
                       << "\tLN:" << get<0>(ref_info);

                if (!get<1>(ref_info).empty())
                    stream << "\t" << get<1>(ref_info);

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write read group (@RG) lines if specified.
            for (auto const & read_group : header_ptr->read_groups)
            {
                stream << "@RG"
                       << "\tID:" << get<0>(read_group);

                if (!get<1>(read_group).empty())
                    stream << "\t" << get<1>(read_group);

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write program (@PG) lines if specified.
            for (auto const & program : header_ptr->program_infos)
            {
                stream << "@PG"
                       << "\tID:" << program.id;

                if (!program.name.empty())
                    stream << "\tPN:" << program.name;

                if (!program.command_line_call.empty())
                    stream << "\tCL:" << program.command_line_call;

                if (!program.previous.empty())
                    stream << "\tPP:" << program.previous;

                if (!program.description.empty())
                    stream << "\tDS:" << program.description;

                if (!program.version.empty())
                    stream << "\tVN:" << program.version;

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write comment (@CO) lines if specified.
            for (auto const & comment : header_ptr->comments)
            {
                stream << "@CO\t" << comment;
                detail::write_eol(stream_it, options.add_carriage_return);
            }
        }
    }
};

} // namespace seqan3
