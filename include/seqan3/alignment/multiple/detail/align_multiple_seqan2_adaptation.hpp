// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides some detail functionality to the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \author Simon Sasse <simon.sasse AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/modifier.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/test_accessor.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/chunk.hpp>

#include <seqan/basic.h>
#include <seqan/reduced_aminoacid.h>

namespace seqan3::detail
{

template <typename seqan3_alphabet_type>
class align_multiple_seqan2_adaptation
{
private:
    // same order
    using seqan3_types = type_list<dna4, dna5, dna15, rna4, rna5, aa27, aa10murphy, aa10li>;
    using seqan2_types = type_list<seqan::Dna,
                                   seqan::Dna5,
                                   seqan::Iupac,
                                   seqan::Rna,
                                   seqan::Rna5,
                                   seqan::AminoAcid,
                                   seqan::ReducedAminoAcid<seqan::Murphy10>,
                                   seqan::ReducedAminoAcid<seqan::Li10>>;
    static constexpr auto index = list_traits::find<seqan3_alphabet_type, seqan3_types>;

public:
    using alphabet_type = list_traits::at<index, seqan2_types>;

    using sequence_type = seqan::String<alphabet_type>;

    using graph_type = seqan::Graph<seqan::Alignment<seqan::StringSet<sequence_type, seqan::Dependent<>>,
                                                     void,
                                                     seqan::WithoutEdgeId>>;

    template <typename seqan3_configuration_t>
    auto create_msa_configuration(seqan3_configuration_t const & config)
    {
        validate_configuration(config);

        auto msaOpt = initialise_scoring_scheme(config);

        seqan::appendValue(msaOpt.method, 0); // global alignment
        seqan::appendValue(msaOpt.method, 1); // local alignment
        msaOpt.build = 0; // neighbour joining to build the guide tree

        if constexpr (config.template exists<seqan3::align_cfg::band_fixed_size>())
        {
            msaOpt.pairwiseAlignmentMethod = 2; // banded
            auto const & band = get<seqan3::align_cfg::band_fixed_size>(config);
            msaOpt.bandWidth = band.upper_diagonal.get() - band.lower_diagonal.get();
        }
        else
        {
            msaOpt.pairwiseAlignmentMethod = 1; // unbanded
        }

        // seqan2 tcoffee app default: gap -1, gap open -13
        auto const & gaps = config.get_or(align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-13}}}).value;
        // convert to seqan2 gap score convention.
        // See: https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1gap__scheme.html
        msaOpt.sc.data_gap_open = gaps.get_gap_open_score() + gaps.get_gap_score();
        msaOpt.sc.data_gap_extend = gaps.get_gap_score();

        return msaOpt;
    }

    template <std::ranges::forward_range range_t>
    auto convert_sequences(std::vector<range_t> const & input)
    {
        seqan::StringSet<sequence_type, seqan::Owner<>> sequenceSet;
        seqan::StringSet<seqan::String<char>> sequenceNames;

        seqan::String<char> dummy_name = "dummy_name";
        for (auto const & seq : input)
        {
            sequence_type tmp;
            for (auto chr : seq)
                seqan::appendValue(tmp, alphabet_type{seqan3::to_char(chr)});

            seqan::appendValue(sequenceSet, tmp);
            seqan::appendValue(sequenceNames, dummy_name);
        }

        return std::make_pair(std::move(sequenceSet), std::move(sequenceNames));
    }

    template <typename graph_type>
    auto create_output(graph_type & gAlign)
    {
        std::string mat;
        seqan::convertAlignment(gAlign, mat);

        using gapped_alphabet_type = seqan3::gapped<seqan3_alphabet_type>;
        std::vector<std::vector<gapped_alphabet_type>> output{seqan::length(seqan::stringSet(gAlign))};
        size_t ali_len = mat.size() / output.size();

        auto iter = output.begin();
        for (auto && gapped_seq : mat | seqan3::views::char_to<gapped_alphabet_type> | seqan3::views::chunk(ali_len))
        {
            iter->reserve(ali_len);
            std::ranges::copy(gapped_seq, std::cpp20::back_inserter(*iter++));
        }

        return output;
    }

private:
    /*!\brief Predicate that represents whether `candidate_t` is a type that is not allowed as MSA configuration.
     * \tparam candidate_t The type to test.
     *
     * \details
     * The predicate's value is set to: NOT (band OR gap OR scoring).
     */
    template <typename candidate_t>
    struct is_invalid_msa_config :
        std::negation<std::disjunction<std::is_same<candidate_t, seqan3::align_cfg::band_fixed_size>,
                                       is_type_specialisation_of<candidate_t, seqan3::align_cfg::gap>,
                                       is_type_specialisation_of<candidate_t, seqan3::align_cfg::scoring>>>{};

    /*!\brief Validate the given MSA configuration.
     * \tparam configs_t The specified configuration elements.
     */
    template <typename ... configs_t>
    void validate_configuration(seqan3::configuration<configs_t...> const &)
    {
        static_assert(seqan3::pack_traits::find_if<is_invalid_msa_config, configs_t...> == -1,
                      "The given MSA configuration is not valid.");
    }

    template <typename seqan3_configuration_t>
        requires seqan3_configuration_t::template exists<seqan3::align_cfg::scoring>()
    auto initialise_scoring_scheme(seqan3_configuration_t const & config)
    {
        auto scoring_scheme = get<seqan3::align_cfg::scoring>(config).value;

        // A seqan3 scoring scheme might be of alphabet type dna15 although the aligned alphabet is dna4
        // seqan3 alphabet type used:
        using scoring_scheme_alphabet_type = typename decltype(scoring_scheme)::alphabet_type;
        constexpr auto scoring_scheme_alphabet_index = list_traits::find<scoring_scheme_alphabet_type, seqan3_types>;
        // seqan2 alphabet type to use accordingly:
        using score_matrix_alphabet_type = list_traits::at<scoring_scheme_alphabet_index, seqan2_types>;
        using score_matrix_type = seqan::Score<int, seqan::ScoreMatrix<score_matrix_alphabet_type>>;

        seqan::MsaOptions<alphabet_type, score_matrix_type> msaOpt{};
        score_matrix_type scMat;

        for (size_t i = 0; i < seqan3::alphabet_size<scoring_scheme_alphabet_type>; ++i)
        {
            for (size_t j = 0; j < seqan3::alphabet_size<scoring_scheme_alphabet_type>; ++j)
            {
                auto seqan3_i = seqan3::assign_rank_to(i, scoring_scheme_alphabet_type{});
                auto seqan3_j = seqan3::assign_rank_to(j, scoring_scheme_alphabet_type{});
                score_matrix_alphabet_type seqan2_i{seqan3::to_char(seqan3_i)};
                score_matrix_alphabet_type seqan2_j{seqan3::to_char(seqan3_j)};

                setScore(scMat, seqan2_i, seqan2_j, scoring_scheme.score(seqan3_i, seqan3_j));
            }
        }

        msaOpt.sc = scMat;
        return msaOpt;
    }

    template <typename seqan3_configuration_t>
        requires !seqan3_configuration_t::template exists<seqan3::align_cfg::scoring>()
    auto initialise_scoring_scheme(seqan3_configuration_t const &)
    {
        using score_type = std::conditional_t<std::same_as<alphabet_type, seqan::AminoAcid> ||
                                              std::same_as<alphabet_type, seqan::ReducedAminoAcid<seqan::Murphy10>> ||
                                              std::same_as<alphabet_type, seqan::ReducedAminoAcid<seqan::Li10>>,
                                              seqan::Blosum62,
                                              seqan::Score<int>>;

        seqan::MsaOptions<alphabet_type, score_type> msaOpt{};

        if constexpr (std::is_same_v<score_type, seqan::Score<int>>)
        {
            msaOpt.sc.data_match = 5;       // tcoffee app default
            msaOpt.sc.data_mismatch = -4;   // tcoffee app default
        }

        return msaOpt;
    }

    //!\brief Befriend seqan3::detail::test_accessor to grant access to private member functions.
    friend struct ::seqan3::detail::test_accessor;
};

} //seqan3::detail
