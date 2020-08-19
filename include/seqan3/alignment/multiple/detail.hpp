// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides some detail functionality to the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

#include <seqan/basic.h>
#include <seqan/reduced_aminoacid.h>

namespace seqan3::detail
{

// unfortunately it doesn't work with our alphabets..
// Nucleotide conversion------------------------------------------------------------------------------------------------
auto convert_alph_3_to_2(seqan3::dna4 chr)
{
    return seqan::Dna{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::dna5 chr)
{
    return seqan::Dna5{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::dna15 chr)
{
    return seqan::Iupac{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::rna4 chr)
{
    return seqan::Rna{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::rna5 chr)
{
    return seqan::Rna5{seqan3::to_char(chr)};
}
// Amino Acid conversion------------------------------------------------------------------------------------------------
auto convert_alph_3_to_2(seqan3::aa27 chr)
{
    return seqan::AminoAcid{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::aa10murphy chr)
{
    return seqan::ReducedAminoAcid<seqan::Murphy10>{seqan3::to_char(chr)};
}

auto convert_alph_3_to_2(seqan3::aa10li chr)
{
    return seqan::ReducedAminoAcid<seqan::Li10>{seqan3::to_char(chr)};
}

/*!
 * \brief Predicate that represents whether `candidate_t` is a type that is not allowed as MSA configuration.
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

/*!
 * \brief Validate the given MSA configuration.
 * \tparam configs_t The specified configuration elements.
 */
template <detail::config_element_specialisation ... configs_t>
void validate_configuration(seqan3::configuration<configs_t...> const &)
{
    static_assert(seqan3::pack_traits::find_if<is_invalid_msa_config, configs_t...> == -1,
                  "The given MSA configuration is not valid.");
}

template <typename alphabet_type, typename seqan3_configuration_t>
    requires seqan3_configuration_t::template exists<seqan3::align_cfg::scoring>()
auto initialise_scoring_scheme(seqan3_configuration_t const & config)
{
    auto scoring_scheme = get<seqan3::align_cfg::scoring>(config).value;
    using scoring_scheme_alphabet_type = typename decltype(scoring_scheme)::alphabet_type; // seqan3 alphabet
    using score_matrix_alphabet_type = decltype(convert_alph_3_to_2(scoring_scheme_alphabet_type{})); // seqan2 alphabet
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

template <typename alphabet_type, typename seqan3_configuration_t>
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

template <typename alphabet_type, typename seqan3_configuration_t>
auto seqan2_msa_configuration(seqan3_configuration_t const & config)
{
    auto msaOpt = initialise_scoring_scheme<alphabet_type>(config);

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

} //seqan3::detail