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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::take_until and seqan3::view::take_until_or_throw.
 */

#pragma once

#include <range/v3/view/take_while.hpp>
#include <range/v3/algorithm/copy.hpp>

#include <seqan3/core/metafunction/iterator.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/type_traits>
#include <seqan3/std/view/view_all.hpp>

namespace seqan3::detail
{

// ============================================================================
//  view_take_until
// ============================================================================

/*!\brief The type returned by seqan3::view::take_until and seqan3::view::take_until_or_throw.
 * \tparam urng_t    The type of the underlying range, must model std::ranges::View.
 * \tparam fun_t     Type of the callable that will be evaluated on every member; must model
 *                   std::Invocable with seqan3::reference_t<urng_t> as argument and return `bool`.
 * \tparam or_throw  Whether to throw an exception when the input is exhausted before the end of line is reached.
 * \implements std::ranges::View
 * \implements std::ranges::RandomAccessRange
 * \ingroup view
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::View urng_t, typename fun_t, bool or_throw>
class view_take_until : public ranges::view_interface<view_take_until<urng_t, fun_t, or_throw>>
{
private:

    static_assert(std::Invocable<fun_t, reference_t<urng_t>>,
                  "The functor type for view::take_until must model std::Invocable<fun_t, reference_t<urng_t>>.");
    static_assert(std::Boolean<std::result_of_t<fun_t&&(reference_t<urng_t>)>>,
                  "The functor type for view::take_until must return std::Boolean.");

    //!\brief The underlying range.
    urng_t urange;

    //!\brief The functor.
    ranges::semiregular_t<fun_t> fun;

    //!\brief Whether this view is const_iterable or not.
    static constexpr bool const_iterable = const_iterable_concept<urng_t> &&
                                           std::RegularInvocable<fun_t, reference_t<urng_t>>;

    //!\brief The sentinel type is identical to that of the underlying range.
    using sentinel_type = std::ranges::sentinel_t<urng_t>;

    //!\brief The iterator type inherits from the underlying type, but overwrites several operators.
    //!\tparam rng_t Should be `urng_t` for defining #iterator and `urng_t const` for defining #const_iterator.
    template <typename rng_t>
    class iterator_type : public inherited_iterator_base<iterator_type<rng_t>, std::ranges::iterator_t<rng_t>>
    {
    private:
        //!\brief The iterator type of the underlying range.
        using base_base_t = std::ranges::iterator_t<rng_t>;
        //!\brief The CRTP wrapper type.
        using base_t      = inherited_iterator_base<iterator_type, std::ranges::iterator_t<rng_t>>;

        //!\brief Auxiliary type.
        using fun_ref_t = std::conditional_t<std::is_const_v<rng_t>,
                                             std::remove_reference_t<fun_t> const &,
                                             std::remove_reference_t<fun_t> &>;
        //!\brief Reference to the functor stored in the view.
        ranges::semiregular_t<fun_ref_t> fun;

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        iterator_type() = default;
        constexpr iterator_type(iterator_type const & rhs) = default;
        constexpr iterator_type(iterator_type && rhs) = default;
        constexpr iterator_type & operator=(iterator_type const & rhs) = default;
        constexpr iterator_type & operator=(iterator_type && rhs) = default;
        ~iterator_type() = default;

        //!\brief Constructor that delegates to the CRTP layer.
        iterator_type(base_base_t const & it) :
            base_t{it}
        {}

        //!\brief Constructor that delegates to the CRTP layer and initialises the callable.
        iterator_type(base_base_t it, fun_ref_t _fun) :
            base_t{it}, fun{_fun} {}
        //!\}

        /*!\name Associated types
         * \brief All are derived from the base_base_t.
         * \{
         */
        using difference_type       = typename std::iterator_traits<base_base_t>::difference_type;
        using value_type            = typename std::iterator_traits<base_base_t>::value_type;
        using reference             = typename std::iterator_traits<base_base_t>::reference;
        using pointer               = typename std::iterator_traits<base_base_t>::pointer;
        using iterator_category     = typename std::iterator_traits<base_base_t>::iterator_category;
        //!\}

        /*!\name Comparison operators
         * \brief We define comparison against self and against the sentinel.
         * \{
         */
        bool operator==(iterator_type const & rhs) const noexcept(!or_throw)
        {
            return static_cast<base_base_t>(*this) == static_cast<base_base_t>(rhs);
        }

        bool operator==(sentinel_type const & rhs) const noexcept(!or_throw)
        {
            if (static_cast<base_base_t>(*this) == rhs)
            {
                if constexpr (or_throw)
                    throw unexpected_end_of_input{"Reached end of input before functor evaluated to true."};
                else
                    return true;
            }

            return fun(**this);
        }

        friend bool operator==(sentinel_type const & lhs, iterator_type const & rhs) noexcept(!or_throw)
        {
            return rhs == lhs;
        }

        bool operator!=(sentinel_type const & rhs) const noexcept(!or_throw)
        {
            return !(*this == rhs);
        }

        bool operator!=(iterator_type const & rhs) const noexcept(!or_throw)
        {
            return static_cast<base_base_t>(*this) != static_cast<base_base_t>(rhs);
        }

        friend bool operator!=(sentinel_type const & lhs, iterator_type const & rhs) noexcept(!or_throw)
        {
            return rhs != lhs;
        }
        //!\}
    }; // class iterator_type

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The reference_type.
    using reference         = reference_t<urng_t>;
    //!\brief The const_reference type is equal to the reference type if the underlying range is const-iterable.
    using const_reference   = detail::transformation_trait_or_t<seqan3::reference<urng_t const>, void>;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = value_type_t<urng_t>;
    //!\brief The size_type is void, because this range is never sized.
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = difference_type_t<urng_t>;
    //!\brief The iterator type of this view (a random access iterator).
    using iterator          = iterator_type<urng_t>;
    //!\brief The const_iterator type is equal to the iterator type if the underlying range is const-iterable.
    using const_iterator    = detail::transformation_trait_or_t<std::type_identity<iterator_type<urng_t const>>, void>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_take_until() = default;
    constexpr view_take_until(view_take_until const & rhs) = default;
    constexpr view_take_until(view_take_until && rhs) = default;
    constexpr view_take_until & operator=(view_take_until const & rhs) = default;
    constexpr view_take_until & operator=(view_take_until && rhs) = default;
    ~view_take_until() = default;

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range.
     * \param[in] _fun    The functor that acts as termination criterium.
     */
    view_take_until(urng_t _urange, fun_t _fun)
        : urange{std::move(_urange)}, fun{std::forward<fun_t>(_fun)}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to seqan3::view_take_until::end().
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
        return {seqan3::begin(urange), static_cast<fun_t &>(fun)};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
        requires const_iterable
    {
        return {seqan3::cbegin(urange), static_cast<fun_t const &>(fun)};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
        requires const_iterable
    {
        return {seqan3::cbegin(urange), static_cast<fun_t const &>(fun)};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
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
    sentinel_type end() noexcept
    {
        return {seqan3::end(urange)};
    }

    //!\copydoc end()
    sentinel_type end() const noexcept
        requires const_iterable
    {
        return {seqan3::cend(urange)};
    }

    //!\copydoc end()
    sentinel_type cend() const noexcept
        requires const_iterable
    {
        return {seqan3::cend(urange)};
    }
    //!\}

    /*!\brief Convert this view into a container implicitly.
     * \tparam container_t Type of the container to convert to; must model seqan3::sequence_container_concept and it's
     *                     seqan3::reference_t must model std::CommonReference with `reference`.
     * \returns This view converted to container_t.
     */
    template <sequence_container_concept container_t>
    operator container_t()
    //!\cond
        requires std::CommonReference<reference_t<container_t>, reference>
    //!\endcond
    {
        container_t ret;
        std::ranges::copy(begin(), end(), std::back_inserter(ret));
        return ret;
    }

    //!\overload
    template <sequence_container_concept container_t>
    operator container_t() const
    //!\cond
        requires std::CommonReference<reference_t<container_t>, reference> && const_iterable_concept<urng_t>
    //!\endcond
    {
        container_t ret;
        std::ranges::copy(cbegin(), cend(), std::back_inserter(ret));
        return ret;
    }
};

//!\brief Type deduction guide that strips references.
//!\relates seqan3::detail::view_take_until
template <typename urng_t, typename fun_t, bool or_throw = false>
view_take_until(urng_t, fun_t) -> view_take_until<std::remove_reference_t<urng_t>, fun_t, or_throw>;

// ============================================================================
//  take_until_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for view::take_until and view::take_until_or_throw.
 * \tparam or_throw Whether to throw an exception when the input is exhausted before the end of line is reached.
 */
template <bool or_throw>
class take_until_fn : public pipable_adaptor_base<take_until_fn<or_throw>>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = pipable_adaptor_base<take_until_fn<or_throw>>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief Call the view's constructor with the given parameters.
     * \tparam    urng_t The underlying range type; must model std::ranges::View.
     * \tparam    fun_t  The type of the callable; concept checks done in class.
     * \param[in] urange The underlying range.
     * \param[in] fun    The callable that will be evaluated on every element.
     * \returns An instance of seqan3::detail::view_take_until.
     */
    template <std::ranges::View urng_t, typename fun_t>
    static auto impl(urng_t urange, fun_t && fun)
    {
        return view_take_until<urng_t, fun_t, or_throw>{std::move(urange), std::forward<fun_t>(fun)};
    }

    /*!\brief Wraps the range argument in seqan3::view::all and forwards to the other overload.
     * \tparam    urng_t The underlying range type; must model std::ranges::ViewableRange (be lvalue reference
     *                   to non-view range).
     * \tparam    fun_t  The type of the callable; concept checks done in class.
     * \param[in] urange The underlying range.
     * \param[in] fun    The callable that will be evaluated on every element.
     * \returns An instance of seqan3::detail::view_take_until.
     */
    template <std::ranges::ViewableRange urng_t, typename fun_t>
    static auto impl(urng_t && urange, fun_t && fun)
    {
        return impl(view::all(std::forward<urng_t>(urange)), std::forward<fun_t>(fun));
    }
};

} // namespace seqan3::detail

// ============================================================================
//  view::take_until (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns elements from the underlying range until the functor evaluates to
 *                      true (or the end of the underlying range is reached).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam fun_t        The type of the functor; must model std::Invocable with seqan3::reference_t<urng_t>
 *                      and return a type convertible to `bool`.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] fun       The functor.
 * \returns             All elements of the underlying range up until (but excluding) the element that evaluates the
 *                      functor to true.
 * \ingroup view
 *
 * \details
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | depends on functor (see below)                     |
 * | std::ranges::BidirectionalRange |                                       | depends on functor (see below)                     |
 * | std::ranges::RandomAccessRange  |                                       | depends on functor (see below)                     |
 * | std::ranges::ContiguousRange    |                                       | depends on functor (see below)                     |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *lost*                                             |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | depends on functor (see below)                     |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * This view only *preserves* certain concepts if the specified functor models
 * seqan3::regular_std::Invocable<fun_t, reference_t<urng_t>, i.e.
 * applying the functor doesn't change the functor. If the functor only models std::Invocable and not
 * seqan3::regular_invocable_concept these concepts are *lost*.
 *
 * ### Example
 *
 * \snippet test/snippet/range/view/take_until.cpp usage
 *
 * \hideinitializer
 */
inline auto constexpr take_until = detail::take_until_fn<false>{};

// ============================================================================
//  view::take_until_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns elements from the underlying range until the functor evaluates to true
 *        (**throws** if the end of the underlying range is reached).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no element that satisfies the functor.
 * \ingroup view
 *
 * \copydetails seqan3::view::take_until
 * \hideinitializer
 */
inline auto constexpr take_until_or_throw = detail::take_until_fn<true>{};

//!\}

} // namespace seqan3::view
