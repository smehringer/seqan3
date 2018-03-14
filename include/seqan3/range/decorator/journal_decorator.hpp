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
 * \brief This file contains the journal_decorator class.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/core.hpp>       // integral_concept
#include <seqan3/core/concept/iterator.hpp>   // input_iterator_concept
#include <seqan3/range/container/concept.hpp> // random_access_sequence_concept
#include <seqan3/range/concept.hpp>           // random_access_range_concept
#include <seqan3/range/decorator/journal_decorator_detail.hpp>

namespace seqan3
{

/*!\interface seqan3::journal_decorator_traits_concept <>
 * \brief Checks the required traits passed to the journal_decorator
 */
//!\cond
template <typename t>
concept bool journal_decorator_traits_concept = requires (t v)
{
    // Journal Node Container
    typename t::template journal_container_type<detail::journal_node<uint32_t, uint32_t>>;
    // the insertion buffer must have random access and modifiers
    requires random_access_sequence_concept<
        typename t::template journal_container_type<detail::journal_node<uint32_t, uint32_t>>>;

    // Insertion Buffer Container
    typename t::template insertion_buffer_type<uint32_t>;
    // the insertion buffer must have random access and modifiers
    requires random_access_sequence_concept<
        typename t::template insertion_buffer_type<uint32_t>>;
};
//!\endCond

/*\brief The default trait_types a custom traits_type can inherit from
 * \tparam position_type The type of a position in the host sequence or insertion buffer.
 * \tparam size_type The type of the length of segments of the host sequence or insertion buffer.
 * \tparam journal_container_type The container type holding the journal entries.
 *
 * Note: The journal_container_type is important for the runtime of many functions
 * of the jourmal_decorator and should be chosen carefully based on the
 * application.
 */
struct journal_decorator_default_traits
{
    template <typename journal_node_type>
    using journal_container_type = std::vector<journal_node_type>;

    template <typename value_type>
    using insertion_buffer_type = std::vector<value_type>;
};

template <random_access_range_concept urng_t,
          journal_decorator_traits_concept traits_type = journal_decorator_default_traits>
class journal_decorator
{
protected:
    //!\brief The container type storing novel inserted sequences.
    using insertion_buffer_type = typename traits_type::template insertion_buffer_type<typename urng_t::value_type>;

public:
    //!\brief Same as the value_type of the underlying range.
    using value_type = typename urng_t::value_type;
    //!\brief Same as the reference type of the underlying range.
    using reference = typename urng_t::reference; //TODO this must be a proxy.
    //!\brief Same as the const_reference type of the underlying range.
    using const_reference = reference;
    //!\brief The maximum size_type of the underlying range or the insertion buffer.
    using size_type = std::conditional_t<
        (sizeof(typename urng_t::size_type) > sizeof(typename insertion_buffer_type::size_type)),
        typename urng_t::size_type,
        typename insertion_buffer_type::size_type >;
    //!\brief Same as the difference_type of the underlying range.
    using difference_type = size_type;
    //!\brief The iterator type that enables to iterate over the modified range.
    // using iterator = journal_decorator_iterator; TODO when iterator class is defined
    //!\brief The const iterator type that enables to iterate over the modified range.
    // using const_iterator = journal_decorator_const_iterator; TODO when iterator class is defined

    /*!\name Constructors / destructor / assignment
     * \{
     */
    /*!\brief Constructs a journal_decorator from a range.
     * \param urange The range to decorate.
     *
     * This is the recommended way to construct a journal_decorator and will
     * allow for runtime efficient modification of the range without actually
     * modifying the underlying range.
     */
    explicit journal_decorator(urng_t const & urange) :
        host_ptr{&urange},
        journal_tree{{{journal_node_type::source::HOST,
                       static_cast<size_type>(urange.size()),
                       0, 0, 0}}},
        length{urange.size()}
    {}

    //!/brief Capturing of rvalues is explicitly prohibited.
    journal_decorator(urng_t const && urange) = delete;

    /*!\brief Constructs a journal_decorator by assigning a sequence of replicated values.
     * \param count The number of replications of \val and after assignment the new length.
     * \param val The value to be replicated as a new range.
     *
     * Delegates to assign(count, value).
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    journal_decorator(size_type const count, value_type const & val)
    {
        assign(count, val);
    }

    /*!\brief Construct a new journal_decorator by assigning from two iterators.
     * \tparam iterator_type Must satisfy the seqan3::input_iterator_concept and dereference to the same value_type.
     * \param begin The begin iterator of the range to be assigned.
     * \param end The end iterator of the range to be assigned.
     *
     * Delegates to assign(begin, end).
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    template <input_iterator_concept iterator_type>
    journal_decorator(iterator_type begin, iterator_type end)
    {
        assign(begin, end);
    }

    /*!\brief Assign a new journal_decorator via an initializer list.
     * \param list The initializer list to be assigned.
     *
     * Delegates to the member function assign(list.begin(), list.end()).
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    explicit journal_decorator(std::initializer_list<value_type> list) :
        journal_decorator{list.begin(), list.end()}
    {}

    /*!\brief Default construction will create a nullptr as the underlying reference.
     *
     * The string will then be constructed by pushing into the insertion buffer.
     * This is highly inefficient and not a use case of this class but required
     * to fulfill the seqan3::container_concept.
     */
    journal_decorator() = default;
    //!\brief Copy construction is defaulted.
    journal_decorator(journal_decorator const &) = default;
    //!\brief Copy assignment is defaulted.
    journal_decorator & operator=(journal_decorator const &) = default;
    //!\brief Move construction is defaulted.
    journal_decorator(journal_decorator &&) = default;
    //!\brief Move assignment is defaulted.
    journal_decorator & operator=(journal_decorator &&) = default;
    //!\brief Destructor is defaulted.
    ~journal_decorator() = default;

    /*!\brief Assign a new range via iterators to the journal_decorator.
     * \tparam iterator_type Must satisfy the seqan3::input_iterator_concept and dereference to the same value_type.
     * \param begin The begin iterator of the range to be assigned.
     * \param end The end iterator of the range to be assigned.
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    template <input_iterator_concept iterator_type>
        requires std::is_same_v<value_type, std::remove_const_t<typename std::iterator_traits<iterator_type>::value_type>>
    void assign(iterator_type begin, iterator_type end)
    {
        host_ptr = nullptr;
        insertion_buffer.assign(begin, end);
        length = insertion_buffer.size();
        journal_tree.assign({{journal_node_type::source::BUFFER,
                              static_cast<size_type>(insertion_buffer.size()),
                              0, 0, 0}});
    }

    /*!\brief Assign a new range via an initializer list to the journal_decorator.
     * \param list The initializer list to be assigned.
     *
     * Merely delegates to the member function assign(list.begin(), list.end()).
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    void assign(std::initializer_list<value_type> list)
    {
        assign(list.begin(), list.end());
    }

    /*!\brief Replace /p count times value /p val as a new range for the journal_decorator.
     * \param count The number of replications of \val and after assignment the new length.
     * \param val The value to be replicated as a new range.
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    void assign(size_type const count, value_type const & val)
    {
        host_ptr = nullptr;
        insertion_buffer.assign(count, val);
        length = insertion_buffer.size();
        journal_tree.assign({{journal_node_type::source::BUFFER,
                              static_cast<size_type>(insertion_buffer.size()),
                              0, 0, 0}});
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns whether the journal_decorator represents an empty range.
    bool empty() const
    {
        return (length == 0);
    }

    //!\brief Returns the length of the current state of the journal_decorator.
    size_type size() const
    {
        return length;
    }

    //!\brief Returns maximal size of buffer of the journal_decorator.
    size_type max_size() const
    {
        return insertion_buffer.max_size();
    }
    //!\}

protected:
    //!\brief The type of journal_node that store the modification information.
    using journal_node_type = detail::journal_node<size_type, size_type>;

    /*!\brief The container type for storing journal_nodes.
     *
     * The type properties of the container are important for runtime of many
     * member/free function of the journal decorator and should be chosen
     * carefully based on the type of application.
     */
    using journal_tree_type = typename traits_type::template journal_container_type<journal_node_type>;

    //!/brief A pointer to the underlying range that is to be decorated.
    urng_t const * host_ptr{nullptr};
    //!\brief A buffer for storing inserted novel subsequences in a concatenated fashion.
    insertion_buffer_type insertion_buffer{};
    //!\brief Stores modifications to the underlying range using journal entries.
    journal_tree_type journal_tree{};
    //!\brief Tracks the length for constant size access.
    size_type length{0};
};

} // namespace seqan3
