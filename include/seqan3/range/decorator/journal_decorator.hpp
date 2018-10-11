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

#include <cassert>
#include <iterator>
#include <memory>

#include <seqan3/std/concepts>       // integral_concept
// #include <seqan3/core/concept/iterator.hpp>   // std::InputIterator
#include <seqan3/std/ranges>                 // std::ranges::RandomAccessRange
#include <seqan3/range/container/concept.hpp>// random_access_sequence_concept
#include <seqan3/range/decorator/journal_decorator_detail.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

namespace seqan3
{

/*!\interface seqan3::journal_decorator_traits_concept <>
 * \brief Checks the required traits passed to the journal_decorator
 */
//!\cond
template <typename t>
concept journal_decorator_traits_concept = requires (t v)
{
    // Journal Node Container
    typename t::template journal_container_type<detail::journal_node<uint32_t, uint32_t>>;
    // the insertion buffer must have random access and modifiers
    requires random_access_container_concept<
        typename t::template journal_container_type<detail::journal_node<uint32_t, uint32_t>>>;

    // Insertion Buffer Container
    typename t::template insertion_buffer_type<uint32_t>;
    // the insertion buffer must have random access and modifiers
    requires random_access_container_concept<
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

//!\cond
template<typename journal_decorator_type>
class journal_decorator_iterator; // Forward declaration
//!\endCond

template <std::ranges::RandomAccessRange urng_t,
          journal_decorator_traits_concept traits_type = journal_decorator_default_traits>
class journal_decorator
{
protected:
    template<typename journal_decorator_type>
    friend class journal_decorator_iterator;

    //!\brief The container type storing novel inserted sequences.
    using insertion_buffer_type = typename traits_type::template insertion_buffer_type<typename urng_t::value_type>;

public:
    //!\brief Same as the value_type of the underlying range.
    using value_type = typename urng_t::value_type;
    //!\brief Same as the reference type of the underlying range.
    using reference = std::remove_const_t<typename urng_t::value_type>; //TODO this must be a proxy.
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
    using iterator = journal_decorator_iterator<journal_decorator>;
    //!\brief The const iterator type that enables to iterate over the modified range.
    using const_iterator = journal_decorator_iterator<journal_decorator const>;

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
     * \tparam iterator_type Must satisfy the std::InputIterator and dereference to the same value_type.
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
    template <std::InputIterator iterator_type>
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
     * \tparam iterator_type Must satisfy the seqan3::std::InputIterator and dereference to the same value_type.
     * \param begin The begin iterator of the range to be assigned.
     * \param end The end iterator of the range to be assigned.
     *
     * IMPORTANT: This will not replace the underlying host sequence but delete
     * it and insert the new range into the insertion buffer. Note that as a
     * consequence every subsequent operation on the journal_decorator is on
     * the insertion buffer and the runtime/memory advantages of the journal_decorator
     * are not existent anymore.
     */
    template <std::InputIterator iterator_type>
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

    /*!\name Element Access
     * \{
     */
    /*!\brief Returns a reference to the element at the specific location \p pos, with bounds checking.
     * \param pos The position of the element to return.
     * \returns The reference to the element at the position \pos.
     * \throws std::out_of_range If the the position \pos is not within the bounds.
     */
    reference at(size_type pos)
    {
        if (pos < 0 || pos > (*this).size())
            std::out_of_range(("position " + std::to_string(pos) +
                               "is greater than range size " +
                               std::to_string((*this).size())));

        return *(iterator{*this} + pos);
    }

    //!\copydoc at
    const_reference at(size_type pos) const
    {
        if (pos < 0 || pos > (*this).size())
            std::out_of_range(("position " + std::to_string(pos) +
                               "is greater than range size " +
                               std::to_string((*this).size())));

        return *(const_iterator{*this} + pos);
    }

    /*!\brief Returns a reference to the element at the specific location \p pos.
     * \param pos The position of the element to return.
     * \returns The reference to the element at the position \pos.
     */
    reference operator[](size_type pos)
    {
        return *(iterator{*this} + pos);
    }

    //!\copydoc operator[]
    const_reference operator[](size_type pos) const
    {
        return *(const_iterator{*this} + pos);
    }

    /*!\brief Returns a reference to the first element.
     * \returns The reference to the element at the position 0.
     */
    reference front()
    {
        return *(iterator{*this});
    }

    //!\copydoc front
    const_reference front() const
    {
        return *(const_iterator{*this});
    }

    /*!\brief Returns a reference to the last element.
     * \returns The reference to the element at the position 0.
     */
    reference back()
    {
        return *(--((*this).end()));
    }

    //!\copydoc back
    const_reference back() const
    {
        return *(--((*this).cend()));
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator pointing to the start of the journal_decorator.
    iterator begin()
    {
        return iterator{*this};
    }

    //!\brief Returns an const iterator pointing to the start of the journal_decorator.
    const_iterator begin() const
    {
        return const_iterator{*this};
    }

    //!\brief Returns an const iterator pointing to the start of the const journal_decorator.
    const_iterator cbegin() const
    {
        return const_iterator{*this};
    }

    //!\brief Returns an iterator pointing to the end of the journal_decorator.
    iterator end()
    {
        return (iterator{*this}).set_to_end();
    }

    //!\brief Returns an const iterator pointing to the end of the journal_decorator.
    const_iterator end() const
    {
        return (const_iterator{*this}).set_to_end();
    }

    //!\brief Returns an const iterator pointing to the end of the const journal_decorator.
    const_iterator cend() const
    {
        return (const_iterator{*this}).set_to_end();
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

    /*!\name Modifiers
     * \{
     */
    /*!\brief Clears the journal_decorator to represent an empty string.
     *
     * Attention: This means that the journal decorator now has no valid
     * underlying range anymore. The type is still the same but the internal
     * pointer is set to a nullptr. The journal_decorator can still be used
     * but subsequent modifications cannot benefit from the initial design
     * anymore and will be less efficient.
     */
    void clear()
    {
        host_ptr = nullptr;
        insertion_buffer.clear();
        journal_tree.clear();
        length = 0;
    }

    /*!\brief Resets the journal_decorator to the original host range.
     *
     * This means that the journal decorator now represents the original host
     * range again. The host range remains valid (the private pointer
     * is still valid) but all members storing modifications are cleared. The
     * size of the journal_decorator is now equal to the underlying host again.
     */
    void reset()
    {
        insertion_buffer.clear();
        length = (*host_ptr).size();
        journal_tree = {{journal_node_type::source::HOST, length, 0, 0, 0}};
    }

    /*!\brief Inserts a range into the journal_decorator at \p pos_it
     * \tparam iterator_type Input range iterator, must satisfy the seqan3::std::InputIterator.
     * \param pos_it The iterator pointing to the position where the range is to be inserted.
     * \pram first The begin iterator of the range to be inserted.
     * \pram last The end iterator of the range to be inserted.
     * \returns An iterator pointing to the start of the inserted sequence or the
     *          same position if `first == last`.
     *
     * Note that by design the underlying host range is not modified by inserting
     * into the journal_decorator. Only the difference will be stored in the
     * internal data structure.
     * ```cpp
     * std::string host{"AAAAAA"};
     * std::string ins{"TT"};
     * journal_decorator<std::string> journal{host};
     * journal.insert(journal.begin() + 2, ins.begin(), ins.end());
     *
     * std::cout << journal << std::endl; // AATTAAAA
     * std::cout << host << std::endl;    // AAAAAA
     * ```
     *
     * ### Exception
     *
     * The exception guarantee depends on the container type storing journal nodes
     * which can be specified over the journal_decorator `traits_type`.
     * When the container type ensures the strong exception guarantee for insert,
     * push_back, swap, move assignment and iterating over it, like it is the
     * case for default container `std::vector` then insert() also guarantees
     * strong exception safety.
     *
     * ### Complexity
`    *
     * Inserting into the journal_decorator at most linear over the size of the
     * input range and linear over the current length of the journal node
     * container.
     *
     * ### Thread safety
     *
     * This container provides no thread-safety beyond the promise given also by the STL that all
     * calls to `const` member function are safe from multiple threads (as long as no thread calls
     * a non-`const` member function at the same time).
     *
     */
    template <std::ForwardIterator iterator_type>
    iterator insert(iterator pos_it, // how can this be const_iterator as documented by cppreference
                    iterator_type first,
                    iterator_type last)
    {
        insertion_buffer.insert(insertion_buffer.end(), first, last);
        return insert_insertion_node(pos_it, static_cast<size_type>(std::distance(first, last)));
    }

    iterator insert(iterator pos_it, std::initializer_list<value_type> ilist)
    {
        insertion_buffer.insert(insertion_buffer.end(), ilist.begin(), ilist.end());
        return insert_insertion_node(pos_it, ilist.size());
    }

    iterator insert(iterator pos_it, value_type const & value)
    {
        insertion_buffer.push_back(value);
        return insert_insertion_node(pos_it, 1);
    }

    iterator insert(iterator pos_it, size_type count, value_type const & value)
    {
        insertion_buffer.insert(insertion_buffer.end(), count, value);
        return insert_insertion_node(pos_it, count);
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

    iterator insert_insertion_node(iterator pos_it, size_type insertion_length)
    {
        /* In general, we need to copy the range into the insertion_buffer
         * and add a new journal node pointing to it.
         * The following cases need to be considered for inserting a node:
         * (1) The journal_tree is empty -> we cannot access the iterator
         *     pos_it and need to insert directly.
         * (2) We want to insert at the very end -> add node right of last node
         * (3) We want to insert in between node -> add node left of pos_it
         * (4) We want to insert inside a node   -> split node into 3
        */
        using source_type = typename journal_node_type::source;
        using tree_iterator_type = typename journal_tree_type::iterator;

        if (insertion_length == 0)
            return pos_it; // nothing to insert

        // New node source - the source type is always the insertion buffer
        source_type const new_node_source{journal_node_type::source::BUFFER};
        // New node length - the length is computed once from the input range
        size_type const new_node_length{insertion_length};

        if (journal_tree.empty()) // case (1)
        {
            // update insertion buffer
            // insertion_buffer.insert(insertion_buffer.end(), first, last);
            // update journal nodes
            journal_tree.push_back({new_node_source, new_node_length, 0, 0, 0});
            // update length
            length += new_node_length;

            return (iterator{*this}).set(journal_tree.begin(), 0);
        }

        // New node virtual position - the vpos is the position pointed to by the iterator pos_it
        size_type const new_node_virtual_position{pos_it.current_node_it->virtual_position + pos_it.offset};
        // New node physical position - the ppos is at end of the insertion buffer where the new range is inserted
        size_type const new_node_physical_position{insertion_buffer.size() - insertion_length};
        // New node physical origin position - will be set differently in each case
        size_type new_node_physical_origin_position{0};
        // current node iterator where sequence is to be inserted.
        tree_iterator_type const current_node_it = pos_it.current_node_it;
        // cash position of current_node_it
        size_type current_node_pos = current_node_it - journal_tree.begin();
        // initialize the position of nodes to be updated later on
        size_type shift_right_of = current_node_pos;

        assert(new_node_virtual_position <= length); // pos_it should be in range [0, length]
        assert(current_node_it->virtual_position <= new_node_virtual_position);
        assert(current_node_it->virtual_position + current_node_it->length >= new_node_virtual_position);

        if (new_node_virtual_position == length) // case (2)
        {
            // update journal nodes
            journal_tree.push_back({new_node_source,
                                    new_node_length,
                                    new_node_virtual_position,
                                    new_node_physical_position,
                                    journal_tree.end()->physical_origin_position});
            ++current_node_pos;
            ++shift_right_of;
        }
        else if (current_node_it->virtual_position == new_node_virtual_position) // case (3)
        {
            if (current_node_it != journal_tree.begin())
                new_node_physical_origin_position = (current_node_it - 1)->physical_origin_position;

            // update journal nodes
            journal_tree.insert(current_node_it, {new_node_source,
                                                  new_node_length,
                                                  new_node_virtual_position,
                                                  new_node_physical_position,
                                                  new_node_physical_origin_position});
            ++shift_right_of;
        }
        else // case (4)
        {
            // split current node into two and place the new insertion node in
            // the middle ---> lhs_node ~ new_node ~ rhs_node

            std::array<journal_node_type, 3> buffer{}; // The buffer to store the 3 nodes

            // Compute the size of the left and right flank of the two nodes
            // flanking the new insertion node.
            size_type const left_flank_len{new_node_virtual_position - current_node_it->virtual_position};
            size_type const right_flank_len{current_node_it->length - left_flank_len};

            buffer[0] = {current_node_it->src,
                         left_flank_len,
                         current_node_it->virtual_position,
                         current_node_it->physical_position,
                         current_node_it->physical_origin_position};

            buffer[1] = {new_node_source,
                         new_node_length,
                         new_node_virtual_position,
                         new_node_physical_position,
                         current_node_it->physical_origin_position};

            buffer[2] = {current_node_it->src,
                         right_flank_len,
                         new_node_virtual_position + new_node_length,
                         current_node_it->physical_position + left_flank_len,
                         ((current_node_it->src == journal_node_type::source::BUFFER) ?
                            current_node_it->physical_origin_position :
                            current_node_it->physical_origin_position + left_flank_len)};

            ++current_node_pos;
            std::swap(*current_node_it, buffer[2]);
            journal_tree.insert(current_node_it, buffer.begin(), buffer.begin() + 2);
            shift_right_of += 3;
        }

        // Update journal entries right of pos.
        assert(shift_right_of <= journal_tree.size());
        for (auto it = journal_tree.begin() + shift_right_of; it != journal_tree.end(); ++it)
        {
            it->virtual_position += new_node_length;
            if (it != journal_tree.begin() && it->src == journal_node_type::source::BUFFER)
                it->physical_origin_position = (it - 1)->physical_origin_position;
        }

        // update the insertion buffer
        // insertion_buffer.insert(insertion_buffer.end(), first, last);
        // update the length
        length += new_node_length;

        return (iterator{*this}).set((journal_tree.begin() + current_node_pos), 0);
    }
};

/*!\name Comparison operators.
 * \brief Enables comparison between a journal_decorator and it's underlying range.
 * \tparam range_type The range to compare to a decorated range of the same type.
 * \tparam traits_type The traits type of the decorated range.
 * \param rhs The first range/decorator to compare.
 * \param rhs The second range/decorator to compare.
 *
 * This way we can check whether e.g. `host == journal_decorated_host` or
 * `journal_decorated_host == host`. Note that one of the comparees must be a
 * journal_decorator.
 * \{
 */
//!brief Return `true` if \p rhs == \p lhs, `false` otherwise.
template <typename range_type, typename traits_type>
constexpr bool operator==(range_type const & rhs, journal_decorator<range_type, traits_type> const & lhs)
{
    if (rhs.size() != lhs.size())
        return false;

    auto rhs_it = rhs.begin();
    auto lhs_it = lhs.begin();

    for (; rhs_it != rhs.end(); ++rhs_it, ++lhs_it)
        if (*rhs_it != *lhs_it)
            return false;

    return true;
}

//!brief Return `true` if \p rhs == \p lhs, `false` otherwise.
template <typename range_type, typename traits_type>
constexpr bool operator==(journal_decorator<range_type, traits_type> const & rhs, range_type const & lhs)
{
    return (lhs == rhs);
}

//!brief Return `true` if \p rhs != \p lhs, `false` otherwise.
template <typename range_type, typename traits_type>
constexpr bool operator!=(range_type const & rhs, journal_decorator<range_type, traits_type> const & lhs)
{
    return !(rhs == lhs);
}

//!brief Return `true` if \p rhs != \p lhs, `false` otherwise.
template <typename range_type, typename traits_type>
constexpr bool operator!=(journal_decorator<range_type, traits_type> const & rhs, range_type const & lhs)
{
    return !(rhs == lhs);
}
//!\}

/*!\brief The iterator for the journal_decorator.
 * \tparam The journal_decorator type to itereate over.
 *
 * You can initialize the iterator on an journal_decorator and then access
 * the decorated sequence just like any other iterator.
 *
 * ```cpp
 * std::string host{"ACTG"};
 * journal_decorator<std::string> journal{host};
 *
 * // construct the iterator
 * journal_decorator_iterator it{journal}; // iterator now points to the start.
 *
 * std::cout << *it << std::endl; // "A"
 * ```
 */
template<typename journal_decorator_type>
class journal_decorator_iterator : public
    detail::random_access_iterator_base<journal_decorator_type, journal_decorator_iterator>
{
private:
    using base = detail::random_access_iterator_base<journal_decorator_type, journal_decorator_iterator>;
public:
    //!\brief Same as the value_type of the journal_decorator.
    using value_type = typename journal_decorator_type::value_type;
    //!\brief Same as the size_type of the journal_decorator.
    using size_type = typename journal_decorator_type::size_type;
    //!\brief Same as the difference_type of the journal_decorator.
    using difference_type = typename journal_decorator_type::difference_type;
    //!\brief Same as the reference type of the journal_decorator.
    using reference = typename journal_decorator_type::reference;
    //!\brief Same as the const_reference type of the journal_decorator.
    using const_reference = typename journal_decorator_type::const_reference;
    //!\brief The pointer type of the value_type.
    using pointer = value_type *;
    //!\brief Import the parent's constructors of the random_access_iterator_base class.
    using base::base;

    /*!\name Constructors / destructor / assignment
     * \{
     */
    //!\brief Construction from the journal_decorator.
    journal_decorator_iterator(journal_decorator_type & decorator) :
        decorator_ptr{&decorator},
        current_node_it{decorator.journal_tree.begin()}
    {}

    //!\brief Default construction is defaulted.
    journal_decorator_iterator() = default;
    //!\brief Copy construction is defaulted.
    journal_decorator_iterator(journal_decorator_iterator const &) = default;
    //!\brief Copy assignment is defaulted.
    journal_decorator_iterator & operator=(journal_decorator_iterator const &) = default;
    //!\brief Move construction is defaulted.
    journal_decorator_iterator(journal_decorator_iterator &&) = default;
    //!\brief Move assignment is defaulted.
    journal_decorator_iterator & operator=(journal_decorator_iterator &&) = default;
    //!\brief Destructor is defaulted.
    ~journal_decorator_iterator() = default;
    //!\}

    /*!\name Comparison operators
     * \brief Compare iterators by position.
     * \{
     */
    bool operator==(journal_decorator_iterator const & rhs) const noexcept
    {
        return (current_node_it == rhs.current_node_it && offset == rhs.offset);
    }

    bool operator!=(journal_decorator_iterator const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    bool operator<(journal_decorator_iterator const & rhs) const noexcept
    {
        if (current_node_it == rhs.current_node_it)
            return offset < rhs.offset;
        return current_node_it < rhs.current_node_it;
    }

    bool operator>(journal_decorator_iterator const & rhs) const noexcept
    {
        if (current_node_it == rhs.current_node_it)
            return offset > rhs.offset;
        return current_node_it > rhs.current_node_it;
    }

    bool operator<=(journal_decorator_iterator const & rhs) const noexcept
    {
        if (current_node_it == rhs.current_node_it)
            return offset <= rhs.offset;
        return current_node_it <= rhs.current_node_it;
    }

    bool operator>=(journal_decorator_iterator const & rhs) const noexcept
    {
        if (current_node_it == rhs.current_node_it)
            return offset >= rhs.offset;
        return current_node_it >= rhs.current_node_it;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Pre-increment, return updated iterator.
    journal_decorator_iterator & operator++() noexcept
    {
        ++offset;
        if (offset == (*current_node_it).length) // I am at the end of the segment
        {
            if (current_node_it != (decorator_ptr->journal_tree.end() - 1))
            {
                ++current_node_it;
                offset = 0;
            }
        }
        return *this;
    }

    //!\brief Post-increment, return previous iterator state.
    journal_decorator_iterator operator++(int) noexcept
    {
        journal_decorator_iterator cpy{*this};
        ++(*this);
        return cpy;
    }

    //!\brief Pre-decrement, return updated iterator.
    journal_decorator_iterator & operator--() noexcept
    {
        if (offset == 0) // I am at the end of the segment
        {
            --current_node_it;
            offset = (*current_node_it).length; // points to the end but will be decremented at the end of this if clause
            assert(offset != 0);
        }
        --offset;
        return *this;
    }

    //!\brief Post-decrement, return previous iterator state.
    journal_decorator_iterator operator--(int) noexcept
    {
        journal_decorator_iterator cpy{*this};
        --(*this);
        return cpy;
    }

    //!\brief Forward this iterator.
    journal_decorator_iterator & operator+=(difference_type const skip) noexcept
    {
        offset += skip;
        // TODO can/shall this be done by binary seach in logarithmic time?
        while (offset > (*current_node_it).length) // I am at the end of the segment
        {
            offset -= (*current_node_it).length;
            ++current_node_it;
        }

        // last move must be handled separately in case the iterator will point to the very end
        if (offset == (*current_node_it).length)
        {
            if (current_node_it != (decorator_ptr->journal_tree.end() - 1))
            {
                ++current_node_it;
                offset = 0;
            }
        }

        return *this;
    }

    //!\brief Forward copy of this iterator.
    journal_decorator_iterator operator+(difference_type const skip) const noexcept
    {
        journal_decorator_iterator cpy{*this};
        return cpy += skip;
    }

    //!\brief Decrement iterator by skip.
    journal_decorator_iterator & operator-=(difference_type const skip) noexcept
    {
        difference_type skip_cpy{skip};
        while (offset < skip_cpy) // otherwise I would surpass the begin of the segment
        {
            skip_cpy -= offset;
            --current_node_it;
            offset = (*current_node_it).length; // points to the end but will be decremented at the end of this if clause
        }
        offset -= skip_cpy;
        return *this;
    }

    //!\brief Return decremented copy of this iterator.
    journal_decorator_iterator operator-(difference_type const skip) const noexcept
    {
        journal_decorator_iterator cpy{*this};
        return cpy -= skip;
    }

    //!\brief Return offset between this and remote iterator's position.
    difference_type operator-(journal_decorator_iterator const lhs) const noexcept
    {
        return static_cast<difference_type>(
            ((*current_node_it).virtual_position + offset) - ((*lhs.current_node_it).virtual_position + lhs.offset));
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \{
    */
    //!\brief Dereference operator returns element currently pointed at.
    reference operator*() const
    {
        size_type src_pos = (*current_node_it).physical_position + offset;

        switch ((*current_node_it).src)
        {
            case journal_node_type::source::HOST:
                return (*((*decorator_ptr).host_ptr))[src_pos];
            case journal_node_type::source::BUFFER:
                return (*decorator_ptr).insertion_buffer[src_pos];
            default:
                assert(false && "Invalid segment source!");
        }
    }

//    TODO when proxy iterator is implemented
//    //!\brief Return pointer to this iterator.
//    std::unique_ptr<value_type> operator->() const
//    {
//        return std::unique_ptr<value_type>{new value_type{*(*this)}};
//    }

//    //!\brief Return underlying container value currently pointed at.
//    reference operator[](position_type const n) const noexcept(noexcept((*host)[pos+n]))
//    {
//        return (*host)[pos + n];
//    }
    //!\}
private:
    template <typename urng_t, typename traits_type>
    friend class journal_decorator;

    //!\brief Same as the journal_node_type of the journal_decorator
    using journal_node_type = typename journal_decorator_type::journal_node_type;
    //!\brief Same as the journal_tree_type of the journal_decorator
    using journal_tree_type = typename journal_decorator_type::journal_tree_type;
    //!\brief The (const) iterator type of the journal container type
    using journal_tree_iterator_type = std::conditional_t<
        std::is_const_v<journal_decorator_type>,
        typename journal_tree_type::const_iterator,
        typename journal_tree_type::iterator>;

    //!\brief sets iterator to the end in constant time.
    journal_decorator_iterator set(journal_tree_iterator_type node_it, size_type off)
    {
        current_node_it = node_it;
        offset = off;
        return *this;
    }

    //!\brief sets iterator to the end in constant time.
    journal_decorator_iterator set_to_end()
    {
        journal_tree_iterator_type tree_end{decorator_ptr->journal_tree.end() - 1};
        return set(tree_end, tree_end->length);
    }

    //!\brief A pointer to the journal_decorator to access or update the data members.
    journal_decorator_type * decorator_ptr{nullptr};

    //!\brief An iterator to the current journal_node of the journal_tree of the journal_decorator.
    journal_tree_iterator_type current_node_it{};

    //!/brief The relative offset to the current node start position to infer the current position.
    size_type offset{0};
};

} // namespace seqan3
