// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================

#include <gtest/gtest.h>

#include <string>

#include <seqan3/range/decorator/journal_decorator_detail.hpp>
#include <seqan3/range/decorator/journal_decorator.hpp>

using namespace seqan3;
using namespace std::literals;

TEST(journal_node, construction)
{
    // default construct journal node
    detail::journal_node<int, int> node{};
    EXPECT_EQ(node.src, (detail::journal_node<int, int>::source::NONE));
    EXPECT_EQ(node.length, 0);
    EXPECT_EQ(node.virtual_position, 0);
    EXPECT_EQ(node.physical_position, 0);
    EXPECT_EQ(node.physical_origin_position, 0);

    // aggregate initialization
    detail::journal_node<int, int> node_agg{detail::journal_node<int, int>::source::HOST, 1, 2, 3, 4};
    EXPECT_EQ(node_agg.src, (detail::journal_node<int, int>::source::HOST));
    EXPECT_EQ(node_agg.length, 1);
    EXPECT_EQ(node_agg.virtual_position, 2);
    EXPECT_EQ(node_agg.physical_position, 3);
    EXPECT_EQ(node_agg.physical_origin_position, 4);
}

TEST(journal_node, comparison_operator)
{
    detail::journal_node<int, int> node_1{};
    detail::journal_node<int, int> node_2{detail::journal_node<int, int>::source::HOST, 1, 2, 3, 4};
    detail::journal_node<int, int> node_3{detail::journal_node<int, int>::source::HOST, 1, 2, 3, 4};

    EXPECT_EQ(node_2, node_3);
    EXPECT_NE(node_1, node_2);
    EXPECT_NE(node_1, node_3);
}

TEST(journal_decorator, concept_checks)
{
    EXPECT_TRUE(journal_decorator_traits_concept<journal_decorator_default_traits>);
    // EXPECT_TRUE(random_access_sequence_concept<journal_decorator<std::string>>);
}

TEST(journal_decorator, construction)
{
    // Constructing with an underlying range.
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};

    // Default construction.
    journal_decorator<std::string> journal_default{};

    // Copy construction.
    journal_decorator<std::string> journal_copied{journal};

    // Copy assignment.
    journal_decorator<std::string> journal_copy_assigned = journal;

    // Move construction.
    journal_decorator<std::string> journal_moved{std::move(journal)};

    // Move assignment.
    journal_decorator<std::string> journal_move_assigned = std::move(journal_copied);

    // Construction from count and value
    journal_decorator<std::string> journal_count_and_value{4, 'A'};

    // Construction from two iterators
    journal_decorator<std::string> journal_two_it{host.begin(), host.end()};

    // Construction from initializer list
    //journal_decorator<std::string> journal_ilist{{'A', 'C', 'T', 'G', 'A', 'C', 'G', 'T'}};
}

TEST(journal_decorator, function_size)
{
    // default has length 0
    journal_decorator<std::string> journal_default{};
    EXPECT_EQ(journal_default.size(), (std::string::size_type)0);

    // construction from a host sequence
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};
    EXPECT_EQ(host.size(), journal.size());
}

TEST(journal_decorator, function_assign)
{
    journal_decorator<std::string> journal{};

    // assign from iterators
    std::string host{"ACTG"};
    journal.assign(host.begin(), host.end());
    EXPECT_EQ(host.size(), journal.size());
    // TODO test for equality when iterator is implemented

    // assign from initializer list
    journal.assign({'A', 'C', 'T', 'G', 'A', 'C', 'G', 'T'});
    EXPECT_EQ(8u, journal.size());
    // TODO test for equality when iterator is implemented

    // assign by replicating a value
    journal.assign(5, 'A');
    EXPECT_EQ(5u, journal.size());
    // TODO test for equality when iterator is implemented
}

TEST(journal_decorator, function_empty)
{
    journal_decorator<std::string> journal{};
    EXPECT_TRUE(journal.empty());
}

TEST(journal_decorator, function_clear)
{
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};
    journal.clear();
    EXPECT_EQ(0u, journal.size());
}

TEST(journal_decorator, function_reset)
{
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};
    // TODO insert, delete or do some other modification here when available
    journal.reset();
    EXPECT_EQ(host.size(), journal.size());
}

TEST(journal_decorator_iterator, construction)
{
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};

    // Construction from journal decorator
    journal_decorator_iterator journal_it{journal};

    // Default construction.
    journal_decorator_iterator<journal_decorator<std::string>> journal_default_iterator{};

    // Copy construction.
    journal_decorator_iterator journal_it_copied{journal_it};

    // Copy assignment.
    journal_decorator_iterator journal_it_copy_assigned = journal_it;

    // Move construction.
    journal_decorator_iterator journal_it_moved{std::move(journal_it)};

    // Move assignment.
    journal_decorator_iterator journal_it_move_assigned = std::move(journal_it_copied);
}

TEST(journal_decorator_iterator, comparison_operator)
{
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};
    journal_decorator_iterator journal_it1{journal};
    journal_decorator_iterator journal_it2{journal};
    journal_decorator_iterator<journal_decorator<std::string>> journal_default_iterator{};

    EXPECT_EQ(journal_it1, journal_it2);
    EXPECT_NE(journal_it1, journal_default_iterator);
    EXPECT_NE(journal_it2, journal_default_iterator);
}

TEST(journal_decorator_iterator, dereference_operator)
{
    std::vector<std::string> host{{"ACTG", "TGCA"}};
    journal_decorator<std::vector<std::string>> journal{host};
    journal_decorator_iterator journal_it{journal};

    // dereference operator
    EXPECT_EQ(*journal_it, host[0]);

    // -> operator TODO
    // EXPECT_EQ(journal_it->size(), 4u);
}

TEST(journal_decorator_iterator, arithmetic_operators_single_node)
{
    // since so far no modification (insert, erase) are implemented for the
    // journal_decorator this first test only checks weather movement inside
    // a single journal node works.

    std::string host{"ACTGAGCTAGAGTCAGAGATCT"};
    journal_decorator<std::string> journal{host};
    journal_decorator_iterator journal_it{journal};

    // pre-increment
    journal_decorator_iterator pre_inc = ++journal_it;
    EXPECT_EQ(*pre_inc, host[1]);
    EXPECT_EQ(*journal_it, host[1]);

    // post-increment
    journal_decorator_iterator post_inc = journal_it++;
    EXPECT_EQ(*post_inc, host[1]);
    EXPECT_EQ(*journal_it, host[2]); //TODO ?

    // pre-decrement
    journal_decorator_iterator pre_dec = --journal_it;
    EXPECT_EQ(*pre_dec, host[1]);
    EXPECT_EQ(*journal_it, host[1]);

    // post-decrement
    journal_decorator_iterator post_dec = journal_it--;
    EXPECT_EQ(*post_dec, host[1]);
    EXPECT_EQ(*journal_it, host[0]); //TODO ?

    // += operator
    journal_it += 10;
    EXPECT_EQ(*journal_it, host[10]);

    // + operator
    EXPECT_EQ(*(3 + journal_it), *(journal_it + 3));
    EXPECT_EQ(*(3 + journal_it), host[13]);
    EXPECT_EQ(*(journal_it + 3), host[13]);

    // -= operator
    journal_it -= 2;
    EXPECT_EQ(*journal_it, host[8]);

    // - operator
    EXPECT_EQ(*(3 - journal_it), *(journal_it - 3));
    EXPECT_EQ(*(3 - journal_it), host[5]);
    EXPECT_EQ(*(journal_it - 3), host[5]);

    // iterator subtraction
    journal_decorator_iterator journal_it2{journal};
    EXPECT_EQ((journal_it - journal_it2), 8u);
}

TEST(journal_decorator, functions_begin_and_end)
{
    using jd_type = journal_decorator<std::string>;
    std::string host{"AAAAAAAT"};
    jd_type journal{host};
    const jd_type const_journal{host};

    // (c)begin
    auto it_begin = journal.begin();
    auto it_cbegin = journal.cbegin();
    auto it_const_begin = const_journal.begin();

    EXPECT_TRUE((std::is_same_v<decltype(it_begin), jd_type::iterator>));
    EXPECT_TRUE((std::is_same_v<decltype(it_cbegin), jd_type::const_iterator>));
    EXPECT_TRUE((std::is_same_v<decltype(it_const_begin), jd_type::const_iterator>));

    EXPECT_EQ(*it_begin, 'A');
    EXPECT_EQ(*it_cbegin, 'A');
    EXPECT_EQ(*it_const_begin, 'A');

    // (c)end
    auto it_end = journal.end();
    auto it_cend = journal.cend();
    auto it_const_end = const_journal.end();

    EXPECT_TRUE((std::is_same_v<decltype(it_end), jd_type::iterator>));
    EXPECT_TRUE((std::is_same_v<decltype(it_cend), jd_type::const_iterator>));
    EXPECT_TRUE((std::is_same_v<decltype(it_const_end), jd_type::const_iterator>));

    EXPECT_EQ(*(--it_end), 'T');
    EXPECT_EQ(*(--it_cend), 'T');
    EXPECT_EQ(*(--it_const_end), 'T');
}

TEST(journal_decorator, hetero_comparison_operator)
{
    std::string host{"ACTG"};
    std::string host2{"GGG"};

    journal_decorator<std::string> journal{host};
    journal_decorator<std::string> journal2{host2};

    EXPECT_TRUE(host == journal);
    EXPECT_TRUE(journal == host);

    EXPECT_FALSE(host2 == journal);
    EXPECT_FALSE(journal == host2);

    EXPECT_FALSE(host != journal);
    EXPECT_FALSE(journal != host);

    EXPECT_TRUE(host2 != journal);
    EXPECT_TRUE(journal != host2);
}

TEST(journal_decorator, function_at)
{
    // reference at
    std::string host{"ACTG"};
    journal_decorator<std::string> journal{host};
    EXPECT_EQ(journal.at(0), 'A');
    EXPECT_EQ(journal.at(1), 'C');
    EXPECT_EQ(journal.at(2), 'T');
    EXPECT_EQ(journal.at(3), 'G');

    // const_reference at
    const journal_decorator<std::string> const_journal{host};
    EXPECT_EQ(const_journal.at(0), 'A');
    EXPECT_EQ(const_journal.at(1), 'C');
    EXPECT_EQ(const_journal.at(2), 'T');
    EXPECT_EQ(const_journal.at(3), 'G');
}

TEST(journal_decorator, function_random_access_operator)
{
    using jd_type = journal_decorator<std::string>;

    // reference operator[]
    std::string host{"ACTG"};
    jd_type journal{host};

    EXPECT_TRUE((std::is_same_v<decltype(journal[0]), jd_type::reference>));

    EXPECT_EQ(journal[0], 'A');
    EXPECT_EQ(journal[1], 'C');
    EXPECT_EQ(journal[2], 'T');
    EXPECT_EQ(journal[3], 'G');

    // const_reference operator[]
    const jd_type const_journal{host};
    EXPECT_TRUE((std::is_same_v<decltype(const_journal[0]), jd_type::const_reference>));

    EXPECT_EQ(const_journal[0], 'A');
    EXPECT_EQ(const_journal[1], 'C');
    EXPECT_EQ(const_journal[2], 'T');
    EXPECT_EQ(const_journal[3], 'G');
}

TEST(journal_decorator, function_front)
{
    using jd_type = journal_decorator<std::string>;

    // reference front
    std::string host{"ACTG"};
    jd_type journal{host};

    EXPECT_TRUE((std::is_same_v<decltype(journal[0]), jd_type::reference>));

    EXPECT_EQ(journal.front(), 'A');

    // const_reference front
    const jd_type const_journal{host};
    EXPECT_TRUE((std::is_same_v<decltype(const_journal[0]), jd_type::const_reference>));

    EXPECT_EQ(const_journal.front(), 'A');
}

TEST(journal_decorator, function_back)
{
    using jd_type = journal_decorator<std::string>;

    // reference back
    std::string host{"ACTG"};
    jd_type journal{host};

    EXPECT_TRUE((std::is_same_v<decltype(journal[0]), jd_type::reference>));

    EXPECT_EQ(journal.back(), 'G');

    // const_reference back
    const jd_type const_journal{host};
    EXPECT_TRUE((std::is_same_v<decltype(const_journal[0]), jd_type::const_reference>));

    EXPECT_EQ(const_journal.back(), 'G');
}

TEST(journal_decorator_iterator, function_insert)
{
    std::string host{"AAAAAAAAAAAAGG"};
    std::string ins{"TT"};
    std::string ins2{"CCC"};
    journal_decorator<std::string> journal{host};
    journal_decorator<std::string> empty_journal{};
    journal_decorator_iterator journal_it{journal};

    // insert into an empty journal
    empty_journal.insert(empty_journal.begin(), host.begin(), host.end());
    EXPECT_EQ(empty_journal, host);

    // insert at the very first node
    journal_it = journal.insert(journal_it, ins.begin(), ins.end());
    EXPECT_EQ(journal, std::string("TTAAAAAAAAAAAAGG"));
    EXPECT_EQ(*journal_it, 'T');
    EXPECT_EQ(*(journal_it + 1), 'T');

    // insert in some middle HOST node that has to be split
    journal_it = journal.insert(journal_it+5, ins.begin(), ins.end());
    EXPECT_EQ(journal, std::string("TTAAATTAAAAAAAAAGG"));
    EXPECT_EQ(*journal_it, 'T');
    EXPECT_EQ(*(journal_it + 1), 'T');

    // insert in some middle BUFFER node that has to be split
    journal_it = journal.insert(journal_it-4, ins2.begin(), ins2.end());
    EXPECT_EQ(journal, std::string("TCCCTAAATTAAAAAAAAAGG"));
    EXPECT_EQ(*journal_it, 'C');
    EXPECT_EQ(*(journal_it + 1), 'C');
    EXPECT_EQ(*(journal_it + 2), 'C');

    // insert at the very end
    journal_it = journal.insert(journal.end(), ins2.begin(), ins2.end());
    EXPECT_EQ(journal, std::string("TCCCTAAATTAAAAAAAAAGGCCC"));
    EXPECT_EQ(*journal_it, 'C');
    EXPECT_EQ(*(journal_it + 1), 'C');
    EXPECT_EQ(*(journal_it + 2), 'C');
}
