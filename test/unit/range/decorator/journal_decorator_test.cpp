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
