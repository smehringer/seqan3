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
 * \brief This file contains helper functionality for the journal_decorator.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

namespace seqan3::detail
{

/*!\brief A journal node represents either a host or insertion segment.
 * \tparam position_type The type of a position in the host sequence or insertion buffer.
 * \tparam size_type The type of the length of segments of the host sequence or insertion buffer.
 *
 * The sequence of the segment described by a journal node can be retrieved
 * via the physical position and the length. The virtual position refers to the
 * position in the "modified" sequence. See the seqan3::journal_decorator for
 * a detailed description of how a container of journal_entries can be used
 * to represent a modified range without actually modifying the underlying range.
 * Note: Always initialize this struct via journal_node{} to ensure proper
 * default construction.
 */
template <typename position_type, typename size_type>
struct journal_node
{
    //!\brief Specifies whether a seqan3::journal_node refers to the original host or an insertion in the buffer.
    enum struct source
    {
        NONE,  //!< Default bit for default construction.
        HOST,  //!< Refers to the underlying range (host) of the seqan3::journal_decorator owning the node.
        BUFFER //!< Refers to the insertion buffer of the seqan3::journal_decorator owning the node.
    };

    source src;                      //!< Flag specifying where to find the segment sequence.
    size_type length;                //!< Length of the segment.
    position_type virtual_position;  //!< Position in the virtual "modified" range representation.
    position_type physical_position; //!< Position in the underlying range (host) or insertion buffer.

    /*!\brief Physical position in the host string.
     *
     * For SOURCE_PATCH entries, this is the position of the closest
     * SOURCE_ORIGINAL node left of it.  If there is no such node then this
     * is 0.  Unused and always 0 for unbalanced tree.
     * This variable was introduced to provide efficient mapping between virtual
     * and physical position.
     */
    position_type physical_origin_position;

    /*!\brief Compares if all member variables are the same.
     * \param rhs The journal node to compare to.
     */
    bool operator==(journal_node const & rhs) const
    {
        return (src == rhs.src &&
                length == rhs.length &&
                virtual_position == rhs.virtual_position &&
                physical_position == rhs.physical_position &&
                physical_origin_position == rhs.physical_origin_position);
    }

    /*!\brief Compares if any of the member variables is not equal.
     * \param rhs The journal node to compare to.
     */
    bool operator!=(journal_node const & rhs) const
    {
        return (!(*this == rhs));
    }
};

}  // namespace seqan3::detail
