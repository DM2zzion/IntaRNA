
#ifndef INTARNA_INDEXRANGE_H_
#define INTARNA_INDEXRANGE_H_


#include "IntaRNA/general.h"

#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

/**
 * Defines an index region with functionality similar to pair
 *
 * @author Martin Mann
 */
class IndexRange {

public:

	//! placeholder to define the whole range (without explicit naming the last
	//! index) is defined
	static const size_t LAST_INDEX;

	//! placeholder for not-defined values
	static const size_t NA_INDEX;

	//! the start of the index range
	size_t from;

	//! the end of the index range
	size_t to;

	//! regular expression that matches valid IndexRange string encodings
	static const boost::regex regex;

public:

	/**
	 * Creates a range
	 * @param from the start index (default 0)
	 * @param to the end index (default NA_INDEX)
	 */
	IndexRange(const size_t from = 0,
			const size_t to = NA_INDEX)
		: from(from), to(to)
	{
	}

	/**
	 * Creates a range from a string encoding
	 * @param stringEncoding string encoding of the range as produced by the
	 *  ostream operator
	 */
	IndexRange(const std::string & stringEncoding)
		: from(0), to(NA_INDEX)
	{
		fromString(stringEncoding);
	}

	/**
	 * destruction
	 */
	virtual ~IndexRange() {}

	/**
	 * Checks whether or not the range encoding is ascending.
	 * @return from <= to
	 */
	bool isAscending() const
	{
		return from <= to;
	}

	/**
	 * Checks whether or not the range encoding is descending.
	 * @return from >= to
	 */
	bool isDescending() const
	{
		return from >= to;
	}

	/**
	 * Adds the given shift to the index range but ensure that the minimal value
	 * is 0.
	 * @param r the range to be added to
	 * @param shift the shift to be added to r
	 * @return an altered range or (NA_INDEX,NA_INDEX) if the range falls
	 *    completely below zero or changing from to the lower bound of 0.
	 */
	IndexRange operator + ( const int shift ) const {
		if (shift == 0) {
			return *this;
		}
		if (shift > 0) {
			return IndexRange(from+shift, to+shift);
		}
		if (to < std::abs(shift)) {
			return IndexRange( NA_INDEX, NA_INDEX);
		}
		return IndexRange( from - std::min(from, (size_t)std::abs(shift)), to + shift );
	}

	/**
	 * Substracts the given shift to the index range but ensure that the minimal
	 * value is 0.
	 * @param r the range to be altered
	 * @param shift the shift to be substracted from r
	 * @return an altered range or (NA_INDEX,NA_INDEX) if the range falls
	 *    completely below zero or changing from to the lower bound of 0.
	 */
	IndexRange operator - ( const int shift ) const {
		// forward implementation using inverted shift
		return this->operator+(-shift);
	}

	/**
	 * Checks whether or not the start of this range preceeds the given range.
	 * @param r the range to compare to
	 * @return (from<r.from) || (from==r.from && to<r.to)
	 */
	const bool operator < ( const IndexRange &r ) const {
		return ( from < r.from || (from==r.from && to<r.to) );
	}

	/**
	 * Checks whether or not two ranges are equivalent
	 * @param r the range to compare to
	 * @return ( from == r.from && to == r.to )
	 */
	const bool operator == ( const IndexRange &r ) const {
		return ( from == r.from && to == r.to );
	}

	/**
	 * Checks whether or not two ranges are different
	 * @param r the range to compare to
	 * @return !( this == r)
	 */
	const bool operator != ( const IndexRange &r ) const {
		return !( this->operator ==(r) );
	}


	/**
	 * Prints the range's boundaries to stream
	 * @param out the ostream to write to
	 * @param range the IndexRange object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const IndexRange& range)
	{
		// has to be in accordance to this.regex
		return (out <<range.from<<"-"<<range.to);
	}

	/**
	 * updates the range data from a valid string encoding (matching regex)
	 * @param stringEncoding the interval string encoding
	 * @throws std::runtime_error if stringEncoding does not match regex
	 */
	void
	fromString( const std::string & stringEncoding )
	{
		if( ! boost::regex_match(stringEncoding, regex, boost::match_perl) ) {
			throw std::runtime_error("IndexRange::fromString("+stringEncoding+") uses no valid index range string encoding");
		}
		// find split position
		const size_t splitPos = stringEncoding.find('-');
		// parse interval boundaries
		from = boost::lexical_cast<size_t>(stringEncoding.substr(0,splitPos));
		to = boost::lexical_cast<size_t>(stringEncoding.substr(splitPos+1));
	}
	
	/**
	 * Cuts an IndexRange into overlapping windows
	 * @param windowWidth the width of a window
	 * @param windowsOverlap the amount of overlap between two windows
	 * @return a vector of IndexRanges which represent windows
	 */
	std::vector<IndexRange> overlappingWindows(const size_t& windowWidth, const size_t& windowsOverlap) const 
	{
		if (windowWidth <= windowsOverlap) 
		{
			throw std::runtime_error("The window width must be larger than the overlap width");
		}
		if ((to - from + 1) <= windowsOverlap)
		{
			throw std::runtime_error("The IndexRange width must be larger than the overlap width");
		}
		
		// size_t numberOfWindows = ceil(double(to - from - windowsOverlap + 1) / (windowWidth - windowsOverlap));
		size_t x = to - from - windowsOverlap + 1;
		size_t y = windowWidth - windowsOverlap;

		if (((std::numeric_limits<size_t>::max)() - x) < y)
		{
			throw std::runtime_error("An overflow occured when calculating the number of windows");
		}

		size_t numberOfWindows = (x + y - 1) / y;
		std::vector<IndexRange> windows = std::vector<IndexRange>(numberOfWindows);

		int i = -1;
		size_t currentIndex = from;
		while (currentIndex + windowsOverlap - 1 < to)
		{
			i++;
			windows[i] = IndexRange(currentIndex, currentIndex + windowWidth - 1);
			currentIndex = windows[i].to - windowsOverlap + 1;
		}

		if (windows[i].to > to)
		{
			windows[i].to = to;
		}

		return windows;
	}
	
	/**
	 * Calculates all combinations of two windows from two different ranges
	 * @param query the first range, usually representing the query sequence
	 * @param target the second range, usually representing the target sequence
	 * @param windowWidth the width of a window (default 20)
	 * @param windowsOverlap the amount of overlap between two windows (default 10)
	 * @return a vector of pairs of two IndexRanges respecitvely windows
	 */
	static
	std::vector<std::pair<IndexRange, IndexRange>> getRangePairs(const IndexRange& query, const IndexRange& target, const size_t& windowWidth = 20, const size_t& windowsOverlap = 10)
	{
		std::vector<IndexRange> queryWindows = query.overlappingWindows(windowWidth, windowsOverlap);
		std::vector<IndexRange> targetWindows = target.overlappingWindows(windowWidth, windowsOverlap);
		std::vector<std::pair<IndexRange, IndexRange>> pairs = std::vector<std::pair<IndexRange, IndexRange>>(queryWindows.size() * targetWindows.size());

		int i = 0;

		BOOST_FOREACH(IndexRange q, queryWindows)
		{
			BOOST_FOREACH(IndexRange t, targetWindows)
			{
				pairs[i] = std::pair<IndexRange, IndexRange>(q, t);
				i++;
			}
		}

		return pairs;
	}
	
};

} // namespace

#endif /* INDEXRANGE_H_ */
