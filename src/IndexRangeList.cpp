
#include "IndexRangeList.h"

#include <algorithm>

//////////////////////////////////////////////////////////////////////

IndexRangeList::IndexRangeList()
:
list()
{
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::IndexRangeList( const IndexRangeList & toCopy )
:
list(toCopy.list)
{
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::~IndexRangeList()
{
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const size_t index ) const
{
	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > i
	const_iterator r = std::upper_bound( list.begin(), list.end(), IndexRange(index,std::numeric_limits<size_t>::max()) );
	if ( r == list.begin() ) {
		return false;
	} else {
		// go to preceding range and check if <= the end of the blocked range
		return index <= (--r)->to;
	}
}

//////////////////////////////////////////////////////////////////////

void
IndexRangeList::
push_back( const IndexRange& range )
{
#if IN_DEBUG_MODE
	if (range.isDescending())  {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") range is descending");
	}
	if (!list.empty() && list.rbegin()->to > range.from) {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") violates order given last range = "+toString(*(list.rbegin())));
	}
#endif
	// sorting should be OK (in debug mode.. ;) )
	list.push_back( range );
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator
IndexRangeList::
insert( const IndexRange& range )
{
#if IN_DEBUG_MODE
	if (range.isDescending())  {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") range is descending");
	}
#endif
	// add first member to list
	if (list.empty()) {
		list.push_back( range );
		return begin();
	} else
	// insert accordingly
	{
		// find first range that with begin > i
		List::iterator r = std::upper_bound( list.begin(), list.end(), range );
		// insert accordingly preserving sorting
		return list.insert( r, range );
	}
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::erase( IndexRangeList::iterator i ) { return list.erase( i ); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::begin() { return list.begin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::end() { return list.end(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_iterator IndexRangeList::begin() const { return list.begin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_iterator IndexRangeList::end() const { return list.end(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::reverse_iterator IndexRangeList::rbegin() { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::reverse_iterator IndexRangeList::rend() { return list.rend(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_reverse_iterator IndexRangeList::rbegin() const { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_reverse_iterator IndexRangeList::rend() const { return list.rend(); }

//////////////////////////////////////////////////////////////////////

bool IndexRangeList::empty() const { return list.empty(); }

//////////////////////////////////////////////////////////////////////

size_t IndexRangeList::size() const { return list.size(); }

//////////////////////////////////////////////////////////////////////

void IndexRangeList::clear() { return list.clear(); }

//////////////////////////////////////////////////////////////////////


