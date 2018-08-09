
#include "OutputHandlerInteractionList.h"

#include <algorithm>

namespace IntaRNA
{

/////////////////////////////////////////////////////////////////////////////

OutputHandlerInteractionList::
OutputHandlerInteractionList(const size_t maxToStore)
 :	storage()
	, maxToStore(maxToStore)
{
}

/////////////////////////////////////////////////////////////////////////////

OutputHandlerInteractionList::
~OutputHandlerInteractionList()
{
	// cleanup stored interactions
	for( auto it = storage.begin(); it != storage.end(); it++ ) {
		delete (*it);
	}
	storage.clear();
}

/////////////////////////////////////////////////////////////////////////////

void
OutputHandlerInteractionList::
add( const Interaction & interaction )
{
	if (interaction.isEmpty()) {
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_OutputHandlerInteractionListUpdate)
#endif
		{
		// count interaction
		reportedInteractions++;
		}
		return;
	}

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_OutputHandlerInteractionListUpdate)
#endif
	{
		// count interaction
		reportedInteractions++;
		if (storage.size() < maxToStore || lessThan_StorageContainer( &interaction, *(storage.rbegin()) )) {

			// find where to insert this interaction
			StorageContainer::iterator insertPos = std::lower_bound( storage.begin(), storage.end(), &interaction, lessThan_StorageContainer );

			// check for duplicates
			// lessThan_StorageContainer(*insertPos, &interaction) obsolete because of lower_bound
			if (insertPos == storage.end() || lessThan_StorageContainer(&interaction, *insertPos)) {				
				
				bool isLastElement = false;
				
				// remove last element if needed
				if (storage.size() >= maxToStore) {
					
					// check if insertPos points to the last element in the list which is about to get deleted -
					// if this is the case, move the iterator one element backwards
					if (insertPos == std::prev(storage.end())){
						insertPos = std::prev(insertPos);
						isLastElement = true;
					}
					
					// delete object
					delete (*(storage.rbegin()));
					// remove pointer
					storage.resize(storage.size()-1);
				}
				
				// if insertPos was moved backwards, move it back to it's former position (== storage.end())
				if (isLastElement){
					insertPos = std::next(insertPos);
				}
				
				storage.insert( insertPos, new Interaction(interaction) );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////

void
OutputHandlerInteractionList::
add( const InteractionRange & range )
{
	INTARNA_NOT_IMPLEMENTED("OutputHandlerInteractionList::add( const InteractionRange & range )");
}

/////////////////////////////////////////////////////////////////////////////

} /* namespace IntaRNA */
