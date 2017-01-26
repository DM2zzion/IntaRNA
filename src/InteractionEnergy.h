
#ifndef INTERACTIONENERGY_H_
#define INTERACTIONENERGY_H_


#include "general.h"
#include "Interaction.h"
#include "Accessibility.h"
#include "ReverseAccessibility.h"
#include <emmintrin.h>
/**
 * Abstract utility class that covers necessary energy related functionalities
 * for the interaction energy computation given two RNAs.
 *
 * @author Martin Mann 2014
 *
 */
class InteractionEnergy {

public:

	/**
	 * Container that provides the different energy contributions for an interaction
	 */
	struct EnergyContributions {
	public:
		//! the energy for all intermolecular loops
		E_type loops;
		//! the energy penalty for initiating the interaction
		E_type init;
		//! the energy penalty for making the interaction site accessible in seq1
		E_type ED1;
		//! the energy penalty for making the interaction site accessible in seq2
		E_type ED2;
		//! the energy for the dangling ends at the left end of the interaction
		E_type dangleLeft;
		//! the energy for the dangling ends at the right end of the interaction
		E_type dangleRight;
		//! the energy penalty for the left end of the interaction
		E_type endLeft;
		//! the energy penalty for the right end of the interaction
		E_type endRight;
	};


	/**
	 * Construct energy utility object given the accessibility ED values for
	 * two RNA sequences.
	 *
	 * @param accS1 accessibility of the first sequence
	 * @param accS2 accessibility of the second sequence (reversed to 3'-5'
	 *          index reading)
	 * @param maxInternalLoopSize1 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 1, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j1-i1) <= (1+maxInternalLoopSize1)
	 * @param maxInternalLoopSize2 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 2, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j2-i2) <= (1+maxInternalLoopSize2)
	 *
	 */
	InteractionEnergy( const Accessibility & accS1
			, const ReverseAccessibility & accS2
			, const size_t maxInternalLoopSize1
			, const size_t maxInternalLoopSize2
			);

	/**
	 * destruction
	 */
	virtual ~InteractionEnergy();


	/**
	 * Provides the overall energy for an interaction from [i1,j1] in the first
	 * sequence and [i2,j2] in the second sequence given the hybridization
	 * energy contribution.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 * @param hybridE the hybridization energy for the interaction
	 *
	 * @return E = hybridE
	 * 				+ ED1(i1,j1) + ED2(i2,j2)
	 * 				+ Edangle(i1,i2) + Edangle(j1,j2)
	 * 				+ Eend(i1,i2) + Eend(j1,j2)
	 */
	virtual
	E_type
	getE( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type hybridE ) const;

	virtual
	__m128
	getE_SSE( const __m128i i1, const __m128i j1
				, const __m128i i2, const __m128i j2
				, const __m128 hybridE ) const;
		
	/**
	 * Provides details about the energy contributions for the given interaction
	 *
	 * @param interaction the interaction of interest
	 *
	 * @return the individual energy contributions
	 */
	virtual
	EnergyContributions
	getE_contributions( const Interaction & interaction ) const;

	/**
	 * Checks whether or not two positions can form a base pair
	 * @param i1 index in first sequence
	 * @param i2 index in second sequence
	 * @return true if seq1(i1) can form a base pair with seq2(i2)
	 */
	virtual
	bool
	areComplementary( const size_t i1, const size_t i2 ) const;

	/**
	 * Length of sequence 1
	 * @return length of sequence 1
	 */
	virtual
	size_t
	size1() const;

	/**
	 * Length of sequence 2
	 * @return length of sequence 2
	 */
	virtual
	size_t
	size2() const;

	/**
	 * Provides the ED penalty for making a region with sequence 1 accessible
	 *
	 * @param i1 the start of the accessible region
	 * @param j1 the end of the accessible region
	 * @return the ED value for [i1,j1]
	 */
	virtual
	E_type
	getED1( const size_t i1, const size_t j1 ) const;

	/**
	 * Provides the ED penalty for making a region with (the reversed)
	 * sequence 2 accessible
	 *
	 * @param i2 the start of the accessible region
	 * @param j2 the end of the accessible region
	 * @return the ED value for [i2,j2]
	 */
	virtual
	E_type
	getED2( const size_t i2, const size_t j2 ) const;

	/**
	 * Whether or not position i is accessible for interaction in sequence 1
	 * @param i the position of interest in sequence 1
	 * @return true if the position can partake in an interaction; false otherwise
	 */
	virtual
	bool
	isAccessible1( const size_t i ) const;

	/**
	 * Whether or not position i is accessible for interaction in sequence 2
	 * @param i the position of interest in sequence 2
	 * @return true if the position can partake in an interaction; false otherwise
	 */
	virtual
	bool
	isAccessible2( const size_t i ) const;


	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 1 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i1 the start of the structured region of seq1
	 * @param j1 the end of the structured region of seq1
	 * @return the ES value for [i1,j1] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES1( const size_t i1, const size_t j1 ) const;

	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 2 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i2 the start of the structured region of seq2
	 * @param j2 the end of the structured region of seq2
	 * @return the ES value for [i2,j2] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES2( const size_t i2, const size_t j2 ) const;

	/**
	 * Provides the energy contribution for a given number of unpaired
	 * nucleotides under the
	 * assumption that the region is part of an (intermolecular) multiloop.
	 *
	 * @param numUnpaired the number of unpaired bases
	 * @return the energy contribution of the given number of unpaired bases
	 *         within an intramolecular multiloop
	 */
	virtual
	E_type
	getEU( const size_t numUnpaired ) const = 0;

	/**
	 * Provides the duplex initiation energy.
	 *
	 * @return the energy for duplex initiation
	 */
	virtual
	E_type
	getE_init( ) const = 0;

	/**
	 * Computes the energy estimate for the 'left side' interaction loop region
	 * closed by the intermolecular base pairs (i1,i2) and enclosing (j1,j2)
	 * where the regions [i1,j1] and [i2,j2] are considered unpaired or E_INF
	 * is the internal loop size exceeds the allowed maximum (see constructor).
	 *
	 * Note, the right interaction base pair (j1,j2) is not included in the
	 * returned energy value.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the energy for the loop
	 *         or E_INF if the allowed loop size is exceeded or no valid internal loop boundaries
	 */
	virtual
	E_type
	getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const = 0;


	/**
	 * Computes the dangling end energy penalties for the left side
	 * (i1-1 and i2-1) of the interaction closed by the intermolecular
	 * base pair (i1,i2).
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the dangling end penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_danglingLeft( const size_t i1, const size_t i2 ) const = 0;


	/**
	 * Computes the dangling end energy penalties for the right side
	 * (j1+1 and j2+1) of the interaction closed by the intermolecular
	 * base pair (j1,j2).
	 *
	 * @param j1 the index of the first sequence interacting with j2
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 * @return the dangling end penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_danglingRight( const size_t j1, const size_t j2 ) const = 0;

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "left side" (i1 = 5' end of seq1 of the interaction)
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the loop closure penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_endLeft( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "right side" (j1 = 3' end of seq1 of the interaction)
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the loop closure penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_endRight( const size_t j1, const size_t j2 ) const = 0;

	/**
	 * Computes the probability of the dangling ends for the left side
	 * (i1-1 and i2-1) of the interaction closed by the intermolecular
	 * base pair (i1,i2) for an interaction of [i1,j1] with [i2,j2].
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the dangling end probability for the left side of the interaction
	 */
	virtual
	E_type
	getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

	/**
	 * Computes the probability of the dangling ends for the right side
	 * (j1+1 and j2+1) of the interaction closed by the intermolecular
	 * base pair (j1,j2) for an interaction of [i1,j1] with [i2,j2].
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the dangling end probability for the right side of the interaction
	 */
	virtual
	E_type
	getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

	/**
	 * Access to the accessibility object of the first sequence
	 * (including sequence access)
	 * @return the accessibility object for the first sequence
	 */
	virtual
	const Accessibility &
	getAccessibility1() const;

	/**
	 * Access to the accessibility object of the second sequence
	 * (including sequence access)
	 * @return the reverse accessibility object for the second sequence
	 */
	virtual
	const ReverseAccessibility &
	getAccessibility2() const;

	/**
	 * Access to the maximal size of an unpaired stretch within seq1 within
	 * an interaction.
	 * @return the maximal internal loop size of an interaction for seq1
	 */
	const size_t getMaxInternalLoopSize1() const {
		return maxInternalLoopSize1;
	}

	/**
	 * Access to the maximal size of an unpaired stretch within seq2 within
	 * an interaction.
	 * @return the maximal internal loop size of an interaction for seq2
	 */
	const size_t getMaxInternalLoopSize2() const {
		return maxInternalLoopSize2;
	}

	/**
	 * Access to the normalized temperature for Boltzmann weight computation
	 */
	virtual
	E_type
	getRT() const = 0;


	/**
	 * Provides the best energy gain via stacking possible for this energy
	 * model
	 * @return the best stacking energy gain produced by getE_interLoop()
	 */
	virtual
	E_type
	getBestE_interLoop() const = 0;

	/**
	 * Provides the best energy gain possible for left/right dangle
	 * for this energy model
	 * @return the best initiation energy gain produced by getE_danglingLeft() or
	 *          getE_danglingRight()
	 */
	virtual
	E_type
	getBestE_dangling() const = 0;

	/**
	 * Provides the best energy gain possible for left/right interaction ends
	 * for this energy model
	 * @return the best end energy gain produced by getE_endLeft() or
	 *          getE_endRight()
	 */
	virtual
	E_type
	getBestE_end() const = 0;


	/**
	 * Provides the Boltzmann weight for a given energy.
	 * @param energ the energy the Boltzmann weight is to be computed for
	 * @return the Boltzmann weight, i.e. exp( - energy / RT );
	 */
	virtual
	E_type
	getBoltzmannWeight( const E_type energy ) const ;


	/**
	 * Provides the base pair encoding for the given indices.
	 * @param i1 the index in the first sequence
	 * @param i2 the index in the (reversed) second sequence
	 * @return the according base pair (i1,reverseIdx(i2))
	 */
	virtual
	Interaction::BasePair
	getBasePair( const size_t i1, const size_t i2 ) const;


	/**
	 * Provides the index within the first sequence of the given base pair.
	 * @return the index of the first sequence within the base pair encoding
	 */
	virtual
	size_t
	getIndex1( const Interaction::BasePair & bp ) const;


	/**
	 * Provides the index within the second sequence of the given base pair.
	 * @return the index of the second sequence within the base pair encoding
	 */
	virtual
	size_t
	getIndex2( const Interaction::BasePair & bp ) const;


protected:

	//! accessibility values for sequence S1
	const Accessibility & accS1;

	//! accessibility values for sequence S2 (reversed index order)
	const ReverseAccessibility & accS2;

	//! maximally allowed unpaired range between two base pairs in sequence S1
	//! forming an intermolecular internal loop
	const size_t maxInternalLoopSize1;

	//! maximally allowed unpaired range between two base pairs in sequence S2
	//! forming an intermolecular internal loop
	const size_t maxInternalLoopSize2;

	/**
	 * Checks whether or not the given indices are valid index region within the
	 * sequence for an intermolecular loop and do not violate the maximal
	 * internal loop size.
	 * @param seq the sequence the indices correspond to
	 * @param i begin index of the region in the sequence
	 * @param j end index of the region in the sequence
	 * @param maxInternalLoopSize the maximally allowed distance of i and j, ie.
	 *        (j-i+1) <= maxInternalLoopSize
	 *
	 * @return true if the indices are fulfilling 0 <= i <= j < seq.length,
	 *           both sequence positions denote non-ambiguous nucleotides (!= N)
	 *           and (j-i+1) <= maxInternalLoopSize; false otherwise
	 */
	static
	bool
	isAllowedLoopRegion( const RnaSequence& seq, const size_t i, const size_t j, const size_t maxInternalLoopSize );

	/**
	 * Checks whether or not the given indices mark valid internal loop
	 * boundaries, i.e.
	 *  - (i1,i2) and (j1,j2) are complementary
	 *  - i1..j1 and i2..j2 are allowed loop regions
	 *  - no boundary overlap ( (j1-i1==0 && j2-i2==0) || (j1-i1>0 && j2-i2>0) )
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return true if the boundaries are sound for internal loop calculation;
	 *         false otherwise
	 */
	bool
	isValidInternalLoop( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

};




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
InteractionEnergy::InteractionEnergy( const Accessibility & accS1
				, const ReverseAccessibility & accS2
				, const size_t maxInternalLoopSize1
				, const size_t maxInternalLoopSize2
		)
  :
	accS1(accS1)
	, accS2(accS2)
	, maxInternalLoopSize1(maxInternalLoopSize1)
	, maxInternalLoopSize2(maxInternalLoopSize2)

{
}

////////////////////////////////////////////////////////////////////////////

inline
InteractionEnergy::~InteractionEnergy()
{
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergy::
isAllowedLoopRegion( const RnaSequence& seq, const size_t i, const size_t j, const size_t maxInternalLoopSize )
{
	// ensure index and loop size validity
	return	   i < seq.size()
			&& j < seq.size()
			&& seq.asString().at(i) != 'N'
			&& seq.asString().at(j) != 'N'
			&& i <= j
			&& (j-i) <= (1+maxInternalLoopSize);

}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergy::
areComplementary( const size_t i1, const size_t i2 ) const
{
	return RnaSequence::areComplementary( accS1.getSequence(), accS2.getSequence(), i1, i2);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergy::
size1() const
{
	return getAccessibility1().getSequence().size();
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergy::
size2() const
{
	return getAccessibility2().getSequence().size();
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergy::
isValidInternalLoop( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return
		   (j1-i1>0 && j2-i2>0)
		&& areComplementary( i1, i2)
		&& areComplementary( j1, j2)
		&& InteractionEnergy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& InteractionEnergy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
		;
}

////////////////////////////////////////////////////////////////////////////

inline
const Accessibility &
InteractionEnergy::
getAccessibility1() const
{
	return accS1;
}

////////////////////////////////////////////////////////////////////////////

inline
const ReverseAccessibility &
InteractionEnergy::
getAccessibility2() const
{
	return accS2;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getBoltzmannWeight( const E_type e ) const
{
	// TODO can be optimized when using exp-energies from VRNA
	return std::exp( - e / getRT() );
}

////////////////////////////////////////////////////////////////////////////

inline
Interaction::BasePair
InteractionEnergy::
getBasePair( const size_t i1, const size_t i2 ) const
{
	return Interaction::BasePair( i1, getAccessibility2().getReversedIndex(i2) );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergy::
getIndex1( const Interaction::BasePair & bp ) const
{
	return bp.first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergy::
getIndex2( const Interaction::BasePair & bp ) const
{
	return getAccessibility2().getReversedIndex( bp.second );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getED1( const size_t i1, const size_t j1 ) const
{
	return getAccessibility1().getED( i1, j1 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getED2( const size_t i2, const size_t j2 ) const
{
	return getAccessibility2().getED( i2, j2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getES1( const size_t i1, const size_t j1 ) const
{
#if IN_DEBUG_MODE
	// sanity check
	if (i1>j1) throw std::runtime_error("InteractionEnergy::getES1(i1="+toString(i1)+" > j1="+toString(j1));
	if (j1>=size1()) throw std::runtime_error("InteractionEnergy::getES1() : j1="+toString(j1)+" >= size1()="+toString(size1()));
#endif

	// return computed value
	return accS1.getES(i1,j1);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getES2( const size_t i2, const size_t j2 ) const
{
#if IN_DEBUG_MODE
	// sanity check
	if (i2>j2) throw std::runtime_error("InteractionEnergy::getES2(i2="+toString(i2)+" > j2="+toString(j2));
	if (j2>=size2()) throw std::runtime_error("InteractionEnergy::getES2() : j2="+toString(j2)+" >= size2()="+toString(size2()));
#endif

	// return computed value
	return accS2.getES(i2,j2);

}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergy::
isAccessible1( const size_t i ) const
{
	return
			(!getAccessibility1().getSequence().isAmbiguous(i))
			&& getAccessibility1().getAccConstraint().isAccessible(i)
			;
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergy::
isAccessible2( const size_t i ) const
{
	return
			(!getAccessibility2().getSequence().isAmbiguous(i))
			&& getAccessibility2().getAccConstraint().isAccessible(i);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// initial probabilities
	E_type probDangle1 = 1.0, probDangle2 = 1.0;

	// if dangle1 possible
	if (i1>0)  {
		// Pr( i1-1 is unpaired | i1..j1 unpaired )
		probDangle1 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED1(i1-1,j1)-getED1(i1,j1) )
							)
					)
			;
	}
	// if dangle2 possible
	if (i2>0)  {
		// Pr( i2-1 is unpaired | i2..j2 unpaired )
		probDangle2 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED2(i2-1,j2)-getED2(i2,j2) )
							)
					)
			;
	}

	// get overall probability
	return probDangle1 * probDangle2;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// initial probabilities
	E_type probDangle1 = 1.0, probDangle2 = 1.0;

	// if dangle1 possible
	if (j1+1<size1())  {
		// Pr( j1+1 is unpaired | i1..j1 unpaired )
		probDangle1 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED1(i1,j1+1)-getED1(i1,j1) )
							)
					)
			;
	}
	// if dangle2 possible
	if (j2+1<size2())  {
		// Pr( j2+1 is unpaired | i2..j2 unpaired )
		probDangle2 =
			std::max( (E_type)0.0
					, std::min( (E_type)1.0
							, getBoltzmannWeight( getED2(i2,j2+1)-getED2(i2,j2) )
							)
					)
			;
	}

	// get overall probability
	return probDangle1 * probDangle2;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergy::
getE( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE ) const
{
	// check if hybridization energy is not infinite
	if ( E_isNotINF(hybridE) ) {
		// compute overall interaction energy
		// std::cout << "foo" << std::endl;
		return hybridE
				// accessibility penalty
				+ getED1( i1, j1 )
				+ getED2( i2, j2 )
				// dangling end penalty
				// weighted by the probability that ends are unpaired
				+ (getE_danglingLeft( i1, i2 )*getPr_danglingLeft(i1,j1,i2,j2))
				+ (getE_danglingRight( j1, j2 )*getPr_danglingRight(i1,j1,i2,j2))
				// helix closure penalty
				+ getE_endLeft( i1, i2 )
				+ getE_endRight( j1, j2 )
				;
	} else {
		// hybridE is infinite, thus overall energy is infinity as well
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
__m128
InteractionEnergy::
getE_SSE( const __m128i i1, const __m128i j1
		, const __m128i i2, const __m128i j2
		, const __m128 hybridE ) const
{

	__m128 result = {0.0, 0.0, 0.0, 0.0};
	__m128 resultMultiplication = {0.0, 0.0, 0.0, 0.0};
	__m128 inf = {E_INF, E_INF, E_INF, E_INF};
	__m128 ed1 = {getED1( i1[0], j1[0] ), getED1( i1[1], j1[1] ), getED1( i1[2], j1[2] ), getED1( i1[3], j1[3] )};
	__m128 ed2 = {getED2( i2[0], j2[0] ), getED2( i2[1], j2[1] ), getED2( i2[2], j2[2] ), getED2( i2[3], j2[3] )};
	__m128 E_danglingLeft = {getE_danglingLeft( i1[0], i2[0] ), getE_danglingLeft( i1[1], i2[1] ), getE_danglingLeft( i1[2], i2[2] ), getE_danglingLeft( i1[3], i2[3] )};
	__m128 Pr_danglingLeft = {getPr_danglingLeft( i1[0], j1[0], i2[0], j2[0]), getPr_danglingLeft( i1[1], j1[1], i2[1], j2[1] ), getPr_danglingLeft( i1[2], j1[2] , i2[2], j2[2]), getPr_danglingLeft( i1[3], j1[3], i2[3], j2[3] )};
	__m128 E_danglingRight = {getE_danglingRight( j1[0], j2[0] ), getE_danglingRight( j1[1], j2[1] ), getE_danglingRight( j1[2], j2[2] ), getE_danglingRight( j1[3], j2[3] )};
	__m128 Pr_danglingRight = {getPr_danglingRight(  i1[0], j1[0], i2[0], j2[0] ), getPr_danglingRight(  i1[1], j1[1], i2[1], j2[1]), getPr_danglingRight(  i1[2], j1[2], i2[2], j2[2]), getPr_danglingRight(  i1[3], j1[3], i2[3], j2[3])};
	__m128 E_endLeft = {getE_endLeft( i1[0], i2[0] ), getE_endLeft( i1[1], i2[1] ), getE_endLeft( i1[2], i2[2] ), getE_endLeft( i1[3], i2[3] )};
	__m128 E_endRight = {getE_endRight( j1[0], j2[0] ), getE_endRight( j1[1], j2[1] ), getE_endRight( j1[2], j2[2] ), getE_endRight( j1[3], j2[3] )};
	
	result = _mm_add_ps(_mm_add_ps(_mm_add_ps(ed1, ed2), E_endLeft), E_endRight);
	resultMultiplication = _mm_add_ps(_mm_mul_ps(E_danglingLeft, Pr_danglingLeft), _mm_mul_ps(E_danglingRight, Pr_danglingRight));
	result = _mm_add_ps(result, resultMultiplication);
	__m128 infCompare = _mm_cmpeq_ps(hybridE, inf);
	return _mm_or_ps(_mm_and_ps(infCompare,inf), _mm_andnot_ps(infCompare,result));
	// return result;
	// _mm_mul_ps(boltzmannWeight_SSE, hybridZ_SSE));
	// // check if hybridization energy is not infinite
	// if ( E_isNotINF(hybridE) ) {
	// 	// compute overall interaction energy
	// 	// std::cout << "foo" << std::endl;
	// 	return hybridE
	// 			// accessibility penalty
	// 			+ getED1( i1, j1 )
	// 			+ getED2( i2, j2 )
	// 			// dangling end penalty
	// 			// weighted by the probability that ends are unpaired
	// 			+ (getE_danglingLeft( i1, i2 )*getPr_danglingLeft(i1,j1,i2,j2))
	// 			+ (getE_danglingRight( j1, j2 )*getPr_danglingRight(i1,j1,i2,j2))
	// 			// helix closure penalty
	// 			+ getE_endLeft( i1, i2 )
	// 			+ getE_endRight( j1, j2 )
	// 			;
	// } else {
	// 	// hybridE is infinite, thus overall energy is infinity as well
	// 	return E_INF;
	// }
}



#endif /* INTERACTIONENERGY_H_ */
