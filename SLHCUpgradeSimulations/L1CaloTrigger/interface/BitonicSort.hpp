#ifndef BitonicSort_h
#define BitonicSort_h

#include <vector>

// <TEMP>
    #include "TString.h"
// </TEMP>

enum sort_direction {up, down};

// ***********************************************************************************************
// *                               Functions for the bitonic sort                                *
// ***********************************************************************************************

template <typename T>
void    BitonicSort( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd);

template <typename T>
void   BitonicMerge( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd);

template <typename T>
void BitonicCompare( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd);


// Bitonic sort routine
template <typename T>
void BitonicSort( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd){
  
	uint32_t lSize( aDataEnd - aDataStart );

	if( lSize != 1 ){
		typename std::vector<T>::iterator lMidpoint( aDataStart+(lSize>>1) );
		BitonicSort<T> ( up, aDataStart , lMidpoint );
		BitonicSort<T> ( down, lMidpoint , aDataEnd );
		BitonicMerge<T> (aDir, aDataStart , aDataEnd );
	}

}

// Bitonic merge
template <typename T>
void BitonicMerge( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd){

	uint32_t lSize( aDataEnd - aDataStart );

	if( lSize != 1 ){
		BitonicCompare<T> ( aDir, aDataStart , aDataEnd );
		typename std::vector<T>::iterator lMidpoint( aDataStart+(lSize>>1) );
		BitonicMerge<T> ( aDir, aDataStart , lMidpoint );
		BitonicMerge<T> ( aDir, lMidpoint , aDataEnd );
	}

}


// Comparison routine
template <typename T>
void BitonicCompare( sort_direction aDir, typename std::vector<T>::iterator & aDataStart, typename std::vector<T>::iterator & aDataEnd){

	uint32_t lSize( aDataEnd - aDataStart );
	
	typename std::vector<T>::iterator lFirst( aDataStart );
	typename std::vector<T>::iterator lSecond( aDataStart+(lSize>>1) );

	for( ; lSecond != aDataEnd ; ++lFirst , ++lSecond ){


		bool lComp( *lFirst > *lSecond );
		if( ( lComp && (aDir == up) ) || ( !lComp && (aDir == down) ) ) {
			std::swap( *lFirst , *lSecond );
		}
	}

}




#endif
