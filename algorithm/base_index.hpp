#ifndef EFANNA_BASE_INDEX_H_
#define EFANNA_BASE_INDEX_H_

#include "general/params.hpp"
#include "general/distance.hpp"
#include "general/matrix.hpp"
#include <fstream>
#include <iostream>


namespace efanna{

template <typename DataType>
class InitIndex{
public:
	InitIndex(const Matrix<DataType>& features, const Distance* d, const IndexParams& params):
		features_(features),
		distance_(d),
		params_(params)
{
}
	virtual ~InitIndex() {};

	virtual void buildIndex(float* data, size_t points_num, int dim)
	{
		buildIndexImpl(data, points_num, dim);
	}
	virtual void buildIndexImpl(float* data, size_t points_num, int dim) = 0;

	virtual void loadIndex(char* filename) = 0;
	virtual void saveIndex(char* filename) = 0;

	void saveResults(char* filename){
		std::ofstream out(filename,std::ios::binary);
		std::vector<std::vector<int>>::iterator i;
		for(i = nn_results.begin(); i!= nn_results.end(); i++){
			std::vector<int>::iterator j;
			int dim = i->size();
			out.write((char*)&dim, sizeof(int));
			for(j = i->begin(); j != i->end(); j++){
				int id = *j;
				out.write((char*)&id, sizeof(int));
			}
		}
		out.close();
	}
	SearchParams SP;
	void setSearchParams(int nGroup, int init_num){
		SP.search_groups = nGroup;
		SP.search_init_num = init_num;
	}


	virtual void knnSearch(int K, const Matrix<DataType>& query, float* query_){
		getNeighbors(K,query,query_);
	}

	virtual void getNeighbors(size_t K, const Matrix<DataType>& query, float* query_) = 0;



protected:
	const Matrix<DataType> features_;
	const Distance* distance_;
	const IndexParams params_;
	std::vector<std::vector<int> > nn_results;
};
#define USING_BASECLASS_SYMBOLS \
		using InitIndex<DataType>::features_;\
		using InitIndex<DataType>::distance_;\
		using InitIndex<DataType>::params_;\
		using InitIndex<DataType>::nn_results;\
		using InitIndex<DataType>::buildIndex;\
		using InitIndex<DataType>::saveResults;\
		using InitIndex<DataType>::SP;\

}
#endif
