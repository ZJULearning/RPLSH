// Copyright (C) 2018 Deng Cai <dengcai@gmail.com>. All Rights Reserved.
#ifndef EFANNA
#define EFANNA
#include "general/distance.hpp"
#include "general/matrix.hpp"
#include "general/params.hpp"
#include "algorithm/init_indices.hpp"
namespace efanna{
template <typename DataType>
class FIndex{
public:
	typedef InitIndex<DataType> IndexType;

	FIndex(const Matrix<DataType>& features, Distance* d, const IndexParams& params)
	: index_params_(params)
	{
		init_algorithm init_index_type= params.init_index_type;
		index_params_ = params;
		initIndex_ = create_index_by_type(init_index_type, features, params, d);
	}

	virtual ~FIndex () {
	}
	void buildIndex(float* data, size_t points_num, int dim)
	{
		initIndex_->buildIndex(data, points_num, dim);
	}
	void saveIndex(char* filename){
		initIndex_->saveIndex(filename);
	}
	void loadIndex(char* filename){
		initIndex_->loadIndex(filename);
	}
	void setSearchParams(int tables, int init_num){
		initIndex_->setSearchParams(tables, init_num);
	}
	void knnSearch(int k, const Matrix<DataType>& query, float* query_){
		initIndex_->knnSearch(k, query, query_);
	}
	void saveResults(char* filename){
		initIndex_->saveResults(filename);
	}
private:
	/** Pointer to actual index class */
	IndexType* initIndex_;
	/** Parameters passed to the index */
	IndexParams index_params_;

};
}
#endif
